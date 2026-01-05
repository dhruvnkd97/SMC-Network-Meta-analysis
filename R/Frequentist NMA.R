# ============================================================
# SMC Network Meta-Analysis (Frequentist) — netmeta
# ============================================================

library(here)
library(readr)
library(dplyr)
library(meta)
library(netmeta)
library(dmetar)

set.seed(123)

# Paths
nma_path      <- here("data", "smc_nma_cleaned.csv")
out_dir       <- here("outputs")

# ------------------------------------------------------------
# Frequentist NMA (contrast-based; exploratory)
# ------------------------------------------------------------

smc_data_freq <- read_csv(nma_path, show_col_types = FALSE) %>%
  mutate(
    trt1 = recode(trt1, "nodrug" = "placebo"),
    trt2 = recode(trt2, "nodrug" = "placebo"),
    trt1 = recode(trt1, "placebo" = "placebo/nodrug"),
    trt2 = recode(trt2, "placebo" = "placebo/nodrug"),
    trt1 = recode(trt1, "spas" = "spas/sp_m", "sp_m" = "spas/sp_m"),
    trt2 = recode(trt2, "spas" = "spas/sp_m", "sp_m" = "spas/sp_m")
  ) %>%
  rename(logTE = logte, seTE = sete)

# Optional: check multi-arm structure
table(smc_data_freq$studyid)


# Fit NMA model (common-effect)
smc_network_freq <- netmeta(
  TE = logTE,
  seTE = seTE,
  treat1 = trt1,
  treat2 = trt2,
  studlab = studyid,
  data = smc_data_freq,
  sm = "IRR",
  common = TRUE,
  reference.group = "placebo/nodrug",
  details.chkmultiarm = TRUE,
  tol.multiarm = 0.1,
  tol.multiarm.se = 0.0001,
  sep.trts = " vs "
)

summary(smc_network_freq)

# Design-by-treatment decomposition (heterogeneity/inconsistency)
smc_nma_freqhet <- decomp.design(smc_network_freq)
smc_nma_freqhet

# ------------------------------------------------------------
# Plots
# ------------------------------------------------------------

long_labels <- c(
  "ASAQ",
  "ASAQ (bimonthly)",
  "DHAPQ",
  "Placebo/No drug",
  "Seasonal RTS,S",
  "SP (bimonthly)",
  "SPAQ",
  "Seasonal SPAQ+RTS,S",
  "SPAS/SP (monthly)",
  "SPPQ"
)

# Network graph
netgraph(smc_network_freq, labels = long_labels)

# Direct evidence plot
smc_devidence <- direct.evidence.plot(smc_network_freq)
plot(smc_devidence)

# Heat plot
netheat(smc_network_freq)

# ------------------------------------------------------------
# Effect tables
# ------------------------------------------------------------

# League table (random effects results printed by netleague; export random table)
netleague_out <- netleague(smc_network_freq, bracket = "(", digits = 2)
write.csv(netleague_out$random, file.path(out_dir, "nma_effecttable.csv"), row.names = FALSE)

# ------------------------------------------------------------
# Ranking
# ------------------------------------------------------------

netrank(smc_network_freq, small.values = "good")

# ------------------------------------------------------------
# Forest plot
# ------------------------------------------------------------

forest(
  smc_network_freq,
  reference.group = "placebo/nodrug",
  pooled = "common",
  sortvar = TE,
  smlab = "SMC vs placebo",
  test.overall.common = TRUE,
  overall = TRUE,
  xlim = c(0.05, 10),
  hetstat = "common",
  label.left = "Favors SMC",
  label.right = "Favors placebo",
  labels = long_labels,
  fontsize = 11
)

# ------------------------------------------------------------
# Node-splitting / inconsistency checks
# ------------------------------------------------------------

smc_netsplit <- netsplit(smc_network_freq)
smc_netsplit

# Forest of split estimates (may be verbose depending on network)
netsplit(smc_network_freq) %>% forest()

netsplit_placebo <- netsplit(
  smc_network_freq,
  method = "back-calculation",
  reference.group = "placebo/nodrug",
  baseline.reference = TRUE,
  common = TRUE,
  backtransf = TRUE,
  ci = TRUE
)

netsplit_spaq <- netsplit(
  smc_network_freq,
  method = "back-calculation",
  reference.group = "spaq",
  baseline.reference = TRUE,
  backtransf = TRUE,
  ci = FALSE
)

netsplit_placebo
netsplit_spaq

# ------------------------------------------------------------
# Publication bias (funnel plot)
# ------------------------------------------------------------

pch18 <- c(16, 17, 15, 18, 3, 4, 8, 1, 2, 0, 5, 6, 7, 9, 10, 11, 12, 13)

funnel(
  smc_network_freq,
  order = c("placebo/nodrug", "sp", "asaq_bi", "asaq", "spas", "rtss",
            "spaq", "dp", "sppq", "spaq+rtss"),
  method.bias = "Egger",
  pooled = "random",
  lump.comparator = FALSE,
  legend = TRUE,
  pch = pch18,
  level = 0.5,
  studlab = TRUE,
  cex.studlab = 0.7,
  pos.studlab = 3
)
