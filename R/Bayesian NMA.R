# ============================================================
# Bayesian Network Meta-Analysis (NMA) for SMC
# Author: Dhruv Darji
# ============================================================

# ---------------------------
# Setup
# ---------------------------


library(here)
library(gemtc)
library(rjags)
library(dplyr)
library(tidyr)
library(tibble)
library(magrittr)
library(ggplot2)
library(scales)
library(svglite)
library(dmetar)

set.seed(1234)


# ---------------------------
# Load data
# ---------------------------

smc_path <- here("data", "smc_nma_cleaned.csv")
dot_path <- here("data", "smc_nma_adherence.csv")


smc_data <- read.csv(smc_path, stringsAsFactors = FALSE)


# ---------------------------
# Harmonize treatment labels
# ---------------------------

# Combine no drug and placebo into a single category
# Combine SPAS and SP (monthly) into a single category
smc_data <- smc_data %>%
  mutate(
    trt1 = recode(trt1, "nodrug" = "placebo"),
    trt2 = recode(trt2, "nodrug" = "placebo"),
    
    trt1 = recode(trt1, "placebo" = "placebo/nodrug"),
    trt2 = recode(trt2, "placebo" = "placebo/nodrug"),
    
    trt1 = recode(trt1, "spas" = "spas/sp_m", "sp_m" = "spas/sp_m"),
    trt2 = recode(trt2, "spas" = "spas/sp_m", "sp_m" = "spas/sp_m")
  )

# ---------------------------
# Reshape to GEMTC "relative effects" format
# ---------------------------
# gemtc contrast-based format (data.re) expects:
#   study, treatment, diff, std.err
# where baseline/reference arm rows have diff = NA and std.err = baseline SE


smc_data_long <- smc_data %>% 
  pivot_longer( #Reshape from wide to long
    cols = c("trt1", "trt2", "logte", "sete"), # Select columns to reshape
    names_to = c(".value"), # Use column name parts as output column names 
    names_pattern = "(..)"  # Group column names by first two columns (trt1 & trt2 == grouped into treatments)
  ) %>% 
  rename( #Rename columns
    diff = lo,
    std.err = se,
    treatment = tr
  ) %>% 
  arrange(study)
    # smc_data_long has data in long format, for each study:
    # - All possible pairwise comparisons are listed (for a 2-arm study = 1 comparison; 3-arm study = 3 comparisons; 4-arm study = 6 comparisons)
    # - Each pairwise comparison is represented in 2 rows, the comparator's row values of the logTE (diff) and standard error of TE, the reference arm below has 'NA'
    # - For GEMTC, each study must have a common reference arm (see handling of multi-arm studies below)

# ---------------------------
# Multi-arm cleaning 
# ---------------------------
# Multi-arm studies are internally consistent i.e. In a study of A vs B vs C, if we know the effect of A vs B and B vs C, then we can estimate A vs C.
# - In the current format, smc_data_long all possible pairwise comparisons for multi-arm studies - this is not required!
# - For multi-arm studies, GEMTC requires pairwise comparisons against a single reference treatment (arbitrarily selected)
# - In a study of A vs B vs C, if we select C as the reference, we only need to supply A vs C and B vs C
# - Hence for multi-arm studies below, I select a reference treatment and drop all pairwise comparisons that do not contain this reference
# - For most multi-arm studies with placebo/nodrug, I (arbitrarily) selected placebo/nodrug and dropped pairwise comparison rows that do not compare against placebo/nodrug in that study.
# - E.g. for Nuwa 2025, I drop the dp vs spaq comparison (rows 3 and 4) and keep spaq vs placebo/nodrug and dp vs placebo/nodrug
# - If a study did not have placebo/nodrug (e.g. Kweku 2008), I (arbitrarily) selected another reference arm, and dropped accordingly



# Note - indices may no longer correspond to the same comparisons, if data changes so always check dataframe using the logic above.

smc_data_long$row_id <- seq_len(nrow(smc_data_long))

drop_rows <- c(3, 4, 9, 10, 13, 14, 37, 38,
               45, 46, 51, 52, 53, 54, 57, 58)

smc_data_long <- smc_data_long %>%
  filter(!(row_id %in% drop_rows)) %>%
  select(-row_id)

# Remove duplicated baseline rows within study x treatment
smc_data_long <- smc_data_long %>%
  group_by(studyid, treatment) %>%
  filter(!(is.na(diff) & duplicated(treatment))) %>%
  ungroup()

# ---------------------------
# Insert baseline SE for reference arms for multi-arm studies
# ---------------------------
baseline_SE <- tibble(
  studyid = c("Nuwa2025",
              "Traore2024",
              "Chandramohan2021",
              "Bojang2010",
              "Sokhna2008",
              "Kweku2008"),
  baseline_se = c(0.055727821,
                  0.113960576,
                  0.056888012,
                  0.5,
                  0.136082763,
                  0.073922127)
)         #Calculated as SE(λ) = 1/sqrt(y); where y = No. events (counts) in arm 

smc_data_long <- smc_data_long %>%
  left_join(baseline_SE, by = "studyid") %>%
  mutate(
    std.err = if_else(is.na(diff) & !is.na(baseline_se), baseline_se, std.err)
  ) %>%
  select(-baseline_se)

# ---------------------------
# Uniform treatment IDs
# ---------------------------

smc_data_long <- smc_data_long %>%
  mutate(
    treatment = recode(treatment,
                       "placebo/nodrug" = "placebo_nodrug",
                       "spas/sp_m"      = "spas_sp_m",
                       "spaq+rtss"      = "spaq_rtss",
                       .default = treatment
    )
  )

#Treatment codes for plots
treat.codes <- c(
  "spaq"           = "SPAQ",
  "placebo_nodrug" = "Placebo/No drug control",
  "dp"             = "DHAPQ",
  "rtss"           = "RTS,S vaccine",
  "spaq_rtss"      = "SPAQ + RTS,S Vaccine",
  "asaq"           = "ASAQ",
  "sppq"           = "SPPQ",
  "spas_sp_m"      = "SPAS/SP (monthly)",
  "sp"             = "SP (bimonthly)",
  "asaq_bi"        = "ASAQ (bimonthly)"
) %>%
  data.frame(description = ., stringsAsFactors = FALSE) %>%
  rownames_to_column("id")

# ---------------------------
# Build network + plot
# ---------------------------

smc_network <- mtc.network(
  data.re    = smc_data_long,
  treatments = treat.codes
)

summary(smc_network, use.description = TRUE)

svglite(filename = file.path(out_dir, "network_plot.svg"), width = 11, height = 5.825)
plot(
  smc_network,
  use.description = TRUE,
  vertex.color = "yellow",
  vertex.label.color = "black",
  vertex.label.family = "Arial",
  vertex.label.dist = 2.1,
  vertex.label.cex = 0.8
)
dev.off()

# ---------------------------
# Fit NMA models
# ---------------------------

# Common effect (fixed)
smc_model_ce <- mtc.model(
  smc_network,
  likelihood  = "normal",
  link        = "identity",
  linearModel = "fixed",
  n.chain     = 4
)

mcmc_ce_quick <- mtc.run(smc_model_ce, n.adapt = 50,   n.iter = 1000,   thin = 10)
mcmc_ce_main  <- mtc.run(smc_model_ce, n.adapt = 5000, n.iter = 100000, thin = 10)

# Random effects
smc_model_re <- mtc.model(
  smc_network,
  likelihood  = "normal",
  link        = "identity",
  linearModel = "random",
  n.chain     = 4
)

mcmc_re_quick <- mtc.run(smc_model_re, n.adapt = 50,   n.iter = 1000,   thin = 10)
mcmc_re_main  <- mtc.run(smc_model_re, n.adapt = 5000, n.iter = 100000, thin = 10)

# Diagnostics
plot(mcmc_ce_quick)
plot(mcmc_ce_main)
gelman.plot(mcmc_ce_main)
gelman.diag(mcmc_ce_main)$mpsrf

summary(mcmc_ce_main)
summary(mcmc_re_main)

# ---------------------------
# Relative effects + forest plot
# ---------------------------

results_vs_placebo <- relative.effect(mcmc_ce_main, t1 = "placebo_nodrug")
results_vs_spaq    <- relative.effect(mcmc_ce_main, t1 = "spaq")

summary(results_vs_placebo)
summary(results_vs_spaq)

forest(results_vs_placebo, use.description = TRUE)

# ---------------------------
# Heterogeneity
# ---------------------------

smc_anohe <- mtc.anohe(
  smc_network,
  factor      = 2.5,
  n.chain     = 4,
  likelihood  = "normal",
  link        = "identity",
  linearModel = "fixed",
  n.adapt     = 5000,
  n.iter      = 100000,
  thin        = 10
)
summary(smc_anohe)

# ---------------------------
# Inconsistency
# ---------------------------

# Node-splitting candidates
mtc.nodesplit.comparisons(smc_network)

# Node-splitting (common effect)
nodesplit_ce <- mtc.nodesplit(
  smc_network,
  linearModel = "fixed",
  likelihood  = "normal",
  link        = "identity",
  n.adapt     = 5000,
  n.iter      = 100000,
  thin        = 10
)
nodesplit_ce_sum <- summary(nodesplit_ce)
plot(nodesplit_ce_sum)

# Node-splitting (random effects)
nodesplit_re <- mtc.nodesplit(
  smc_network,
  linearModel = "random",
  likelihood  = "normal",
  link        = "identity",
  n.adapt     = 5000,
  n.iter      = 100000,
  thin        = 10
)
plot(summary(nodesplit_re))

# Export node-splitting estimates (common effect)
df_nodesplit <- cbind(
  as.data.frame(nodesplit_ce_sum$dir.effect,  stringsAsFactors = FALSE),
  as.data.frame(nodesplit_ce_sum$ind.effect,  stringsAsFactors = FALSE),
  as.data.frame(nodesplit_ce_sum$cons.effect, stringsAsFactors = FALSE),
  as.data.frame(nodesplit_ce_sum$p.value,     stringsAsFactors = FALSE)
)

write.csv(df_nodesplit, file.path(out_dir, "nodesplit_output.csv"), row.names = FALSE)

# ---------------------------
# Treatment Ranking
# ---------------------------

rank_probs <- rank.probability(mcmc_ce_main, preferredDirection = -1)
plot(rank_probs, beside = TRUE)

rank_probs_all <- rank.probability(mcmc_ce_main)
sucra_obj <- dmetar::sucra(rank_probs_all, lower.is.better = TRUE)
print(sucra_obj)

svglite(filename = file.path(out_dir, "sucra.svg"), width = 10, height = 5.625)
plot(sucra_obj, ylab = "Cumulative rank probability")
dev.off()

# ---------------------------
# League table (relative effects for all treatments)
# ---------------------------

league_table <- relative.effect.table(mcmc_ce_main)
write.csv(league_table, file.path(out_dir, "smc_ltable_uncompmalaria.csv"), row.names = FALSE)

# ---------------------------
# Meta-regression (Adherence / DOT)
# ---------------------------

  smc_nma_dot <- read.csv(dot_path, header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(dot = ifelse(dot == "Fully supervised", 1, 0))
  
  network_ad <- mtc.network(
    data.re    = smc_data_long,
    studies    = smc_nma_dot,
    treatments = treat.codes
  )
  
  adherence <- list(
    coefficient = "shared",
    variable    = "dot",
    control     = "placebo_nodrug"
  )
  
  mr_ad <- mtc.model(
    network_ad,
    likelihood  = "normal",
    link        = "identity",
    type        = "regression",
    linearModel = "random",
    regressor   = adherence
  )
  
  mcmc_ad <- mtc.run(mr_ad, n.adapt = 5000, n.iter = 10000, thin = 10)
  summary(mcmc_ad)
  
  forest(relative.effect(mcmc_ad, t1 = "placebo_nodrug", covariate = 1),
         use.description = TRUE)
  title("Fully supervised")
  
  forest(relative.effect(mcmc_ad, t1 = "placebo_nodrug", covariate = 0),
         use.description = TRUE)
  title("Partially supervised")
  

