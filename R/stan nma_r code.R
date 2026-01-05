# ============================================================
# Contrast-based Bayesian NMA in Stan (cmdstanr)
#
# Repo assumptions:
# - data/ contains:
#     1) smc_nma_cleaned2.csv
#     2) smc_nma_adherence.csv (not used here; this is the Stan NMA runner)
# - outputs/ exists (this script will write results there)
# - Stan model file exists (default: stan/nma_contrast.stan)
#
# IMPORTANT NOTE ON MULTI-ARM TRIALS:
# - The provided Stan likelihood uses a *diagonal* within-study covariance.
# - This treats contrasts as conditionally independent given their SEs.
# - It is NOT a full multi-arm correlation correction unless you add off-diagonals.
# ============================================================

# ---------------------------
# 0) Setup
# ---------------------------

library(here)
library(dplyr)
library(cmdstanr)
library(posterior)

set.seed(123)

# Paths 
data_path <- here("data", "smc_nma_cleaned.csv")
stan_path <- here("stan", "nma_contrast.stan")   # recommended location: stan/
out_dir   <- here("outputs")


# Reference treatment label MUST match your CSV exactly
REF <- "placebo/nodrug"

# Drop specific rows to keep baseline comparisons only
# Always check data indices depend on input CSV ordering!)
drop_rows <- c(2, 5, 7, 19, 23, 26, 28, 29)

# ---------------------------
# 1) Read + basic checks
# ---------------------------

smc_data <- read.csv(data_path, stringsAsFactors = FALSE)

smc_data <- smc_data %>%
  mutate(
    logte = as.numeric(logte),
    sete  = as.numeric(sete)
  )

# ---------------------------
# 2) Keep only baseline comparisons (as per your original approach)
# ---------------------------
smc_data <- smc_data %>%
  mutate(.row_id = row_number()) %>%
  filter(!(.row_id %in% drop_rows)) %>%
  select(-.row_id)

# ---------------------------
# 3) Index studies
# ---------------------------
studies <- smc_data %>%
  distinct(studyid) %>%
  arrange(studyid) %>%
  mutate(s_id = row_number())

S <- nrow(studies)
s_index <- setNames(studies$s_id, studies$studyid)

smc_data <- smc_data %>%
  mutate(s_id = unname(s_index[studyid])) %>%
  arrange(s_id)

# ---------------------------
# 4) Index treatments (ensure reference is first)
# ---------------------------
treatments_all <- sort(unique(c(smc_data$trt1, smc_data$trt2)))

treatments <- c(REF, setdiff(treatments_all, REF))
Tn <- length(treatments)
t_index <- setNames(seq_len(Tn), treatments)

# ---------------------------
# 5) Build design matrix C (K x (T-1))
# ---------------------------
# Parameterization:
# - d has length (T-1) and corresponds to treatments[-1] vs REF
# - For a contrast A vs B (y = effect A relative to B), the design row is:
#     +1 for A (if not REF)
#     -1 for B (if not REF)
#
# Your original code uses:
#   a <- trt1; b <- trt2; set +1 for a, -1 for b
# so interpret y accordingly (logIRR trt1 vs trt2). Make sure this matches your extraction.

K <- nrow(smc_data)
C <- matrix(0, nrow = K, ncol = Tn - 1)

col_of <- function(trt) {
  j <- t_index[[trt]]
  if (is.null(j) || j == 1) return(0L) else return(j - 1L)  # REF has no column
}

for (k in seq_len(K)) {
  a <- smc_data$trt1[k]
  b <- smc_data$trt2[k]
  ja <- col_of(a)
  jb <- col_of(b)
  
  if (ja > 0) C[k, ja] <- C[k, ja] + 1
  if (jb > 0) C[k, jb] <- C[k, jb] - 1
}

# ---------------------------
# 6) Build ragged (per-study) block structure
# ---------------------------
y        <- smc_data$logte
se       <- smc_data$sete
study_id <- smc_data$s_id

m <- as.integer(tabulate(study_id, nbins = S))
start <- c(1L, 1L + cumsum(m)[1:(S - 1)])

# ---------------------------
# 7) Covariates (none for now)
# ---------------------------

M <- 0L
X <- matrix(numeric(0), nrow = S, ncol = 0)

# ---------------------------
# 8) Assemble Stan data
# ---------------------------

stan_data <- list(
  T = Tn,
  S = S,
  K = K,
  M = M,
  y = as.vector(y),
  se = as.vector(se),
  study_id = as.array(as.integer(study_id)),
  C = C,
  X = X,
  m = as.array(m),
  start = as.array(start)
)

# Optional: inspect
str(stan_data, max.level = 1)

# ---------------------------
# 9) Compile + sample
# ---------------------------

mod <- cmdstan_model(stan_path)

# Quick run (sanity check)
fit_quick <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 1000,
  iter_sampling = 2000,
  adapt_delta = 0.99,
  max_treedepth = 15,
  refresh = 200
)

fit_quick$cmdstan_diagnose()

# Main run
fit_main <- mod$sample(
  data = stan_data,
  chains = 4,
  parallel_chains = 4,
  iter_warmup = 5000,
  iter_sampling = 100000,
  adapt_delta = 0.9,
  max_treedepth = 12,
  refresh = 500
)

fit_main$cmdstan_diagnose()

# ---------------------------
# 10) Diagnostics
# ---------------------------

sum_df <- fit_main$summary()
flags <- subset(sum_df, rhat > 1.01 | ess_bulk < 400)
flags

# Divergences + treedepth
sd <- fit_main$sampler_diagnostics(format = "array")
div_rate      <- mean(sd[, , "divergent__"])
div_per_chain <- colMeans(sd[, , "divergent__"])
td_hits_rate  <- mean(sd[, , "treedepth__"] >= 12)

div_rate
div_per_chain
td_hits_rate

# ---------------------------
# 11) Summaries: treatment effects vs reference
# ---------------------------

# d[1]..d[T-1] correspond to treatments[-1]
lab <- treatments[-1]

d_draws <- as_draws_matrix(fit_main$draws(variables = "d"))
# Ensure columns correspond to d[1], d[2], ...
colnames(d_draws) <- lab

qs <- function(x) c(
  med = median(x),
  l95 = unname(quantile(x, 0.025)),
  u95 = unname(quantile(x, 0.975))
)

summ_log <- t(apply(d_draws, 2, qs))
summ_irr <- exp(summ_log)

res_vs_ref <- data.frame(
  treatment = lab,
  log_med = summ_log[, "med"],
  log_l95 = summ_log[, "l95"],
  log_u95 = summ_log[, "u95"],
  IRR_med = summ_irr[, "med"],
  IRR_l95 = summ_irr[, "l95"],
  IRR_u95 = summ_irr[, "u95"],
  row.names = NULL
)

res_vs_ref

write.csv(res_vs_ref, file.path(out_dir, "stan_results_vs_ref.csv"), row.names = FALSE)

# ---------------------------
# 12) Ranking + SUCRA (computed in R from posterior draws)
# ---------------------------

# Include reference with effect 0
d_full <- cbind(REF = 0, d_draws)
# Ensure correct column order
d_full <- d_full[, treatments, drop = FALSE]

# If lower is better (e.g., log IRR where negative is good), use score = d_full
# If higher is better, use score = -d_full
score <- d_full

# ranks: 1 = best
ranks <- t(apply(score, 1, function(row) rank(row, ties.method = "average")))
colnames(ranks) <- treatments

# P(Rank = r)
rank_probs <- sapply(1:Tn, function(r) colMeans(ranks == r))
colnames(rank_probs) <- paste0("Rank", 1:Tn)
rownames(rank_probs) <- treatments

# Mean rank
mean_rank <- as.numeric(rank_probs %*% (1:Tn))
names(mean_rank) <- treatments

# SUCRA
cum_rank_probs <- t(apply(rank_probs, 1, cumsum))
SUCRA <- rowSums(1 - cum_rank_probs[, -ncol(cum_rank_probs), drop = FALSE]) / (Tn - 1)
names(SUCRA) <- treatments

summary_ranks <- data.frame(
  treatment = treatments,
  prob_best = rank_probs[, "Rank1"],
  mean_rank = mean_rank,
  SUCRA = SUCRA,
  row.names = NULL
) %>%
  arrange(desc(SUCRA))

summary_ranks

write.csv(summary_ranks, file.path(out_dir, "stan_ranking_summary.csv"), row.names = FALSE)

# ---------------------------
# 13) Save session info
# ---------------------------

