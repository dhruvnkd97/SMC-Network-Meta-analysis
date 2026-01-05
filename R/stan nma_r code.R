# Contrast-based Bayesian NMA in Stan (cmdstanr)

library(here)
library(dplyr)
library(cmdstanr)
library(posterior)

set.seed(123)

# Paths
data_path <- here("data", "smc_nma_cleaned.csv")
stan_path <- here("stan", "nma_contrast.stan")
out_dir   <- here("outputs")

# Reference treatment (must match CSV)
REF <- "placebo/nodrug"

# Keep baseline comparisons only (row indices depend on CSV ordering)
drop_rows <- c(2, 5, 7, 19, 23, 26, 28, 29)

# ---- Load + filter ----
smc_data <- read.csv(data_path, stringsAsFactors = FALSE) %>%
  mutate(
    logte = as.numeric(logte),
    sete  = as.numeric(sete)
  ) %>%
  mutate(.row_id = row_number()) %>%
  filter(!(.row_id %in% drop_rows)) %>%
  select(-.row_id)

# ---- Index studies ----
studies <- smc_data %>%
  distinct(studyid) %>%
  arrange(studyid) %>%
  mutate(s_id = row_number())

S <- nrow(studies)
s_index <- setNames(studies$s_id, studies$studyid)

smc_data <- smc_data %>%
  mutate(s_id = unname(s_index[studyid])) %>%
  arrange(s_id)

# ---- Index treatments (REF first) ----
treatments_all <- sort(unique(c(smc_data$trt1, smc_data$trt2)))
treatments <- c(REF, setdiff(treatments_all, REF))
Tn <- length(treatments)
t_index <- setNames(seq_len(Tn), treatments)

# ---- Design matrix C ----
K <- nrow(smc_data)
C <- matrix(0, nrow = K, ncol = Tn - 1)

col_of <- function(trt) {
  j <- t_index[[trt]]
  if (is.null(j) || j == 1) return(0L) else return(j - 1L)
}

for (k in seq_len(K)) {
  a <- smc_data$trt1[k]
  b <- smc_data$trt2[k]
  ja <- col_of(a)
  jb <- col_of(b)
  if (ja > 0) C[k, ja] <- C[k, ja] + 1
  if (jb > 0) C[k, jb] <- C[k, jb] - 1
}

# ---- Ragged blocks ----
y        <- smc_data$logte
se       <- smc_data$sete
study_id <- smc_data$s_id

m <- as.integer(tabulate(study_id, nbins = S))
start <- c(1L, 1L + cumsum(m)[1:(S - 1)])

# ---- No covariates ----
M <- 0L
X <- matrix(numeric(0), nrow = S, ncol = 0)

stan_data <- list(
  T = Tn, S = S, K = K, M = M,
  y = as.vector(y),
  se = as.vector(se),
  study_id = as.array(as.integer(study_id)),
  C = C,
  X = X,
  m = as.array(m),
  start = as.array(start)
)

# ---- Fit ----
mod <- cmdstan_model(stan_path)

fit_quick <- mod$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4,
  iter_warmup = 1000, iter_sampling = 2000,
  adapt_delta = 0.99, max_treedepth = 15
)

fit_main <- mod$sample(
  data = stan_data,
  chains = 4, parallel_chains = 4,
  iter_warmup = 5000, iter_sampling = 100000,
  adapt_delta = 0.9, max_treedepth = 12
)

# ---- Effects vs reference ----
lab <- treatments[-1]
d_draws <- as_draws_matrix(fit_main$draws(variables = "d"))
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

write.csv(res_vs_ref, file.path(out_dir, "stan_results_vs_ref.csv"), row.names = FALSE)

# ---- Ranking + SUCRA ----
d_full <- cbind(REF = 0, d_draws)
d_full <- d_full[, treatments, drop = FALSE]

score <- d_full  # use -d_full if higher is better

ranks <- t(apply(score, 1, function(row) rank(row, ties.method = "average")))
colnames(ranks) <- treatments

rank_probs <- sapply(1:Tn, function(r) colMeans(ranks == r))
colnames(rank_probs) <- paste0("Rank", 1:Tn)
rownames(rank_probs) <- treatments

mean_rank <- as.numeric(rank_probs %*% (1:Tn))
names(mean_rank) <- treatments

cum_rank_probs <- t(apply(rank_probs, 1, cumsum))
SUCRA <- rowSums(1 - cum_rank_probs[, -ncol(cum_rank_probs), drop = FALSE]) / (Tn - 1)

summary_ranks <- data.frame(
  treatment = treatments,
  prob_best = rank_probs[, "Rank1"],
  mean_rank = mean_rank,
  SUCRA = SUCRA,
  row.names = NULL
) %>%
  arrange(desc(SUCRA))

write.csv(summary_ranks, file.path(out_dir, "stan_ranking_summary.csv"), row.names = FALSE)

# ---- Session info ----
writeLines(capture.output(sessionInfo()),
           con = file.path(out_dir, "sessionInfo_stan.txt"))
