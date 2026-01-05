// ============================================================
// Contrast-based Bayesian NMA (Common-Effect) with
// study-level meta-regression (shared slope across contrasts)
//
// NOTES
// - This model treats contrasts as conditionally independent given
//   their reported SEs, even when contrasts come from multi-arm
//   studies. (Sigma is diagonal within study blocks.)
// - Treatment effects are coded vs an implicit reference treatment
//   with effect fixed to 0 (so parameters are d[1:(T-1)]).
//
// DATA 
// - y[k] is the observed contrast (log IRR) for contrast k
// - se[k] is the reported standard error for y[k]
// - C[k,] maps treatment effects d into contrast k
// - X[s,] are study-level covariates 
// - study_id[k] identifies the study for contrast k
// - m[s], start[s] define contiguous blocks of contrasts per study
//   in y/se/C/study_id
// ============================================================

data {
  // ----------------------------
  // Dimensions
  // ----------------------------
  int<lower=2> T;                        // number of treatments (including reference)
  int<lower=1> S;                        // number of studies
  int<lower=1> K;                        // total number of contrasts across all studies
  int<lower=0> M;                        // number of study-level covariates (can be 0)

  // ----------------------------
  // Contrast-level data
  // ----------------------------
  vector[K] y;                           // observed contrast (logIRR (A vs B))
  vector<lower=0>[K] se;                 // SE of observed contrast
  array[K] int<lower=1, upper=S> study_id;

  // Design matrix mapping treatment effects to contrasts:
  // C is K x (T-1), d is length (T-1)
  matrix[K, T - 1] C;

  // Study-level covariates (S x M); if M=0 this is S x 0
  matrix[S, M] X;

  // ----------------------------
  // Ragged structure (per-study blocks)
  // ----------------------------
  array[S] int<lower=1> m;               // number of contrasts in study s
  array[S] int<lower=1> start;           // 1-based start index for study s in y/se/etc
}

parameters {
  vector[T - 1] d;                       // treatment effects vs reference
  vector[M] beta;                        // shared covariate slopes (only if M > 0)
}

transformed parameters {
  vector[K] mu;                          // linear predictor for each contrast

  for (k in 1:K) {
    real trt_part = dot_product(C[k], d);

    // Guard for M == 0 (empty covariate matrix)
    real cov_part = (M > 0)
      ? dot_product(to_row_vector(X[study_id[k]]), beta)
      : 0.0;

    mu[k] = trt_part + cov_part;
  }
}

model {
  // ----------------------------
  // Priors
  // ----------------------------
  d    ~ normal(0, 1.5);
  beta ~ normal(0, 1);

  // ----------------------------
  // Likelihood
  // ----------------------------
  // Within each study, the contrasts are modeled as multivariate normal with
  // diagonal covariance (conditional independence given SEs).
  for (s in 1:S) {
    int ms = m[s];
    int st = start[s];

    vector[ms] y_s  = segment(y,  st, ms);
    vector[ms] mu_s = segment(mu, st, ms);
    vector[ms] se_s = segment(se, st, ms);

    // Diagonal covariance: Sigma = diag(se^2)
    // => Cholesky factor is simply diag(se)
    matrix[ms, ms] L = diag_matrix(se_s);

    y_s ~ multi_normal_cholesky(mu_s, L);
  }
}

generated quantities {
  // Future extensions:
  // - posterior predictive checks
  // - pairwise treatment effects
  // - ranking probabilities / SUCRA 
}
