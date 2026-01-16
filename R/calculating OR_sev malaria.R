# Data 
y_tx   <- 0 # Number of events in treatment arm ('successes')
y_ctrl <- 5 # Number of events in control arm ('failures')
dat <- data.frame(y_tx = y_tx, y_ctrl = y_ctrl)


fit <- stan_glm(
  cbind(y_tx, y_ctrl) ~ 1,
  family = binomial(link = "logit"), # Severe malaria treated as a binary outcome (assuming rare recurrence)
  data = dat,
  prior_intercept = normal(0, 1),  # θ = log(OR) ~ Normal(0, 1), implying a prior median OR of 1
  chains = 4, 
  iter = 5000, 
  seed = 1, 
  refresh = 0
              )

# Posterior draws of θ = log(OR)
th <- as.vector(as.matrix(fit, pars = "(Intercept)"))


OR <- exp(th)
    output <-  c(
      #Log-odds ratios
      theta_mean = mean(th),
      theta_median = median(th),
      theta_95CrI_lo = quantile(th, .025),
      theta_95CrI_hi = quantile(th, .975),
      #Odds ratios
      OR_mean = mean(OR),
      OR_median = median(OR),
      OR_95CrI_lo = quantile(OR, .025),
      OR_95CrI_hi = quantile(OR, .975),
      # Posterior probability that the odds ratio is less than 1
      Pr_OR_less_1 = mean(OR < 1)
    )
    
    output


