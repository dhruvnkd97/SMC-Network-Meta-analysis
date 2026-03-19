library(rstanarm)
library(meta)



    #Bayesian logistic regression method
    tine <- data.frame(
      events = c(0,5),
      n = c(500,500),
      trt = c(1,0)
    )
    
    
    dicko <- data.frame(
      events = c(0,5),
      n = c(131,131),
      trt = c(1,0)
    )
    
    
    fit <- stan_glm(
      cbind(events, n - events) ~ trt,
      family = binomial(link = "logit"),
      data = tine, #replace with dicko
      prior = normal(0,1),
      prior_intercept = normal(0,5),
      chains = 4,
      iter = 5000,
      seed = 1
      )

    
    beta <- as.vector(as.matrix(fit, pars="trt"))
    
    OR <- exp(beta)
    sd(beta) #Posterior SD of log OR
    
    output <- c(
      theta_mean = mean(beta),
      theta_median = median(beta),
      theta_95CrI_lo = quantile(beta,.025),
      theta_95CrI_hi = quantile(beta,.975),
      OR_mean = mean(OR),
      OR_median = median(OR),
      OR_95CrI_lo = quantile(OR,.025),
      OR_95CrI_hi = quantile(OR,.975),
     Pr_OR_less_1 = mean(OR < 1)
    )

    output
    
    
    
    #Peto method
    t1 <- 0; n1 <- 500
    t2 <- 5; n2 <- 500
    
    m <- metabin(event.e = t1, n.e = n1,
                 event.c = t2, n.c = n2,
                 sm = "OR", method = "Peto",
                 common = TRUE)
    
    m
    
    
    

    

    
    