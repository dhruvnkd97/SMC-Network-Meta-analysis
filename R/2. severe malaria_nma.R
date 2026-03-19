#Title: SMC Network meta-analysis - Severe malaria
#Author: Dhruv Darji

# ============================================================
# SMC Network Meta-Analysis : SEVERE MALARIA 
# ============================================================

library(readr)
library(dplyr)
library(tidyr)
library(tibble)
library(meta)
library(netmeta)
library(dmetar)


#---- Frequentist (contrast-based) NMA ----



sm_data <- read_csv("data/2. smc_nma_severe malaria_analysis dataset.csv") 


#---- Fit NMA model (common-effect) ----
smc_sm_network <- netmeta(
  TE = logTE,
  seTE = seTE,
  treat1 = treat1,
  treat2 = treat2,
  studlab = studyid,
  data = sm_data,
  sm = "OR",
  common = TRUE,
  random = TRUE,
  reference.group = "placebo/nodrug",
  details.chkmultiarm = TRUE,
  tol.multiarm = 0.5,
  tol.multiarm.se = 1,
  sep.trts = " vs "
)

summary(smc_sm_network)


# Plots
# ------------------------------------------------------------


# ---- Network graph ----
# Treatment labels

  long_labels <- c(
    "ASAQ",
    "ASAQ (bimonthly)",
    "Placebo/No drug",
    "Seasonal RTS,S",
    "SP (bimonthly)",
    "SPAQ",
    "Seasonal SPAQ+RTS,S"
  )

  netgraph(smc_sm_network, labels = long_labels)

# ---- Direct evidence plot ---- 
sm_devidence <- direct.evidence.plot(smc_sm_network)
plot(sm_devidence)
      #For all comparisons, the total evidence is either direct or indirect.

# ---- Forest plot of relative efficacy ----
  forest(
    smc_sm_network,
    reference.group = "placebo/nodrug",
    pooled = "common",
    sortvar = TE,
    smlab = "SMC vs placebo\n(Severe malaria)",
    test.overall.common = TRUE,
    #test.overall.random = TRUE,
    #overall = TRUE,
    xlim = c(0.05, 10),
    hetstat = "common",
    #label.left = "Favors SMC",
    #label.right = "Favors placebo",
    labels = long_labels,
    fontsize = 11,
    #print.I2    = FALSE,
    #print.tau2  = FALSE,
    #print.pval  = FALSE,
    cex = 0.85
  )



# ---- Node-splitting / inconsistency checks ----
  sm_netsplit <- netsplit(smc_sm_network)
  sm_netsplit                 #For all comparisons, the total evidence is either direct or indirect.
  

  
  netsplit_placebo <- netsplit(
    smc_sm_network,
    method = "SIDDE",
    reference.group = "placebo/nodrug",
    show = "all",
    baseline.reference = TRUE,
    common = TRUE,
    backtransf = TRUE,
    ci = TRUE
  )
  

  forest(netsplit_placebo)    #For all comparisons, the total evidence is either direct or indirect.



  # ---- Bayesian NMA ----
  library(gemtc)
  library(rjags)
  
  # Reshape to GEMTC "relative effects" format

  sm_data_long <- sm_data %>% 
    pivot_longer( #Reshape from wide to long
      cols = c("treat1", "treat2", "logTE", "seTE"), # Select columns to reshape
      names_to = c(".value"), # Use column name parts as output column names 
      names_pattern = "(..)"  # Group column names by first two columns (trt1 & trt2 == grouped into treatments)
    ) %>% 
    rename( #Rename columns
      study = studyid,
      diff = lo,
      std.err = se,
      treatment = tr
    ) %>% 
    arrange(study)
  
  
  
  # Multi-arm cleaning 
  sm_data_long$row_id <- seq_len(nrow(sm_data_long))
  drop_rows <- c(1, 2, 
                 21, 22, 
                 23, 24, 
                 25, 26) #Drop redundant multi-arm comparisons (not needed by GEMTC)
  
  sm_data_long <- sm_data_long %>%
    filter(!(row_id %in% drop_rows)) %>%
    select(-row_id)
  
  sm_data_long <- sm_data_long %>%   # Remove duplicated baseline rows within study x treatment
    group_by(study, treatment) %>%
    filter(!(is.na(diff) & duplicated(treatment))) %>%
    ungroup()
  
  
  baseline_SE <- tibble(
    study = c("Chandramohan2021",
                "Kweku2008"),
    baseline_se = c(0.3015113,
                    0.2294157)     
  )                 #standard error in baseline arm for multiarm studies, required by GEMTC
                    #Calculated as SE(λ) = 1/sqrt(y); where y = No. events (counts) in arm 
  
  
  sm_data_long <- sm_data_long %>%
    left_join(baseline_SE, by = "study") %>%
    mutate(
      std.err = if_else(is.na(diff) & !is.na(baseline_se), baseline_se, std.err)
    ) %>%
    select(-baseline_se)
  
  
  sm_data_long <- sm_data_long %>%
    mutate(
      treatment = recode(treatment,
                         "placebo/nodrug" = "placebo_nodrug",
                         "spaq+rtss"      = "spaq_rtss"
                          )
                        )
  
  #Treatment codes for plots
  treat.codes <- c(
    "spaq"           = "SPAQ",
    "placebo_nodrug" = "Placebo/No drug control",
    "rtss"           = "RTS,S vaccine",
    "spaq_rtss"      = "SPAQ + RTS,S Vaccine",
    "asaq"           = "ASAQ",
    "sp"             = "SP (bimonthly)",
    "asaq_bi"        = "ASAQ (bimonthly)"
  ) %>%
    data.frame(description = ., stringsAsFactors = FALSE) %>%
    rownames_to_column("id")
  
  
  
  
  # ---- NMA model ----
  sm_network <- mtc.network(
    data.re    = sm_data_long,
    treatments = treat.codes
  )
  
  summary(sm_network, use.description = TRUE)

  # Common effect (fixed)
  sm_model <- mtc.model(
    sm_network,
    likelihood  = "normal",
    link        = "identity",
    linearModel = "fixed",
    n.chain     = 4
  )
  
  mcmc  <- mtc.run(sm_model, n.adapt = 5000, n.iter = 100000, thin = 10)
  
  plot(mcmc)
  gelman.plot(mcmc)
  gelman.diag(mcmc)$mpsrf
  
  summary(mcmc)
  
  
  
  # ---- network plot ----
  
  plot(
    sm_network,
    use.description = TRUE,
    vertex.color = "yellow",
    vertex.label.color = "black",
    vertex.label.family = "Arial",
    vertex.label.dist = 2.1,
    vertex.label.cex = 0.8
  )


  # --- Relative effects + forest plot ----
  
  #Summary of relative effects vs placebo
  results_vs_placebo <- relative.effect(mcmc, t1 = "placebo_nodrug")
  summary(results_vs_placebo)
  
  
  #Extract results
  output <- summary(results_vs_placebo) 
  
  irr_re_plac <- as.data.frame(output[["summaries"]][["statistics"]]) #Extract logIRRs
  cri_re_plac <- as.data.frame(output[["summaries"]][["quantiles"]]) #Extract credible intervals
  
  re_plac <- cbind(irr_re_plac, cri_re_plac)
  re_plac$comparison <- rownames(re_plac)
  
  
  #Forest plot (quick)
  forest(results_vs_placebo, use.description = TRUE)
  

  
  
  #Forest plot (clean)
    #Prepare data
  #Rename treatments
  re_plac <- re_plac %>% 
      mutate(label = case_when(
        comparison == "d.placebo_nodrug.spaq" ~ "SPAQ",
        comparison == "d.placebo_nodrug.rtss" ~ "Seasonal RTS,S",
        comparison == "d.placebo_nodrug.spaq_rtss" ~ "Seasonal RTS,S + SPAQ",
        comparison == "d.placebo_nodrug.asaq" ~ "ASAQ",
        comparison == "d.placebo_nodrug.sp" ~ "SP (bimonthly)",
        comparison == "d.placebo_nodrug.asaq_bi" ~ "ASAQ (bimonthly)",
      )
             )
  #Back-transform estimates
  re_plac <- re_plac %>%
    mutate(est = exp(`50%`))%>%
    mutate(lo = exp(`2.5%`))%>%
    mutate(hi = exp(`97.5%`))

  #Weights
  re_plac <- re_plac%>%
    mutate(weight = ((hi - lo) / (2 * 1.96)^2))
  
  #Final dataset for plot
  forest_data_sevmal <- re_plac %>%
    select(label, est, lo, hi, weight)
  
  
  # Make forest plot:
  # - Insert forest_data_sevmal object into forest plot_main.R 
  
  
  
  
  #Ranking
  rank_sm <- rank.probability(mcmc_ce_main, preferredDirection = -1)
  
  plot(rank_sm, beside = TRUE)
  
  
  
  
  
  # ---- Pairwise meta-analysis ----

  #Only SPAQ vs placebo
  spaq_placebo <- subset(sm_data, 
                         treat1 == "spaq" & treat2 == "placebo/nodrug") 
  
  #Pairwise meta-analysis
  spaq_plac_pairwise <- metagen(logTE, 
                                seTE, 
                                data = spaq_placebo,
                                sm = "OR",
                                test.subgroup = TRUE,
                                backtransform = TRUE)
     spaq_plac_pairwise
  
  
  #Forest plot
     #Labels
  spaq_studlab <- c("Dicko2011",
                    "Konate2011",
                    "Cisse2016",
                    "Tine2011")
    #Plot
  forest(spaq_plac_pairwise,
               studlab = spaq_studlab,
               random = FALSE,
               header.line = "both",
               leftcols = "studlab",
               label.eff = "",
               rightlabs = c("OR", "95% CI", "Weight"),
               xlab = "",
               smlab = "",
               cex         = 0.6)
  





