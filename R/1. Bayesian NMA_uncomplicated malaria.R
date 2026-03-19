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

smc_data <- read.csv("data/smc_nma_cleaned.csv", stringsAsFactors = FALSE)
smc_dot <- read.csv("data/smc_nma_adherence.csv", stringsAsFactors = FALSE)  




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
                       "spaq+rtss"      = "spaq_rtss"
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

# - # - # - # - # - # - # 
# Main NMA model====
# - # - # - # - # - # - # 


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

forest(relative.effect(mcmc_ce_main, t1 = "placebo_nodrug"),
      use.description = TRUE,
      xlim = c(-1.5, 0.5))



      output_main <- summary(results_vs_placebo) #Summary of relative effects vs placebo
      
      s_re_plac <- as.data.frame(output_main[["summaries"]][["statistics"]]) #Extract logIRRs
      q_re_plac <- as.data.frame(output_main[["summaries"]][["quantiles"]]) #Extract logIRRs
      
      df_re_plac <- cbind(s_re_plac, q_re_plac)
      df_re_plac$comparison <- rownames(df_re_plac)
      
      write.csv(df_re_plac, "re_plac_main.csv", row.names = FALSE) #Dataframe with relative effects vs placebo
      # Make forest plot:
      # - Exponentiate posterior means for 50%, 2.5% and 97.% 
      # - Calculate weight as inverse variance
      # - Run plot using forest plot_main. R


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
write.csv(league_table, "smc_ltable_uncompmalaria.csv", row.names = FALSE)


ltable_exp <- exp(league_table)
league_exp_round <- round(ltable_exp, 2)


write.csv(league_exp_round, "smc_ltable_uncompmalaria_exp.csv", row.names = TRUE)


# ---------------------------
# Meta-regression 1 (Adherence / DOT)
# ---------------------------

  smc_nma_dot <- read.csv("smc_nma_adherence.csv", header = TRUE) %>%
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
    linearModel = "fixed",
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
  

  # ---------------------------
  # Meta-regression 2 (median a437g)
  # ---------------------------
  
        # --- West africa-only data
        smc_data_long_wa <- smc_data_long %>% filter(studyid != "Nuwa2025")
  
        # --- West Africa-only model
        smc_network_wa <- mtc.network(
          data.re    = smc_data_long_wa,
          treatments = treat.codes
        )
        
        
        # --- NMA model (Common effect)
        smc_model_wa <- mtc.model(
                        smc_network_wa,
                        likelihood  = "normal",
                        link        = "identity",
                        linearModel = "fixed",
                        n.chain     = 4
                        )
        
        mcmc_wa_quick <- mtc.run(smc_model_wa, n.adapt = 50,   n.iter = 1000,   thin = 10)
        mcmc_wa_main  <- mtc.run(smc_model_wa, n.adapt = 5000, n.iter = 100000, thin = 10)
        
        
        summary(mcmc_wa_main)
  
  
  # A437G dataset
  smc_nma_437 <- read.csv("a437g_smc_covariate.csv", header = TRUE, stringsAsFactors = FALSE) %>%
    mutate(
      median_100_bin_label = ifelse(
        median_100_bin == 2, "High", "Low")
    ) 
  
  smc_nma_437 <- smc_nma_437 %>% select(study, studyid, median_100_bin)
  
  
  # --- meta-regression network
  network_437 <- mtc.network(
    data.re    = smc_data_long_wa,
    studies    = smc_nma_437,
    treatments = treat.codes
  )
  
  # --- regressor variable
  a437g <- list(
    coefficient = "shared",
    variable    = "median_100_bin",
    control     = "placebo_nodrug"
  )
  
  # ---- meta regression NMA model
  mr_437 <- mtc.model(
    network_437,
    likelihood  = "normal",
    link        = "identity",
    type        = "regression",
    linearModel = "fixed",
    regressor   = a437g
  )
  
  mcmc_437 <- mtc.run(mr_437, n.adapt = 5000, n.iter = 10000, thin = 10)
  summary(mcmc_437)
  
  
  
  
  # ---------------------------
  # Meta-regression 4 (cRCT vs iRCT)
  # ---------------------------
  smc_nma_rct <- read.csv("smc_nma_design.csv", header = TRUE, stringsAsFactors = FALSE) 
  
  
  network_rct <- mtc.network(
    data.re    = smc_data_long,
    studies    = smc_nma_rct,
    treatments = treat.codes
  )
  
  rct <- list(
    coefficient = "shared",
    variable    = "design",
    control     = "placebo_nodrug"
  )
  
  mr_rct <- mtc.model(
    network_rct,
    likelihood  = "normal",
    link        = "identity",
    type        = "regression",
    linearModel = "fixed",
    regressor   = rct
  )
  
  mcmc_rct <- mtc.run(mr_rct, n.adapt = 5000, n.iter = 10000, thin = 10)
  summary(mcmc_rct)
  
  
  # ---------------------------
  # Meta-regression 5 (median pfpr2-10)
  # ---------------------------
  
  
  
  
  
  
  
  
  
  # # - # # - # # - # # - # # - # # - #  
  # Meta-regression 6: Risk of bias ====
  # # - # # - # # - # # - # # - # # - #  
  
  smc_nma_rob <- read.csv("smc_nma_rob.csv", header = TRUE, stringsAsFactors = FALSE) #Dataset
  
  
  network_rob <- mtc.network(
    data.re    = smc_data_long,
    studies    = smc_nma_rob,
    treatments = treat.codes
  ) #Network
  
  rob <- list(
    coefficient = "shared",
    variable    = "rob",
    control     = "placebo_nodrug"
  ) #Regressor specification
  
  mr_rob <- mtc.model(
    network_rob,
    likelihood  = "normal",
    link        = "identity",
    type        = "regression",
    linearModel = "fixed",
    regressor   = rob
  ) #Network meta-regression model
  
  mcmc_rob <- mtc.run(mr_rob, n.adapt = 5000, 
                      n.iter = 10000, 
                      thin = 10) #Run model
  
  summary(mcmc_rob) #Summary
  
  
  
  re_rob <- relative.effect(mcmc_rob, t1 = "placebo_nodrug")
  
  
  rank_rob <- rank.probability(mcmc_rob, preferredDirection = -1)
  plot(rank_rob, beside = TRUE)
  
  sucra_rob <- dmetar::sucra(rank_rob, lower.is.better = FALSE)
  print(sucra_rob)
  
  plot(mcmc_rob)
  gelman.plot(mcmc_rob)
  gelman.diag(mcmc_rob)$mpsrf
  
  
  dev_main <- mtc.deviance(mcmc_ce_main)
  dev_rob <- mtc.deviance(mcmc_rob)
  
  mtc.devplot(dev_main, auto.layout = TRUE)
  mtc.levplot(dev_main)
  
  
  forest(relative.effect(mcmc_rob, t1 = "placebo_nodrug", covariate = 1),
         use.description = TRUE, xlim = c(-5, 5))
  title("High Risk of Bias")
  
  
  forest(relative.effect(mcmc_rob, t1 = "placebo_nodrug", covariate = 0),
         use.description = TRUE, xlim = c(-5, 5))
  title("Low Risk of Bias")
  
  
  
  barplot(
       dev_main$dev.re,
       las = 2,
       main = "mean posterior residual deviance of each study (main model)",
       ylab = "Deviance"
     )
  
  
  
  study_map <- unique(smc_data[, c("study", "studyid")])
  
  study_map
  study_labels <- study_map$studyid
  
  
  # # - # # - # # - # # - # # - # # - #  
  # Sensitivity analysis 1: West African trials only ====
  # # - # # - # # - # # - # # - # # - #  
  
  smc_nma_long_wa <- smc_data_long %>%
    filter(!studyid %in% c("Nuwa2025"))
  
  
  smc_bnma_wa <- mtc.network(
    data.re = smc_nma_long_wa,
    treatments = treat.codes
  )
  
  summary(smc_bnma_wa,
          use.description = TRUE)
  
  
  bnma_model_wa <- mtc.model(
    smc_bnma_wa,
    likelihood  = "normal",
    link        = "identity",
    linearModel = "fixed",
    n.chain     = 4
  )
  
  mcmc_wa_quick <- mtc.run(bnma_model_wa, n.adapt = 50,   n.iter = 1000,   thin = 10)
  mcmc_wa_main  <- mtc.run(bnma_model_wa, n.adapt = 5000, n.iter = 100000, thin = 10)
  
  
  # Node-splitting
  nodesplit_wa <- mtc.nodesplit(
    smc_bnma_wa,
    linearModel = "fixed",
    likelihood  = "normal",
    link        = "identity",
    n.adapt     = 5000,
    n.iter      = 100000,
    thin        = 10
  )
  nodesplit_wa_sum <- summary(nodesplit_wa)
  plot(nodesplit_wa_sum)
  
  
  
  
  results_vs_placebo_wa <- relative.effect(mcmc_wa_main, t1 = "placebo_nodrug")
  summary(results_vs_placebo_wa)
  
  
  waresults_vs_placebo <- relative.effect(mcmc_wa_main, t1 = "placebo_nodrug") #mtc.result object with relative effects vs placebo
  output_wa <- summary(waresults_vs_placebo) #Summary of relative effects vs placebo
  
  s_wa_re_plac <- as.data.frame(output_wa[["summaries"]][["statistics"]]) #Extract logIRRs
  q_wa_re_plac <- as.data.frame(output_wa[["summaries"]][["quantiles"]]) #Extract logIRRs
         
  df_wa_re_plac <- cbind(s_wa_re_plac, q_wa_re_plac)
  df_wa_re_plac$comparison <- rownames(df_wa_re_plac)
  
  write.csv(df_wa_re_plac, "wa_re_plac.csv", row.names = FALSE) #Dataframe with relative effects vs placebo
            # Make forest plot:
            # - Exponentiate posterior means for 50%, 2.5% and 97.% 
            # - Calculate weight as inverse variance
            # - Run plot using forest plot_main. R
  
  
  
  
  rank_prob_wa <- rank.probability(mcmc_wa_main, preferredDirection = -1)
  summary(rank_prob_wa)
  
  plot(rank_prob_wa, beside = TRUE)
  
  sucra_wa <- dmetar::sucra(rank_prob_wa, lower.is.better = FALSE)
  
  sucra_wa
  
  

  
  # - # - # - # - # - # - # 
  # Sensitivity analysis 3: Effects of all DP containing studies on DP estimates vs main model ====
  # - # - # - # - # - # - # 
  
  # Exclude Traore 2024 ====
  
  # Dataset
  smc_nma_extraore24 <- smc_data_long %>%
    filter(!studyid %in% c("Traore2024"))
  
  
  # Network
  smc_bnma_extraore24 <- mtc.network(
    data.re = smc_nma_extraore24,
    treatments = treat.codes
  )
  
    summary(smc_bnma_extraore24) #Network summary
    plot(smc_bnma_extraore24)
      
  #NMA model
  bnma_model_extraore24 <- mtc.model(
    smc_bnma_extraore24,
    likelihood  = "normal",
    link        = "identity",
    linearModel = "fixed",
    n.chain     = 4
  )
  
  mcmc_extraore24_main  <- mtc.run(bnma_model_extraore24, 
                                   n.adapt = 5000, 
                                   n.iter = 100000, 
                                   thin = 10)
  
  summary(mcmc_extraore24_main) # NMA summary
  
  
  
  # Relative effects vs placebo/no drug
  re_extraore24 <- relative.effect(mcmc_extraore24_main, t1 = "placebo_nodrug")
  
        #Quick forest plot (on log-scale)
        forest(relative.effect(mcmc_extraore24_main, t1 = "placebo_nodrug"),
               use.description = TRUE)
  
        # Extract summary statistics for forest plot
        output_extraore24 <- summary(re_extraore24)
        s_extraore24 <- as.data.frame(output_extraore24[["summaries"]][["statistics"]])
        q_extraore24 <- as.data.frame(output_extraore24[["summaries"]][["quantiles"]])
        df_extraore24 <- cbind(s_extraore24, q_extraore24)
        df_extraore24$comparison <- rownames(df_extraore24)
        write.csv(df_extraore24, "nma_extraore2024.csv", row.names = FALSE)
                #Prepare csv file for forest plot in excel as follows:
                  # - Exponentiate the median (50%), 2.5% and 97.5% estiamtes to get Median IRR and confidence intervals respectively
                  # - Calculate weight by inverse variance of exponentiated confidence intervals
                  # - Align with example csv files suited for forest plot code
                  # - produce forest plot in forest plot_main.R code
        
  
  
  # Treatment ranking
      #Rank probabilities  
      rankprob_extraore24 <- rank.probability(mcmc_extraore24_main, preferredDirection = -1)
      rankprob_extraore24 #Copy-paste rank probabilities into rankheat plot.R (use chatgpt for cleaning)
      plot(rankprob_extraore24, beside = TRUE) #Quick plot of rank probabilities
  
      #SUCRA
      sucra_extraore24 <- dmetar::sucra(rankprob_extraore24, lower.is.better = FALSE)
      sucra_extraore24 
  
  
  # Node-splitting
  nodesplit_extraore24 <- mtc.nodesplit(
    smc_bnma_extraore24,
    linearModel = "fixed",
    likelihood  = "normal",
    link        = "identity",
    n.adapt     = 5000,
    n.iter      = 100000,
    thin        = 10
  )
  nodesplit_extraore24_sum <- summary(nodesplit_extraore24)
  plot(nodesplit_extraore24)
  
  
  summary(nodesplit_extraore24)
  
  
  
  # Exclude Bojang 2010 ====
  
  # Dataset
  smc_nma_exbojang10 <- smc_data_long %>%
    filter(!studyid %in% c("Bojang2010"))
  
  
  # Network
  
  treat.codes_exbojang2010 <- c(
    "spaq"           = "SPAQ",
    "placebo_nodrug" = "Placebo/No drug control",
    "dp"             = "DHAPQ",
    "rtss"           = "RTS,S vaccine",
    "spaq_rtss"      = "SPAQ + RTS,S Vaccine",
    "asaq"           = "ASAQ",
    #"sppq"           = "SPPQ",
    "spas_sp_m"      = "SPAS/SP (monthly)",
    "sp"             = "SP (bimonthly)",
    "asaq_bi"        = "ASAQ (bimonthly)"
  ) %>%
    data.frame(description = ., stringsAsFactors = FALSE) %>%
    rownames_to_column("id")
  
  
  smc_bnma_exbojang10 <- mtc.network(
    data.re = smc_nma_exbojang10,
    treatments = treat.codes_exbojang2010
  )
  
  summary(smc_bnma_exbojang10) # Network summary
  plot(smc_bnma_exbojang10)
  
  
      plot(smc_bnma_exbojang10,
           use.description = TRUE,
           vertex.color = "yellow",
           vertex.label.color = "black",
           vertex.label.family = "Arial",
           vertex.label.dist = 2.1,
           vertex.label.cex = 0.8) # Network plot
  
  # NMA model
  bnma_model_exbojang10 <- mtc.model(
    smc_bnma_exbojang10,
    likelihood  = "normal",
    link        = "identity",
    linearModel = "fixed",
    n.chain     = 4
  )
  
  mcmc_exbojang10_main <- mtc.run(
    bnma_model_exbojang10, 
    n.adapt = 5000, 
    n.iter  = 100000, 
    thin    = 10
  )
  
  summary(mcmc_exbojang10_main) # NMA summary
  
  
  # Relative effects vs placebo/no drug
  re_exbojang10 <- relative.effect(
    mcmc_exbojang10_main,
    t1 = "placebo_nodrug"
  )
  
  # Quick forest plot (on log-scale)
  forest(
    relative.effect(mcmc_exbojang10_main, t1 = "placebo_nodrug"),
    use.description = TRUE
  )
  
  # Extract summary statistics for forest plot
  output_exbojang10 <- summary(re_exbojang10)
  s_exbojang10 <- as.data.frame(output_exbojang10[["summaries"]][["statistics"]])
  q_exbojang10 <- as.data.frame(output_exbojang10[["summaries"]][["quantiles"]])
  
  df_exbojang10 <- cbind(s_exbojang10, q_exbojang10)
  df_exbojang10$comparison <- rownames(df_exbojang10)
  
  write.csv(df_exbojang10, "nma_exbojang10.csv", row.names = FALSE)
      #Prepare csv file for forest plot in excel as follows:
      # - Exponentiate the median (50%), 2.5% and 97.5% estimates to get Median IRR and confidence intervals respectively
      # - Calculate weight by inverse variance of exponentiated confidence intervals
      # - Align with example csv files suited for forest plot code
      # - produce forest plot in forest plot_main.R code
  
  
  # Treatment ranking
  # Rank probabilities
  rankprob_exbojang10 <- rank.probability(
    mcmc_exbojang10_main,
    preferredDirection = -1
  )
  
  rankprob_exbojang10  # Copy-paste rank probabilities into rankheat plot.R
  plot(rankprob_exbojang10, beside = TRUE) # Quick plot of rank probabilities
  
  
  # SUCRA
  sucra_exbojang10 <- dmetar::sucra(
    rankprob_exbojang10,
    lower.is.better = FALSE
  )
  
  sucra_exbojang10
  
  
  # Node-splitting
  nodesplit_exbojang10 <- mtc.nodesplit(
    smc_bnma_exbojang10,
    linearModel = "fixed",
    likelihood  = "normal",
    link        = "identity",
    n.adapt     = 5000,
    n.iter      = 100000,
    thin        = 10
  )
  
  nodesplit_exbojang10_sum <- summary(nodesplit_exbojang10)
  
  summary(nodesplit_exbojang10)
  
  
  
  
  # Exclude Zongo 2015 ====
  
  # Dataset
  smc_nma_exzongo15 <- smc_data_long %>%
    filter(!studyid %in% c("Zongo2015"))
  
  
  # Network
  smc_bnma_exzongo15 <- mtc.network(
    data.re = smc_nma_exzongo15,
    treatments = treat.codes
  )
  
  summary(smc_bnma_exzongo15) # Network summary
  
  
  # NMA model
  bnma_model_exzongo15 <- mtc.model(
    smc_bnma_exzongo15,
    likelihood  = "normal",
    link        = "identity",
    linearModel = "fixed",
    n.chain     = 4
  )
  
  mcmc_exzongo15_main <- mtc.run(
    bnma_model_exzongo15, 
    n.adapt = 5000, 
    n.iter  = 100000, 
    thin    = 10
  )
  
  summary(mcmc_exzongo15_main) # NMA summary
  
  
  # Relative effects vs placebo/no drug
  re_exzongo15 <- relative.effect(
    mcmc_exzongo15_main,
    t1 = "placebo_nodrug"
  )
  
  # Quick forest plot (on log-scale)
  forest(
    relative.effect(mcmc_exzongo15_main, t1 = "placebo_nodrug"),
    use.description = TRUE
  )
  
  
  # Extract summary statistics for forest plot
  output_exzongo15 <- summary(re_exzongo15)
  s_exzongo15 <- as.data.frame(output_exzongo15[["summaries"]][["statistics"]])
  q_exzongo15 <- as.data.frame(output_exzongo15[["summaries"]][["quantiles"]])
  
  df_exzongo15 <- cbind(s_exzongo15, q_exzongo15)
  df_exzongo15$comparison <- rownames(df_exzongo15)
  
  write.csv(df_exzongo15, "nma_exzongo15.csv", row.names = FALSE)
  # Prepare csv file for forest plot in excel as follows:
  # - Exponentiate the median (50%), 2.5% and 97.5% estimates to get Median IRR and confidence intervals respectively
  # - Calculate weight by inverse variance of exponentiated confidence intervals
  # - Align with example csv files suited for forest plot code
  # - produce forest plot in forest plot_main.R code
  
  
  # Treatment ranking
  # Rank probabilities
  rankprob_exzongo15 <- rank.probability(
    mcmc_exzongo15_main,
    preferredDirection = -1
  )
  
  rankprob_exzongo15  # Copy-paste rank probabilities into rankheat plot.R
  plot(rankprob_exzongo15, beside = TRUE) # Quick plot of rank probabilities
  
  
  # SUCRA
  sucra_exzongo15 <- dmetar::sucra(
    rankprob_exzongo15,
    lower.is.better = FALSE
  )
  
  sucra_exzongo15
  
  
  # Node-splitting
  nodesplit_exzongo15 <- mtc.nodesplit(
    smc_bnma_exzongo15,
    linearModel = "fixed",
    likelihood  = "normal",
    link        = "identity",
    n.adapt     = 5000,
    n.iter      = 100000,
    thin        = 10
  )
  
  nodesplit_exzongo15_sum <- summary(nodesplit_exzongo15)
  
  summary(nodesplit_exzongo15)
  
  
  
  
  
  
