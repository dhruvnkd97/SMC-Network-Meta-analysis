library(tidyr)
library(tibble)

# ------------------------------------------------------------
# Preparing analysis dataset: Severe malaria NMA
# ------------------------------------------------------------


sm_pairwise_input <- read_csv("data/2. smc_nma_severe malaria_raw data.csv")
  View(sm_pairwise_input)

  #Pairwise contrast-based TE / seTE estimation
  for_nma <- pairwise(
    treat = list(trt1, trt2, trt3, trt4),
    event = list(n1, n2, n3, n4),
    n = list(t1, t2, t3, t4),
    data = sm_pairwise_input, 
    studlab = studyid,
    sm = "OR",
  )    
  
  
  
  
  for_nma #LogORs for Chandramohan 2021, Dicko 2011 and Konate 2011

  write.csv(for_nma, file = "output/sm_pairwise_output.csv", row.names = FALSE) #store; outputs are on log scale



sm_data <- read.csv("output/sm_pairwise_output.csv") %>%
  select (studlab, treat1, treat2, TE, seTE) %>% 
  mutate(
    treat1 = recode(treat1, "nodrug" = "placebo"),
    treat2 = recode(treat2, "nodrug" = "placebo"),
    treat1 = recode(treat1, "placebo" = "placebo/nodrug"),
    treat2 = recode(treat2, "placebo" = "placebo/nodrug"),
  ) %>%
  rename(logTE = TE, studyid = studlab) 



rows <- tribble(
  ~studyid,      ~treat1,    ~treat2,           ~logTE,        ~seTE,
  "Cisse2016",   "spaq",     "placebo/nodrug",  -0.597837,     0.274887469, #see 2. kweku cisse_treatment effects_severe malaria.xlsx
  "Kweku2008",   "sp",       "placebo/nodrug",  -0.076881044,  0.383127375, #see 2. kweku cisse_treatment effects_severe malaria.xlsx
  "Kweku2008",   "asaq_bi",  "placebo/nodrug",   0.099845335,  0.445954381,
  "Kweku2008",   "asaq",     "placebo/nodrug",  -0.163696093,  0.445745462,
  "Kweku2008",   "asaq_bi",  "sp",               0.176726379,  0.445954381,
  "Kweku2008",   "asaq",     "sp",              -0.086815048,  0.445745462,
  "Kweku2008",   "asaq",     "asaq_bi",         -0.263541428,  0.445745462,
  "Tine2011",    "spaq",     "placebo/nodrug",  -1.17669054,   0.72935, #Bayesian logistic regression model estimates
  "Dicko2008",   "sp",       "placebo/nodrug",  -1.21623965,  0.7385527 #Bayesian logistic regression model estimates
)


sm_data <- bind_rows(sm_data, rows)

write.csv(sm_data, file = "data/SMC_NMA_severe malaria_analysis dataset.csv", row.names = FALSE) #Main analysis dataset
