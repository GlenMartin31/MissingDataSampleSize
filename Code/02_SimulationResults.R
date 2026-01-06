# ##############################################################################

# Author of code: Glen P. Martin.

# This is code for a simulation study presented in a manuscript entitled: 
# Impact of Missing Data on Sample Size Requirements for Developing Clinical 
# Prediction Models
# Authors:
#   Glen P. Martin
#   Matthew Sperrin
#   Richard D. Riley

# ##############################################################################

library(tidyverse)

#Read in the result files from the CSF run
True_mod_performance <- read_rds(file = here::here("Outputs", 
                                                   "True_mod_performance.RDS")) %>%
  arrange(Scenario)
prediction_results <- read_rds(file = here::here("Outputs",
                                                 "prediction_results.RDS")) %>%
  mutate(Scenario = as.numeric(str_remove(Scenario, pattern = "scenario_"))) %>%
  arrange(Scenario)
EVPI <- read_rds(file = here::here("Outputs", 
                                   "EVPI.RDS")) %>%
  mutate(Scenario = as.numeric(str_remove(Scenario, pattern = "scenario_"))) %>%
  arrange(Scenario)
SampleSizeInfo <- read_rds(file = here::here("Outputs", 
                                             "SampleSizeInfo.RDS")) %>%
  arrange(Scenario)

sims_parameters <- read_rds(file = here::here("Outputs",
                                              "sims_parameters.RDS")) %>%
  mutate(Scenario = 1:n(), 
         .before = "P") %>%
  mutate("VarType" = map_chr(beta_X, 
                             function(X){
                               ifelse(any(X == 0),
                                      "SomeNoiseVariables",
                                      "NoNoiseVariables")
                               }
                             ))

####----------------------------------------------------------------------------
## Plot Calibration Slope results to compare with RR sample size overfitting control
####----------------------------------------------------------------------------

################################################################################
#### Plot median results across scenarios for Logistic model, with 
#### prevalence of Y as 20%, across both MAR and MNAR mechanisms
################################################################################
ggpubr::ggarrange(
  prediction_results %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise(dplyr::across(tidyr::everything(), ~median(., na.rm = T)),
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "Logistic") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, CalSlope) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = CalSlope, group = Imputation, colour = Imputation)) +
    geom_line() +
    geom_hline(aes(yintercept = 0.9), linetype = "dashed", colour = "black") +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Calibration Slope") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise(dplyr::across(tidyr::everything(), ~median(., na.rm = T)),
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MNAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "Logistic") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, CalSlope) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = CalSlope, group = Imputation, colour = Imputation)) +
    geom_line() +
    geom_hline(aes(yintercept = 0.9), linetype = "dashed", colour = "black") +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Calibration Slope") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)

################################################################################
#### Plot median results across scenarios for AIC model, with 
#### prevalence of Y as 20%, across both MAR and MNAR mechanisms
################################################################################
ggpubr::ggarrange(
  prediction_results %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise(dplyr::across(tidyr::everything(), ~median(., na.rm = T)),
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "AIC") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, CalSlope) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = CalSlope, group = Imputation, colour = Imputation)) +
    geom_line() +
    geom_hline(aes(yintercept = 0.9), linetype = "dashed", colour = "black") +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Calibration Slope") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise(dplyr::across(tidyr::everything(), ~median(., na.rm = T)),
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MNAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "AIC") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, CalSlope) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = CalSlope, group = Imputation, colour = Imputation)) +
    geom_line() +
    geom_hline(aes(yintercept = 0.9), linetype = "dashed", colour = "black") +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Calibration Slope") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)

################################################################################
#### Plot median results across scenarios for LASSO model, with 
#### prevalence of Y as 20%, across both MAR and MNAR mechanisms
################################################################################
ggpubr::ggarrange(
  prediction_results %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise(dplyr::across(tidyr::everything(), ~median(., na.rm = T)),
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "LASSO") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, CalSlope) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = CalSlope, group = Imputation, colour = Imputation)) +
    geom_line() +
    geom_hline(aes(yintercept = 0.9), linetype = "dashed", colour = "black") +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Calibration Slope") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise(dplyr::across(tidyr::everything(), ~median(., na.rm = T)),
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MNAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "LASSO") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, CalSlope) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = CalSlope, group = Imputation, colour = Imputation)) +
    geom_line() +
    geom_hline(aes(yintercept = 0.9), linetype = "dashed", colour = "black") +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Calibration Slope") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)



################################################################################
#### Plot instability results across scenarios for Logistic model, with 
#### prevalence of Y as 20%, across both MAR and MNAR mechanisms
################################################################################
ggpubr::ggarrange(
  prediction_results %>%
    dplyr::select(Scenario, iter, Model, Imputation, CalSlope) %>%
    dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                             CalSlope <= 1.1,
                                           1,
                                           0)) %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise("AssuranceProp" = (sum(CalSlopeInRange)/n())*100,
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "Logistic") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, AssuranceProp) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = AssuranceProp, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    dplyr::select(Scenario, iter, Model, Imputation, CalSlope) %>%
    dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                             CalSlope <= 1.1,
                                           1,
                                           0)) %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise("AssuranceProp" = (sum(CalSlopeInRange)/n())*100,
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MNAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "Logistic") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, AssuranceProp) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = AssuranceProp, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)

################################################################################
#### Plot instability results across scenarios for AIC model, with 
#### prevalence of Y as 20%, across both MAR and MNAR mechanisms
################################################################################
ggpubr::ggarrange(
  prediction_results %>%
    dplyr::select(Scenario, iter, Model, Imputation, CalSlope) %>%
    dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                             CalSlope <= 1.1,
                                           1,
                                           0)) %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise("AssuranceProp" = (sum(CalSlopeInRange)/n())*100,
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "AIC") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, AssuranceProp) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = AssuranceProp, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    dplyr::select(Scenario, iter, Model, Imputation, CalSlope) %>%
    dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                             CalSlope <= 1.1,
                                           1,
                                           0)) %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise("AssuranceProp" = (sum(CalSlopeInRange)/n())*100,
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MNAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "AIC") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, AssuranceProp) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = AssuranceProp, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)

################################################################################
#### Plot instability results across scenarios for LASSO model, with 
#### prevalence of Y as 20%, across both MAR and MNAR mechanisms
################################################################################
ggpubr::ggarrange(
  prediction_results %>%
    dplyr::select(Scenario, iter, Model, Imputation, CalSlope) %>%
    dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                             CalSlope <= 1.1,
                                           1,
                                           0)) %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise("AssuranceProp" = (sum(CalSlopeInRange)/n())*100,
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "LASSO") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, AssuranceProp) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = AssuranceProp, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    dplyr::select(Scenario, iter, Model, Imputation, CalSlope) %>%
    dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                             CalSlope <= 1.1,
                                           1,
                                           0)) %>%
    dplyr::group_by(Scenario, Model, Imputation) %>%
    dplyr::summarise("AssuranceProp" = (sum(CalSlopeInRange)/n())*100,
                     .groups = "drop") %>%
    dplyr::left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  missing_mech == "MNAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "LASSO") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, AssuranceProp) %>%
    dplyr::mutate(rho_X = factor(paste("rho=", rho_X, sep ="")),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = AssuranceProp, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(rho_X ~ prop_missing, scales = "fixed") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)


####----------------------------------------------------------------------------
## Plot EVPI results
####----------------------------------------------------------------------------

ggpubr::ggarrange(
  EVPI %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    filter(P == 10,
           VarType == "SomeNoiseVariables",
           rho_X == 0,
           Prev_Y == 0.2,
           Model == "Logistic",
           missing_mech == "MAR") %>%
    select(multipliers, prop_missing, Z, ends_with("_Mean")) %>%
    pivot_longer(cols = ends_with("_Mean")) %>%
    mutate(name = str_remove(name, pattern = "EVPI_NB_model_"),
           name = str_remove(name, pattern = "_Mean")) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  name = fct_recode(factor(name,
                                           c("FullyObserved",
                                             "CCA",
                                             "MeanImpute",
                                             "RI",
                                             "RF",
                                             "MI")),
                                    "Complete Case Analysis" = "CCA",
                                    "Fully Observed" = "FullyObserved",
                                    "Mean Imputation" = "MeanImpute",
                                    "Multiple Imputation" = "MI",
                                    "Random Forest Imputation" = "RF",
                                    "Single Regression Imputation" = "RI"),
                  multipliers = paste("delta=", multipliers, sep="")) %>%
    ggplot(aes(x = Z, y = value, group = name, colour = name)) +
    geom_line() +
    facet_grid(multipliers~prop_missing) +
    ylab("EVPI") +
    xlab("Threshold") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  EVPI %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    filter(P == 10,
           VarType == "SomeNoiseVariables",
           rho_X == 0,
           Prev_Y == 0.2,
           Model == "Logistic",
           missing_mech == "MNAR") %>%
    select(multipliers, prop_missing, Z, ends_with("_Mean")) %>%
    pivot_longer(cols = ends_with("_Mean")) %>%
    mutate(name = str_remove(name, pattern = "EVPI_NB_model_"),
           name = str_remove(name, pattern = "_Mean")) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  name = fct_recode(factor(name,
                                           c("FullyObserved",
                                             "CCA",
                                             "MeanImpute",
                                             "RI",
                                             "RF",
                                             "MI")),
                                    "Complete Case Analysis" = "CCA",
                                    "Fully Observed" = "FullyObserved",
                                    "Mean Imputation" = "MeanImpute",
                                    "Multiple Imputation" = "MI",
                                    "Random Forest Imputation" = "RF",
                                    "Single Regression Imputation" = "RI"),
                  multipliers = paste("delta=", multipliers, sep="")) %>%
    ggplot(aes(x = Z, y = value, group = name, colour = name)) +
    geom_line() +
    facet_grid(multipliers~prop_missing) +
    ylab("EVPI") +
    xlab("Threshold") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)


ggpubr::ggarrange(
  EVPI %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    filter(P == 10,
           VarType == "SomeNoiseVariables",
           rho_X == 0,
           Prev_Y == 0.2,
           Model == "AIC",
           missing_mech == "MAR") %>%
    select(multipliers, prop_missing, Z, ends_with("_Mean")) %>%
    pivot_longer(cols = ends_with("_Mean")) %>%
    mutate(name = str_remove(name, pattern = "EVPI_NB_model_"),
           name = str_remove(name, pattern = "_Mean")) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  name = fct_recode(factor(name,
                                           c("FullyObserved",
                                             "CCA",
                                             "MeanImpute",
                                             "RI",
                                             "RF",
                                             "MI")),
                                    "Complete Case Analysis" = "CCA",
                                    "Fully Observed" = "FullyObserved",
                                    "Mean Imputation" = "MeanImpute",
                                    "Multiple Imputation" = "MI",
                                    "Random Forest Imputation" = "RF",
                                    "Single Regression Imputation" = "RI"),
                  multipliers = paste("delta=", multipliers, sep="")) %>%
    ggplot(aes(x = Z, y = value, group = name, colour = name)) +
    geom_line() +
    facet_grid(multipliers~prop_missing) +
    ylab("EVPI") +
    xlab("Threshold") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  EVPI %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    filter(P == 10,
           VarType == "SomeNoiseVariables",
           rho_X == 0,
           Prev_Y == 0.2,
           Model == "AIC",
           missing_mech == "MNAR") %>%
    select(multipliers, prop_missing, Z, ends_with("_Mean")) %>%
    pivot_longer(cols = ends_with("_Mean")) %>%
    mutate(name = str_remove(name, pattern = "EVPI_NB_model_"),
           name = str_remove(name, pattern = "_Mean")) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  name = fct_recode(factor(name,
                                           c("FullyObserved",
                                             "CCA",
                                             "MeanImpute",
                                             "RI",
                                             "RF",
                                             "MI")),
                                    "Complete Case Analysis" = "CCA",
                                    "Fully Observed" = "FullyObserved",
                                    "Mean Imputation" = "MeanImpute",
                                    "Multiple Imputation" = "MI",
                                    "Random Forest Imputation" = "RF",
                                    "Single Regression Imputation" = "RI"),
                  multipliers = paste("delta=", multipliers, sep="")) %>%
    ggplot(aes(x = Z, y = value, group = name, colour = name)) +
    geom_line() +
    facet_grid(multipliers~prop_missing) +
    ylab("EVPI") +
    xlab("Threshold") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)

ggpubr::ggarrange(
  EVPI %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    filter(P == 10,
           VarType == "SomeNoiseVariables",
           rho_X == 0,
           Prev_Y == 0.2,
           Model == "LASSO",
           missing_mech == "MAR") %>%
    select(multipliers, prop_missing, Z, ends_with("_Mean")) %>%
    pivot_longer(cols = ends_with("_Mean")) %>%
    mutate(name = str_remove(name, pattern = "EVPI_NB_model_"),
           name = str_remove(name, pattern = "_Mean")) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  name = fct_recode(factor(name,
                                           c("FullyObserved",
                                             "CCA",
                                             "MeanImpute",
                                             "RI",
                                             "RF",
                                             "MI")),
                                    "Complete Case Analysis" = "CCA",
                                    "Fully Observed" = "FullyObserved",
                                    "Mean Imputation" = "MeanImpute",
                                    "Multiple Imputation" = "MI",
                                    "Random Forest Imputation" = "RF",
                                    "Single Regression Imputation" = "RI"),
                  multipliers = paste("delta=", multipliers, sep="")) %>%
    ggplot(aes(x = Z, y = value, group = name, colour = name)) +
    geom_line() +
    facet_grid(multipliers~prop_missing) +
    ylab("EVPI") +
    xlab("Threshold") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  EVPI %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    filter(P == 10,
           VarType == "SomeNoiseVariables",
           rho_X == 0,
           Prev_Y == 0.2,
           Model == "LASSO",
           missing_mech == "MNAR") %>%
    select(multipliers, prop_missing, Z, ends_with("_Mean")) %>%
    pivot_longer(cols = ends_with("_Mean")) %>%
    mutate(name = str_remove(name, pattern = "EVPI_NB_model_"),
           name = str_remove(name, pattern = "_Mean")) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  name = fct_recode(factor(name,
                                           c("FullyObserved",
                                             "CCA",
                                             "MeanImpute",
                                             "RI",
                                             "RF",
                                             "MI")),
                                    "Complete Case Analysis" = "CCA",
                                    "Fully Observed" = "FullyObserved",
                                    "Mean Imputation" = "MeanImpute",
                                    "Multiple Imputation" = "MI",
                                    "Random Forest Imputation" = "RF",
                                    "Single Regression Imputation" = "RI"),
                  multipliers = paste("delta=", multipliers, sep="")) %>%
    ggplot(aes(x = Z, y = value, group = name, colour = name)) +
    geom_line() +
    facet_grid(multipliers~prop_missing) +
    ylab("EVPI") +
    xlab("Threshold") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)


####----------------------------------------------------------------------------
## Plot Predictive Performance results
####----------------------------------------------------------------------------
ggpubr::ggarrange(
  prediction_results %>%
    group_by(Scenario, Model, Imputation) %>%
    summarise(across(everything(), ~median(., na.rm = T)),
              .groups = "drop") %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  rho_X == 0.5,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "Logistic") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, ends_with("_degradation")) %>%
    pivot_longer(cols = ends_with("_degradation")) %>%
    dplyr::mutate(name = str_remove(name, pattern = "_degradation"),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = value, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(prop_missing ~ name, scales = "fixed") +
    ylab("Model Degradation") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  prediction_results %>%
    group_by(Scenario, Model, Imputation) %>%
    summarise(across(everything(), ~median(., na.rm = T)),
              .groups = "drop") %>%
    left_join(sims_parameters, by = "Scenario")  %>%
    dplyr::filter(P == 10,
                  Prev_Y == 0.2,
                  rho_X == 0.5,
                  missing_mech == "MAR",
                  VarType == "SomeNoiseVariables", 
                  Model == "LASSO") %>%
    dplyr::select(Scenario, multipliers, prop_missing, 
                  rho_X, Model, Imputation, ends_with("_degradation")) %>%
    pivot_longer(cols = ends_with("_degradation")) %>%
    dplyr::mutate(name = str_remove(name, pattern = "_degradation"),
                  prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  Imputation = fct_recode(factor(Imputation,
                                                 c("FullyObserved",
                                                   "CCA",
                                                   "MeanImpute",
                                                   "RI",
                                                   "RF",
                                                   "MI")),
                                          "Complete Case Analysis" = "CCA",
                                          "Fully Observed" = "FullyObserved",
                                          "Mean Imputation" = "MeanImpute",
                                          "Multiple Imputation" = "MI",
                                          "Random Forest Imputation" = "RF",
                                          "Single Regression Imputation" = "RI")) %>%
    ggplot(aes(x = multipliers, y = value, group = Imputation, colour = Imputation)) +
    geom_line() +
    facet_grid(prop_missing ~ name, scales = "fixed") +
    ylab("Model Degradation") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank())
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)


