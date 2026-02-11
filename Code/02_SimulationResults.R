# ##############################################################################

# Author of code: Glen P. Martin.

# This is code for a simulation study presented in a manuscript entitled: 
# Impact of Missing Data on Sample Size Requirements for Developing Clinical 
# Prediction Models
# Authors:
#   Glen P. Martin
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

table(SampleSizeInfo$required_SS)
table(SampleSizeInfo$N_dev)

####----------------------------------------------------------------------------
## Function to plot Calibration Slope/ Overfitting results
####----------------------------------------------------------------------------
overfitting_fnc <- function(predictiveperformance,
                            simulation_paramaters,
                            missingnessmech = c("MAR", "MNAR"),
                            VariableScenario = c("SomeNoiseVariables",
                                                 "NoNoiseVariables"),
                            rhovalue,
                            metric = c("median", 
                                       "PosteriorProb0.9to1.1")) {
  missingnessmech <- as.character(match.arg(missingnessmech))
  VariableScenario <- as.character(match.arg(VariableScenario))
  metric <- as.character(match.arg(metric))
  
  if(metric == "median") {
    predictiveperformance_subset <- predictiveperformance %>%
      dplyr::left_join(simulation_paramaters, by = "Scenario")  %>%
      dplyr::filter(missing_mech == missingnessmech,
                    VarType == VariableScenario, 
                    rho_X == rhovalue) %>%
      dplyr::group_by(Scenario, ModellingMethod, DevelopmentImputationMethod, IPDImputationMethod,
                      P, beta_X, rho_X, Prev_Y, prop_missing, missing_mech,
                      multipliers, VarType) %>%
      dplyr::summarise("CalSlope" = median(CalSlope),
                       .groups = "drop") %>%
      dplyr::mutate(IPDImputationMethod = ifelse(IPDImputationMethod == "FullyObserved",
                                                 "Fully Observed Target Population",
                                                 "Imputed Target Population")) %>%
      dplyr::select(Scenario, ModellingMethod,
                    DevelopmentImputationMethod, IPDImputationMethod,
                    CalSlope,
                    rho_X,
                    prop_missing,
                    multipliers)
  } else {
    predictiveperformance_subset <- predictiveperformance %>%
      dplyr::left_join(simulation_paramaters, by = "Scenario")  %>%
      dplyr::filter(missing_mech == missingnessmech,
                    VarType == VariableScenario, 
                    rho_X == rhovalue) %>%
      dplyr::mutate(CalSlopeInRange = ifelse(CalSlope >= 0.9 & 
                                               CalSlope <= 1.1,
                                             1,
                                             0)) %>%
      dplyr::group_by(Scenario, ModellingMethod, DevelopmentImputationMethod, IPDImputationMethod,
                      P, beta_X, rho_X, Prev_Y, prop_missing, missing_mech,
                      multipliers, VarType) %>%
      dplyr::summarise("CalSlope" = (sum(CalSlopeInRange)/n())*100,
                       .groups = "drop") %>%
      dplyr::mutate(IPDImputationMethod = ifelse(IPDImputationMethod == "FullyObserved",
                                                 "Fully Observed Target Population",
                                                 "Imputed Target Population")) %>%
      dplyr::select(Scenario, ModellingMethod,
                    DevelopmentImputationMethod, IPDImputationMethod,
                    CalSlope,
                    rho_X,
                    prop_missing,
                    multipliers)
  }
  
  FullyObservedPlot <- predictiveperformance_subset %>%
    dplyr::filter(IPDImputationMethod == "Fully Observed Target Population") %>%
    dplyr::mutate(ModellingMethod = fct_recode(factor(ModellingMethod,
                                                      levels = c("logistic",
                                                                 "AIC",
                                                                 "lasso")),
                                               "Unpenalised Model\n no variable selection" = "logistic",
                                               "Unpenalised Model\n AIC variable Selection" = "AIC",
                                               "LASSO\n Penalised Model" = "lasso"),
                  
                  prop_missing = factor(paste((prop_missing*100), 
                                              "%", " Missing", 
                                              sep = "")),
                  
                  DevelopmentImputationMethod = fct_recode(factor(DevelopmentImputationMethod,
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
    ggplot(aes(x = multipliers, y = CalSlope, 
               group = DevelopmentImputationMethod,
               colour = DevelopmentImputationMethod)) +
    facet_grid(prop_missing~ModellingMethod) +
    geom_line() +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank(),
          legend.position = "top")
  
  ImputedTargetPlot <- predictiveperformance_subset %>%
    dplyr::filter(IPDImputationMethod == "Imputed Target Population") %>%
    dplyr::bind_rows(predictiveperformance_subset %>%
                       dplyr::filter(DevelopmentImputationMethod %in% c("FullyObserved",
                                                                        "CCA"),
                                     IPDImputationMethod == "Fully Observed Target Population") %>%
                       dplyr::mutate(IPDImputationMethod = "Imputed Target Population")) %>%
    dplyr::arrange(Scenario) %>%
    dplyr::mutate(ModellingMethod = fct_recode(factor(ModellingMethod,
                                                      levels = c("logistic",
                                                                 "AIC",
                                                                 "lasso")),
                                               "Unpenalised Model\n no variable selection" = "logistic",
                                               "Unpenalised Model\n AIC variable Selection" = "AIC",
                                               "LASSO\n Penalised Model" = "lasso"),
                  
                  prop_missing = factor(paste((prop_missing*100), 
                                              "%", " Missing", 
                                              sep = "")),
                  
                  DevelopmentImputationMethod = fct_recode(factor(DevelopmentImputationMethod,
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
    ggplot(aes(x = multipliers, y = CalSlope, 
               group = DevelopmentImputationMethod,
               colour = DevelopmentImputationMethod)) +
    facet_grid(prop_missing~ModellingMethod) +
    geom_line() +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank(),
          legend.position = "top")
  
  if(metric == "median") {
    FullyObservedPlot <- FullyObservedPlot +
      geom_hline(yintercept = 0.9, linetype = "dotted") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      ylab("Calibration Slope")
    
    ImputedTargetPlot <- ImputedTargetPlot +
      geom_hline(yintercept = 0.9, linetype = "dotted") +
      geom_hline(yintercept = 1, linetype = "dashed") +
      ylab("Calibration Slope")
  } else{
    FullyObservedPlot <- FullyObservedPlot +
      ylab("Prob(Calibration Slope between 0.9 and 1.1)")
    
    ImputedTargetPlot <- ImputedTargetPlot +
      ylab("Prob(Calibration Slope between 0.9 and 1.1)")
  }
  
  ggpubr::ggarrange(FullyObservedPlot,
                    ImputedTargetPlot,
                    common.legend = TRUE,
                    labels = c("A", "B"))
}



####----------------------------------------------------------------------------
## Function to plot degradation in performance due to missingness (step 5.5)
####----------------------------------------------------------------------------
performancedegradation_fnc <- function(predictiveperformance,
                                       simulation_paramaters,
                                       DGM_performance,
                                       missingnessmech = c("MAR", "MNAR"),
                                       VariableScenario = c("SomeNoiseVariables",
                                                            "NoNoiseVariables"),
                                       modelname = c("logistic",
                                                     "AIC",
                                                     "lasso"),
                                       rhovalue) {
  missingnessmech <- as.character(match.arg(missingnessmech))
  VariableScenario <- as.character(match.arg(VariableScenario))
  modelname <- as.character(match.arg(modelname))
  
  predictiveperformance_subset <- predictiveperformance %>%
    dplyr::left_join(simulation_paramaters, by = "Scenario")  %>%
    dplyr::filter(missing_mech == missingnessmech,
                  VarType == VariableScenario, 
                  ModellingMethod == modelname,
                  rho_X == rhovalue) %>%
    dplyr::group_by(Scenario, ModellingMethod, DevelopmentImputationMethod, IPDImputationMethod,
                    P, beta_X, rho_X, Prev_Y, prop_missing, missing_mech,
                    multipliers, VarType) %>%
    dplyr::summarise(across(c("OE", "CalInt", "CalSlope", "AUC"),
                            ~median(.)),
                     .groups = "drop") %>%
    dplyr::mutate(IPDImputationMethod = ifelse(IPDImputationMethod == "FullyObserved",
                                               "Fully Observed Target Population",
                                               "Imputed Target Population"))
  
  
  ## Summarise posterior distribution of fully observed CPM within the fully
  ## observed target population
  FullyObservedCPM_FullyObservedTarget <- predictiveperformance_subset %>%
    dplyr::filter(DevelopmentImputationMethod == "FullyObserved",
                  IPDImputationMethod == "Fully Observed Target Population") %>%
    dplyr::select(Scenario, ModellingMethod,
                  DevelopmentImputationMethod, IPDImputationMethod,
                  OE:AUC,
                  rho_X,
                  prop_missing,
                  multipliers)
  
  ## Performance of the reference model within the fully observed target
  ## population Note: this is also the performance of the reference model within
  ## the imputed target population, given the design of the simulation (since
  ## true risks are the same and there is no missing outcomes so imputed target
  ## population outcomes are the same as the fully observed target outcomes)
  RefModel_FullyObservedTarget <- DGM_performance %>%
    dplyr::mutate(Scenario = as.numeric(str_remove(Scenario, "scenario_"))) %>%
    dplyr::arrange(Scenario)
  
  RefModel_ImputedTarget <- RefModel_FullyObservedTarget
  
  ## Summarise posterior distribution of developed CPM within the imputed
  ## target population
  CPM_ImputedTarget <- predictiveperformance_subset %>%
    dplyr::filter(DevelopmentImputationMethod != "FullyObserved",
                  IPDImputationMethod == "Imputed Target Population") %>%
    dplyr::select(Scenario, ModellingMethod,
                  DevelopmentImputationMethod, IPDImputationMethod,
                  OE:AUC,
                  rho_X,
                  prop_missing,
                  multipliers) %>%
    dplyr::bind_rows(predictiveperformance_subset %>%
                       dplyr::filter(DevelopmentImputationMethod == "CCA",
                                     IPDImputationMethod == "Fully Observed Target Population") %>%
                       dplyr::mutate(IPDImputationMethod = "Imputed Target Population") %>%
                       dplyr::select(Scenario, ModellingMethod,
                                     DevelopmentImputationMethod, IPDImputationMethod,
                                     OE:AUC,
                                     rho_X,
                                     prop_missing,
                                     multipliers)) %>%
    dplyr::arrange(Scenario)
  
  FullyObservedCPM_Degradation <- FullyObservedCPM_FullyObservedTarget %>%
    dplyr::left_join(RefModel_FullyObservedTarget %>%
                       dplyr::rename("True_OE" = "OE",
                                     "True_CalInt" = "CalInt",
                                     "True_CalSlope" = "CalSlope",
                                     "True_AUC" = "AUC"),
                     by = "Scenario") %>%
    dplyr::mutate("FullyObserved_OE_Degradation" = True_OE - OE,
                  "FullyObserved_CalInt_Degradation" = True_CalInt - CalInt,
                  "FullyObserved_CalSlope_Degradation" = True_CalSlope - CalSlope,
                  "FullyObserved_AUC_Degradation" = True_AUC - AUC) %>%
    dplyr::select(Scenario,
                  prop_missing, multipliers, rho_X,
                  tidyr::ends_with("_Degradation"))
  
  CPM_Degradation <- CPM_ImputedTarget %>%
    dplyr::left_join(RefModel_ImputedTarget %>%
                       dplyr::rename("True_OE" = "OE",
                                     "True_CalInt" = "CalInt",
                                     "True_CalSlope" = "CalSlope",
                                     "True_AUC" = "AUC"),
                     by = "Scenario") %>%
    dplyr::mutate("OE_Degradation" = True_OE - OE,
                  "CalInt_Degradation" = True_CalInt - CalInt,
                  "CalSlope_Degradation" = True_CalSlope - CalSlope,
                  "AUC_Degradation" = True_AUC - AUC) %>%
    dplyr::select(Scenario, DevelopmentImputationMethod,
                  prop_missing, multipliers, rho_X,
                  tidyr::ends_with("_Degradation"))
  
  
  DegredationDueMissingness <- CPM_Degradation %>%
    dplyr::left_join(FullyObservedCPM_Degradation,
                     by = c("Scenario", "prop_missing", "multipliers", "rho_X")) %>%
    dplyr::mutate("OE_MissingnessDegradation" = FullyObserved_OE_Degradation - OE_Degradation, 
                  "CalInt_MissingnessDegradation" = FullyObserved_CalInt_Degradation - CalInt_Degradation,
                  "CalSlope_MissingnessDegradation" = FullyObserved_CalSlope_Degradation - CalSlope_Degradation,
                  "AUC_MissingnessDegradation" = FullyObserved_AUC_Degradation - AUC_Degradation) %>%
    dplyr::select(DevelopmentImputationMethod, 
                  prop_missing,
                  multipliers,
                  rho_X,
                  tidyr::ends_with("_MissingnessDegradation")) %>%
    tidyr::pivot_longer(cols = tidyr::ends_with("_MissingnessDegradation")) %>%
    tidyr::separate_wider_delim(name, delim = "_", names = c("Metric", "name")) %>%
    dplyr::select(-name)
  
  DegredationDueMissingness %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), 
                                              "%", " Missing", 
                                              sep = "")),
                  
                  DevelopmentImputationMethod = fct_recode(factor(DevelopmentImputationMethod,
                                                                  c("CCA",
                                                                    "MeanImpute",
                                                                    "RI",
                                                                    "RF",
                                                                    "MI")),
                                                           "Complete Case Analysis" = "CCA",
                                                           "Mean Imputation" = "MeanImpute",
                                                           "Multiple Imputation" = "MI",
                                                           "Random Forest Imputation" = "RF",
                                                           "Single Regression Imputation" = "RI"),
                  
                  Metric = fct_recode(factor(Metric,
                                             c("OE",
                                               "CalInt",
                                               "CalSlope",
                                               "AUC")),
                                      "AUC" = "AUC",
                                      "Calibration Intercept" = "CalInt",
                                      "Calibration Slope" = "CalSlope",
                                      "Observed:Expected Ratio" = "OE")) %>%
    ggplot(aes(x = multipliers, y = value, 
               group = DevelopmentImputationMethod,
               colour = DevelopmentImputationMethod)) +
    facet_grid(prop_missing~Metric) +
    geom_line() +
    geom_hline(yintercept = 0, linetype = "dotted") +
    ylab("Model Degradation due to Missing Data") +
    xlab(expression(paste("Development Sample Size Relative to ", N[min], sep=""))) +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank(),
          legend.position = "top")
}


####----------------------------------------------------------------------------
## Function to plot EVPI results
####----------------------------------------------------------------------------
EVPI_fnc <- function(EVPI_df,
                     simulation_paramaters,
                     missingnessmech = c("MAR", "MNAR"),
                     VariableScenario = c("SomeNoiseVariables",
                                          "NoNoiseVariables"),
                     modelname = c("logistic",
                                   "AIC",
                                   "lasso"),
                     rhovalue,
                     metric = c("EVPI", "REVPI")) {
  missingnessmech <- as.character(match.arg(missingnessmech))
  VariableScenario <- as.character(match.arg(VariableScenario))
  modelname <- as.character(match.arg(modelname))
  metric <- as.character(match.arg(metric))
  
  EVPI_subset <- EVPI %>%
    dplyr::left_join(simulation_paramaters, by = "Scenario")  %>%
    dplyr::filter(missing_mech == missingnessmech,
                  VarType == VariableScenario, 
                  rho_X == rhovalue) %>%
    dplyr::select(Scenario, multipliers, prop_missing, Z, DevelopmentImputationMethod,
                  IPDImputationMethod, ends_with("_Mean")) %>%
    tidyr::pivot_longer(cols = ends_with("_Mean"), names_to = "ModellingMethod") %>%
    dplyr::mutate(ModellingMethod = str_remove(ModellingMethod, pattern = "_Mean")) %>%
    tidyr::separate_wider_delim(cols = ModellingMethod, names = c("Metric", 
                                                                  "Model"),
                                delim = "_") %>%
    dplyr::filter(Model == modelname,
                  Metric == metric) %>%
    dplyr::mutate(prop_missing = factor(paste((prop_missing*100), "%", " Missing", sep = "")),
                  DevelopmentImputationMethod = fct_recode(factor(DevelopmentImputationMethod,
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
                  IPDImputationMethod = ifelse(IPDImputationMethod == "FullyObserved",
                                               "Target Fully Observed",
                                               "Missingness in Target"),
                  multipliers = paste("delta=", multipliers, sep=""))
  
  if(metric == "EVPI") {
    ylabel <- "EVPI"
  } else {
    ylabel <- "REVPI"
  }
  
  EVPI_subset %>%
    dplyr::filter(IPDImputationMethod == "Target Fully Observed") %>%
      ggplot(aes(x = Z, y = value, 
                 group = DevelopmentImputationMethod, 
                 colour = DevelopmentImputationMethod)) +
      geom_line() +
      facet_grid(multipliers~prop_missing) +
      ylab(ylabel) +
      xlab("Threshold") +
      theme_bw(base_size = 12) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.title=element_blank(),
            legend.position = "top")
}



####----------------------------------------------------------------------------
## Create the plots for the manuscript
####----------------------------------------------------------------------------

#Figure 1
overfitting_fnc(predictiveperformance = prediction_results,
                simulation_paramaters = sims_parameters,
                missingnessmech = "MAR",
                VariableScenario = "NoNoiseVariables",
                rhovalue =  0.5,
                metric = "median")
ggsave(filename = here::here("Manuscript", "Fig1.tiff"), dpi = 300)

#Supplementary Figure 1
overfitting_fnc(predictiveperformance = prediction_results,
                simulation_paramaters = sims_parameters,
                missingnessmech = "MNAR",
                VariableScenario = "NoNoiseVariables",
                rhovalue =  0.5,
                metric = "median")

#Figure 2
overfitting_fnc(predictiveperformance = prediction_results,
                simulation_paramaters = sims_parameters,
                missingnessmech = "MAR",
                VariableScenario = "NoNoiseVariables",
                rhovalue =  0.5,
                metric = "PosteriorProb0.9to1.1")
ggsave(filename = here::here("Manuscript", "Fig2.tiff"), dpi = 300)

#Supplementary Figure 2
ggpubr::ggarrange(performancedegradation_fnc(predictiveperformance = prediction_results,
                                             simulation_paramaters = sims_parameters,
                                             DGM_performance = True_mod_performance,
                                             missingnessmech = "MAR",
                                             VariableScenario = "NoNoiseVariables",
                                             modelname = "logistic",
                                             rhovalue = 0.5),
                  performancedegradation_fnc(predictiveperformance = prediction_results,
                                             simulation_paramaters = sims_parameters,
                                             DGM_performance = True_mod_performance,
                                             missingnessmech = "MNAR",
                                             VariableScenario = "NoNoiseVariables",
                                             modelname = "logistic",
                                             rhovalue = 0.5),
                  common.legend = TRUE,
                  labels = c("A", "B"))

#Figure 3
ggpubr::ggarrange(EVPI_fnc(EVPI_df = EVPI,
                           simulation_paramaters = sims_parameters,
                           missingnessmech = "MAR",
                           VariableScenario = "NoNoiseVariables",
                           modelname = "logistic",
                           rhovalue = 0.5,
                           metric = "EVPI"),
                  EVPI_fnc(EVPI_df = EVPI,
                           simulation_paramaters = sims_parameters,
                           missingnessmech = "MNAR",
                           VariableScenario = "NoNoiseVariables",
                           modelname = "logistic",
                           rhovalue = 0.5,
                           metric = "EVPI"),
                  common.legend = TRUE,
                  labels = c("A", "B"))
ggsave(filename = here::here("Manuscript", "Fig3.tiff"), dpi = 300)


#Supplementary Figure 3
overfitting_fnc(predictiveperformance = prediction_results,
                simulation_paramaters = sims_parameters,
                missingnessmech = "MNAR",
                VariableScenario = "SomeNoiseVariables",
                rhovalue =  0.5,
                metric = "median")

#Supplementary Figure 4
overfitting_fnc(predictiveperformance = prediction_results,
                simulation_paramaters = sims_parameters,
                missingnessmech = "MNAR",
                VariableScenario = "SomeNoiseVariables",
                rhovalue =  0.5,
                metric = "PosteriorProb0.9to1.1")

#Supplementary Figure 5
ggpubr::ggarrange(EVPI_fnc(EVPI_df = EVPI,
                           simulation_paramaters = sims_parameters,
                           missingnessmech = "MAR",
                           VariableScenario = "SomeNoiseVariables",
                           modelname = "logistic",
                           rhovalue = 0.5,
                           metric = "EVPI"),
                  EVPI_fnc(EVPI_df = EVPI,
                           simulation_paramaters = sims_parameters,
                           missingnessmech = "MNAR",
                           VariableScenario = "SomeNoiseVariables",
                           modelname = "logistic",
                           rhovalue = 0.5,
                           metric = "EVPI"),
                  common.legend = TRUE,
                  labels = c("A", "B"))




