# ##############################################################################

# Author of code: Glen P. Martin.

# This is code for the example 2 presented in a manuscript entitled: 
# Impact of Missing Data on Sample Size Requirements for Developing Clinical 
# Prediction Models
# Authors:
#   Glen P. Martin
#   Richard D. Riley

# ##############################################################################

library(tidyverse)
library(pmsampsize)
library(mice)
library(glmnet)
library(predRupdate)
library(missForestPredict)

#Load the data used in Martin et al. (2021) DOI: 10.1177/09622802211046388
MIMIC_df <- read_csv(file = here::here("Data", "mimic_penalised_cpms_cohort.csv"), 
                     col_names = TRUE,
                     # Adjust data structures, as needed:
                     col_types = cols(
                       gender = col_factor(levels = c("M", "F")),
                       admission_type = col_factor(levels = c("ELECTIVE", "URGENT", "EMERGENCY")),
                       ethnicity_grouped = col_factor(levels = c("white", "black", "asian",
                                                                 "hispanic", "other", "unknown"))
                     ))

####----------------------------------------------------------------------------
## Apply the same data cleaning steps as in Marin et al. 
#https://github.com/GlenMartin31/Penalised-CPMs-In-Minimum-Sample-Sizes
####----------------------------------------------------------------------------

MIMIC_df <- MIMIC_df %>%
  mutate(age_grouped = factor(ifelse(age < 30, "<30",
                                     ifelse(age < 40, "<40",
                                            ifelse(age < 50, "<50",
                                                   ifelse(age < 60, "<60",
                                                          ifelse(age < 70, "<70",
                                                                 ifelse(age < 80, "<80", 
                                                                        ">80")))))),
                              levels = c("<30", "<40", "<50", "<60", "<70", "<80", ">80"))
  ) %>%
  mutate(admission_type = fct_recode(admission_type,
                                     "Emergency" = "EMERGENCY", 
                                     "Non-Emergency" = "ELECTIVE", 
                                     "Non-Emergency" = "URGENT"),
         ethnicity_grouped = fct_recode(ethnicity_grouped,
                                        "white" = "white", 
                                        "black" = "black", 
                                        "other" = "asian",
                                        "other" = "hispanic", 
                                        "other" = "other", 
                                        NULL = "unknown"))

CandidatePredictorVariables <- c(names(MIMIC_df %>%
                                         select(age_grouped,
                                                gender,
                                                admission_type,
                                                ethnicity_grouped,
                                                ends_with("_mean"))))

Martin_et_al_analysis_cohort <-  MIMIC_df %>%
  select(contains(CandidatePredictorVariables),
         hospital_expire_flag) 

#Martin et al. didnt report the coefficients of the model, so we re-fit 
#the model here, and take that model as the reference model
reference_model <- glm(hospital_expire_flag ~ .,
                       data = Martin_et_al_analysis_cohort,
                       family = binomial(link = "logit"))
reference_model_coefs <- data.frame("Variable" = names(coef(reference_model)),
                                    "Estimate" = as.numeric(coef(reference_model)))

#Sample Size based on closed form solutions:
RR_sample_size <- pmsampsize(type = "b",
                             cstatistic = 0.73, #based on Martin et al. 
                             seed = 1234,
                             parameters = ncol(model.matrix(hospital_expire_flag ~ .,
                                                            data = Martin_et_al_analysis_cohort)),
                             shrinkage = 0.9,
                             prevalence = mean(Martin_et_al_analysis_cohort$hospital_expire_flag))


####----------------------------------------------------------------------------
##Practical Example 2 - assume we do not have access to the above data to inform 
##the sample size procedure
####----------------------------------------------------------------------------

############# Step 1.1 - define the reference model ############################

# assume these are obtained/requested: reference_model_coefs above

############# Step 2.1a - generate a large target dataset ######################

set.seed(983257)
continuous_Predictors <- MASS::mvrnorm(n = 500000,
                                       mu = Martin_et_al_analysis_cohort %>%
                                         select(ends_with("_mean")) %>%
                                         summarise(across(everything(), ~mean(., na.rm=T))) %>%
                                         as.numeric(),
                                       Sigma = diag(diag(var(Martin_et_al_analysis_cohort %>% 
                                                               select(ends_with("_mean")), na.rm =T))))
colnames(continuous_Predictors) <- names(Martin_et_al_analysis_cohort %>%
                                           select(ends_with("_mean")))

target_population <- data.frame(ID = 1:500000,
                                "age_grouped" = factor(sample(c("<30", "<40", "<50", "<60", "<70", "<80", ">80"),
                                                              size = 500000,
                                                              replace = TRUE,
                                                              prob = table(Martin_et_al_analysis_cohort$age_grouped)/nrow(Martin_et_al_analysis_cohort)),
                                                       levels = c("<30", "<40", "<50", "<60", "<70", "<80", ">80")),
                                
                                "gender" = factor(sample(c("M", "F"),
                                                         size = 500000,
                                                         replace = TRUE,
                                                         prob = table(Martin_et_al_analysis_cohort$gender)/nrow(Martin_et_al_analysis_cohort)),
                                                  levels = c("M", "F")),
                                
                                "admission_type" = factor(sample(c("Non-Emergency", "Emergency"),
                                                                 size = 500000,
                                                                 replace = TRUE,
                                                                 prob = table(Martin_et_al_analysis_cohort$admission_type)/nrow(Martin_et_al_analysis_cohort)),
                                                          levels = c("Non-Emergency", "Emergency")),
                                
                                "ethnicity_grouped" = factor(sample(c("white", "black", "other"),
                                                                    size = 500000,
                                                                    replace = TRUE,
                                                                    prob = table(Martin_et_al_analysis_cohort$ethnicity_grouped)/nrow(Martin_et_al_analysis_cohort)),
                                                             levels = c("white", "black", "other"))) %>%
  bind_cols(continuous_Predictors)

#Generate outcomes
DM <- model.matrix(~.,
                   target_population %>% 
                     dplyr::select(tidyr::contains(CandidatePredictorVariables)))
Pi <- 1 / (1 + exp(-as.numeric(DM %*% reference_model_coefs$Estimate)))
target_population <- target_population %>%
  mutate(Pi = Pi,
         hospital_expire_flag = rbinom(n = n(), 1, prob = Pi))

rm(DM, Pi)

############ Step 2.1b - generate a large target dataset with missingness ######
#for the practical example, we assume targeting ideal performance (see 
#Box 2 of paper), so skip this step.

############ Define a function to run steps 2.2a - 5 ##############
Example2_Analysis_fnc <- function(target_IPD,
                                  raw_IPD_data,
                                  ref_model_coefs,
                                  CandidatePredictorVars,
                                  N_dev,
                                  bootstrap_iterations = 100,
                                  model_type = c("logistic", "lasso")) {
  
  model_type <- as.character(match.arg(model_type))
  N_dev <- ceiling(N_dev)
  
  #Define an internal function to fit the CPM
  if (model_type == "logistic") {
    mod_fnc <- function(df) {
      logistic_mod <- glm(hospital_expire_flag ~.,
                          data = df %>%
                            dplyr::select(age_grouped,
                                          gender,
                                          admission_type,
                                          ethnicity_grouped,
                                          ends_with("_mean"),
                                          hospital_expire_flag),
                          family = binomial(link = "logit"))
      
      data.frame("Variable" = names(coef(logistic_mod)),
                 "Estimate" = as.numeric(coef(logistic_mod)))
    }
  } else if (model_type == "lasso") {
    mod_fnc <- function(df) {
      lasso_mod <- cv.glmnet(x = model.matrix(hospital_expire_flag ~ .,
                                              data = df %>%
                                                dplyr::select(age_grouped,
                                                              gender,
                                                              admission_type,
                                                              ethnicity_grouped,
                                                              ends_with("_mean"),
                                                              hospital_expire_flag))[,-1],
                             y = df %>% dplyr::pull("hospital_expire_flag"),
                             family = "binomial")
      
      data.frame("Variable" = rownames(coef(lasso_mod, 
                                            s = "lambda.min")),
                 "Estimate" = as.numeric(coef(lasso_mod, 
                                              s = "lambda.min")))
    } 
  }
  
  FullyObervedIPD_Performance <- NULL
  for(boot in 1:bootstrap_iterations) {
    #### Step 2.2a & 2.2b - generate development data
    # For practical example 2, use normal distributions based on univariate info
    continuous_Predictors <- MASS::mvrnorm(n = N_dev,
                                           mu = raw_IPD_data %>%
                                             select(ends_with("_mean")) %>%
                                             summarise(across(everything(), ~mean(., na.rm=T))) %>%
                                             as.numeric(),
                                           Sigma = diag(diag(var(Martin_et_al_analysis_cohort %>% 
                                                                   select(ends_with("_mean")), na.rm =T))))
    colnames(continuous_Predictors) <- names(raw_IPD_data %>%
                                               select(ends_with("_mean")))
    
    dev_data <- data.frame(ID = 1:N_dev,
                           "age_grouped" = factor(sample(c("<30", "<40", "<50", "<60", "<70", "<80", ">80"),
                                                         size = N_dev,
                                                         replace = TRUE,
                                                         prob = table(raw_IPD_data$age_grouped)/nrow(raw_IPD_data)),
                                                  levels = c("<30", "<40", "<50", "<60", "<70", "<80", ">80")),
                           
                           "gender" = factor(sample(c("M", "F"),
                                                    size = N_dev,
                                                    replace = TRUE,
                                                    prob = table(raw_IPD_data$gender)/nrow(raw_IPD_data)),
                                             levels = c("M", "F")),
                           
                           "admission_type" = factor(sample(c("Non-Emergency", "Emergency"),
                                                            size = N_dev,
                                                            replace = TRUE,
                                                            prob = table(raw_IPD_data$admission_type)/nrow(raw_IPD_data)),
                                                     levels = c("Non-Emergency", "Emergency")),
                           
                           "ethnicity_grouped" = factor(sample(c("white", "black", "other"),
                                                               size = N_dev,
                                                               replace = TRUE,
                                                               prob = table(raw_IPD_data$ethnicity_grouped)/nrow(raw_IPD_data)),
                                                        levels = c("white", "black", "other"))) %>%
      bind_cols(continuous_Predictors)
    
    #Generate outcomes
    DM <- model.matrix(~.,
                       dev_data %>% 
                         dplyr::select(tidyr::contains(CandidatePredictorVars)))
    Pi <- 1 / (1 + exp(-as.numeric(DM %*% ref_model_coefs$Estimate)))
    dev_data <- dev_data %>%
      mutate(Pi = Pi,
             hospital_expire_flag = rbinom(n = n(), 1, prob = Pi))
    rm(DM, Pi, continuous_Predictors)
    
    
    #### Step 2.2b - induce missingness into development data
    #in this example, assume an overall 20% missing data (noting its actually 22% in the real data):
    prop_missings <- 0.2
    #also assume no information on patterns of missingness, so take all 
    #possible patterns
    #assume MNAR in the absence of other information
    #ampute missing data under these assumptions
    amp <- mice::ampute(data = dev_data %>% 
                          dplyr::select(tidyr::contains(CandidatePredictorVars)),
                        prop = prop_missings,
                        mech = "MNAR")
    dev_data_withmissing <- dplyr::bind_cols("ID" = dev_data$ID,
                                             amp$amp,
                                             "Pi" = dev_data$Pi,
                                             "hospital_expire_flag" = dev_data$hospital_expire_flag)
    rm(amp)
    
    
    #### Step 3.1 - handle missingness in the development data. In practice one
    #would choose one imputation method based on the strategy they would use to
    #develop the CPM. For illustration here, we consider multiple approaches.
    
    #complete case analysis
    CCA_imputed_dev_df <- dev_data_withmissing %>%
      na.omit()
    
    #random forest imputation
    predmat <- missForestPredict::create_predictor_matrix(dev_data_withmissing)
    predmat["ID", ] <- predmat[, "ID"] <- 0
    predmat["Pi", ] <- predmat[, "Pi"] <- 0
    predmat["hospital_expire_flag", ] <- predmat[, "hospital_expire_flag"] <- 0
    RF_imputation <- missForestPredict::missForest(dev_data_withmissing,
                                                   save_models = TRUE,
                                                   predictor_matrix = predmat,
                                                   num.trees = 200,
                                                   num.threads = 2,
                                                   verbose = FALSE)
    RF_imputed_dev_df <- RF_imputation$ximp
    rm(predmat)
    
    #single regression imputation
    RI_imp <- mice::mice(dev_data_withmissing, m = 1, maxit = 0)
    predmat <- RI_imp$predictorMatrix
    predmat["ID", ] <- predmat[, "ID"] <- 0
    predmat["Pi", ] <- predmat[, "Pi"] <- 0
    predmat["hospital_expire_flag", ] <- predmat[, "hospital_expire_flag"] <- 0
    impmethods <- RI_imp$method
    impmethods[which(impmethods == "pmm")] <- "norm.predict"
    RI_imp <- mice::mice(dev_data_withmissing, 
                         m = 1, 
                         method = impmethods, 
                         predictorMatrix = predmat,
                         maxit = 25,
                         printFlag = 0)
    RI_imputed_dev_df <- mice::complete(RI_imp) 
    rm(predmat, impmethods)
    
    #multiple imputation
    MI_imp <- mice::mice(dev_data_withmissing, m = 1, maxit = 0)
    predmat <- MI_imp$predictorMatrix
    predmat["ID", ] <- predmat[, "ID"] <- 0
    predmat["Pi", ] <- predmat[, "Pi"] <- 0
    impmethods <- MI_imp$method
    impmethods[which(impmethods == "pmm")] <- "norm"
    MI_imp <- mice::mice(dev_data_withmissing,
                         m = 10, 
                         method = impmethods, 
                         predictorMatrix = predmat,
                         maxit = 25,
                         printFlag = 0)
    MI_imputed_dev_df <- mice::complete(MI_imp, action = "long") 
    rm(predmat, impmethods)
    
    
    #### Step 3.2 - fit the CPM 
    #define an internal function to allow repeating
    #of the model fitting on each imputed development dataset:
    model_fitting_fnc <- function(development_dataset,
                                  imputation_method) {
      if(imputation_method == "MI") {
        mod_coefs <- development_dataset %>%
          dplyr::group_by(.imp) %>%
          tidyr::nest() %>%
          dplyr::mutate("MI_coefs" = map(data, mod_fnc)) %>%
          dplyr::select(-data) %>%
          tidyr::unnest(cols = "MI_coefs") %>%
          dplyr::mutate(Variable = factor(Variable,
                                          levels = colnames(model.matrix(hospital_expire_flag ~ .,
                                                                         data = development_dataset %>%
                                                                           dplyr::select(age_grouped,
                                                                                         gender,
                                                                                         admission_type,
                                                                                         ethnicity_grouped,
                                                                                         ends_with("_mean"),
                                                                                         hospital_expire_flag))))) %>%
          dplyr::group_by(Variable) %>%
          dplyr::summarise(Estimate = mean(Estimate)) %>%
          as.data.frame()
      } else {
        mod_coefs <- mod_fnc(df = development_dataset) %>%
          dplyr::mutate(Variable = factor(Variable,
                                          levels = colnames(model.matrix(hospital_expire_flag ~ .,
                                                                         data = development_dataset %>%
                                                                           dplyr::select(age_grouped,
                                                                                         gender,
                                                                                         admission_type,
                                                                                         ethnicity_grouped,
                                                                                         ends_with("_mean"),
                                                                                         hospital_expire_flag)))))
      }
      mod_coefs
    }
    
    #fit the CPM to each version of the imputed development dataset
    dev_df_list <- list("FullyObserved" = dev_data,
                        "CCA" = CCA_imputed_dev_df,
                        "MissForest" = RF_imputed_dev_df,
                        "RI" = RI_imputed_dev_df,
                        "MI" = MI_imputed_dev_df)
    fitted_models <- purrr::map2(.x = dev_df_list,
                                 .y = names(dev_df_list),
                                 .f = model_fitting_fnc)
    
    #### Step 4.1 - validate the CPM within the fully observed target population
    #Define an internal function that tests performance
    predictive_performance_fnc <- function(model_info, 
                                           df,
                                           imputation_method) {
      if(imputation_method == "MI") {
        
        MI_long <- mice::complete(df, action = "long")
        
        coef_vals <- model_info$Estimate
        predictions <- MI_long %>% 
          dplyr::group_by(.imp) %>% 
          tidyr::nest() %>%
          dplyr::mutate("DM" = map(data, 
                                   function(X){
                                     model.matrix(hospital_expire_flag ~ .,
                                                  data = X %>%
                                                    dplyr::select(age_grouped,
                                                                  gender,
                                                                  admission_type,
                                                                  ethnicity_grouped,
                                                                  ends_with("_mean"),
                                                                  hospital_expire_flag))
                                   })) %>%
          dplyr::mutate("PR" = map(DM,
                                   function(X){
                                     1 / (1 + exp(-as.numeric(X %*% coef_vals)))
                                   })) %>% 
          tidyr::unnest(cols = c(".imp", "data", "PR")) %>% 
          dplyr::select(.imp, ID, hospital_expire_flag, PR)
        
        performance <- predictions %>%
          dplyr::group_by(.imp) %>%
          tidyr::nest() %>%
          dplyr::mutate("Performance" = purrr::map(data, 
                                                   function(X) {
                                                     performance <- predRupdate::pred_val_probs(binary_outcome = X$hospital_expire_flag,
                                                                                                Prob = X$PR,
                                                                                                cal_plot = FALSE)
                                                     performance <- data.frame("OE" = performance$OE_ratio,
                                                                               "CalInt" = performance$CalInt,
                                                                               "CalSlope" = performance$CalSlope,
                                                                               "AUC" = performance$AUC)
                                                   })) %>%
          dplyr::select(-data) %>%
          tidyr::unnest(cols = c(".imp", "Performance")) %>%
          dplyr::ungroup() %>%
          dplyr::summarise(dplyr::across(OE:AUC,
                                         ~mean(.)))
        performance
        
      } else {
        coef_vals <- model_info$Estimate
        DM <- model.matrix(hospital_expire_flag ~ .,
                           data = df %>%
                             dplyr::select(age_grouped,
                                           gender,
                                           admission_type,
                                           ethnicity_grouped,
                                           ends_with("_mean"),
                                           hospital_expire_flag))
        predictions <- 1 / (1 + exp(-as.numeric(DM %*% coef_vals)))
        
        performance <- predRupdate::pred_val_probs(binary_outcome = df$hospital_expire_flag,
                                                   Prob = predictions,
                                                   cal_plot = FALSE)
        performance <- data.frame("OE" = performance$OE_ratio,
                                  "CalInt" = performance$CalInt,
                                  "CalSlope" = performance$CalSlope,
                                  "AUC" = performance$AUC)
        performance
      }
    }
    
    predictive_perfm <- purrr::map(.x = fitted_models,
                                   .f = predictive_performance_fnc,
                                   df = target_IPD,
                                   imputation_method = "none") 
    FullyObervedIPD_Performance <- FullyObervedIPD_Performance %>%
      bind_rows(predictive_perfm %>%
                  dplyr::bind_rows(.id = "DevDataImputationMethod") %>%
                  dplyr::mutate("iter" = boot)) 
    rm(predictive_perfm)
    
    #### Step 4.2 - skip as targeting ideal performance
  }
  ### Step 4.1b - validate the reference model within the fully observed target population
  true_performance <- predRupdate::pred_val_probs(binary_outcome = target_IPD$hospital_expire_flag,
                                                  Prob = target_IPD$Pi,
                                                  cal_plot = FALSE)
  RefModelPerformance <- data.frame("OE" = true_performance$OE_ratio,
                                    "CalInt" = true_performance$CalInt,
                                    "CalSlope" = true_performance$CalSlope,
                                    "AUC" = true_performance$AUC)
  
  ### Return results across the repeats
  return(list("RefModelPerformance" = RefModelPerformance,
              "FullyObervedIPD_Performance" = FullyObervedIPD_Performance))
}


# Use the above function to run steps 2.2a - 5 ##############
library(furrr)
plan(multisession, workers = 3)
PracticalExample2Results <- crossing("N_dev" = c(RR_sample_size$sample_size,
                                                 RR_sample_size$sample_size * 1.5,
                                                 RR_sample_size$sample_size * 2),
                                     "model_type" = c("logistic")) %>%
  mutate("results" = future_pmap(.l = list(N_dev = N_dev,
                                           model_type = model_type),
                                 .f = Example2_Analysis_fnc,
                                 target_IPD = target_population,
                                 raw_IPD_data = Martin_et_al_analysis_cohort,
                                 ref_model_coefs = reference_model_coefs,
                                 CandidatePredictorVars = CandidatePredictorVariables,
                                 bootstrap_iterations = 100,
                                 .options = furrr_options(seed=TRUE),
                                 .progress = TRUE))

write_rds(PracticalExample2Results,
          file = here::here("Outputs", "PracticalExample2Results.RDS"))
# PracticalExample2Results <- read_rds(file = here::here("Outputs",
#                                                        "PracticalExample2Results.RDS"))

### Plot Results
Summarised_PracticalExample2Results <- PracticalExample2Results %>%
  dplyr::mutate("FullyObserved_Performance" = purrr::map(results, 
                                                         function(X) {
                                                           X$FullyObervedIPD_Performance %>%
                                                             dplyr::group_by(DevDataImputationMethod) %>%
                                                             dplyr::summarise(dplyr::across(OE:AUC,
                                                                                            ~median(.)))
                                                         }),
                
                "FullyObserved_Degradation" = purrr::map(results, 
                                                         function(X) {
                                                           
                                                           X$FullyObervedIPD_Performance %>% 
                                                             dplyr::bind_cols(X$RefModelPerformance %>% 
                                                                                dplyr::rename("OE_true" = "OE", 
                                                                                              "CalInt_True" = "CalInt", 
                                                                                              "CalSlope_True" = "CalSlope", 
                                                                                              "AUC_True" = "AUC")) %>% 
                                                             dplyr::mutate("OE_Degradation" = OE_true - OE, 
                                                                           "CalInt_Degradataion" = CalInt_True - CalInt, 
                                                                           "CalSlope_Degradation" = CalSlope_True - CalSlope, 
                                                                           "AUC_Degradation" = AUC_True - AUC) %>%
                                                             dplyr::group_by(DevDataImputationMethod) %>%
                                                             dplyr::summarise(dplyr::across(ends_with("_Degradation"),
                                                                                            ~median(.)))
                                                         }),
                
                "FullyObservedCalSlope_Posterior_Prob" = purrr::map(results, 
                                                                    function(X) {
                                                                      X$FullyObervedIPD_Performance %>%
                                                                        dplyr::group_by(DevDataImputationMethod) %>%
                                                                        dplyr::mutate("CalSlopeInRange" = ifelse(CalSlope >= 0.9 &
                                                                                                                   CalSlope <= 1.1,
                                                                                                                 1,
                                                                                                                 0)) %>%
                                                                        dplyr::summarise("ProbCalSlope" = (sum(CalSlopeInRange)/n())*100)
                                                                    }))


MedianPerformanceResults <- Summarised_PracticalExample2Results %>%
  dplyr::select(N_dev, model_type, FullyObserved_Performance) %>%
  tidyr::unnest(cols = c(FullyObserved_Performance)) %>%
  tidyr::pivot_longer(cols = OE:AUC) %>%
  dplyr::mutate("Data" = "Fully Observed Performance") %>%
  dplyr::mutate(DevDataImputationMethod = fct_recode(factor(DevDataImputationMethod,
                                                      c("FullyObserved",
                                                        "CCA",
                                                        "MissForest",
                                                        "RI",
                                                        "MI")),
                                               "Fully Observed (0% missing)" = "FullyObserved",
                                               "Complete Case Analysis" = "CCA",
                                               "Random Forest Imputation" = "MissForest",
                                               "Multiple Imputation" = "MI",
                                               "Regression Imputation" = "RI"),
                
                N_dev = paste("N=", ceiling(N_dev)))

DegradationResults <- Summarised_PracticalExample2Results %>%
  dplyr::select(N_dev, model_type, FullyObserved_Degradation) %>%
  tidyr::unnest(cols = c(FullyObserved_Degradation)) %>%
  tidyr::pivot_longer(cols = OE_Degradation:AUC_Degradation) %>%
  dplyr::mutate("Data" = "Fully Observed Performance") %>%
  dplyr::mutate(DevDataImputationMethod = fct_recode(factor(DevDataImputationMethod,
                                                      c("FullyObserved",
                                                        "CCA",
                                                        "MissForest",
                                                        "RI",
                                                        "MI")),
                                               "Fully Observed (0% missing)" = "FullyObserved",
                                               "Complete Case Analysis" = "CCA",
                                               "Random Forest Imputation" = "MissForest",
                                               "Multiple Imputation" = "MI",
                                               "Regression Imputation" = "RI"),
                
                N_dev = paste("N=", ceiling(N_dev)))


AssuranceResults <- Summarised_PracticalExample2Results %>%
  dplyr::select(N_dev, model_type, FullyObservedCalSlope_Posterior_Prob) %>%
  tidyr::unnest(cols = c(FullyObservedCalSlope_Posterior_Prob)) %>%
  dplyr::mutate("Data" = "Fully Observed Performance") %>%
  dplyr::mutate(DevDataImputationMethod = fct_recode(factor(DevDataImputationMethod,
                                                      c("FullyObserved",
                                                        "CCA",
                                                        "MissForest",
                                                        "RI",
                                                        "MI")),
                                               "Fully Observed (0% missing)" = "FullyObserved",
                                               "Complete Case Analysis" = "CCA",
                                               "Random Forest Imputation" = "MissForest",
                                               "Multiple Imputation" = "MI",
                                               "Regression Imputation" = "RI"),
                
                N_dev = paste("N=", ceiling(N_dev)))




ggpubr::ggarrange(
  MedianPerformanceResults %>%
    dplyr::filter(name %in% c("CalSlope")) %>%
    dplyr::mutate(name = fct_recode(name, 
                                    "Calibration Slope" = "CalSlope")) %>%
    ggplot(aes(y = value, x = N_dev,
               group = DevDataImputationMethod,
               colour = DevDataImputationMethod)) +
    geom_line() +
    geom_hline(data = data.frame(name = c("Calibration Slope"),
                                 value = c(0.9)),
               aes(yintercept = value),
               linetype = "dotted") +
    geom_hline(data = data.frame(name = c("Calibration Slope"),
                                 value = c(1)),
               aes(yintercept = value),
               linetype = "dashed") +
    xlab("Development Dataset Sample Size") +
    ylab("Median Calibration Slope") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank(),
          legend.position = "top") +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                  "#0072B2", "#D55E00", "#CC79A7"))
  ,
  AssuranceResults %>%
    ggplot(aes(y = ProbCalSlope, x = N_dev,
               group = DevDataImputationMethod,
               colour = DevDataImputationMethod)) +
    geom_line() +
    xlab("Development Dataset Sample Size") +
    ylab("Prob(Calibration Slope between 0.9 and 1.1)") +
    theme_bw(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
          legend.title=element_blank(),
          legend.position = "top") +
    scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                  "#0072B2", "#D55E00", "#CC79A7"))
  ,
  common.legend = TRUE,
  labels = c("A", "B")
)


MedianPerformanceResults %>%
  dplyr::filter(name %in% c("AUC")) %>%
  ggplot(aes(y = value, x = N_dev,
             group = DevDataImputationMethod,
             colour = DevDataImputationMethod)) +
  geom_line() +
  xlab("Development Dataset Sample Size") +
  ylab("Median C-statistic") +
  theme_bw(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        legend.title=element_blank(),
        legend.position = "right") +
  scale_color_manual(values = c("#E69F00", "#56B4E9", "#009E73", "#F0E442", 
                                "#0072B2", "#D55E00", "#CC79A7"))
