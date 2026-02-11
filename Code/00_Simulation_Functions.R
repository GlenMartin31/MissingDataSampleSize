# ##############################################################################

# Author of code: Glen P. Martin.

# This is code for a simulation study presented in a manuscript entitled: 
# Impact of Missing Data on Sample Size Requirements for Developing Clinical 
# Prediction Models
# Authors:
#   Glen P. Martin
#   Richard D. Riley

# ##############################################################################

####----------------------------------------------------------------------------
## Define a function to repeat the simulation across all iterations 
####----------------------------------------------------------------------------
simulation_nrun_fnc <- function(n_iter,
                                P, 
                                beta_X, 
                                rho_X,
                                Prev_Y, 
                                prop_missing,
                                missing_mech = c("MCAR", "MAR", "MNAR"),
                                RR_multiplier) {
  #Input: n_iter = the number of iterations to repeat the simulation over
  #       P = the number of predictor variables
  #       beta_X = association between Xs and prob of outcome
  #       rho_X = correlation between pairwise covariates
  #       Prev_Y = the overall event rate in the population
  #       prop_missing = Overall missingness proportion 
  #       missing_mech = type of missingness mechanism - MCAR, MAR or MNAR
  #       RR_multiplier = multiplier to apply relative to Riley's baseline N
  
  library(MASS)
  library(tidyverse)
  library(mice)
  library(glmnet)
  library(pROC)
  library(pmsampsize)
  library(missForestPredict)
  library(predRupdate)
  
  missing_mech <- as.character(match.arg(missing_mech))
  
  ##### Input Checking
  if(n_iter <= 0) {stop("n_iter should be positive")}
  if(P <= 0) {stop("P should be positive")}
  if(length(beta_X) != P) {stop("Length of beta_X is not P")}
  if(rho_X < 0 | rho_X > 1) {stop("rho_X should be between 0 and 1")}
  if(Prev_Y <= 0 |
     Prev_Y>=1) {stop("Prev_Y not between 0 and 1")}
  if(prop_missing <= 0 |
     prop_missing>=1) {stop("prop_missing not between 0 and 1")}
  if(RR_multiplier <= 0) {stop("RR_multiplier should be positive")}
  
  ##############################################################################
  ##### Generate the overarching large target population
  ##############################################################################
  IPD <- data_generating_fnc(P = P,
                             N = 500000,
                             beta_X = beta_X,
                             rho_X = rho_X,
                             Prev_Y = Prev_Y)
  
  #estimate the predictive performance of the 'true' model - useful for 
  #later calculations:
  true_model_performance <- predRupdate::pred_val_probs(binary_outcome = IPD$Y,
                                                        Prob = IPD$True_Pi,
                                                        cal_plot = FALSE)
  #note: this is also equivalent to the performance of the DGM on imputed
  #versions of the target population (regardless of missingness pattern or
  #imputation approach) because the 'true' risks wouldnt change and we dont
  #consider missing outcomes here
  
  ##############################################################################
  ##### Calculate the minimum required sample size for developing a CPM within
  ##### this target population:
  ##############################################################################
  required_SS <- pmsampsize::pmsampsize(type = "b", 
                                        parameters = P, 
                                        nagrsquared = 0.15, 
                                        shrinkage = 0.9,
                                        prevalence = Prev_Y)
  #apply the multiplier
  N_dev <- ceiling(required_SS$sample_size*RR_multiplier)
  
  ##############################################################################
  ##### Run the simulation scenario across all iterations
  ##############################################################################
  #intialise some variables to store results across iterations:
  prediction_results <- NULL
  model_coefs <- NULL
  EVPI <- NULL
  
  #calculate some useful values for EVPI calculations later in code:
  z <- seq(from = 0, to = 0.99, length.out = 100)
  NB_z <- IPD$True_Pi - outer((1-IPD$True_Pi), 
                              z/(1-z), '*') 
  NB_all <- colMeans(NB_z)
  NB_max <- colMeans((IPD$True_Pi > matrix(rep(z, each = length(IPD$True_Pi)), 
                                           nrow = length(IPD$True_Pi), 
                                           ncol = length(z)))*NB_z)
  for(iter in 1:n_iter) {
    # Generate a fully observed development cohort using the same DGM
    fully_observed_dev_data <- data_generating_fnc(P = P,
                                                   N = N_dev,
                                                   beta_X = beta_X,
                                                   rho_X = rho_X,
                                                   Prev_Y = Prev_Y)
    #run the simulation for this iteration:
    iteration_results <- simulation_singlerun_fnc(IPD = IPD,
                                                  fully_observed_dev_data = fully_observed_dev_data,
                                                  P = P,
                                                  prop_missing = prop_missing,
                                                  missing_mech = missing_mech,
                                                  true_model_performance = true_model_performance,
                                                  z = z,
                                                  NB_z = NB_z,
                                                  NB_all = NB_all,
                                                  NB_max = NB_max)
    
    #store the results from this iteration:
    prediction_results <- prediction_results %>% 
      dplyr::bind_rows(iteration_results$predictive_performance %>%
                         dplyr::mutate("iter" = iter,
                                       .before = "ModellingMethod"))
    model_coefs <- model_coefs %>%
      dplyr::bind_rows(iteration_results$model_coef_mat %>%
                         dplyr::mutate("iter" = iter,
                                       .before = "Variable"))
    EVPI <- EVPI %>%
      dplyr::bind_rows(iteration_results$EVPI %>%
                         dplyr::mutate("iter" = iter,
                                       .before = "IPDImputationMethod"))
  }
  
  return(list("True_mod_perfm" = data.frame("OE" = true_model_performance$OE_ratio,
                                            "CalInt" = true_model_performance$CalInt,
                                            "CalSlope" = true_model_performance$CalSlope,
                                            "AUC" = true_model_performance$AUC),
              "prediction_results" = prediction_results,
              "model_coefs" = model_coefs,
              "EVPI" = EVPI,
              "N_dev" = N_dev,
              "required_SS" = required_SS))
}

####----------------------------------------------------------------------------
## Function that runs a single run of the simulation
####----------------------------------------------------------------------------
simulation_singlerun_fnc <- function(IPD,
                                     fully_observed_dev_data,
                                     P,
                                     prop_missing,
                                     missing_mech,
                                     true_model_performance,
                                     z,
                                     NB_z,
                                     NB_all,
                                     NB_max) {
  
  ##############################################################################
  ##### Ampute missing data into the target population data (IPD)
  ##############################################################################
  miss_patterns <- expand.grid(rep(list(0:1), P))
  miss_patterns <- miss_patterns[-c(1,nrow(miss_patterns)),]
  #we dont want all possible combinations of missingness patterns (not realistic), 
  #so randomly sample 100 of them in each iteration:
  miss_patterns <- miss_patterns %>% 
    dplyr::slice_sample(n = 100) 
  
  amp <- mice::ampute(data = IPD %>%
                        dplyr::select(starts_with("X")),
                      prop = prop_missing,
                      mech = missing_mech,
                      patterns = miss_patterns)
  
  IPD_with_missing <- dplyr::bind_cols("ID" = IPD$ID,
                                       amp$amp,
                                       "True_Pi" = IPD$True_Pi,
                                       "Y" = IPD$Y)
  rm(amp)
  
  ##############################################################################
  ##### Ampute missing data into the development dataset
  ##############################################################################
  amp <- mice::ampute(data = fully_observed_dev_data %>%
                        dplyr::select(starts_with("X")),
                      prop = prop_missing,
                      mech = missing_mech,
                      patterns = miss_patterns)

  dev_data <- dplyr::bind_cols("ID" = fully_observed_dev_data$ID,
                               amp$amp,
                               "True_Pi" = fully_observed_dev_data$True_Pi,
                               "Y" = fully_observed_dev_data$Y)
  rm(amp)
  
  ##############################################################################
  ##### Impute the dev_data using different imputation strategies, and 
  ##### transport these imputation models to apply to impute IPD_with_missing
  ##############################################################################
  imputation_objects <- imputation_fnc(df = dev_data)
  
  imputed_dev_data <- list("CCA" = imputation_objects$CCA,
                           "MeanImpute" = imputation_objects$MeanImpute,
                           "RI" = imputation_objects$RI,
                           "MI" = imputation_objects$MI,
                           "RF" = imputation_objects$RF)
  
  #apply the learned imputation models onto the target population:
  imputed_IPD_data <- list(
    "MeanImpute" = IPD_with_missing %>%
      dplyr::mutate(dplyr::across(.cols = tidyr::starts_with("X"),
                                  .fns = ~ifelse(is.na(.), mean(., na.rm = T), .))),
    "RI" = mice::complete(mice::mice.mids(imputation_objects$RI_impmodels, 
                                          newdata = IPD_with_missing,
                                          print = F)),
    "MI" = mice::mice.mids(imputation_objects$MI, 
                           newdata = IPD_with_missing,
                           print = F),
    "RF" = missForestPredict::missForestPredict(imputation_objects$RF_impmodels,
                                                newdata = IPD_with_missing)
  )
  
  ##############################################################################
  ##### Fit the CPMs under MLE, AIC and LASSO likelihoods to each imputed 
  ##### version of the development dataset, and the fully observed version
  ##############################################################################
  fitted_models <- purrr::set_names(c("logistic", "AIC", "lasso")) %>%
    purrr::map(.x = .,
               .f = function(X) {
                 model_fit <- lapply(imputed_dev_data, 
                                     model_fit_fnc, 
                                     P = P,
                                     model_type = X)
                 model_fit$FullyObserved <- model_fit_fnc(df = fully_observed_dev_data,
                                                          P = P,
                                                          model_type = X)
                 model_fit
                 }) %>% 
    purrr::list_flatten()
  
  ## Extract the coefficients of each model
  model_coef_mat <- fitted_models %>% 
    dplyr::bind_rows(.id = "imputation") %>% 
    tidyr::separate_wider_delim(cols = "imputation",
                                delim = "_",
                                names = c("ModellingMethod", "imputation")) %>%
    tidyr::pivot_wider(id_cols = c("Variable", "ModellingMethod"), 
                       names_from = c("imputation"), 
                       values_from = "Estimate") 
  
  
  ##############################################################################
  ##### Quantify every CPM's predictive performance, EVPI and degradation in
  ##### the target population dataset (both fully observed and imputed versions)
  ##############################################################################
  FullyObervedIPD_Performance <- purrr::map(.x = fitted_models,
                                            .f = predictive_performance_fnc,
                                            df = IPD,
                                            true_model_performance = true_model_performance,
                                            z = z,
                                            NB_z = NB_z,
                                            NB_all = NB_all,
                                            NB_max = NB_max)
  
  MeanImputeIPD_Performance <- purrr::map(.x = fitted_models[c("logistic_MeanImpute",
                                                               "AIC_MeanImpute",
                                                               "lasso_MeanImpute")],
                                          .f = predictive_performance_fnc,
                                          df = imputed_IPD_data$MeanImpute,
                                          true_model_performance = true_model_performance,
                                          z = z,
                                          NB_z = NB_z,
                                          NB_all = NB_all,
                                          NB_max = NB_max)
  
  RIIPD_Performance <- purrr::map(.x = fitted_models[c("logistic_RI",
                                                       "AIC_RI",
                                                       "lasso_RI")],
                                  .f = predictive_performance_fnc,
                                  df = imputed_IPD_data$RI,
                                  true_model_performance = true_model_performance,
                                  z = z,
                                  NB_z = NB_z,
                                  NB_all = NB_all,
                                  NB_max = NB_max)
  
  MIIPD_Performance <- purrr::map(.x = fitted_models[c("logistic_MI",
                                                       "AIC_MI",
                                                       "lasso_MI")],
                                  .f = predictive_performance_fnc,
                                  df = imputed_IPD_data$MI,
                                  true_model_performance = true_model_performance,
                                  z = z,
                                  NB_z = NB_z,
                                  NB_all = NB_all,
                                  NB_max = NB_max)
  
  RFIPD_Performance <- purrr::map(.x = fitted_models[c("logistic_RF",
                                                       "AIC_RF",
                                                       "lasso_RF")],
                                  .f = predictive_performance_fnc,
                                  df = imputed_IPD_data$RF,
                                  true_model_performance = true_model_performance,
                                  z = z,
                                  NB_z = NB_z,
                                  NB_all = NB_all,
                                  NB_max = NB_max)
  
  ##############################################################################
  ##### Extract and data wrangle the results ready for function return
  ##############################################################################
  predictive_performance <- lapply(FullyObervedIPD_Performance, 
                                   function(X)X["performance"]) %>% 
    purrr::list_flatten() %>% 
    dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
    tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                delim = "_",
                                names = c("ModellingMethod", 
                                          "DevelopmentImputationMethod",
                                          "Metric")) %>%
    dplyr::select(-Metric) %>%
    dplyr::mutate("IPDImputationMethod" = "FullyObserved",
                  .after = "DevelopmentImputationMethod") %>%
    dplyr::bind_rows(lapply(MeanImputeIPD_Performance, 
                            function(X)X["performance"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "MeanImpute",
                                     .after = "DevelopmentImputationMethod")) %>%
    dplyr::bind_rows(lapply(RIIPD_Performance, 
                            function(X)X["performance"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "RI",
                                     .after = "DevelopmentImputationMethod")) %>%
    dplyr::bind_rows(lapply(MIIPD_Performance, 
                            function(X)X["performance"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "MI",
                                     .after = "DevelopmentImputationMethod")) %>%
    dplyr::bind_rows(lapply(RFIPD_Performance, 
                            function(X)X["performance"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "RF",
                                     .after = "DevelopmentImputationMethod"))
  
  
  
  EVPI <- lapply(FullyObervedIPD_Performance, 
         function(X)X["EVPI"]) %>% 
    purrr::list_flatten() %>% 
    dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
    tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                delim = "_",
                                names = c("ModellingMethod", 
                                          "DevelopmentImputationMethod",
                                          "Metric")) %>%
    dplyr::select(-Metric) %>%
    dplyr::mutate("IPDImputationMethod" = "FullyObserved",
                  .after = "DevelopmentImputationMethod") %>%
    pivot_wider(id_cols = c("DevelopmentImputationMethod", "IPDImputationMethod", "Z"),
                names_from = c("ModellingMethod"),
                values_from = c("EVPI", "REVPI")) %>%
    dplyr::bind_rows(lapply(MeanImputeIPD_Performance, 
                            function(X)X["EVPI"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "MeanImpute",
                                     .after = "DevelopmentImputationMethod") %>%
                       pivot_wider(id_cols = c("DevelopmentImputationMethod", "IPDImputationMethod", "Z"),
                                   names_from = c("ModellingMethod"),
                                   values_from = c("EVPI", "REVPI"))) %>%
    dplyr::bind_rows(lapply(RIIPD_Performance, 
                            function(X)X["EVPI"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "RI",
                                     .after = "DevelopmentImputationMethod") %>%
                       pivot_wider(id_cols = c("DevelopmentImputationMethod", "IPDImputationMethod", "Z"),
                                   names_from = c("ModellingMethod"),
                                   values_from = c("EVPI", "REVPI"))) %>%
    dplyr::bind_rows(lapply(MIIPD_Performance, 
                            function(X)X["EVPI"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "MI",
                                     .after = "DevelopmentImputationMethod") %>%
                       pivot_wider(id_cols = c("DevelopmentImputationMethod", "IPDImputationMethod", "Z"),
                                   names_from = c("ModellingMethod"),
                                   values_from = c("EVPI", "REVPI"))) %>%
    dplyr::bind_rows(lapply(RFIPD_Performance, 
                            function(X)X["EVPI"]) %>% 
                       purrr::list_flatten() %>% 
                       dplyr::bind_rows(.id = "Model_DevelopmentImputationMethod_Metric") %>%
                       tidyr::separate_wider_delim(cols = "Model_DevelopmentImputationMethod_Metric",
                                                   delim = "_",
                                                   names = c("ModellingMethod", 
                                                             "DevelopmentImputationMethod",
                                                             "Metric")) %>%
                       dplyr::select(-Metric) %>%
                       dplyr::mutate("IPDImputationMethod" = "RF",
                                     .after = "DevelopmentImputationMethod") %>%
                       pivot_wider(id_cols = c("DevelopmentImputationMethod", "IPDImputationMethod", "Z"),
                                   names_from = c("ModellingMethod"),
                                   values_from = c("EVPI", "REVPI")))
  
  
  ##############################################################################
  ##### Function return
  ##############################################################################
  return(list("model_coef_mat" = model_coef_mat,
              "predictive_performance" = predictive_performance,
              "EVPI" = EVPI))
}

####----------------------------------------------------------------------------
## Function that generates data according to the data-generating models
####----------------------------------------------------------------------------
data_generating_fnc <- function(P,
                                N,
                                beta_X,
                                rho_X,
                                Prev_Y) {
  # Inputs: P = number of covaraites to simulate
  #         N = number of observations to simulate
  #         beta_X = 'true' association between the covariates and outcome
  #         rho_X = correlation between pairwise covariates
  #         Prev_Y = prevalence of the binary outcome to simulate
  
  if(length(beta_X) != P) {stop("Length of beta_X is not P")}
  if(rho_X < 0 | rho_X > 1) {stop("rho_X should be between 0 and 1")}
  
  Sigma_X <- diag(1, ncol = P, nrow = P)
  #Create pair-wise correlations:
  Sigma_X[(row(Sigma_X) - col(Sigma_X)) == 1 | 
            (row(Sigma_X) - col(Sigma_X)) == -1][c(TRUE, TRUE, 
                                                   FALSE, FALSE)] <- rho_X
  
  X <- MASS::mvrnorm(n = N, 
                     mu = rep(0, P),
                     Sigma = Sigma_X)
  LP <- as.numeric(X %*% beta_X)
  beta_0 <- as.numeric(coef(glm(rbinom(N, 1, prob = Prev_Y) ~ offset(LP), 
                                family = binomial(link = "logit")))[1])
  LP <- beta_0 + LP
  Pi <- exp(LP) / (1+exp(LP))
  Y <- rbinom(N, 1, Pi)
  
  IPD <- tibble::tibble(dplyr::bind_cols(data.frame(X))) %>%
    dplyr::mutate(True_Pi = Pi,
                  Y = Y)  %>%
    dplyr::mutate("ID" = 1:n(),
                  .before = "X1")
  
  IPD
}

####----------------------------------------------------------------------------
## Data imputation function
####----------------------------------------------------------------------------
imputation_fnc <- function(df) {
  # Inputs: df = data on which we wish to impute
  
  CCA_df <- df %>%
    dplyr::select("Y",
                  tidyr::starts_with("X")) %>%
    na.omit()
  
  MeanImpute_df <- df %>%
    dplyr::select("Y",
                  tidyr::starts_with("X")) %>%
    dplyr::mutate(dplyr::across(.cols = tidyr::starts_with("X"),
                                .fns = ~ifelse(is.na(.), mean(., na.rm = T), .)))
  
  #Regression imputation of missing data:
  RI_imp <- mice::mice(df, m = 1, maxit = 0)
  predmat <- RI_imp$predictorMatrix
  predmat["ID", ] <- predmat[, "ID"] <- 0
  predmat["True_Pi", ] <- predmat[, "True_Pi"] <- 0
  predmat["Y", ] <- predmat[, "Y"] <- 0
  RI_imp <- mice::mice(df, 
                       m = 1, 
                       method = "norm.predict", 
                       predictorMatrix = predmat,
                       maxit = 15,
                       printFlag = 0)
  RI_df <- mice::complete(RI_imp) %>%
    dplyr::select("Y", tidyr::starts_with("X"))
  
  #Multiple imputation of missing data:
  MI_imp <- mice::mice(df, m = 20, maxit = 0)
  predmat <- MI_imp$predictorMatrix
  predmat["ID", ] <- predmat[, "ID"] <- 0
  predmat["True_Pi", ] <- predmat[, "True_Pi"] <- 0
  MI_imp <- mice::mice(df, 
                       m = 20, 
                       method = "norm", 
                       predictorMatrix = predmat,
                       maxit = 15,
                       printFlag = 0)
  
  # Random forest imputation
  predmat <- missForestPredict::create_predictor_matrix(df)
  predmat["ID", ] <- predmat[, "ID"] <- 0
  predmat["True_Pi", ] <- predmat[, "True_Pi"] <- 0
  predmat["Y", ] <- predmat[, "Y"] <- 0
  RF_imputation <- missForestPredict::missForest(df,
                                                 save_models = TRUE,
                                                 predictor_matrix = predmat,
                                                 num.trees = 200,
                                                 num.threads = 2,
                                                 verbose = FALSE)
  RF_df <- RF_imputation$ximp %>%
    dplyr::select("Y", tidyr::starts_with("X"))
                         
  return(list("CCA" = CCA_df,
              "MeanImpute" = MeanImpute_df,
              "RI" = RI_df,
              "MI" = MI_imp,
              "RF" = RF_df,
              
              "RI_impmodels" = RI_imp,
              "RF_impmodels" = RF_imputation))
}

####----------------------------------------------------------------------------
## Model fitting function
####----------------------------------------------------------------------------
model_fit_fnc <- function(df,
                          P,
                          model_type = c("logistic", "AIC", "lasso")) {
  # Inputs: df = data on which we wish to fit the CPM to
  #         P = number of predictor variables
  
  # Internal functions to fit the specified model to a given dataset
  model_type <- as.character(match.arg(model_type))
  if (model_type == "logistic") {
    mod_fnc <- function(df) {
      logistic_mod <- glm(Y ~.,
                          data = df %>%
                            dplyr::select("Y",
                                          tidyr::starts_with("X")),
                          family = binomial(link = "logit"))
      
      data.frame("Variable" = names(coef(logistic_mod)),
                 "Estimate" = as.numeric(coef(logistic_mod)))
    }
  } else if (model_type == "AIC") {
    mod_fnc <- function(df) {
      logistic_mod <- glm(Y ~.,
                          data = df %>%
                            dplyr::select("Y",
                                          tidyr::starts_with("X")),
                          family = binomial(link = "logit"))
      
      AIC_mod <- step(logistic_mod, direction = "backward", trace = FALSE)
      
      data.frame("Variable" = names(coef(logistic_mod))) %>%
        dplyr::left_join(data.frame("Variable" = names(coef(AIC_mod)),
                                    "Estimate" = as.numeric(coef(AIC_mod))),
                         by = "Variable")%>%
        tidyr::replace_na(list("Estimate" = 0))
    }
  } else if (model_type == "lasso") {
    mod_fnc <- function(df) {
      lasso_mod <- cv.glmnet(x = df %>%
                               dplyr::select(tidyr::starts_with("X")) %>%
                               data.matrix(),
                             y = df %>% dplyr::pull("Y"),
                             family = "binomial")
      
      data.frame("Variable" = rownames(coef(lasso_mod, 
                                            s = "lambda.min")),
                 "Estimate" = as.numeric(coef(lasso_mod, 
                                              s = "lambda.min")))
    } 
  }
  
  if(is.mids(df)) {
    MI_df <- mice::complete(df, action = "long")
    
    mod_coefs <- MI_df %>%
      dplyr::group_by(.imp) %>%
      tidyr::nest() %>%
      dplyr::mutate("MI_coefs" = map(data, mod_fnc)) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = "MI_coefs") %>%
      dplyr::mutate(Variable = factor(Variable,
                                      levels = c("(Intercept)",
                                                 paste("X", 1:P, sep = "")))) %>%
      dplyr::group_by(Variable) %>%
      dplyr::summarise(Estimate = mean(Estimate)) %>%
      as.data.frame()
  } else {
    
    mod_coefs <- mod_fnc(df) %>%
      dplyr::mutate(Variable = factor(Variable,
                                      levels = c("(Intercept)",
                                                 paste("X", 1:P, sep = ""))))
  }
  mod_coefs
}

####----------------------------------------------------------------------------
## Predictive performance function
####----------------------------------------------------------------------------
predictive_performance_fnc <- function(model_info, 
                                       df,
                                       true_model_performance,
                                       z,
                                       NB_z,
                                       NB_all,
                                       NB_max) {
  if(mice::is.mids(df)) {
    
    MI_long <- mice::complete(df, action = "long")
    
    coef_vals <- model_info$Estimate
    predictions <- MI_long %>% 
      dplyr::group_by(.imp) %>% 
      tidyr::nest() %>%
      dplyr::mutate("DM" = map(data, 
                               function(X){
                                 cbind(1, X %>% 
                                         dplyr::select(tidyr::starts_with("X")) %>% 
                                         data.matrix())
                               })) %>%
      dplyr::mutate("Pi" = map(DM,
                               function(X){
                                 1 / (1 + exp(-as.numeric(X %*% coef_vals)))
                               })) %>% 
      tidyr::unnest(cols = c(".imp", "data", "Pi")) %>% 
      dplyr::select(.imp, ID, Y, Pi)
    
    performance <- predictions %>%
      dplyr::group_by(.imp) %>%
      tidyr::nest() %>%
      dplyr::mutate("Performance" = purrr::map(data, 
                                               function(X) {
                                                 performance <- predRupdate::pred_val_probs(binary_outcome = X$Y,
                                                                                            Prob = X$Pi,
                                                                                            cal_plot = FALSE)
                                                 performance <- data.frame("OE" = performance$OE_ratio,
                                                                           "CalInt" = performance$CalInt,
                                                                           "CalSlope" = performance$CalSlope,
                                                                           "AUC" = performance$AUC) %>%
                                                   dplyr::mutate("OE_degradation" = OE - true_model_performance$OE_ratio,
                                                                 "CalInt_degradation" = CalInt - true_model_performance$CalInt,
                                                                 "CalSlope_degradation" = CalSlope - true_model_performance$CalSlope,
                                                                 "AUC_degradation" = AUC - true_model_performance$AUC)
                                               })) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = c(".imp", "Performance")) %>%
      dplyr::ungroup() %>%
      dplyr::summarise(dplyr::across(OE:AUC_degradation,
                                     ~mean(.)))
    
    EVPI <- predictions %>%
      dplyr::group_by(.imp) %>%
      tidyr::nest() %>%
      dplyr::mutate("EVPI" = purrr::map(data, 
                                        function(X,
                                                 NB_max,
                                                 NB_all,
                                                 z) {
                                          NB_model <- colMeans((X$Pi > matrix(rep(z, each = length(X$Pi)),
                                                                              nrow = length(X$Pi), 
                                                                              ncol = length(z)))*NB_z)
                                          EVPI <- dplyr::bind_cols("NB_model" = NB_model,
                                                                   "NB_max" = NB_max,
                                                                   "NB_all" = NB_all,
                                                                   "Z" = z) %>%
                                            dplyr::mutate("EVPI" = NB_max-pmax(0,NB_model,NB_all),
                                                          "REVPI" = 100*((pmax(0,NB_model,NB_all) - pmax(0,NB_all))/
                                                                           (NB_max - pmax(0,NB_all))),
                                                          "Z" = Z,
                                                          .keep = "none") 
                                        },
                                        NB_max = NB_max,
                                        NB_all = NB_all,
                                        z = z)) %>%
      dplyr::select(-data) %>%
      tidyr::unnest(cols = c(".imp", "EVPI")) %>%
      dplyr::group_by(Z) %>%
      dplyr::summarise(dplyr::across(c(EVPI, REVPI),
                                     ~mean(.)))
    
  } else {
    coef_vals <- model_info$Estimate
    DM <- cbind(1, df %>% 
                  dplyr::select(tidyr::starts_with("X")) %>% 
                  data.matrix())
    predictions <- 1 / (1 + exp(-as.numeric(DM %*% coef_vals)))
    
    performance <- predRupdate::pred_val_probs(binary_outcome = df$Y,
                                               Prob = predictions,
                                               cal_plot = FALSE)
    performance <- data.frame("OE" = performance$OE_ratio,
                              "CalInt" = performance$CalInt,
                              "CalSlope" = performance$CalSlope,
                              "AUC" = performance$AUC) %>%
      dplyr::mutate("OE_degradation" = OE - true_model_performance$OE_ratio,
                    "CalInt_degradation" = CalInt - true_model_performance$CalInt,
                    "CalSlope_degradation" = CalSlope - true_model_performance$CalSlope,
                    "AUC_degradation" = AUC - true_model_performance$AUC)
    
    NB_model <- colMeans((predictions > matrix(rep(z, each = length(predictions)),
                                               nrow = length(predictions), 
                                               ncol = length(z)))*NB_z)
    EVPI <- dplyr::bind_cols("NB_model" = NB_model,
                             "NB_max" = NB_max,
                             "NB_all" = NB_all,
                             "Z" = z) %>%
      dplyr::mutate("EVPI" = NB_max-pmax(0,NB_model,NB_all),
                    "REVPI" = 100*((pmax(0,NB_model,NB_all) - pmax(0,NB_all))/
                                     (NB_max - pmax(0,NB_all))),
                    "Z" = Z,
                    .keep = "none") 
    
  }
  return(list("performance" = performance,
              "EVPI" = EVPI))
}



