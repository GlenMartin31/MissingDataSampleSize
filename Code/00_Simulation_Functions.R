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
    iteration_results <- simulation_singlerun_fnc(IPD = IPD,
                                                  N_dev = N_dev,
                                                  P = P,
                                                  prop_missing = prop_missing,
                                                  missing_mech = missing_mech,
                                                  true_model_performance = true_model_performance,
                                                  z = z,
                                                  NB_z = NB_z,
                                                  NB_all = NB_all,
                                                  NB_max = NB_max)
    
    prediction_results <- prediction_results %>% 
      dplyr::bind_rows(iteration_results$predictive_performance %>%
                         dplyr::mutate("iter" = iter,
                                       .before = "Model"))
    
    model_coefs <- model_coefs %>%
      dplyr::bind_rows(iteration_results$model_coef_mat %>%
                         dplyr::mutate("iter" = iter,
                                       .before = "Model"))
    
    EVPI <- EVPI %>%
      dplyr::bind_rows(iteration_results$EVPI %>%
                         dplyr::mutate("iter" = iter,
                                       .before = "Model"))
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
                                     N_dev,
                                     P,
                                     prop_missing,
                                     missing_mech,
                                     true_model_performance,
                                     z,
                                     NB_z,
                                     NB_all,
                                     NB_max) {
  
  ##############################################################################
  ##### Sample a fully observed development cohort from the IPD
  ##############################################################################
  fully_observed_dev_data <- IPD %>%
    dplyr::slice_sample(n = N_dev, replace = FALSE) %>%
    dplyr::arrange(ID)
  
  ##############################################################################
  ##### Create a version of the development data with missingness amputed
  ##############################################################################
  miss_patterns <- expand.grid(rep(list(0:1), P))
  miss_patterns <- miss_patterns[-c(1,nrow(miss_patterns)),]
  #given the sizes of the development data, we dont want all possible 
  #combinations of missingness patterns, so randomly sample 100 of them in 
  #each iteration:
  miss_patterns <- miss_patterns %>% 
    dplyr::slice_sample(n = 100) 
  amp <- mice::ampute(data = fully_observed_dev_data %>%
                        dplyr::select(starts_with("X")),
                      prop = prop_missing,
                      mech = missing_mech,
                      patterns = miss_patterns)

  dev_data <- dplyr::bind_cols("ID" = fully_observed_dev_data$ID,
                               amp$amp,
                               "True_Pi" = fully_observed_dev_data$True_Pi,
                               "Y" = fully_observed_dev_data$Y)
  
  
  ##############################################################################
  ##### Impute the dev_data using different imputation strategies
  ##############################################################################
  imputed_dfs <- imputation_fnc(df = dev_data)
  
  ##############################################################################
  ##### Fit the CPMs under MLE, AIC and LASSO likelihoods to each imputed 
  ##### version of the development dataset, and the fully observed version
  ##############################################################################
  logistic_models <- lapply(imputed_dfs, model_fit_fnc, 
                            P = P,
                            model_type = "logistic")
  logistic_models$FullyObserved <- model_fit_fnc(df = fully_observed_dev_data,
                                                     P = P,
                                                     model_type = "logistic")
  
  AIC_models <- lapply(imputed_dfs, model_fit_fnc, 
                       P = P,
                       model_type = "AIC")
  AIC_models$FullyObserved <- model_fit_fnc(df = fully_observed_dev_data,
                                            P = P,
                                            model_type = "AIC")
  
  lasso_models <- lapply(imputed_dfs, model_fit_fnc, 
                         P = P,
                         model_type = "lasso")
  lasso_models$FullyObserved <- model_fit_fnc(df = fully_observed_dev_data,
                                                  P = P,
                                                  model_type = "lasso")
  
  ## Extract the coefficients of each model
  logistic_coefs_mat <- dplyr::bind_rows(logistic_models, .id = "imputation") %>% 
    tidyr::pivot_wider(id_cols = "Variable", 
                       names_from = "imputation", values_from = "Estimate")
  AIC_coefs_mat <- dplyr::bind_rows(AIC_models, .id = "imputation") %>% 
    tidyr::pivot_wider(id_cols = "Variable", 
                       names_from = "imputation", values_from = "Estimate")
  lasso_coefs_mat <- dplyr::bind_rows(lasso_models, .id = "imputation") %>% 
    tidyr::pivot_wider(id_cols = "Variable", 
                       names_from = "imputation", values_from = "Estimate")
  
  
  ##############################################################################
  ##### Apply every developed CPM to make predictions for each individual in 
  ##### the target population dataset
  ##############################################################################
  IPD_logistic_predictions <- data.frame(sapply(logistic_models, 
                                                expit_fnc, df = IPD))
  IPD_AIC_predictions <- data.frame(sapply(AIC_models, 
                                           expit_fnc, df = IPD))
  IPD_lasso_predictions <- data.frame(sapply(lasso_models, 
                                             expit_fnc, df = IPD))
  
  ##############################################################################
  ##### Quantify every CPM's predictive performance and degradation in
  ##### the target population dataset
  ##############################################################################
  predictive_performance <- apply(IPD_logistic_predictions, 2, 
                                  performance_fnc, 
                                  outcome = IPD$Y,
                                  true_model_performance = true_model_performance) %>%
    dplyr::bind_rows(.id = "Imputation") %>%
    dplyr::mutate("Model" = "Logistic",
                  .before = "Imputation") %>%
    dplyr::bind_rows(dplyr::bind_rows(apply(IPD_AIC_predictions, 2, 
                                            performance_fnc, 
                                            outcome = IPD$Y,
                                            true_model_performance = true_model_performance),
                                      .id = "Imputation") %>%
                       dplyr::mutate("Model" = "AIC",
                                     .before = "Imputation")) %>%
    dplyr::bind_rows(dplyr::bind_rows(apply(IPD_lasso_predictions, 2, 
                                            performance_fnc, 
                                            outcome = IPD$Y,
                                            true_model_performance = true_model_performance),
                                      .id = "Imputation") %>%
                       dplyr::mutate("Model" = "LASSO",
                                     .before = "Imputation"))
  
  ##############################################################################
  ##### Expected Value of Perfect Information calculations
  ##############################################################################
  NB_models_logistic <- apply(IPD_logistic_predictions, 2, 
                              function(X) {
                                colMeans((X > matrix(rep(z, each = length(X)),
                                                     nrow = length(X), 
                                                     ncol = length(z)))*NB_z)
                                })
  colnames(NB_models_logistic) <- paste("NB_model_", colnames(NB_models_logistic),
                                        sep = "")
  NB_logistic <- dplyr::bind_cols(NB_models_logistic,
                                  "NB_max" = NB_max,
                                  "NB_all" = NB_all,
                                  "Z" = z) %>%
    dplyr::mutate(dplyr::across(starts_with("NB_model"),
                                function(X){NB_max-pmax(0,X,NB_all)},
                                .names = "{paste('EVPI_', col, sep = '')}"),
                  "Z" = Z,
                  .keep = "none")
  
  NB_models_AIC <- apply(IPD_AIC_predictions, 2, 
                         function(X) {
                           colMeans((X > matrix(rep(z, each = length(X)),
                                                nrow = length(X), 
                                                ncol = length(z)))*NB_z)
                         })
  colnames(NB_models_logistic) <- paste("NB_model_", colnames(NB_models_AIC),
                                        sep = "")
  NB_AIC <- dplyr::bind_cols(NB_models_logistic,
                             "NB_max" = NB_max,
                             "NB_all" = NB_all,
                             "Z" = z) %>%
    dplyr::mutate(dplyr::across(starts_with("NB_model"),
                                function(X){NB_max-pmax(0,X,NB_all)},
                                .names = "{paste('EVPI_', col, sep = '')}"),
                  "Z" = Z,
                  .keep = "none")
  
  NB_models_lasso <- apply(IPD_lasso_predictions, 2,
                           function(X) {
                             colMeans((X > matrix(rep(z, each = length(X)),
                                                  nrow = length(X), 
                                                  ncol = length(z)))*NB_z)
                             })
  colnames(NB_models_lasso) <- paste("NB_model_", colnames(NB_models_lasso),
                                     sep = "")
  NB_lasso <- dplyr::bind_cols(NB_models_lasso,
                               "NB_max" = NB_max,
                               "NB_all" = NB_all,
                               "Z" = z) %>%
    dplyr::mutate(dplyr::across(starts_with("NB_model"),
                                function(X){NB_max-pmax(0,X,NB_all)},
                                .names = "{paste('EVPI_', col, sep = '')}"),
                  "Z" = Z,
                  .keep = "none")
  
  ##############################################################################
  ##### Process results for function return
  ##############################################################################
  model_coef_mat <- logistic_coefs_mat %>%
    dplyr::mutate("Model" = "Logistic",
                  .before = "Variable") %>%
    dplyr::bind_rows(AIC_coefs_mat %>%
                       dplyr::mutate("Model" = "AIC",
                                     .before = "Variable")) %>%
    dplyr::bind_rows(lasso_coefs_mat %>%
                       dplyr::mutate("Model" = "LASSO",
                                     .before = "Variable"))
  
  EVPI <- NB_logistic %>%
    dplyr::mutate("Model" = "Logistic",
                  .before = "Z") %>%
    dplyr::bind_rows(NB_AIC %>%
                       dplyr::mutate("Model" = "AIC",
                                     .before = "Z")) %>%
    dplyr::bind_rows(NB_lasso %>%
                       dplyr::mutate("Model" = "LASSO",
                                     .before = "Z"))
  
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
  imp <- mice::mice(df, m = 1, maxit = 0)
  predmat <- imp$predictorMatrix
  predmat["True_Pi", ] <- predmat[, "True_Pi"] <- 0
  predmat["Y", ] <- predmat[, "Y"] <- 0
  imp <- mice::mice(df, 
                    m = 1, 
                    method = "norm.predict", 
                    predictorMatrix = predmat,
                    maxit = 15,
                    printFlag = 0)
  RI_df <- mice::complete(imp) %>%
    dplyr::select("Y", tidyr::starts_with("X"))
  
  #Multiple imputation of missing data:
  imp <- mice::mice(df, m = 20, maxit = 0)
  predmat <- imp$predictorMatrix
  predmat["True_Pi", ] <- predmat[, "True_Pi"] <- 0
  MI <- mice::mice(df, 
                   m = 20, 
                   method = "norm", 
                   predictorMatrix = predmat,
                   maxit = 15,
                   printFlag = 0)
  
  # Random forest imputation
  predmat <- missForestPredict::create_predictor_matrix(df)
  predmat["True_Pi", ] <- predmat[, "True_Pi"] <- 0
  predmat["Y", ] <- predmat[, "Y"] <- 0
  RF_imputation <- missForestPredict::missForest(df,
                                                 save_models = FALSE,
                                                 predictor_matrix = predmat,
                                                 num.trees = 200,
                                                 num.threads = 2,
                                                 verbose = FALSE)
  RF_df <- RF_imputation$ximp
                         
  return(list("CCA" = CCA_df,
              "MeanImpute" = MeanImpute_df,
              "RI" = RI_df,
              "MI" = MI,
              "RF" = RF_df))
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
performance_fnc <- function(outcome, prediction,
                            true_model_performance) {
  library(predRupdate)
  performance <- predRupdate::pred_val_probs(binary_outcome = outcome,
                                             Prob = prediction,
                                             cal_plot = FALSE)
  data.frame("OE" = performance$OE_ratio,
             "CalInt" = performance$CalInt,
             "CalSlope" = performance$CalSlope,
             "AUC" = performance$AUC) %>%
    dplyr::mutate("OE_degradation" = OE - true_model_performance$OE_ratio,
                  "CalInt_degradation" = CalInt - true_model_performance$CalInt,
                  "CalSlope_degradation" = CalSlope - true_model_performance$CalSlope,
                  "AUC_degradation" = AUC - true_model_performance$AUC)
}

####----------------------------------------------------------------------------
## util functions
####----------------------------------------------------------------------------
expit_fnc <- function(coefs, df) {
  DM <- cbind(1, df %>% 
                dplyr::select(tidyr::starts_with("X")) %>% 
                data.matrix())
  
  Estimate <- coefs$Estimate
  
  1 / (1 + exp(-as.numeric(DM %*% Estimate)))
}