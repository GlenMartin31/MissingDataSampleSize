#Script to combine results from separate files on CSF
library(tidyverse)

result_list <- list.files(pattern = "simulation_results_") 

results_files <- result_list %>%
  purrr::map(read_rds)
names(results_files) <- paste("scenario", readr::parse_number(result_list),
                              sep = "_")

True_mod_performance <- purrr::map_dfr(results_files,
                                       ~dplyr::bind_rows(.x$True_mod_perfm),
                                       .id = "Scenario") 

prediction_results <- purrr::map_dfr(results_files,
                                     ~dplyr::bind_rows(.x$prediction_results),
                                     .id = "Scenario") 

model_coefs <- purrr::map_dfr(results_files,
                              ~dplyr::bind_rows(.x$model_coefs),
                              .id = "Scenario") %>%
  dplyr::group_by(Scenario, ModellingMethod, Variable) %>%
  dplyr::summarise(dplyr::across(CCA:FullyObserved,
                                 list("Mean" = ~mean(.),
                                      "SD" = ~sqrt(var(.)),
                                      "Min" = ~min(.),
                                      "Max" = ~max(.))),
                   .groups = "drop")

EVPI <- purrr::map_dfr(results_files,
                       ~dplyr::bind_rows(.x$EVPI),
                       .id = "Scenario") %>%
  dplyr::group_by(Scenario, IPDImputationMethod, DevelopmentImputationMethod, Z) %>%
  dplyr::summarise(dplyr::across(EVPI_logistic:REVPI_lasso,
                                 list("Mean" = ~mean(.),
                                      "SD" = ~sqrt(var(.)),
                                      "Min" = ~min(.),
                                      "Max" = ~max(.))),
                   .groups = "drop")

SampleSizeInfo <- purrr::map_dfr(results_files,
                                 .f = function(X) {
                                   data.frame("N_dev" = X$N_dev,
                                              "required_SS" = X$required_SS$sample_size)
                                 }) %>%
  dplyr::mutate("Scenario" = readr::parse_number(result_list),
                .before = "N_dev")


write_rds(True_mod_performance, 
          file = "./True_mod_performance.RDS")
write_rds(prediction_results, 
          file = "./prediction_results.RDS")
write_rds(model_coefs, 
          file = "./model_coefs.RDS")
write_rds(EVPI, 
          file = "./EVPI.RDS")
write_rds(SampleSizeInfo, 
          file = "./SampleSizeInfo.RDS")

warnings()