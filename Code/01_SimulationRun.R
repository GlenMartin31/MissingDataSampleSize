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
## This script runs the simulations across all scenarios
####----------------------------------------------------------------------------

#Load the simulation functions
# source(here::here("Code", "00_Simulation_Functions.R")) #local run
source("./00_Simulation_Functions.R") #csf run

library(tidyverse)

# Define a dataset that includes all combinations of simulation parameters:
sims_parameters <- tidyr::crossing(
  P = c(10),
  beta_X = list(c(0.5, 0, 0.3, 0, 0.5, 0, 0.3, 0, 0.1, 0),
                c(0.5, 0.2, 0.3, 0.1, 0.5, 0.2, 0.3, 0.05, 0.1, 0.15)),
  rho_X = c(0, 0.5, 0.75),
  Prev_Y = c(0.2),
  prop_missing = c(0.1, 0.2, 0.4, 0.6),
  missing_mech = c("MAR", "MNAR"),
  multipliers = c(1, 1.25, 1.5, 1.75, 2)
)

####----------------------------------------------------------------------------
## Run the scenarios - uses the computational shared facility at 
## University of Manchester
####----------------------------------------------------------------------------

taskid <- commandArgs(trailingOnly = T) 
taskid <- as.numeric(taskid)
set.seed(69 * taskid)

if(taskid == 1) {
  write_rds(sims_parameters, 
            file = "./sims_parameters.RDS")
}

sim_results <- simulation_nrun_fnc(n_iter = 100,
                                   P = sims_parameters$P[taskid],
                                   beta_X = sims_parameters$beta_X[[taskid]],
                                   rho_X = sims_parameters$rho_X[taskid],
                                   Prev_Y = sims_parameters$Prev_Y[taskid],
                                   prop_missing = sims_parameters$prop_missing[taskid],
                                   missing_mech = sims_parameters$missing_mech[taskid],
                                   RR_multiplier = sims_parameters$multipliers[taskid])

write_rds(sim_results, 
          file = paste("./simulation_results_", taskid, ".RDS", sep = ""))

warnings()