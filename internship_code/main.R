# install if not installed
list.of.packages <- c(
  "tidyverse", "graphicalVAR", "rstan", "reshape2", "parSim",
  "psychonetrics", "lqmm", "rjson", "parallel", "Matrix"
)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if (length(new.packages)) install.packages(new.packages)

# load all necessary packages and scripts/objects
library(parSim)
library(tidyverse)
library(psychonetrics)
library(graphicalVAR)
library(lqmm)
library(rstan)
library(reshape2)
library(rjson)
library(parallel)
library(Matrix)


# script with functions to calculate outcome measures (sensitivity etc.)
source("compare_networks.R")

obj <- new.env()
source("compare_networks.R", local = obj)

# script to estimate network with a prior

source("construct_prior.R")
obj2 <- new.env()
source("construct_prior.R", local = obj2) # this one actually needs the stan file in the zip

# the data simulating model ('true model')
model_gvar <- readRDS("true_network_MD_phq.RDS")

# contains one large function that simulates data and estimates networks
source("simulator_networks.R")
obj3 <- new.env()
source("simulator_networks.R", local = obj3)

# everything in the script can stay the same if you change the file here:
config <- fromJSON(file = "config/config_20_80_realistic_100.json")
# config <- fromJSON(file = "/tmp/R_simulations_config.json")

timepoints <- config$number_of_timepoints
df <- config$df
quality_prior <- config$prior_cov_matrix
repetitions <- config$repetitions

results <- parSim(

  # Conditions:
  timepoints = timepoints,
  df = df,
  quality_prior = quality_prior,


  # Setup:
  name = paste0("sumScoreSims2_v1_", args[1]),
  write = FALSE,
  reps = repetitions, # maybe first change to 2?
  debug = FALSE,
  progressbar = FALSE,
  nCores = 8,
  export = c("model_gvar", ls(obj2), ls(obj), ls(obj3)),

  # The simulation code:
  expression = {
    library(parSim)
    library(tidyverse)
    library(psychonetrics)
    library(graphicalVAR)
    library(lqmm)
    library(rstan)
    library(reshape2)
    library(rjson)
    library(parallel)
    library(Matrix)
    # if number of cores is known maybe change to that an remove parallel from packages
    options(mc.cores = detectCores())
    estimate_model(N = timepoints, model_gvar, prior = "both", 
    prior_matrix = quality_prior, mass_prior = df)$quality
  }
)

saveRDS(results, file = paste0(
  "/tmp/sim_result_", timepoints, "_", df, "_",
  quality_prior, "_", config$repetitions, ".RDS"
))

