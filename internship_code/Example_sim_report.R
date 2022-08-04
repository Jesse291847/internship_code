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

# script to estimate network with a prior
source("construct_prior.R")

# the data simulating model ('true model')
model_gvar <- randomGVARmodel(9)

# contains one large function that simulates data and estimates networks
source("simulator_networks.R")

results <- parSim(
  
  # Conditions:
  time_points = c(70, 150),
  
  
  # Setup:
  name = paste0("sumScoreSims2_v1_", args[1]),
  write = FALSE,
  reps = 30,  
  debug = FALSE,
  progressbar = TRUE,

  
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
    estimate_model(N = time_points, model_gvar)$quality
  }
)


results <- results %>% select(-c(error, errorMessage))
results <- melt(results, id = c("rep", "id", "time_points"))
results$type <- ifelse(grepl("temp", results$variable), "Temporal", ifelse(grepl('prior', results$variable), "Prior",
                                                                           "Contemporaneous"))
results$outcome <- ifelse(grepl("Sen", results$variable), "Sensitivity", 
                          ifelse(grepl("Spec", results$variable), "Specificity", 
                                 ifelse(grepl("ME", results$variable), "Estimation Error", 
                                        ifelse(grepl('Det', results$variable), 'Detected', 
                                               ifelse(grepl('Mis', results$variable), 'Missed', 
                                                      ifelse(grepl('Neg', results$variable), 'Spurious Negative Edges', 'Correlation'))))))

#change the filter if you want to see different outcome measures
results %>% filter(outcome %in% c("Estimation Error")) %>% 
  ggplot(results, mapping =aes(x= type, y = value, fill = as.factor(time_points))) +
  facet_wrap(~outcome) + 
  geom_boxplot() +
  ylim(0,1) +
  labs(title = "Priors vs no priors aggregated over all time points", y = "", fill = "Edge", x = "") +
  theme_minimal()







model <- randomGVARmodel(5)
colnames(model$PCC) <- paste0("V", 1:5)
colnames(model$PDC) <- paste0("V", 1:5)


pdf("plot_network_report.pdf")
par(mfrow = c(1,2))
qgraph(model$PCC, theme = "colorblind")
qgraph(model$PDC, theme = "colorblind")
dev.off()
