library(psychonetrics)
source("compare_networks.R")
library(graphicalVAR)
library(tidyverse)

estimate_model <- function(N, true_model, make_ordinal = FALSE, prior = "no",
                           prior_matrix = "perfect", mass_prior = 30, missing = "no", prop_missing = 0.2) {
  #options(warn = -1)
  if (class(true_model) %in% "psychonetrics") {
    # get information required for simulating data and comparing models
    kappa <- try(getmatrix(true_model, "kappa_zeta_within"), silent = TRUE)
    if (!is.matrix(kappa)) {
      kappa <- getmatrix(true_model, "kappa_zeta")
    }
    beta <- getmatrix(true_model, "beta")
    true_temp <- getmatrix(true_model, "PDC")
    true_con <- try(getmatrix(true_model, "omega_zeta_within"), silent = TRUE)
    if (!is.matrix(true_con)) {
      true_con <- getmatrix(true_model, "omega_zeta")
    }
  } else if(class(true_model) %in% "gVARmodel"){
    true_con <- true_model$PCC
    true_temp <- true_model$PDC
    kappa <- true_model$kappa
    beta <- true_model$beta
  }
    
    else{
    stop("Only supports  models of class 'psychonetrics' and 'gVARmodel'")
  }

  if (prior_matrix == "perfect") {
    prior_matrix <- true_con
  } else if (prior_matrix == "test") {
    prior_matrix <- matrix(-1, nrow(true_con), ncol(true_con))
  } else if (prior_matrix == "random") {
    prior_matrix <- matrix(runif(ncol(beta)^2) * 2 - 1, ncol = ncol(beta))
    prior_matrix <- t(prior_matrix) %*% prior_matrix
  } else if (prior_matrix == "realistic") {
    test <- true_con
    test[test > 0] <- sample(-1:1, length(test[test > 0]), prob = c(0, 0.35, 0.65), replace = TRUE)
    test[test < 0] <- sample(-1:1, length(test[test < 0]), prob = c(0.65, 0.3, 0.05), replace = TRUE)
    test[test == 0] <- sample(-1:1, length(test[test == 0]), prob = c(0.1, 0.8, 0.1), replace = TRUE)
    # make symmetrical again
    diag(test) <- 0
    test[lower.tri(test)] <- t(test)[lower.tri(test)]
    prior_matrix <- test * 0.5
  } else {
    (stop("You can only use a perfect, reaslistic or random prior matrix"))
  }

  # sim data
  sim_data <- graphicalVARsim(N, beta = beta, kappa = kappa)
  # if( apply(sim_data,2, function(sim_data) is.na(sim_data)) ==  TRUE)
  # {print("There are NA's")} else {
  #   print("There are no NA's")
  # }

  # keep continuous or make ordinal
  if (make_ordinal == FALSE) {
    sim_data <- sim_data
  } else if (make_ordinal > 1) {
    if (make_ordinal < 4) {
      print("The grahicalVAR model should not be used for data with too few levels in ordinal data")
    }
    vectordata <- as.numeric(unlist(sim_data))

    Q <- c(as.numeric(quantile(vectordata, probs = seq(0, 1, 1 / make_ordinal))))

    sim_data_cat <- data.frame(matrix(ncol = ncol(sim_data), nrow = nrow(sim_data)))
    for (i in seq_len(nrow(sim_data))) {
      y <- cut(as.numeric(sim_data[i, ]),
        breaks = Q, labels = 0:(make_ordinal - 1),
        include.lowest = TRUE
      )
      sim_data_cat[i, ] <- as.integer(y) - 1
    }
    sim_data <- sim_data_cat
  } else {
    stop("Make ordinal needs to be FALSE or a number greater than 1")
  }


  if(missing %in% "no"){
  sim_data <- sim_data
  }
    else if (missing == "items") {
      sim_data[sample(c(TRUE, FALSE), prod(dim(sim_data)), replace = TRUE, prob = c(prop_missing, 1 - prop_missing))] <- NA
    } else if (missing %in% "timepoints") {
      sim_data[sample(seq_len(nrow(sim_data)), prop_missing * nrow(sim_data)), ] <- NA
    } else if (missing %in% "realistic") {
      sim_data[sample(seq_len(nrow(sim_data)), prop_missing * nrow(sim_data), prob = seq(1, nrow(sim_data), 1)), ] <- NA
    }
   else {
    stop("Missing can be changed to 'items', 'timepoints' or 'realistic'")
  }


  est_model <- gvar(sim_data, estimator = "FIML", standardize = "z") %>% runmodel()
  raw_contemptoreneous <- getmatrix(est_model, "omega_zeta")
  raw_temporal <- getmatrix(est_model, "PDC")

  est_model <- est_model %>% prune()
  est_temp <- getmatrix(est_model, "PDC")

  if (prior == "no") {
    est_con <- getmatrix(est_model, "omega_zeta")
  } else if (prior == "yes") {
    est_con <- get_prior_con(input = sim_data, prior = prior_matrix, temporal = raw_temporal, df = mass_prior)
  } else if (prior == "both") {
    est_con <- getmatrix(est_model, "omega_zeta")
    est_con_prior <- get_prior_con(input = sim_data, prior = prior_matrix, temporal = raw_temporal, df = mass_prior)
  }

  if (prior != "both") {
    results <- data.frame(
      ME_temp = calculate_estimation_error(est_temp, true_temp),
      Spec_temp = calculate_specificity(est_temp, true_temp),
      Sen_temp = calculate_sensitivity(est_temp, true_temp),
      Cor_temp = calculate_correlation(est_temp, true_temp),
      Det_temp = calculate_strength_detected(est_temp, true_temp),
      Mis_temp = calculate_strength_missed(est_temp, true_temp),
      SpurNeg_temp = calculate_spur_neg_edges(est_temp, true_temp),
      ME_con = calculate_estimation_error(est_con, true_con),
      Sen_con = calculate_sensitivity(est_con, true_con),
      Spec_con = calculate_specificity(est_con, true_con),
      Cor_con = calculate_correlation(est_con, true_con),
      Det_con = calculate_strength_detected(est_con, true_con),
      Mis_con = calculate_strength_missed(est_con, true_con),
      SpurNeg_con = calculate_spur_neg_edges(est_con, true_con)
    )

    output <- list(results, est_temp, est_con, raw_contemptoreneous, raw_contemptoreneous)
    names(output) <- c(
      "quality", "temporal",
      "contemporaneous", "raw_con", "raw_temp"
    )
  } else {
    results <- data.frame(
      ME_temp = calculate_estimation_error(est_temp, true_temp),
      Spec_temp = calculate_specificity(est_temp, true_temp),
      Sen_temp = calculate_sensitivity(est_temp, true_temp),
      Cor_temp = calculate_correlation(est_temp, true_temp),
      Det_temp = calculate_strength_detected(est_temp, true_temp),
      Mis_temp = calculate_strength_missed(est_temp, true_temp),
      SpurNeg_temp = calculate_spur_neg_edges(est_temp, true_temp),
      ME_con = calculate_estimation_error(est_con, true_con),
      Sen_con = calculate_sensitivity(est_con, true_con),
      Spec_con = calculate_specificity(est_con, true_con),
      Cor_con = calculate_correlation(est_con, true_con),
      Det_con = calculate_strength_detected(est_con, true_con),
      Mis_con = calculate_strength_missed(est_con, true_con),
      SpurNeg_con = calculate_spur_neg_edges(est_con, true_con),
      ME_con_prior = calculate_estimation_error(est_con_prior, true_con),
      Sen_con_prior = calculate_sensitivity(est_con_prior, true_con),
      Spec_con_prior = calculate_specificity(est_con_prior, true_con),
      Cor_con_prior = calculate_correlation(est_con_prior, true_con)
    )


    output <- list(results, est_temp, est_con, est_con_prior, raw_contemptoreneous, raw_temporal)
    names(output) <- c(
      "quality", "temporal",
      "contemporaneous", "prior", "raw_con", "raw_temp"
    )
  }

  return(output)
}


true_model <- randomGVARmodel(9)

