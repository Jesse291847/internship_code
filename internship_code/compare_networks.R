calculate_correlation <- function(est_network, true_network) {
  if (isSymmetric(true_network)) {
    # Edges of both networks:
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # Calculate correlation:
  correlation <- cor(true_edges, est_edges)
  return(correlation)
}


calculate_sensitivity <- function(est_network, true_network) {

  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }
  # True positives cutof:
  true_pos <- sum(true_edges != 0 & est_edges != 0)

  # False Negatives cutoff:
  false_neg <- sum(true_edges != 0 & est_edges == 0)

  # Calculate Sensitivity:
  return((true_pos) / (true_pos + false_neg))
}




calculate_specificity <- function(est_network, true_network) {

  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # False positives cutoff:
  false_pos <- sum(true_edges == 0 & est_edges != 0)

  # True Negatives cutoff:
  true_neg <- sum(true_edges == 0 & est_edges == 0)

  # Calculate Specificity:
  return((true_neg) / (true_neg + false_pos))
}



calculate_spur_neg_edges <- function(est_network, true_network) {

  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # number of spurious negative edges
  spurious_edges <- sum(est_edges < 0 & true_edges >= 0)

  # number of possible edges
  all_possible_edges <- sum(complete.cases(true_edges))

  # return number of spurious negative edges relative to the number of possible edges
  return(spurious_edges / all_possible_edges)
}

calculate_density <- function(est_network) {

  # Edges:

  est_edges <- est_network[upper.tri(est_network)]

  # Number of edges:
  amount_edges <- sum(est_edges != 0)

  # Number of possible edges
  all_possible_edges <- sum(complete.cases(est_edges))

  # Return number of connections relative to the number of possible connection
  return(amount_edges / all_possible_edges)
}

calculate_estimation_error <- function(est_network, true_network) {
  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # Difference
  diff <- true_edges - est_edges

  # Estimation Error
  me <- sum(abs(diff))

  # Return ME
  return(me)
}



calculate_strength_missed <- function(est_network, true_network) {

  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # Mean strength of edges missed

  average_missed_edge <- mean(abs(true_edges[(true_edges != 0 & est_edges == 0)]))

  return(average_missed_edge)
}

calculate_strength_detected <- function(est_network, true_network) {
  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # Mean strength of edges missed

  average_detected_edge <- mean(abs(true_edges[(true_edges != 0 & est_edges != 0)]))

  
  return(average_detected_edge)
}

calculate_weakest_detected <- function(est_network, true_network) {
  # Edges of both networks:
  if (isSymmetric(true_network)) {
    true_edges <- true_network[upper.tri(true_network)]
    est_edges <- est_network[upper.tri(est_network)]
  } else {
    true_edges <- c(true_network)
    est_edges <- c(est_network)
  }

  # smallest edge detected

  weakest_detected_edge <- min(abs(true_edges[(true_edges != 0 & est_edges != 0)]))

  return(weakest_detected_edge)
}
