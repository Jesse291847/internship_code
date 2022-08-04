library(rstan)


get_prior_con <- function(input, prior, temporal, df, pruning = c(0.005, 0.995)) {
  residuals <- input

  for (i in seq_len(nrow(input))) {
    for (j in seq_len(ncol(input))) {
      residuals[i, j] <- residuals[i, j] - sum(input[i, j] * temporal[1:ncol(input), j])
    }
  }
  S <- (df - ncol(input) - 1) * prior
  S <- as.matrix(nearPD(S)$mat)

  data_stan <- list(
    N = nrow(input), #
    dim = ncol(input),
    x = residuals,
    df = df,
    S = S
  )

  # options(mc.cores = 8 )
  fit.stan <- stan(
    file = "data_driven_pt029.stan",
    data = data_stan,
    chains = 16,
    iter = 2000 # , refresh = 0
  )
  res_w1 <- summary(fit.stan, pars = "pcor", probs = pruning)
  est_con <- matrix(
    ifelse(res_w1$summary[, 4] <= 0 & res_w1$summary[, 5] >= 0,
      0, res_w1$summary[, 1]
    ),
    ncol(input), ncol(input)
  )
  diag(est_con) <- 0
  return(est_con)
}
