K_fun <- function(eta) {
  K <- lgamma((eta + 1) / eta * 0.5) - log(pi / eta) / 2 - lgamma(1 / (2 * eta))

  return(K)
}

Bound_f <- function(y, a, b) {
  x <- rep(NA, length(y))
  for (i in 1:length(y)) {
    x[i] <- a + (b - a) * (1 / (1 + exp(-2 * y[i])))
  }
  return(x)
}

Skt_f <- function(v_and_parameters) {
  v_t <- v_and_parameters[, 1]
  sigma <- v_and_parameters[, 2]
  nu1 <- v_and_parameters[, 3]
  alpha <- v_and_parameters[, 5]

  eta <- 1 / v_and_parameters[, 3]
  f_t <- exp(K_fun(eta) - 0.5 * log(sigma^2) - (1 + eta) / (2 * eta) * log(1 + eta * v_t^2 / ((1 - sign(v_t) * alpha)^2 * sigma^2)))

  return(f_t)
}

Skt_F <- function(v_and_parameters) {
  ## Equation 5 of AGQ (2005)

  v_t <- v_and_parameters[, 1]
  sigma <- v_and_parameters[, 2]
  nu1 <- v_and_parameters[, 3]
  alpha <- v_and_parameters[, 5]

  a <- (1 - alpha)
  b <- (1 + alpha)

  if (v_t < 0) {
    CDF <- 2 * b / (a + b) * pt(v_t / (sigma * b), nu1)
  } else {
    CDF <- (b - a) / (a + b) + 2 * a / (a + b) * pt(v_t / (sigma * a), nu1)
  }

  return(CDF)
}

rSkt <- function(N, parameters) {
  sigma <- parameters[1]
  nu1 <- parameters[2]
  alpha <- parameters[4]

  test <- rt(N * 2 + 0.1 * N, nu1)
  test_p <- (1 - alpha) * sigma * test[test > 0]
  test_n <- (1 + alpha) * sigma * test[test < 0]

  lenght_min <- min(c(length(test_p), length(test_n)))
  pp <- (1 - alpha) / 2
  pn <- (1 + alpha) / 2

  X <- rep(NA, lenght_min)
  X[1:(lenght_min * pp)] <- test_p[1:(lenght_min * pp)]
  X[(lenght_min * pp + 1):lenght_min] <- test_n[1:(lenght_min * pn)]
  X <- X[!is.na(X)]
  X <- sample(X)

  return(X[1:N])
}

Skt_score <- function(v_and_parameters, derivative) {
  v_t <- v_and_parameters[, 1]
  sigma <- v_and_parameters[, 2]
  nu1 <- v_and_parameters[, 3]
  alpha <- v_and_parameters[, 5]

  eta <- 1 / v_and_parameters[, 3]
  zeta <- v_t / sigma
  w <- (1 + eta) / ((1 - sign(v_t) * alpha)^2 + eta * zeta^2)

  if (derivative == "location") {
    score <- 1 / sigma * w * zeta
  }

  if (derivative == "scale") {
    score <- 1 / (2 * sigma^2) * (w * zeta^2 - 1)
  }

  if (derivative == "shape") {
    score <- -sign(v_t) / (1 - sign(v_t) * alpha) * w * zeta^2
  }

  return(score)
}

############# Large Factor Model with Location, Scale and Shape Common Factors
estimation <- function(model, nb_it, full = TRUE, check = FALSE, print_out = TRUE) {
  modeli <- loglik_list(rep(0, 5e3), model)$model
  init_par <- rep(0, modeli$n_pars)
  start <- 1
  min_inv <- NULL
  if (length(model$init_par_0) > 0) {
    init_par[1:length(model$init_par_0)] <- model$init_par_0 - 1e-4
  }

  MLE <- optim(init_par + 1e-4, wrapped_loglik_DFM,
    gr = wrapped_loglik_grad_DFM,
    model = model, method = "BFGS",
    control = list(trace = T, maxit = 1e6, reltol = 1e-20, maxit = 1e3)
  )

  if (print_out) {
    print(paste0("Start with MLE = ", MLE$value))
  }

  dsb <- 0.05

  min_MLE <- MLE

  restart <- TRUE

  while (restart == TRUE) {
    if (nb_it > 0) {
      count <- 0
      continue <- TRUE
      while (continue) {
        count <- count + 1
        MLE <- ucminf(rnorm(length(min_MLE$par), min_MLE$par, dsb),
          wrapped_loglik_DFM,
          gr = wrapped_loglik_grad_DFM,
          model = model
        )

        if (round(min_MLE$value - MLE$value, 2) > 0) {
          min_MLE <- MLE
          count <- 0
          if (print_out) {
            cat(paste0("new min ", round(min_MLE$value, 3)))
          }
        }
        if (count %% 5 == 0) {
          cat(paste0("it ", count, "; "))
        }

        if (count == nb_it) {
          continue <- FALSE
        }
      }
      if (print_out) {
        print(paste0("Step 1 finished with MLE = ", min_MLE$value))
      }
    }

    restart <- FALSE

    if (model$options$stoch_shape & model$options$freq[1] == "M") {
      filter_result <- loglik_list(min_MLE$par, model)
      if (sd(filter_result$a_shape_f[ncol(y) + 1, ]) > 1) {
        restart <- TRUE
        nb_it <- nb_it + 50
        min_MLE$par <- rep(0, length(min_MLE$par))
        min_MLE$par[length(min_MLE$par)] <- 1e2
        dsb <- dsb + 0.05
      }
    }
  }

  filter_result <- loglik_list(min_MLE$par, model)

  return(list(MLE = min_MLE, model = model, filter_result = filter_result))
}

seq_esti <- function(y, init_par, nb_it0, spec, freq = "Q", check = FALSE, mix = F) {
  corr <- rep(0, ncol(y) - 1)
  corr_sq <- rep(0, ncol(y) - 1)
  if (freq == "M") {
    start_s <- 1
  } else {
    start_s <- 2
  }
  for (i in 1:ncol(y)) {
    if (i >= start_s | freq == "M") {
      y[, i] <- y[, i] - mean(y[1:573, i], na.rm = T)
      y[, i] <- y[, i] / sd(y[1:573, i], na.rm = T)
    }

    yi <- na.omit(y[, c(1, i)])

    corr[i - 1] <- cor(yi[, 1], yi[, 2])
    corr_sq[i - 1] <- cor(yi[, 1]^2, yi[, 2]^2)
  }


  if (spec %in% c("L", "LS", "LSS")) {
    stoch_loc <- T
  } else {
    stoch_loc <- F
  }

  if (spec %in% c("LS", "LSS")) {
    stoch_scale <- T
  } else {
    stoch_scale <- F
  }

  if (spec == "LSS") {
    stoch_shape <- T
  } else {
    stoch_shape <- F
  }

  if (freq == "Q" & mix == F) {
    div <- 3
  } else {
    div <- 1
  }

  if (is.null(init_par)) {
    options <- list()
    options$freq <- rep(freq, ncol(y))

    if (mix) {
      options$freq <- c("Q", rep("M", ncol(y) - 1))
    }

    options$stoch_loc <- stoch_loc
    options$stoch_vol <- stoch_scale
    options$stoch_shape <- stoch_shape

    model <- model_components_DFM(y, options)

    model$guess_loc <- rep(NA, ncol(model$y))
    for (i in 1:length(model$guess_loc)) {
      model$guess_loc[i] <- mean(model$y[!is.na(model$y[, i]), i][1:12]) / div
    }

    model$guess_scale <- rep(NA, ncol(model$y))
    for (i in 1:length(model$guess_scale)) {
      model$guess_scale[i] <- sd(model$y[!is.na(model$y[, i]), i][1:12])
    }
    #

    run_0 <- loglik_list(rep(0, 2000), model)
    init_par <- rep(0, run_0$model$n_pars)
    if (spec == "LSS") {
      init_par[run_0$model$index_shape] <- corr
    }
    if (spec %in% c("LSS", "LS")) {
      init_par[run_0$model$index_scale] <- corr_sq
    }

    init_par[run_0$model$index_loc] <- corr

    model$init_par_0 <- init_par

    MLE_model <- estimation(model, nb_it = nb_it0, full = T, check = check)
  } else {
    options <- list()
    options$freq <- rep(freq, ncol(y))

    if (mix) {
      options$freq <- c("Q", rep("M", ncol(y) - 1))
    }

    options$stoch_loc <- stoch_loc
    options$stoch_vol <- stoch_scale
    options$stoch_shape <- stoch_shape

    model <- model_components_DFM(y, options)

    model$guess_loc <- rep(NA, ncol(model$y))
    for (i in 1:length(model$guess_loc)) {
      model$guess_loc[i] <- mean(model$y[!is.na(model$y[, i]), i][1:12]) / div
    }

    model$guess_scale <- rep(NA, ncol(model$y))
    for (i in 1:length(model$guess_scale)) {
      model$guess_scale[i] <- sd(model$y[!is.na(model$y[, i]), i][1:12])
    }
    #
    model$init_par_0 <- init_par

    MLE_model <- estimation(model, nb_it = nb_it0, full = T, check = F)
  }

  return(MLE_model)
}

model_components_DFM <- function(y, options) {
  N <- dim(t(y))[2]
  m <- dim(t(y))[1]

  #
  NAs <- matrix(0L, nrow(y), ncol(y))
  NAs[which(!is.na(y))] <- 1L

  # Location matrics
  cf_id <- m * 5 + 1
  T_loc <- matrix(0, m * 5 + 5, m * 5 + 5)
  Z_loc <- matrix(0, m, nrow(T_loc))

  for (i in 1:m) {
    j <- (i - 1) * 5 + 1

    # Transition matrix
    T_loc[j, j] <- 1 # trend
    T_loc[j + 1, j] <- 1 # trend L1
    T_loc[j + 2, j + 1] <- 1 # trend L2
    T_loc[j + 3, j + 2] <- 1 # trend L3
    T_loc[j + 4, j + 3] <- 1 # trend L4

    # Observation matrix
    if (options$freq[i] == "M") {
      Z_loc[i, j] <- 1
      Z_loc[i, cf_id] <- 1
    } else {
      Z_loc[i, j:(j + 4)] <- c(1 / 3, 2 / 3, 1, 2 / 3, 1 / 3)
      Z_loc[i, cf_id:(cf_id + 4)] <- c(1 / 3, 2 / 3, 1, 2 / 3, 1 / 3)
    }
  }

  T_loc[cf_id, cf_id] <- 0 # factor
  T_loc[cf_id + 1, cf_id] <- 1 # L1
  T_loc[cf_id + 2, cf_id + 1] <- 1 # L2
  T_loc[cf_id + 3, cf_id + 2] <- 1 # L3
  T_loc[cf_id + 4, cf_id + 3] <- 1 # L4

  # Scale matrices
  T_scale <- matrix(0, m + 2, m + 2)
  Z_scale <- matrix(0, m, nrow(T_scale))

  for (i in 1:m) {
    # Transition matrix
    T_scale[i, i] <- 1 # trend

    # Observation matrix
    Z_scale[i, i] <- 1
    Z_scale[i, m + 1] <- 1
  }
  T_scale[m + 2, m + 1] <- 1 # factor

  #### Shape components
  T_shape <- matrix(0, m + 2, m + 2)
  Z_shape <- matrix(0, m, nrow(T_shape))

  for (i in 1:m) {
    # Transition matrix
    T_shape[i, i] <- 1 # trend

    # Observation matrix
    Z_shape[i, i] <- 1
    Z_shape[i, m + 1] <- 1
  }
  T_shape[m + 2, m + 1] <- 1 # factor

  return(list(
    y = matrix(y, nrow = N, ncol = m), Z_loc = Z_loc, T_loc = T_loc, N = N,
    m = m, Z_scale = Z_scale, T_scale = T_scale,
    Z_shape = Z_shape, T_shape = T_shape, options = options,
    NAs = NAs
  ))
}

DGP_DFM <- function(filter_results, Tt = NULL) {
  if (is.null(Tt)) {
    Tt <- dim(filter_results$model$y)[1]
  }

  a <- matrix(NA, nrow(filter_results$a), Tt)
  a[, 1] <- filter_results$a[, ncol(filter_results$a)]
  a_scale <- matrix(NA, nrow(filter_results$a_scale), Tt)
  a_scale[, 1] <- filter_results$a_scale[, ncol(filter_results$a_scale)]
  a_shape <- matrix(NA, nrow(filter_results$a_shape), Tt)
  a_shape[, 1] <- filter_results$a_shape[, ncol(filter_results$a_shape)]
  K_loc <- filter_results$model$K_loc
  K_scale <- filter_results$model$K_scale
  K_shape <- filter_results$model$K_shape
  T_loc <- filter_results$model$T_loc
  T_scale <- filter_results$model$T_scale
  T_shape <- filter_results$model$T_shape
  Z_loc <- filter_results$model$Z_loc
  Z_scale <- filter_results$model$Z_scale
  Z_shape <- filter_results$model$Z_shape

  vp <- cbind(NA, filter_results$parameters)
  signal <- y

  for (t in 1:(Tt - 1)) {
    signal[t, ] <- t(Z_loc %*% a[, t])
    vp[, 1] <- Z_loc %*% a[, t]
    vp[, 2] <- exp(Z_scale %*% a_scale[, t])
    vp[, 5] <- tanh(Z_shape %*% a_shape[, t])

    sim_obs <- rep(NA, ncol(y))
    for (i in 1:ncol(y)) {
      sim_obs[i] <- rSkt(1e2, vp[i, -1])[1]
    }

    vp[, 1] <- sim_obs
    score_loc <- Skt_score(vp, "location")
    score_scale <- Skt_score(vp, "scale") * 2 * exp(2 * Z_scale %*% a_scale[, t])
    score_shape <- Skt_score(vp, "shape") * (1 - tanh(Z_shape %*% a_shape[, t]) * tanh(Z_shape %*% a_shape[, t]))

    a[, t + 1] <- T_loc %*% a[, t] + K_loc %*% t(Z_loc) %*% score_loc
    a_scale[, t + 1] <- T_scale %*% a_scale[, t] + K_scale %*% t(Z_scale) %*% score_scale
    a_shape[, t + 1] <- T_shape %*% a_shape[, t] + K_shape %*% t(Z_shape) %*% score_shape
  }

  return(list(
    a = a,
    a_scale = a_scale,
    a_shape = a_shape
  ))
}

gen_obs_DFM <- function(estimation_results) {
  model <- estimation_results$model
  a_f <- estimation_results$a_f
  a_scale_f <- estimation_results$a_scale_f
  a_shape_f <- estimation_results$a_shape_f

  which_NA <- is.na(model$y)

  parameters <- model$parameters

  ###
  y <- matrix(NA, nrow(model$y), ncol(model$y))
  v <- matrix(NA, nrow(model$y), ncol(model$y))

  for (t in 1:nrow(model$y)) {
    parameters[, 1] <- exp(model$Z_scale %*% a_scale_f[, t])
    parameters[, 4] <- tanh(model$Z_shape %*% a_shape_f[, t])

    for (i in 1:ncol(model$y)) {
      ff <- rSkt(50, parameters[i, ])
      ff <- ff[is.finite(ff)]
      v[t, i] <- ff[1]
    }

    y[t, ] <- model$Z_loc[, ] %*% a_f[, t] + v[t, ]
  }

  y[which_NA] <- NA
  return(list(y = y))
}

wrapped_loglik_DFM <- function(pars, model) {
  loglik <- loglik_rcpp(
    pars,
    model$y,
    model$T_loc,
    model$T_scale,
    model$T_shape,
    model$Z_loc,
    model$Z_scale,
    model$Z_shape,
    model$options$stoch_loc,
    model$options$stoch_vol,
    model$options$stoch_shape,
    model$NAs,
    model$guess_loc,
    model$guess_scale
  )


  return(loglik)
}

wrapped_loglik_DFM2 <- function(pars, model) {
  loglik <- loglik_rcpp(
    pars,
    model$y,
    model$T_loc,
    model$T_scale,
    model$T_shape,
    model$Z_loc,
    model$Z_scale,
    model$Z_shape,
    model$options$stoch_loc,
    model$options$stoch_vol,
    model$options$stoch_shape,
    model$NAs,
    model$guess_loc,
    model$guess_scale
  )


  return(-loglik)
}

wrapped_loglik_grad_DFM <- function(pars, model) {
  gradient <- loglik_deriv(
    pars,
    model$y,
    model$T_loc,
    model$T_scale,
    model$T_shape,
    model$Z_loc,
    model$Z_scale,
    model$Z_shape,
    model$options$stoch_loc,
    model$options$stoch_vol,
    model$options$stoch_shape,
    model$NAs,
    model$guess_loc,
    model$guess_scale
  )

  return(gradient)
}

wrapped_loglik_grad_DFM2 <- function(pars, model) {
  gradient <- loglik_deriv(
    pars,
    model$y,
    model$T_loc,
    model$T_scale,
    model$T_shape,
    model$Z_loc,
    model$Z_scale,
    model$Z_shape,
    model$options$stoch_loc,
    model$options$stoch_vol,
    model$options$stoch_shape,
    model$NAs,
    model$guess_loc,
    model$guess_scale
  )

  return(-gradient)
}

wrapped_loglik_Hessian_DFM <- function(pars, model) {
  Hessian <- loglik_hess(
    pars,
    model$y,
    model$T_loc,
    model$T_scale,
    model$T_shape,
    model$Z_loc,
    model$Z_scale,
    model$Z_shape,
    model$options$stoch_loc,
    model$options$stoch_vol,
    model$options$stoch_shape,
    model$NAs,
    model$guess_loc,
    model$guess_scale
  )
  return(Hessian)
}

############# Benchmark + Factor augmented models
obj_AR2 <- function(pars, model, estimate = TRUE) {
  y <- model$y
  spec <- model$spec

  v <- rep(NA, length(y))
  y_pred <- rep(NA, length(y) + 1)
  loc <- rep(NA, length(y) + 1)
  AR2 <- tanh(pars[1])
  AR1 <- Bound_f(pars[2], AR2 - 1, 1 - AR2)
  scale <- pars[3] + log(sd(y, na.rm = TRUE))

  nu <- 1000
  if (spec$dist == "Student") {
    nu <- Bound_f(pars[4] - 3, 1.1, 1000)
  }

  loc <- pars[5] + mean(y, na.rm = T)

  loglik <- rep(0, length(y))

  for (t in 3:length(y)) {
    y_pred[t] <- loc + AR1 * y[t - 1] + AR2 * y[t - 2]
    v[t] <- y[t] - y_pred[t]
    vp <- t(c(v[t], exp(scale), nu, nu, 0))
    loglik[t] <- log(Skt_f(vp))
  }

  y_pred[t + 1] <- loc + AR1 * y[t] + AR2 * y[t - 1]

  v <- v[-c(1, 2)]

  if (estimate) {
    output <- -sum(loglik)
  } else {
    output <- list(
      y_pred = y_pred,
      scale = exp(scale),
      nu = nu,
      N = length(y) + 1
    )
  }
  return(output)
}

wrapped_loglik_reg <- function(pars, model) {
  options <- model$options

  obj <- loglik_reg_rcpp(pars,
    y = model$y, stoch_loc = options$stoch_loc, stoch_vol = options$stoch_vol, stoch_shape = options$stoch_shape,
    X_loc = model$X_loc, X_scale = model$X_scale, X_shape = model$X_shape, model$NAs
  )

  if (!is.finite(obj)) {
    obj <- 1e9
  }
  return(obj)
}

wrapped_loglik_reg_grad <- function(pars, model) {
  options <- model$options

  gradient <- loglik_reg_deriv(pars,
    y = model$y, stoch_loc = options$stoch_loc, stoch_vol = options$stoch_vol, stoch_shape = options$stoch_shape,
    X_loc = model$X_loc, X_scale = model$X_scale, X_shape = model$X_shape, model$NAs
  )
  return(gradient)
}

wrapped_loglik_reg_hess <- function(pars, model) {
  options <- model$options

  hessian <- loglik_reg_hess(pars,
    y = model$y, stoch_loc = options$stoch_loc, stoch_vol = options$stoch_vol, stoch_shape = options$stoch_shape,
    X_loc = model$X_loc, X_scale = model$X_scale, X_shape = model$X_shape, model$NAs
  )
  return(hessian)
}

estimation_reg <- function(model, nb_it, full = TRUE, check = TRUE, print_out = TRUE) {
  modeli <- loglik_reg_list(rep(0, 1e3), model)$model
  init_par <- rep(0, modeli$n_pars)

  if (length(model$init_par_0) > 0) {
    init_par[1:length(model$init_par_0)] <- model$init_par_0
  }

  MLE <- optim(init_par, wrapped_loglik_reg,
    gr = wrapped_loglik_reg_grad,
    model = model, method = "BFGS", control = list(trace = T, maxit = 1000, reltol = 1e-20)
  )

  if (print_out) {
    print(paste0("Start with MLE = ", MLE$value))
  }

  if (MLE$value == 1e9) {
    MLE$par <- rep(0, length(MLE$par))
    nb_it <- 100
  }

  min_MLE <- MLE

  restart <- TRUE

  while (restart == TRUE) {
    if (nb_it >= 1) {
      count <- 0
      continue <- TRUE
      while (continue) {
        count <- count + 1
        MLE <- ucminf(rnorm(length(min_MLE$par), min_MLE$par, 0.5),
          wrapped_loglik_reg,
          gr = wrapped_loglik_reg_grad,
          model = model
        )

        if (round(min_MLE$value - MLE$value, 2) > 0) {
          min_MLE <- MLE
          count <- 0
          if (print_out) {
            cat(paste0("new min ", round(min_MLE$value, 3)))
          }
        }
        if (count %% 5 == 0) {
          cat(paste0("it ", count, "; "))
        }
        if (count >= nb_it) {
          continue <- FALSE
        }
      }
      if (print_out) {
        print(paste0("Step 0 finished with MLE = ", min_MLE$value))
      }
    }

    restart <- FALSE
  }

  filter_result <- loglik_reg_list(min_MLE$par, model)

  return(list(MLE = min_MLE, model = model, filter_result = filter_result))
}

estimation_AR2 <- function(model, nb_it, print_out = TRUE) {
  init_par <- rep(0, 6)

  if (length(model$init_par_0) > 0) {
    init_par[1:length(model$init_par_0)] <- model$init_par_0
  }

  MLE <- optim(init_par, obj_AR2,
    model = model, method = "BFGS", control = list(
      trace = T,
      maxit = 1000,
      reltol = 1e-20
    )
  )

  if (print_out) {
    print(paste0("Start with MLE = ", MLE$value))
  }

  min_MLE <- MLE

  count <- 0
  continue <- TRUE
  if (nb_it > 0) {
    while (continue) {
      count <- count + 1
      MLE <- ucminf(rnorm(length(min_MLE$par), min_MLE$par, 0.2),
        obj_AR2,
        model = model
      )
      if (round(min_MLE$value - MLE$value, 2) > 0) {
        min_MLE <- MLE
        count <- 0
        if (print_out) {
          cat(paste0("new min ", round(min_MLE$value, 3)))
        }
      }
      if (count %% 5 == 0) {
        cat(paste0("it ", count, "; "))
      }
      if (count == nb_it) {
        continue <- FALSE
      }
    }
  }
  if (print_out) {
    print(paste0("Finished with MLE = ", min_MLE$value))
  }

  pred <- obj_AR2(min_MLE$par, model, FALSE)

  return(list(MLE = min_MLE, model = model, pred = pred))
}

############# Evaluation of the density nowcasts
out_of_sample_evaluation <- function(filter_results, model, target_gdp_figure, steps = NULL) {
  y <- model$y
  N <- dim(y)[1]

  if (abs(target_gdp_figure) < 2) {
    limit <- Inf
  } else {
    limit <- Inf
  }

  if (model$options$freq[1] == "Q") {
    st <- 2
  } else {
    st <- 1
    id <- 1
  }

  if (is.null(steps)) {
    if (!any(is.na(y[dim(y)[1], st:dim(y)[2]]))) {
      steps <- 0
    }
    if (any(is.na(y[dim(y)[1], st:dim(y)[2]]))) {
      steps <- 1
    }
    if (any(is.na(y[dim(y)[1] - 1, st:dim(y)[2]]))) {
      steps <- 2
    }
    if (any(is.na(y[dim(y)[1] - 2, st:dim(y)[2]]))) {
      steps <- 3
    }
  }

  ############# Step 2
  if (steps == 1) {
    print_colour("step 1", "cyan")
    n_draws <- 100
    GDP_draws_prediction <- rep(NA, n_draws^(steps + 1))

    if (model$options$freq[1] == "M") {
      location_CF <- matrix(filter_results$a_f[ncol(y) * 5 + 1, ], dim(filter_results$a_f)[2], n_draws^steps)
      scale_CF <- matrix(filter_results$a_scale_f[ncol(y) + 1, ], dim(filter_results$a_scale_f)[2], n_draws^steps)
      shape_CF <- matrix(filter_results$a_shape_f[ncol(y) + 1, ], dim(filter_results$a_shape_f)[2], n_draws^steps)
    }


    ## drawing missing observations in t
    t_N_draws <- matrix(NA, n_draws, ncol(y))

    which_y_is_NAs_N <- which(is.na(y[N, 1:dim(y)[2]]))
    which_y_is_NAs_N <- which_y_is_NAs_N[-1] # dropping GDP

    for (i in which_y_is_NAs_N) {
      draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
        filter_results$scale_f[N, i],
        filter_results$parameters[i, 2],
        filter_results$parameters[i, 3],
        filter_results$shape_f[N, i]
      )))
      draw_first <- draw_first[abs(draw_first) < limit]

      t_N_draws[, i] <- filter_results$location_f[N, i] + draw_first[1:n_draws]
    }
    ##

    ## getting filtered parameters for each draw
    j <- 1
    y_with_draw_t <- y
    for (draw3 in 1:n_draws) {
      model_m0 <- filter_results$model

      ## allocating the draws
      y_with_draw_t[N, which_y_is_NAs_N] <- t_N_draws[draw3, which_y_is_NAs_N]
      NAs <- matrix(0L, nrow(y_with_draw_t), ncol(y_with_draw_t))
      NAs[which(!is.na(y_with_draw_t))] <- 1L
      model_m0$NAs <- NAs
      model_m0$y <- y_with_draw_t
      ##

      # filtering step
      filter_results_m0 <- filter_list(model_m0)

      if (model$options$freq[1] == "Q") {
        draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
          filter_results_m0$scale_f[N, 1],
          filter_results_m0$parameters[1, 2],
          filter_results_m0$parameters[1, 3],
          filter_results_m0$shape_f[N, 1]
        )))
        draw_first <- draw_first[abs(draw_first) < limit]

        GDP_draws_prediction[j:(j + n_draws - 1)] <- filter_results_m0$location_f[N, 1] + draw_first[1:n_draws]
        j <- j + n_draws
      } else {
        location_CF[, id] <- filter_results_m0$a_f[ncol(y) * 5 + 1, ]
        scale_CF[, id] <- filter_results_m0$a_scale_f[ncol(y) + 1, ]
        shape_CF[, id] <- filter_results_m0$a_shape_f[ncol(y) + 1, ]
        id <- id + 1
      }
    }
  }

  ############# Step 2
  if (steps == 2) {
    print_colour("step 2", "cyan")
    n_draws <- 25
    GDP_draws_prediction <- rep(NA, n_draws^(steps + 1))

    if (model$options$freq[1] == "M") {
      location_CF <- matrix(filter_results$a_f[ncol(y) * 5 + 1, ], dim(filter_results$a_f)[2], n_draws^steps)
      scale_CF <- matrix(filter_results$a_scale_f[ncol(y) + 1, ], dim(filter_results$a_scale_f)[2], n_draws^steps)
      shape_CF <- matrix(filter_results$a_shape_f[ncol(y) + 1, ], dim(filter_results$a_shape_f)[2], n_draws^steps)
    }

    ## drawing missing observations in t-1
    t_Nm1_draws <- matrix(NA, n_draws, ncol(y))
    which_y_is_NAs_Nm1 <- which(is.na(y[N - 1, 1:dim(y)[2]]))
    which_y_is_NAs_Nm1 <- which_y_is_NAs_Nm1[-1] # dropping GDP

    for (i in which_y_is_NAs_Nm1) {
      draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
        filter_results$scale_f[N - 1, i],
        filter_results$parameters[i, 2],
        filter_results$parameters[i, 3],
        filter_results$shape_f[N - 1, i]
      )))
      draw_first <- draw_first[abs(draw_first) < limit]

      t_Nm1_draws[, i] <- filter_results$location_f[N - 1, i] + draw_first[1:n_draws]
    }
    ##

    ## Predicting parameters in t for each draw
    j <- 1
    y_with_draw_t1 <- y

    for (draw2 in 1:n_draws) {
      model_m1 <- filter_results$model

      ## allocating the draws
      y_with_draw_t1[N - 1, which_y_is_NAs_Nm1] <- t_Nm1_draws[draw2, which_y_is_NAs_Nm1]
      NAs <- matrix(0L, nrow(y_with_draw_t1), ncol(y_with_draw_t1))
      NAs[which(!is.na(y_with_draw_t1))] <- 1L
      model_m1$NAs <- NAs
      model_m1$y <- y_with_draw_t1
      ##

      ## prediction step
      filter_results_m1 <- filter_list(model_m1)

      ## drawing missing observations in t
      t_N_draws <- matrix(NA, n_draws, ncol(y))

      which_y_is_NAs_N <- which(is.na(y[N, 1:dim(y)[2]]))
      which_y_is_NAs_N <- which_y_is_NAs_N[-1] # dropping GDP

      for (i in which_y_is_NAs_N) {
        draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
          filter_results_m1$scale_f[N, i],
          filter_results_m1$parameters[i, 2],
          filter_results_m1$parameters[i, 3],
          filter_results_m1$shape_f[N, i]
        )))
        draw_first <- draw_first[abs(draw_first) < limit]

        t_N_draws[, i] <- filter_results_m1$location_f[N, i] + draw_first[1:n_draws]
      }
      ##

      ## getting filtered parameters for each draw
      y_with_draw_t <- y_with_draw_t1
      for (draw3 in 1:n_draws) {
        model_m0 <- model_m1

        ## allocating the draws
        y_with_draw_t[N, which_y_is_NAs_N] <- t_N_draws[draw3, which_y_is_NAs_N]
        NAs <- matrix(0L, nrow(y_with_draw_t), ncol(y_with_draw_t))
        NAs[which(!is.na(y_with_draw_t))] <- 1L
        model_m0$NAs <- NAs
        model_m0$y <- y_with_draw_t
        ##

        # filtering step
        filter_results_m0 <- filter_list(model_m0)

        if (model$options$freq[1] == "Q") {
          draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
            filter_results_m0$scale_f[N, 1],
            filter_results_m0$parameters[1, 2],
            filter_results_m0$parameters[1, 3],
            filter_results_m0$shape_f[N, 1]
          )))
          draw_first <- draw_first[abs(draw_first) < limit]

          GDP_draws_prediction[j:(j + n_draws - 1)] <- filter_results_m0$location_f[N, 1] + draw_first[1:n_draws]
          j <- j + n_draws
        } else {
          location_CF[, id] <- filter_results_m0$a_f[ncol(y) * 5 + 1, ]
          scale_CF[, id] <- filter_results_m0$a_scale_f[ncol(y) + 1, ]
          shape_CF[, id] <- filter_results_m0$a_shape_f[ncol(y) + 1, ]
          id <- id + 1
        }
      }
    }
  }

  ############# Step 3
  if (steps == 3) {
    print_colour("step 3", "cyan")
    n_draws <- 10
    GDP_draws_prediction <- rep(NA, n_draws^(steps + 1))

    if (model$options$freq[1] == "M") {
      n_draws <- 10
      location_CF <- matrix(filter_results$a_f[ncol(y) * 5 + 1, ], dim(filter_results$a_f)[2], n_draws^3)
      scale_CF <- matrix(filter_results$a_scale_f[ncol(y) + 1, ], dim(filter_results$a_scale_f)[2], n_draws^3)
      shape_CF <- matrix(filter_results$a_shape_f[ncol(y) + 1, ], dim(filter_results$a_shape_f)[2], n_draws^3)
    }

    ## drawing missing observations in t-2
    which_y_is_NAs_Nm2 <- which(is.na(y[N - 2, 1:dim(y)[2]]))
    which_y_is_NAs_Nm2 <- which_y_is_NAs_Nm2[-1] # dropping GDP

    t_Nm2_draws <- matrix(NA, nrow = n_draws, ncol = ncol(y))
    for (i in which_y_is_NAs_Nm2) {
      draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
        filter_results$scale_f[N - 2, i],
        filter_results$parameters[i, 2],
        filter_results$parameters[i, 3],
        filter_results$shape_f[N - 2, i]
      )))
      draw_first <- draw_first[abs(draw_first) < limit]

      t_Nm2_draws[, i] <- filter_results$location_f[N - 2, i] + draw_first[1:n_draws]
    }
    ##

    ## Predicting parameters in t-1 for each draw
    j <- 1
    y_with_draw_t2 <- y

    for (draw in 1:n_draws) {
      model_m2 <- filter_results$model

      ## allocating the draws
      y_with_draw_t2[N - 2, which_y_is_NAs_Nm2] <- t_Nm2_draws[draw, which_y_is_NAs_Nm2]
      NAs <- matrix(0L, nrow(y_with_draw_t2), ncol(y_with_draw_t2))
      NAs[which(!is.na(y_with_draw_t2))] <- 1L
      model_m2$NAs <- NAs
      model_m2$y <- y_with_draw_t2
      ##

      ## prediction step
      filter_results_m2 <- filter_list(model_m2)

      ## drawing missing observations in t-1
      t_Nm1_draws <- matrix(NA, n_draws, ncol(y))
      which_y_is_NAs_Nm1 <- which(is.na(y[N - 1, 1:dim(y)[2]]))
      which_y_is_NAs_Nm1 <- which_y_is_NAs_Nm1[-1] # dropping GDP

      for (i in which_y_is_NAs_Nm1) {
        draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
          filter_results_m2$scale_f[N - 1, i],
          filter_results_m2$parameters[i, 2],
          filter_results_m2$parameters[i, 3],
          filter_results_m2$shape_f[N - 1, i]
        )))
        draw_first <- draw_first[abs(draw_first) < limit]

        t_Nm1_draws[, i] <- filter_results_m2$location_f[N - 1, i] + draw_first[1:n_draws]
      }
      ##

      ## Predicting parameters in t for each draw
      y_with_draw_t1 <- y_with_draw_t2

      for (draw2 in 1:n_draws) {
        model_m1 <- model_m2

        ## allocating the draws
        y_with_draw_t1[N - 1, which_y_is_NAs_Nm1] <- t_Nm1_draws[draw2, which_y_is_NAs_Nm1]
        NAs <- matrix(0L, nrow(y_with_draw_t1), ncol(y_with_draw_t1))
        NAs[which(!is.na(y_with_draw_t1))] <- 1L
        model_m1$NAs <- NAs
        model_m1$y <- y_with_draw_t1
        ##

        ## prediction step
        filter_results_m1 <- filter_list(model_m1)

        ## drawing missing observations in t
        t_N_draws <- matrix(NA, n_draws, ncol(y))

        which_y_is_NAs_N <- which(is.na(y[N, 1:dim(y)[2]]))
        which_y_is_NAs_N <- which_y_is_NAs_N[-1] # dropping GDP

        for (i in which_y_is_NAs_N) {
          draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
            filter_results_m1$scale_f[N, i],
            filter_results_m1$parameters[i, 2],
            filter_results_m1$parameters[i, 3],
            filter_results_m1$shape_f[N, i]
          )))
          draw_first <- draw_first[abs(draw_first) < limit]

          t_N_draws[, i] <- filter_results_m1$location_f[N, i] + draw_first[1:n_draws]
        }
        ##

        ## getting filtered parameters for each draw
        y_with_draw_t <- y_with_draw_t1
        for (draw3 in 1:n_draws) {
          model_m0 <- model_m1

          ## allocating the draws
          y_with_draw_t[N, which_y_is_NAs_N] <- t_N_draws[draw3, which_y_is_NAs_N]
          NAs <- matrix(0L, nrow(y_with_draw_t), ncol(y_with_draw_t))
          NAs[which(!is.na(y_with_draw_t))] <- 1L
          model_m0$NAs <- NAs
          model_m0$y <- y_with_draw_t
          ##

          # filtering step
          filter_results_m0 <- filter_list(model_m0)

          if (model$options$freq[1] == "Q") {
            draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
              filter_results_m0$scale_f[N, 1],
              filter_results_m0$parameters[1, 2],
              filter_results_m0$parameters[1, 3],
              filter_results_m0$shape_f[N, 1]
            )))
            draw_first <- draw_first[abs(draw_first) < limit]

            GDP_draws_prediction[j:(j + n_draws - 1)] <- filter_results_m0$location_f[N, 1] + draw_first[1:n_draws]
            j <- j + n_draws
          } else {
            location_CF[, id] <- filter_results_m0$a_f[ncol(y) * 5 + 1, ]
            scale_CF[, id] <- filter_results_m0$a_scale_f[ncol(y) + 1, ]
            shape_CF[, id] <- filter_results_m0$a_shape_f[ncol(y) + 1, ]
            id <- id + 1
          }
        }
      }
    }
  }

  ############# Step 0
  if (steps == 0) {
    print_colour("step 0", "cyan")
    if (model$options$freq[1] == "Q") {
      n_draws <- 1e5
      draw_first <- rSkt(n_draws + n_draws * 0.5, t(c(
        filter_results$scale_f[nrow(y), 1],
        filter_results$parameters[1, 2],
        filter_results$parameters[1, 3],
        filter_results$shape_f[nrow(y), 1]
      )))
      draw_first <- draw_first[abs(draw_first) < limit]

      GDP_draws_prediction <- filter_results$location_f[nrow(y), 1] + draw_first[1:n_draws]
    } else {
      location_CF <- filter_results$a_f[ncol(y) * 5 + 1, ]
      scale_CF <- filter_results$a_scale_f[ncol(y) + 1, ]
      shape_CF <- filter_results$a_shape_f[ncol(y) + 1, ]
    }
  }

  if (model$options$freq[1] == "Q") {
    output <- forecast_evaluation(GDP_draws_prediction, target_gdp_figure, steps, filter_results)
  } else {
    if (steps > 0) {
      location_CF <- location_CF[, sample(min(1e3, ncol(location_CF)))]
      scale_CF <- scale_CF[, sample(min(1e3, ncol(location_CF)))]
      shape_CF <- shape_CF[, sample(min(1e3, ncol(location_CF)))]
    }

    output <- list(
      location_CF = location_CF,
      scale_CF = scale_CF,
      shape_CF = shape_CF
    )
  }

  return(output)
}

forecast_evaluation <- function(GDP_draws_prediction, target_gdp_figure, steps, MLE_filter_result) {
  stats <- list()
  stats$steps <- steps
  stats$target_gdp_figure <- target_gdp_figure

  GDP_draws_prediction <- GDP_draws_prediction[!is.na(GDP_draws_prediction)]

  density_t <- density(GDP_draws_prediction)

  ### Get empirical CDF from empirical PDF
  stats$condi_mean <- mean(GDP_draws_prediction)
  stats$CRPS <- crps_sample(target_gdp_figure, dat = GDP_draws_prediction)

  stats$sq_error <- (stats$condi_mean - target_gdp_figure)^2
  stats$abs_error <- abs(stats$condi_mean - target_gdp_figure)

  HDI <- hdi(GDP_draws_prediction, credMass = 0.68)
  stats$HDI_low <- HDI[1]
  stats$HDI_high <- HDI[2]

  stats$parameters <- c(
    target_gdp_figure - MLE_filter_result$location_f[nrow(y), 1],
    MLE_filter_result$scale_f[nrow(y), 1],
    MLE_filter_result$parameters[1, 2],
    MLE_filter_result$parameters[1, 3],
    MLE_filter_result$shape_f[nrow(y), 1]
  )

  if (steps == 0) {
    stats$log_score <- log(Skt_f(t(stats$parameters)))
    stats$PIT <- Skt_F(t(stats$parameters))
  } else {
    stats$log_score <- logscore_function(density_t, target_gdp_figure)
    stats$PIT <- pit_function(density_t, target_gdp_figure)
  }
  stats$density <- density_t
  stats$n_series <- ncol(MLE_filter_result$model$y)
  return(stats)
}

simulation_FA_LSS <- function(MLEpar, model, simulated_factors) {
  model$y <- c(model$y, NA)
  model$NAs <- c(0, which(is.na(model$y)))

  n_draws <- 100

  data_t <- matrix(NA, ncol(simulated_factors$location_CF), n_draws)

  for (i in 1:ncol(simulated_factors$location_CF)) {
    model_i <- list(
      y = model$y,
      options = model$options,
      NAs = model$NAs
    )
    N1 <- length(y)

    if (length(model$X_loc) > 1) {
      model_i$X_loc <- simulated_factors$location_CF[, i]
    } else {
      model_i$X_loc <- 0
    }

    if (length(model$X_scale) > 1) {
      model_i$X_scale <- simulated_factors$scale_CF[, i]
    } else {
      model_i$X_scale <- 0
    }

    if (length(model$X_shape) > 1) {
      model_i$X_shape <- simulated_factors$shape_CF[, i]
    } else {
      model_i$X_shape <- 0
    }

    Predic <- loglik_reg_list(MLEpar, model_i)
    data_t[i, ] <- as.vector(rSkt(n_draws, t(c(
      Predic$scale[N1],
      Predic$model$parameters[2],
      Predic$model$parameters[2],
      Predic$shape[N1]
    ))) + Predic$location[N1])
  }

  return(as.vector(data_t))
}

pit_function <- function(density, target_gdp_figure) {
  CDF <- list()
  CDF$y <- cumsum(density$y) / sum(density$y)
  CDF$x <- density$x
  PIT <- CDF$y[which(abs(density$x - target_gdp_figure) == min(abs(density$x - target_gdp_figure)))]

  return(PIT)
}

logscore_function <- function(density, target_gdp_figure) {
  LS <- log(density$y[which(abs(density$x - target_gdp_figure) == min(abs(density$x - target_gdp_figure)))])
  if (is.infinite(LS)) {
    LS <- log(min(unique(density$y[density$y != 0])))
  }
  return(LS)
}

############# Tests
Dielbolt_test <- function(pits, step, title_add = "") {
  library(binom)

  p <- c("cyan", "red", "blue")

  alpha <- 0.05
  BC_alpha <- 1 - alpha
  CI <- binom.confint(x = length(pits) / 5, n = length(pits), conf.level = BC_alpha)
  hist1 <- hist(pits, breaks = 5, plot = F)
  hist1$counts <- hist1$counts / sum(hist1$counts)

  test_outcome <- "Accept"
  if (any(hist1$counts > CI$upper[2]) | any(hist1$counts < CI$lower[2])) {
    test_outcome <- "Reject"
  }

  plot(hist1, ylim = c(0, CI$upper[2] + 0.1), col = "grey88", main = "")
  title(paste0("DGT test at step ", step, "; ", test_outcome), line = 0) # line=-1
  abline(h = CI$lower[2], col = "red", lty = 2)
  abline(h = CI$upper[2], col = "red", lty = 2)

  if (step == 1) {
    title(title_add, line = 1)
  }
}

Rossi_test <- function(pits, step) {
  p <- c("cyan", "red", "blue")

  ###### Tests after transforming the pits to a standard normal (using the quantile function, or normal inverse transform)
  ## transforming to normality
  pits[pits == 0] <- 1e-10 # otherwise gives -Inf after normal transformation
  pits[pits == 1] <- 1 - 1e-10 # otherwise gives Inf after normal transformation
  normal_pits <- qnorm(pits)
  upper <- 1.34 / sqrt(length(pits))
  lower <- -1.34 / sqrt(length(pits))

  ## Test of Rossi and Sekhposyan (2014)
  CDF_pits <- ecdf(pits)

  ## Manual empirical CDF to flag the points data cross the critical region of the test
  n <- length(pits)
  sorted_pits <- sort(pits)
  x <- seq(1, n) / n
  m_ecdf <- rep(NA, n)
  for (i in 1:n) {
    m_ecdf[i] <- length(which(pits <= x[i])) / n
  }

  test_outcome <- "Accept"
  critical_v <- 1.34 / sqrt(length(pits))
  if (any(abs(m_ecdf - x) > critical_v)) {
    test_outcome <- "Reject"
  }

  plot(CDF_pits, xlim = c(0, 1), pch = ".", main = paste0("RS test at step ", step, "; ", test_outcome))
  abline(coef = c(0, 1))
  abline(coef = c(0 + 1.34 / sqrt(length(pits)), 1), col = "red", lty = 2)
  abline(coef = c(0 - 1.34 / sqrt(length(pits)), 1), col = "red", lty = 2)
}

DGT_test <- function(pits) {
  library(binom)

  hist1 <- hist(pits, breaks = 5, plot = F)
  hist1$counts <- hist1$counts / sum(hist1$counts)

  alpha <- 0.10
  BC_alpha <- 1 - alpha
  CI <- binom.confint(x = length(pits) / 5, n = length(pits), conf.level = BC_alpha)
  test_outcome <- "Accept"
  if (any(hist1$counts > CI$upper[2]) | any(hist1$counts < CI$lower[2])) {
    test_outcome <- "Reject at 10%"
  }

  alpha <- 0.05
  BC_alpha <- 1 - alpha
  CI <- binom.confint(x = length(pits) / 5, n = length(pits), conf.level = BC_alpha)
  if (any(hist1$counts > CI$upper[2]) | any(hist1$counts < CI$lower[2])) {
    test_outcome <- "Reject at 5%"
  }

  alpha <- 0.01
  BC_alpha <- 1 - alpha
  CI <- binom.confint(x = length(pits) / 5, n = length(pits), conf.level = BC_alpha)
  if (any(hist1$counts > CI$upper[2]) | any(hist1$counts < CI$lower[2])) {
    test_outcome <- "Reject at 1%"
  }

  return(test_outcome)
}

calibration_tests <- function(pits) {
  ###### Tests after transforming the pits to a standard normal (using the quantile function, or normal inverse transform)
  ## transforming to normality
  pitsb <- pits
  pitsb[pits == 0] <- 1e-10 # otherwise gives -Inf after normal transformation
  pitsb[pits == 1] <- 1 - 1e-10 # otherwise gives Inf after normal transformation
  normal_pits <- qnorm(pitsb)

  ###### Tests on the pits directly
  ## AD test
  AD <- ad.test(pitsb)
  #

  ## Independence
  ## Box-Pierce and Ljung-Box Tests (portmanteau tests) for examining the null hypothesis of independence
  LB <- Box.test(pitsb, lag = 4, type = c("Ljung-Box"), fitdf = 0)

  ## LR test of Berkowitz
  BT <- BerkowitzTest(normal_pits)

  outcome <- round(c(AD$statistic, LB$statistic, BT$LR), 2)

  return(outcome)
}

dm_test <- function(e1, e2, alternative = c("two.sided", "less", "greater"),
                    h = 1, power = 0, varestimator = c("acf", "bartlett")) {
  # modified function from the forecast package

  alternative <- match.arg(alternative)
  varestimator <- match.arg(varestimator)

  if (power == 0) {
    d <- e1 - e2
  } else {
    d <- c(abs(e1))^power - c(abs(e2))^power
  }

  d.cov <- acf(d,
    na.action = na.omit, lag.max = h - 1, type = "covariance",
    plot = FALSE
  )$acf[, , 1]
  if (varestimator == "acf") {
    d.var <- sum(c(d.cov[1], 2 * d.cov[-1])) / length(d)
  } else {
    m <- h - 1
    b <- seq.int(from = 1, to = m)
    d.var <- sum(c(d.cov[1], 2 * (1 - b / (m + 1)) * d.cov[-1])) / length(d)
  }
  dv <- d.var
  if (dv > 0) {
    STATISTIC <- mean(d, na.rm = TRUE) / sqrt(dv)
  } else if (h == 1) {
    stop("Variance of DM statistic is zero")
  } else {
    warning("Variance is negative. Try varestimator = bartlett. Proceeding with horizon h=1.")
    return(dm.test(e1, e2, alternative, h = 1, power, varestimator))
  }
  n <- length(d)
  k <- ((n + 1 - 2 * h + (h / n) * (h - 1)) / n)^(1 / 2)
  STATISTIC <- STATISTIC * k
  names(STATISTIC) <- "DM"
  if (alternative == "two.sided") {
    PVAL <- 2 * pt(-abs(STATISTIC), df = Inf) # n - 1
  } else if (alternative == "less") {
    PVAL <- pt(STATISTIC, df = Inf)
  } else if (alternative == "greater") {
    PVAL <- pt(STATISTIC, df = Inf, lower.tail = FALSE)
  }
  PARAMETER <- c(h, power)
  names(PARAMETER) <- c("Forecast horizon", "Loss function power")
  structure(list(
    statistic = STATISTIC, parameter = PARAMETER,
    alternative = alternative, varestimator = varestimator,
    p.value = PVAL, method = "Diebold-Mariano Test", data.name = c(
      deparse(substitute(e1)),
      deparse(substitute(e2))
    )
  ), class = "htest")
}
