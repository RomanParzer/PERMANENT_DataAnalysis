# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # Benchmark Models for Comparison with RFNR # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Three benchmark model families, each with fit / predict / plot / summary S3 methods,
# structured to mirror RFNR so they can be dropped into the existing foreach loop.
#
# Models implemented:
#   1. WienerDeg  – Wiener process with linear drift on log(y), MLE via optim()
#                   Jump-augmented variant uses observed dummy indicators (same as RFNR)
#   2. GammaDeg   – Gamma process on -diff(log(y)) for the base degradation trend;
#                   jumps handled as additive Gaussian shocks on log(y)
#   3. EKF_JD     – Linear Gaussian State-Space model on log(y) estimated via the
#                   Kalman filter (FKF package); jumps enter as deterministic exogenous
#                   input (mean-shift) identified from the observed dummy vector.
#   4. SparseGP   – Sparse / local GP on log(y) via laGP package (handles 2k-10k pts).
#                   Jumps are absorbed into the mean function as step dummies.
#
# Shared evaluation helper (mirrors RFNR eval_fit):
#   eval_fit_bench(log_y_true, log_y_pred) -> c(RMSE, Cor, MAE)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pacman::p_load(optimx, FKF, laGP, dplyr, ggplot2, tidyr)

# ── shared metric (identical to RFNR eval_fit, works on log scale) ──────────────────────
eval_fit_bench <- function(log_y_true, log_y_pred) {
  resid  <- log_y_true - log_y_pred
  rmse   <- sqrt(mean(resid^2))
  cor_v  <- cor(log_y_true, log_y_pred)
  mae    <- mean(abs(resid))
  c(RMSE = rmse, Cor = cor_v, MAE = mae)
}

# ── shared cumulative-dummy builder ─────────────────────────────────────────────────────
# Converts a binary indicator vector into a step function that switches ON at each 1
# and stays ON until the next 1.  Used by EKF_JD and SparseGP mean functions.
make_step_mat <- function(dum, n) {
  jump_times <- which(dum == 1)
  if (length(jump_times) == 0) return(matrix(0, n, 0))
  Z <- matrix(0, n, length(jump_times))
  for (j in seq_along(jump_times)) {
    Z[jump_times[j]:n, j] <- 1
  }
  Z
}


# ═══════════════════════════════════════════════════════════════════════════════════════════
# 1.  WIENER DEGRADATION MODEL  (WienerDeg)
# ═══════════════════════════════════════════════════════════════════════════════════════════
#
# Model for log(y):
#   log y_t = mu + beta * t + sum_k [ delta_k * I(t >= tau_k) ] + eps_t
#   eps_t ~ N(0, sigma^2)   (i.i.d.; Wiener increments on the residual)
#
# Parameters (log-space for positivity where needed):
#   mu        – intercept
#   beta      – drift (expected negative for degradation)
#   log_sigma – log observation noise sd
#   delta_1 … delta_K – jump sizes for each observed event (if jump_direction != "none")
#
# ── constructor ─────────────────────────────────────────────────────────────────────────
WienerDeg <- function(y, x,
                      dum            = NULL,
                      ndx            = 1:length(y),
                      jump_direction = c("up", "down", "none"),
                      init           = NULL,
                      methods        = c("nvm", "subplex", "lbfgsb3c")) {
  
  jump_direction <- match.arg(jump_direction)
  x_use <- x[ndx]
  y_use <- y[ndx]
  ly    <- log(y_use)
  n     <- length(x_use)
  
  # ── jump design matrix ────
  if (jump_direction == "none" || is.null(dum)) {
    Z         <- matrix(0, n, 0)
    dum_use   <- rep(0, n)
    has_jumps <- FALSE
  } else {
    dum_use   <- if (length(dum) == length(y)) dum[ndx] else dum
    Z         <- make_step_mat(dum_use, n)
    has_jumps <- ncol(Z) > 0
  }
  
  n_jump <- ncol(Z)
  
  # ── objective (negative log-likelihood) ────
  nll <- function(par) {
    mu     <- par[1]
    beta   <- par[2]
    sigma  <- exp(par[3])
    deltas <- if (n_jump > 0) par[4:(3 + n_jump)] else numeric(0)
    mu_vec <- mu + beta * x_use
    if (n_jump > 0) mu_vec <- mu_vec + Z %*% deltas
    -sum(dnorm(ly, mu_vec, sigma, log = TRUE))
  }
  
  # ── analytical gradient ────
  gr_nll <- function(par) {
    mu     <- par[1]
    beta   <- par[2]
    sigma  <- exp(par[3])
    deltas <- if (n_jump > 0) par[4:(3 + n_jump)] else numeric(0)
    mu_vec <- mu + beta * x_use
    if (n_jump > 0) mu_vec <- mu_vec + Z %*% deltas
    e      <- ly - mu_vec
    s2     <- sigma^2
    g      <- numeric(3 + n_jump)
    g[1]   <- sum(e) / s2 * (-1)          # d/d_mu
    g[2]   <- sum(e * x_use) / s2 * (-1)  # d/d_beta
    g[3]   <- -n + sum(e^2) / s2          # d/d_log_sigma (chain rule)
    if (n_jump > 0)
      g[4:(3 + n_jump)] <- -t(Z) %*% e / s2
    g
  }
  
  # ── starting values ────
  if (is.null(init)) {
    beta0  <- coef(lm(ly ~ x_use))[2]
    mu0    <- coef(lm(ly ~ x_use))[1]
    start  <- c(mu0, beta0, log(sd(ly - mu0 - beta0 * x_use)))
    if (n_jump > 0) start <- c(start, rep(0.05, n_jump))
  } else {
    start <- init
  }
  
  # ── optimisation ────
  opm_res <- opm(start, nll, gr = gr_nll, method = methods)
  
  best_cand <- which(opm_res$kkt1 & opm_res$kkt2)
  if (length(best_cand) == 0) best_cand <- which(opm_res$convergence == 0)
  if (length(best_cand) == 0) best_cand <- seq_along(methods)
  best_ind  <- best_cand[which.min(opm_res$value[best_cand])]
  pars      <- as.numeric(opm_res[best_ind, 1:(ncol(opm_res) - 8)])
  par_names <- c("mu", "beta", "log_sigma")
  if (n_jump > 0) par_names <- c(par_names, paste0("delta_", seq_len(n_jump)))
  names(pars) <- par_names
  
  # ── fitted values ────
  mu_vec <- pars["mu"] + pars["beta"] * x_use
  fitted_base  <- exp(mu_vec)
  if (n_jump > 0) {
    deltas <- pars[4:(3 + n_jump)]
    fitted_jumps <- exp(mu_vec + Z %*% deltas)
  } else {
    fitted_jumps <- fitted_base
  }
  
  sigma <- exp(pars["log_sigma"])
  
  res <- list(
    pars           = pars,
    fitted_base    = as.numeric(fitted_base),
    fitted_jumps   = as.numeric(fitted_jumps),
    x_use          = x_use,
    y_use          = y_use,
    dum_use        = dum_use,
    Z              = Z,
    opm_res        = opm_res,
    sigma          = sigma,
    jump_direction = jump_direction,
    log_order      = 1L,
    rho_upper      = NA_real_,
    residuals      = y_use - as.numeric(fitted_jumps),
    log_residuals  = ly - as.numeric(log(fitted_jumps))
  )
  class(res) <- "WienerDeg"
  res
}

# ── S3 methods ──────────────────────────────────────────────────────────────────────────
coef.WienerDeg     <- function(res_obj, ...) res_obj$pars
residuals.WienerDeg <- function(res_obj, ...) as.numeric(res_obj$log_residuals)

predict.WienerDeg <- function(res_obj, xnew = NULL, dumnew = NULL, type = c("linear", "response"), ...) {
  type <- match.arg(type)
  if (is.null(xnew)) {
    pred_log <- log(res_obj$fitted_jumps)
  } else {
    n_new  <- length(xnew)
    pred_log <- res_obj$pars["mu"] + res_obj$pars["beta"] * xnew
    if (!is.null(dumnew) && ncol(res_obj$Z) > 0) {
      # extend step matrix with new dummies (conservative: treat as new jumps)
      dum_full <- c(res_obj$dum_use, dumnew)
      Z_full   <- make_step_mat(dum_full, length(res_obj$x_use) + n_new)
      Z_new    <- Z_full[(length(res_obj$x_use) + 1):(length(res_obj$x_use) + n_new), ,
                         drop = FALSE]
      n_old_j  <- ncol(res_obj$Z)
      n_new_j  <- ncol(Z_new) - n_old_j
      deltas   <- res_obj$pars[4:(3 + n_old_j)]
      if (n_new_j > 0) deltas <- c(deltas, rep(mean(deltas), n_new_j))
      pred_log <- pred_log + Z_new %*% deltas
    }
  }
  if (type == "response") exp(pred_log) else as.numeric(pred_log)
}

plot.WienerDeg <- function(res_obj, type = c("scatter", "residuals", "res-QQ", "res-hist", "res-acf"), ...) {
  type   <- match.arg(type)
  tmp_df <- data.frame(x      = res_obj$x_use,
                       y      = res_obj$y_use,
                       base   = res_obj$fitted_base,
                       fitted = res_obj$fitted_jumps,
                       lres   = res_obj$log_residuals)
  .bench_plot(tmp_df, type, res_obj$sigma, length(res_obj$pars), model_label = "WienerDeg")
}

summary.WienerDeg <- function(res_obj, ...) {
  .bench_summary(res_obj, "WienerDeg")
}

print.WienerDeg <- function(res_obj, ...) print(summary(res_obj))


# ═══════════════════════════════════════════════════════════════════════════════════════════
# 2.  GAMMA DEGRADATION MODEL  (GammaDeg)
# ═══════════════════════════════════════════════════════════════════════════════════════════
#
# The Gamma process is defined on POSITIVE, MONOTONE increments.
# Because log(y) is not monotone (reversible shocks), we split the signal:
#   log y_t = trend_t + jump_t + noise_t
#
# TREND: modelled as a Gamma process on cumulative degradation
#   D_t = -log(y_t / y_0) restricted to its increasing part (we use inter-jump segments)
#   Gamma process:  D_t ~ Gamma(alpha * t^beta_g, nu)
#   MLE on increments: dD_k ~ Gamma(alpha * dt_k^beta_g, nu)
#
# JUMPS: additive Gaussian shocks on log(y) at observed dummy times (same as WienerDeg)
#
# Parameters:
#   log_alpha – log shape rate
#   log_beta_g – log time-scale exponent (= 1 → homogeneous GP)
#   log_nu     – log inverse-scale (rate parameter)
#   delta_1 … delta_K – jump magnitudes
#   log_sigma_noise – residual noise after GP + jumps
#
# ── constructor ─────────────────────────────────────────────────────────────────────────
GammaDeg <- function(y, x,
                     dum            = NULL,
                     ndx            = 1:length(y),
                     jump_direction = c("up", "down", "none"),
                     init           = NULL,
                     methods        = c("nvm", "subplex", "lbfgsb3c")) {
  
  jump_direction <- match.arg(jump_direction)
  x_use <- x[ndx]
  y_use <- y[ndx]
  ly    <- log(y_use)
  n     <- length(x_use)
  
  # ── jump design matrix ────
  if (jump_direction == "none" || is.null(dum)) {
    Z         <- matrix(0, n, 0)
    dum_use   <- rep(0, n)
  } else {
    dum_use   <- if (length(dum) == length(y)) dum[ndx] else dum
    Z         <- make_step_mat(dum_use, n)
  }
  n_jump <- ncol(Z)
  
  # ── compute base degradation increments (remove jump effect first) ────
  # We subtract a step-function approximation of the jumps before computing increments.
  # At fit time the jump deltas are unknown; we iterate, but for MLE we treat jumps as
  # nuisance parameters estimated jointly.
  
  dt  <- c(x_use[1], diff(x_use))   # time increments (first point uses x[1] as dt)
  dt  <- pmax(dt, 1e-6)             # guard against duplicates
  
  # ── negative log-likelihood ────
  nll <- function(par) {
    alpha   <- exp(par[1])
    beta_g  <- exp(par[2])
    nu      <- exp(par[3])
    sig_n   <- exp(par[4])
    deltas  <- if (n_jump > 0) par[5:(4 + n_jump)] else numeric(0)
    
    # reconstructed log(y) minus jumps
    ly_base <- ly
    if (n_jump > 0) ly_base <- ly - Z %*% deltas
    
    # base degradation: D_t = y_0 - ly_base (we want positive increments downward)
    D       <- ly_base[1] - ly_base          # cumulative degradation from start
    dD      <- pmax(diff(c(0, D)), 1e-10)   # incremental degradation
    
    # Gamma increment log-likelihoods
    shape_k <- alpha * dt^beta_g
    shape_k <- pmax(shape_k, 1e-8)
    ll_gam  <- sum(dgamma(dD, shape = shape_k, rate = nu, log = TRUE))
    
    # residual noise (after GP trend + jumps removed)
    D_fitted <- cumsum(shape_k / nu)          # E[D_t] under GP
    ly_trend <- ly_base[1] - D_fitted
    ll_noise <- sum(dnorm(ly_base - ly_trend, 0, sig_n, log = TRUE))
    
    -(ll_gam + ll_noise)
  }
  
  # ── starting values ────
  if (is.null(init)) {
    start <- c(log(1e-4), log(1), log(1e-4), log(sd(diff(ly))))
    if (n_jump > 0) start <- c(start, rep(0.05, n_jump))
  } else {
    start <- init
  }
  
  # ── optimisation ────
  opm_res <- opm(start, nll, method = methods,
                 control = list(maxit = 5000))
  best_cand <- which(opm_res$kkt1 & opm_res$kkt2)
  if (length(best_cand) == 0) best_cand <- which(opm_res$convergence == 0)
  if (length(best_cand) == 0) best_cand <- seq_along(methods)
  best_ind  <- best_cand[which.min(opm_res$value[best_cand])]
  pars      <- as.numeric(opm_res[best_ind, 1:(ncol(opm_res) - 8)])
  par_names <- c("log_alpha", "log_beta_g", "log_nu", "log_sigma_noise")
  if (n_jump > 0) par_names <- c(par_names, paste0("delta_", seq_len(n_jump)))
  names(pars) <- par_names
  
  # ── fitted values ────
  alpha   <- exp(pars["log_alpha"])
  beta_g  <- exp(pars["log_beta_g"])
  nu      <- exp(pars["log_nu"])
  deltas  <- if (n_jump > 0) pars[5:(4 + n_jump)] else numeric(0)
  shape_k <- alpha * dt^beta_g
  D_fitted <- cumsum(shape_k / nu)
  fitted_base_log  <- ly[1] - D_fitted
  fitted_base      <- exp(fitted_base_log)
  if (n_jump > 0) {
    fitted_jumps <- exp(fitted_base_log + Z %*% deltas)
  } else {
    fitted_jumps <- fitted_base
  }
  sigma <- exp(pars["log_sigma_noise"])
  
  res <- list(
    pars           = pars,
    fitted_base    = as.numeric(fitted_base),
    fitted_jumps   = as.numeric(fitted_jumps),
    x_use          = x_use,
    y_use          = y_use,
    dum_use        = dum_use,
    Z              = Z,
    dt             = dt,
    opm_res        = opm_res,
    sigma          = sigma,
    jump_direction = jump_direction,
    log_order      = 1L,
    rho_upper      = NA_real_,
    residuals      = y_use - as.numeric(fitted_jumps),
    log_residuals  = ly - log(as.numeric(fitted_jumps))
  )
  class(res) <- "GammaDeg"
  res
}

# ── S3 methods ──────────────────────────────────────────────────────────────────────────
coef.GammaDeg      <- function(res_obj, ...) res_obj$pars
residuals.GammaDeg <- function(res_obj, ...) as.numeric(res_obj$log_residuals)

predict.GammaDeg <- function(res_obj, xnew = NULL, dumnew = NULL, type = c("linear", "response"), ...) {
  type <- match.arg(type)
  if (is.null(xnew)) {
    pred_log <- log(res_obj$fitted_jumps)
  } else {
    # Extrapolate GP trend from last training point
    alpha  <- exp(res_obj$pars["log_alpha"])
    beta_g <- exp(res_obj$pars["log_beta_g"])
    nu     <- exp(res_obj$pars["log_nu"])
    dt_new <- c(xnew[1] - tail(res_obj$x_use, 1), diff(xnew))
    dt_new <- pmax(dt_new, 1e-6)
    D_extra  <- cumsum(alpha * dt_new^beta_g / nu)
    pred_log <- log(tail(res_obj$fitted_base, 1)) - D_extra
    
    if (!is.null(dumnew) && ncol(res_obj$Z) > 0) {
      n_old   <- length(res_obj$x_use)
      n_new   <- length(xnew)
      dum_full <- c(res_obj$dum_use, dumnew)
      Z_full   <- make_step_mat(dum_full, n_old + n_new)
      Z_new    <- Z_full[(n_old + 1):(n_old + n_new), , drop = FALSE]
      n_old_j  <- ncol(res_obj$Z)
      n_new_j  <- ncol(Z_new) - n_old_j
      deltas   <- res_obj$pars[5:(4 + n_old_j)]
      if (n_new_j > 0) deltas <- c(deltas, rep(mean(deltas), n_new_j))
      pred_log <- pred_log + Z_new %*% deltas
    }
  }
  if (type == "response") exp(pred_log) else as.numeric(pred_log)
}

plot.GammaDeg <- function(res_obj, type = c("scatter", "residuals", "res-QQ", "res-hist", "res-acf"), ...) {
  type   <- match.arg(type)
  tmp_df <- data.frame(x      = res_obj$x_use,
                       y      = res_obj$y_use,
                       base   = res_obj$fitted_base,
                       fitted = res_obj$fitted_jumps,
                       lres   = res_obj$log_residuals)
  .bench_plot(tmp_df, type, res_obj$sigma, length(res_obj$pars), model_label = "GammaDeg")
}

summary.GammaDeg <- function(res_obj, ...) .bench_summary(res_obj, "GammaDeg")
print.GammaDeg   <- function(res_obj, ...) print(summary(res_obj))


# ═══════════════════════════════════════════════════════════════════════════════════════════
# 3.  EXTENDED KALMAN FILTER WITH JUMP-DIFFUSION  (EKF_JD)
# ═══════════════════════════════════════════════════════════════════════════════════════════
#
# State-space model on log(y_t):
#   State:   alpha_t  = alpha_{t-1} + beta * dt + sigma_w * w_t    w_t ~ N(0,1)
#   Obs:     log y_t  = alpha_t + jump_t + sigma_v * v_t            v_t ~ N(0,1)
#   Jump:    jump_t   = delta * dum_t  (single shared delta — avoids overfitting)
#
# Speed: uses FKF::fkf() (C implementation) instead of R loop — ~50x faster.
# Parameters: beta, log_sigma_w, log_sigma_v, delta (4 total)
#
# ── constructor ─────────────────────────────────────────────────────────────────────────
EKF_JD <- function(y, x,
                   dum            = NULL,
                   ndx            = 1:length(y),
                   jump_direction = c("up", "down", "none"),
                   init           = NULL,
                   methods        = c("nvm", "lbfgsb3c")) {
  
  if (!requireNamespace("FKF", quietly = TRUE)) stop("Please install the FKF package.")
  
  jump_direction <- match.arg(jump_direction)
  x_use  <- x[ndx]
  y_use  <- y[ndx]
  ly     <- log(y_use)
  n      <- length(x_use)
  dt_vec <- c(x_use[1], diff(x_use))
  
  if (jump_direction == "none" || is.null(dum)) {
    dum_use   <- rep(0, n)
    has_jumps <- FALSE
  } else {
    dum_use   <- if (length(dum) == length(y)) dum[ndx] else dum
    has_jumps <- any(dum_use == 1)
  }
  
  # ── fast KF log-likelihood via FKF (C backend) ───────────────────────────────────────
  # FKF::fkf() expects arrays: all system matrices as [m x m x n] or [p x m x n]
  # Here m=1 (state dim), p=1 (obs dim), so everything is scalar wrapped in arrays.
  #
  # State eq:  a_t = a_{t-1} + beta*dt_t + w_t,   w_t ~ N(0, sig_w^2 * dt_t)
  # Obs eq:    ly_t - delta*dum_t = a_t + v_t,     v_t ~ N(0, sig_v^2)
  #
  nll_kf <- function(par) {
    beta   <- par[1]
    sig_w  <- exp(par[2])
    sig_v  <- exp(par[3])
    delta  <- if (has_jumps) par[4] else 0
    
    # FKF expects dt/ct as matrices (1 x n), Tt/Zt/HHt/GGt as arrays (1 x 1 x n)
    ct  <- matrix(beta * dt_vec,      nrow = 1)   # state intercept
    dt_obs <- matrix(-delta * dum_use, nrow = 1)  # obs intercept (subtract jump)
    Tt  <- array(1,        dim = c(1, 1, n))
    Zt  <- array(1,        dim = c(1, 1, n))
    Qt  <- array(sig_w^2 * dt_vec, dim = c(1, 1, n))
    HHt <- array(sig_v^2, dim = c(1, 1, n))
    
    kf <- tryCatch(
      FKF::fkf(a0  = ly[1],
               P0  = matrix(sig_v^2),
               dt  = ct,
               ct  = dt_obs,
               Tt  = Tt,
               Zt  = Zt,
               HHt = Qt,
               GGt = HHt,
               yt  = matrix(ly, nrow = 1)),
      error = function(e) NULL
    )
    if (is.null(kf) || !all(kf$status == 0)) return(1e12)
    -kf$logLik
  }
  
  # ── starting values ────
  if (is.null(init)) {
    beta0 <- coef(lm(ly ~ x_use))[2]
    start <- c(beta0, log(1e-3), log(1e-2))
    if (has_jumps) start <- c(start, 0.05)
  } else {
    start <- init
  }
  
  opm_res   <- opm(start, nll_kf, method = methods, control = list(maxit = 5000))
  best_cand <- which(opm_res$kkt1 & opm_res$kkt2)
  if (length(best_cand) == 0) best_cand <- which(opm_res$convergence == 0)
  if (length(best_cand) == 0) best_cand <- seq_along(methods)
  best_ind  <- best_cand[which.min(opm_res$value[best_cand])]
  pars      <- as.numeric(opm_res[best_ind, 1:(ncol(opm_res) - 8)])
  par_names <- c("beta", "log_sigma_w", "log_sigma_v")
  if (has_jumps) par_names <- c(par_names, "delta")
  names(pars) <- par_names
  
  # ── run FKF one final time to get filtered state ─────────────────────────────────────
  beta   <- pars["beta"]
  sig_w  <- exp(pars["log_sigma_w"])
  sig_v  <- exp(pars["log_sigma_v"])
  delta  <- if (has_jumps) pars["delta"] else 0
  
  ct_f    <- matrix(beta * dt_vec,       nrow = 1)
  dt_f    <- matrix(-delta * dum_use,    nrow = 1)
  Tt_f    <- array(1,                    dim = c(1, 1, n))
  Qt_f    <- array(sig_w^2 * dt_vec,     dim = c(1, 1, n))
  Zt_f    <- array(1,                    dim = c(1, 1, n))
  HHt_f   <- array(sig_v^2,             dim = c(1, 1, n))
  
  kf_final <- FKF::fkf(a0  = ly[1],
                       P0  = matrix(sig_v^2),
                       dt  = ct_f,
                       ct  = dt_f,
                       Tt  = Tt_f,
                       Zt  = Zt_f,
                       HHt = Qt_f,
                       GGt = HHt_f,
                       yt  = matrix(ly, nrow = 1))
  
  # att[1,1,] = filtered state mean at each time point
  a_filt   <- as.numeric(kf_final$att[1, ])
  jump_vec <- delta * dum_use
  Z        <- matrix(dum_use, ncol = 1)
  
  fitted_base  <- exp(a_filt)
  fitted_jumps <- exp(a_filt + jump_vec)
  sigma        <- exp(pars["log_sigma_v"])
  
  res <- list(
    pars           = pars,
    fitted_base    = fitted_base,
    fitted_jumps   = fitted_jumps,
    x_use          = x_use,
    y_use          = y_use,
    dum_use        = dum_use,
    Z              = Z,
    dt_vec         = dt_vec,
    opm_res        = opm_res,
    sigma          = sigma,
    jump_direction = jump_direction,
    log_order      = 1L,
    rho_upper      = NA_real_,
    residuals      = y_use - fitted_jumps,
    log_residuals  = ly - log(fitted_jumps)
  )
  class(res) <- "EKF_JD"
  res
}

# ── S3 methods ──────────────────────────────────────────────────────────────────────────
coef.EKF_JD      <- function(res_obj, ...) res_obj$pars
residuals.EKF_JD <- function(res_obj, ...) as.numeric(res_obj$log_residuals)

predict.EKF_JD <- function(res_obj, xnew = NULL, dumnew = NULL, type = c("linear", "response"), ...) {
  type <- match.arg(type)
  if (is.null(xnew)) {
    pred_log <- log(res_obj$fitted_jumps)
  } else {
    beta       <- res_obj$pars["beta"]
    n_new      <- length(xnew)
    dt_new     <- c(xnew[1] - tail(res_obj$x_use, 1), diff(xnew))
    dt_new     <- pmax(dt_new, 1e-6)
    a_last     <- log(tail(res_obj$fitted_base, 1))
    a_pred_log <- a_last + beta * cumsum(dt_new)
    
    # single shared delta applied to new dummy indicators
    jump_vec_new <- rep(0, n_new)
    if (!is.null(dumnew) && "delta" %in% names(res_obj$pars)) {
      jump_vec_new <- res_obj$pars["delta"] * dumnew
    }
    pred_log <- a_pred_log + jump_vec_new
  }
  if (type == "response") exp(pred_log) else as.numeric(pred_log)
}

plot.EKF_JD <- function(res_obj, type = c("scatter", "residuals", "res-QQ", "res-hist", "res-acf"), ...) {
  type   <- match.arg(type)
  tmp_df <- data.frame(x      = res_obj$x_use,
                       y      = res_obj$y_use,
                       base   = res_obj$fitted_base,
                       fitted = res_obj$fitted_jumps,
                       lres   = res_obj$log_residuals)
  .bench_plot(tmp_df, type, res_obj$sigma, length(res_obj$pars), model_label = "EKF_JD")
}

summary.EKF_JD <- function(res_obj, ...) .bench_summary(res_obj, "EKF_JD")
print.EKF_JD   <- function(res_obj, ...) print(summary(res_obj))


# ═══════════════════════════════════════════════════════════════════════════════════════════
# 4.  SPARSE / LOCAL GAUSSIAN PROCESS  (SparseGP)
# ═══════════════════════════════════════════════════════════════════════════════════════════
#
# For 2k-10k points, full GP is O(n^3) and infeasible.
# We use laGP::aGP() (local approximate GP), which fits a GP using a neighbourhood of
# m << n nearest points around each prediction location.
#
# Mean function: linear trend + step dummies for jumps (detrend first, GP on residuals)
# Kernel: Matérn 5/2 (default in laGP via "matern5_2" correlation)
#
# Workflow:
#   1. Fit linear mean + jump dummies via OLS → get detrended residuals r_t
#   2. Fit laGP on (x, r) to learn the residual covariance structure
#   3. Predictions = OLS mean + GP correction
#
# Parameters returned: OLS coefficients + GP hyperparameters (d, g) from laGP MLE
#
# ── constructor ─────────────────────────────────────────────────────────────────────────
SparseGP <- function(y, x,
                     dum            = NULL,
                     ndx            = 1:length(y),
                     jump_direction = c("up", "down", "none"),
                     m_local        = 50L) {   # neighbourhood size for aGP
  
  jump_direction <- match.arg(jump_direction)
  x_use  <- x[ndx]
  y_use  <- y[ndx]
  ly     <- log(y_use)
  n      <- length(x_use)
  
  if (jump_direction == "none" || is.null(dum)) {
    Z         <- matrix(0, n, 0)
    dum_use   <- rep(0, n)
  } else {
    dum_use   <- if (length(dum) == length(y)) dum[ndx] else dum
    Z         <- make_step_mat(dum_use, n)
  }
  n_jump <- ncol(Z)
  
  # ── step 1: OLS mean ────
  if (n_jump > 0) {
    X_ols <- cbind(1, x_use, Z)
  } else {
    X_ols <- cbind(1, x_use)
  }
  ols_fit <- lm(ly ~ X_ols - 1)
  beta_ols <- coef(ols_fit)
  r_ols    <- as.numeric(residuals(ols_fit))
  
  # ── step 2: local GP on residuals ────
  # aGP needs matrix input for X; we use 1-D time scaled to [0,1]
  x_sc   <- (x_use - min(x_use)) / diff(range(x_use))
  X_mat  <- matrix(x_sc, ncol = 1)
  
  # Fit GP hyperparameters on a subsample (max 2000 pts) to speed up MLE
  sub_n   <- min(2000L, n)
  sub_idx <- sort(sample(n, sub_n))
  gp_fit  <- laGP::newGPsep(X_mat[sub_idx, , drop = FALSE],
                            r_ols[sub_idx],
                            d = 0.1, g = 0.01 * var(r_ols),
                            dK = TRUE)
  hp_mle <- laGP::mleGPsep(gp_fit, param = "both",
                           tmin = c(1e-6, 1e-8), tmax = c(10, var(r_ols)))
  hp <- list(theta = hp_mle$d, g = hp_mle$g)
  laGP::deleteGPsep(gp_fit)
  
  # ── step 3: aGP prediction at training locations ────
  # We predict in batches to avoid memory overload
  batch_size <- 500L
  r_pred <- numeric(n)
  r_var  <- numeric(n)
  
  for (start_i in seq(1, n, by = batch_size)) {
    end_i   <- min(start_i + batch_size - 1, n)
    idx_b   <- start_i:end_i
    agp_out <- laGP::aGP(X = X_mat, Z = r_ols,
                         XX = X_mat[idx_b, , drop = FALSE],
                         start = 6L, end = m_local,
                         d = list(mle = FALSE, start = hp$theta),
                         g = list(mle = FALSE, start = hp$g),
                         method = "nn", verb = 0)
    r_pred[idx_b] <- agp_out$mean
    r_var[idx_b]  <- agp_out$var
  }
  
  fitted_log_base  <- as.numeric(X_ols[, 1:2] %*% beta_ols[1:2])
  fitted_log_jumps <- as.numeric(X_ols %*% beta_ols)
  fitted_base      <- exp(fitted_log_base  + r_pred)
  fitted_jumps     <- exp(fitted_log_jumps + r_pred)
  sigma            <- sqrt(mean(r_var))
  
  res <- list(
    pars           = beta_ols,
    hp             = hp,
    beta_ols       = beta_ols,
    fitted_base    = fitted_base,
    fitted_jumps   = fitted_jumps,
    x_use          = x_use,
    x_sc           = x_sc,
    X_mat          = X_mat,
    r_ols          = r_ols,
    y_use          = y_use,
    dum_use        = dum_use,
    Z              = Z,
    m_local        = m_local,
    opm_res        = NULL,
    sigma          = sigma,
    jump_direction = jump_direction,
    log_order      = 1L,
    rho_upper      = NA_real_,
    residuals      = y_use - fitted_jumps,
    log_residuals  = ly - log(fitted_jumps)
  )
  class(res) <- "SparseGP"
  res
}

# ── S3 methods ──────────────────────────────────────────────────────────────────────────
coef.SparseGP      <- function(res_obj, ...) res_obj$pars
residuals.SparseGP <- function(res_obj, ...) as.numeric(res_obj$log_residuals)

predict.SparseGP <- function(res_obj, xnew = NULL, dumnew = NULL, type = c("linear", "response"), ...) {
  type <- match.arg(type)
  if (is.null(xnew)) {
    pred_log <- log(res_obj$fitted_jumps)
  } else {
    n_new  <- length(xnew)
    n_old  <- length(res_obj$x_use)
    x_range <- diff(range(res_obj$x_use))
    x_min   <- min(res_obj$x_use)
    xnew_sc <- (xnew - x_min) / x_range
    Xnew_mat <- matrix(xnew_sc, ncol = 1)
    
    # OLS mean at new points
    beta_ols <- res_obj$beta_ols
    if (!is.null(dumnew) && ncol(res_obj$Z) > 0) {
      dum_full <- c(res_obj$dum_use, dumnew)
      Z_full   <- make_step_mat(dum_full, n_old + n_new)
      Z_new    <- Z_full[(n_old + 1):(n_old + n_new), , drop = FALSE]
      n_old_j  <- ncol(res_obj$Z)
      n_new_j  <- ncol(Z_new) - n_old_j
      delta_ols <- beta_ols[-(1:2)]
      if (n_new_j > 0) delta_ols <- c(delta_ols, rep(mean(delta_ols), n_new_j))
      mean_new <- beta_ols[1] + beta_ols[2] * xnew + Z_new %*% delta_ols
    } else {
      mean_new <- beta_ols[1] + beta_ols[2] * xnew
    }
    
    # aGP correction in batches
    batch_size <- 500L
    r_pred_new <- numeric(n_new)
    for (start_i in seq(1, n_new, by = batch_size)) {
      end_i  <- min(start_i + batch_size - 1, n_new)
      idx_b  <- start_i:end_i
      agp_out <- laGP::aGP(X      = res_obj$X_mat,
                           Z      = res_obj$r_ols,
                           XX     = Xnew_mat[idx_b, , drop = FALSE],
                           start  = 6L, end = res_obj$m_local,
                           d      = list(mle = FALSE, start = res_obj$hp$theta),
                           g      = list(mle = FALSE, start = res_obj$hp$g),
                           method = "nn", verb = 0)
      r_pred_new[idx_b] <- agp_out$mean
    }
    pred_log <- as.numeric(mean_new) + r_pred_new
  }
  if (type == "response") exp(pred_log) else as.numeric(pred_log)
}

plot.SparseGP <- function(res_obj, type = c("scatter", "residuals", "res-QQ", "res-hist", "res-acf"), ...) {
  type   <- match.arg(type)
  tmp_df <- data.frame(x      = res_obj$x_use,
                       y      = res_obj$y_use,
                       base   = res_obj$fitted_base,
                       fitted = res_obj$fitted_jumps,
                       lres   = res_obj$log_residuals)
  .bench_plot(tmp_df, type, res_obj$sigma, length(res_obj$pars), model_label = "SparseGP")
}

summary.SparseGP <- function(res_obj, ...) .bench_summary(res_obj, "SparseGP")
print.SparseGP   <- function(res_obj, ...) print(summary(res_obj))


# ═══════════════════════════════════════════════════════════════════════════════════════════
# SHARED INTERNAL HELPERS  (not exported to user)
# ═══════════════════════════════════════════════════════════════════════════════════════════

# ── unified diagnostic plots ─────────────────────────────────────────────────────────────
.bench_plot <- function(tmp_df, type, sigma, n_pars, model_label) {
  if (type == "scatter") {
    tmp_df %>%
      pivot_longer(c(y, base, fitted), names_to = "type", values_to = "value") %>%
      ggplot(aes(x = x, y = value, col = type)) +
      geom_line() +
      labs(y = " ", title = model_label) +
      theme_classic()
  } else if (type == "residuals") {
    ggplot(tmp_df, aes(x = x, y = lres)) +
      geom_point(size = 0.4) +
      geom_hline(yintercept = 0, col = 2) +
      labs(y = "log residuals", title = model_label) +
      theme_classic()
  } else if (type == "res-QQ") {
    mysd <- sqrt(sum(tmp_df$lres^2) / (nrow(tmp_df) - n_pars))
    ggplot(tmp_df, aes(sample = lres / mysd)) +
      stat_qq() + stat_qq_line() +
      labs(x = "Theoretical quantiles", y = "Standardised log-residuals",
           title = model_label) +
      theme_classic()
  } else if (type == "res-hist") {
    ggplot(tmp_df, aes(x = lres)) +
      geom_histogram(aes(y = after_stat(density)), colour = "black", fill = "white") +
      stat_function(fun = dnorm, args = list(mean = 0, sd = sigma)) +
      labs(x = "log residuals", title = model_label) +
      theme_classic()
  } else if (type == "res-acf") {
    forecast::ggAcf(tmp_df$lres) +
      labs(title = model_label, x = "lag", y = "acf")
  }
}

# ── unified summary ──────────────────────────────────────────────────────────────────────
.bench_summary <- function(res_obj, model_label) {
  n    <- length(res_obj$x_use)
  k    <- length(res_obj$pars)
  ss   <- sum(res_obj$log_residuals^2)
  rmse <- sqrt(mean(res_obj$log_residuals^2))
  mae  <- mean(abs(res_obj$log_residuals))
  r2   <- 1 - ss / sum((log(res_obj$y_use) - mean(log(res_obj$y_use)))^2)
  ll   <- -n / 2 * log(2 * pi * res_obj$sigma^2) - ss / (2 * res_obj$sigma^2)
  aic  <- -2 * ll + 2 * k
  bic  <- -2 * ll + log(n) * k
  cat(sprintf("── %s Summary ────────────────────────────────\n", model_label))
  cat(sprintf("  n = %d,  k = %d\n", n, k))
  cat(sprintf("  RMSE (log) = %.5f   MAE (log) = %.5f   R2 = %.4f\n", rmse, mae, r2))
  cat(sprintf("  Approx. log-lik = %.2f   AIC = %.2f   BIC = %.2f\n", ll, aic, bic))
  cat("  Parameters:\n")
  print(round(res_obj$pars, 5))
  invisible(list(pars = res_obj$pars, RMSE = rmse, MAE = mae, R2 = r2,
                 logLik = ll, AIC = aic, BIC = bic))
}
