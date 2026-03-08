# ── carica funzioni ────────────────────────────────────────────────────────
source("./RFNR_model_functions.R")
source("./benchmark_model_functions.R")

pacman::p_load(readxl, dplyr, ggplot2, tidyr, optimx, laGP)

# ── scegli il dataset (cambia i qui) ──────────────────────────────────────
i <- 1

# ── carica e prepara dati ─────────────────────────────────────────────────
data <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx", sheet = i)
colnames(data) <- c("time", "current")
data$time <- (data$time - 1) * 60

jump_direction <- if (i == 10) "down" else "up"

grad        <- mygrad(log(data$current), hs[i])
grad_cutoff <- grad_cutoffs[i]
if (i == 10) grad <- -grad
ind <- grad > grad_cutoff
for (k in seq_len(length(ind) - 1)) {
  if (!is.na(ind[k] * ind[k+1]) && ind[k] && ind[k+1]) ind[k] <- FALSE
}
if (i %in% c(10, 13, 26)) ind[1] <- TRUE

x   <- data$time
y   <- data$current
ndx <- if (i %in% c(10, 23)) seq(1, length(x), 12*4) else seq(1, length(x), 60*4)

dum <- numeric(length(ndx))
for (k in which(ind)) dum[which.min(abs(k - ndx))] <- 1

ind_train <- 1:round(0.7 * length(ndx))
ind_test  <- (round(0.7 * length(ndx)) + 1):length(ndx)

# ── fit modelli ───────────────────────────────────────────────────────────
cat("Fitting RFNR Full...\n")
fit_rfnr <- RFNR(y, x, dum, ndx, log_order=2, jump_direction=jump_direction,
                 init=c(-1,-10,-1,-10,0,1))

cat("Fitting WienerDeg...\n")
fit_wiener <- WienerDeg(y, x, dum, ndx, jump_direction=jump_direction)

cat("Fitting GammaDeg...\n")
fit_gamma <- GammaDeg(y, x, dum, ndx, jump_direction=jump_direction)

cat("Fitting EKF_JD...\n")
fit_ekf <- EKF_JD(y, x, dum, ndx, jump_direction=jump_direction)

cat("Fitting SparseGP...\n")
fit_gp <- SparseGP(y, x, dum, ndx, jump_direction=jump_direction, m_local=50L)

fit_gp <- tryCatch(
  SparseGP(y, x, dum, ndx, jump_direction=jump_direction, m_local=50L),
  error = function(e) { cat("ERRORE SparseGP:", conditionMessage(e), "\n"); NULL }
)

# ── metriche in-sample ────────────────────────────────────────────────────
models_fitted <- list(
  RFNR     = fit_rfnr,
  WienerDeg = fit_wiener,
  GammaDeg  = fit_gamma,
  EKF_JD    = fit_ekf,
  SparseGP  = fit_gp
)

metrics <- lapply(names(models_fitted), function(nm) {
  m <- models_fitted[[nm]]
  data.frame(
    model = nm,
    RMSE  = sqrt(mean(m$log_residuals^2)),
    Cor   = cor(log(m$y_use), log(m$fitted_jumps)),
    MAE   = mean(abs(m$log_residuals))
  )
}) |> bind_rows()

print(metrics)


plot_df <- data.frame(
  x        = x[ndx],
  observed = y[ndx],
  RFNR     = fit_rfnr$fitted_jumps,
  Wiener   = fit_wiener$fitted_jumps,
  Gamma    = fit_gamma$fitted_jumps,
  EKF      = fit_ekf$fitted_jumps,
  SparseGP = fit_gp$fitted_jumps
) |>
  pivot_longer(-c(x, observed), names_to="model", values_to="fitted")

ggplot(plot_df, aes(x=x)) +
  geom_line(aes(y=observed), color="black", linewidth=0.4, alpha=0.5) +
  geom_line(aes(y=fitted, color=model), linewidth=0.6) +
  facet_wrap(~model, ncol=2) +
  labs(x="time (s)", y="current", title=paste("Dataset:", datanames[i])) +
  theme_classic() +
  theme(legend.position="none")



pred_test_ekf <- predict(fit_ekf, 
                         xnew   = x[ndx[ind_test]], 
                         dumnew = dum[ind_test], 
                         type   = "linear")
cat("EKF_JD out-of-sample RMSE:", 
    sqrt(mean((log(y[ndx[ind_test]]) - pred_test_ekf)^2)), "\n")



# parametri di tuinnig  ---------------------------------------------------
auto_tune <- function(log_y, jump_direction = "up",
                      h_grid   = c(8, 20, 50, 80, 100, 150, 200),
                      tau_grid = seq(0.005, 0.15, by = 0.005)) {
  
  results <- expand.grid(h = h_grid, tau = tau_grid)
  results$n_jumps <- NA
  
  for (r in 1:nrow(results)) {
    h   <- results$h[r]
    tau <- results$tau[r]
    grad <- mygrad(log_y, h)
    if (jump_direction == "down") grad <- -grad
    ind  <- !is.na(grad) & grad > tau
    for (k in seq_len(length(ind)-1))
      if (!is.na(ind[k]*ind[k+1]) && ind[k] && ind[k+1]) ind[k] <- FALSE
    results$n_jumps[r] <- sum(ind, na.rm=TRUE)
  }
  
  # plateau più stabile per ogni h
  stability <- results %>%
    group_by(h, n_jumps) %>%
    summarise(tau_range = diff(range(tau)),
              tau_lo    = min(tau),   # bordo basso del plateau
              tau_mid   = mean(range(tau)),
              .groups   = "drop") %>%
    filter(n_jumps > 0) %>%
    arrange(desc(tau_range))
  
  best <- stability[1, ]
  
  # usa tau al 25° percentile del plateau (più sensibile, meno salti persi)
  tau_chosen <- best$tau_lo + 0.25 * best$tau_range
  
  list(h          = best$h,
       tau        = tau_chosen,
       n_jumps    = best$n_jumps,
       tau_range  = best$tau_range,
       stability  = stability)
}

# prova su dataset 1
res_auto <- auto_tune(log_y, jump_direction = "up")
cat(sprintf("h=%d  tau=%.4f  n_jumps=%d  (plateau width=%.3f)\n",
            res_auto$h, res_auto$tau, res_auto$n_jumps, res_auto$tau_range))
cat(sprintf("Manuale: h=%d  tau=%.4f\n", hs[1], grad_cutoffs[1]))