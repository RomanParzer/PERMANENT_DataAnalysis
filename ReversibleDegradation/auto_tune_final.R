# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # #  AUTO-TUNING: quantile-based threshold selection  # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# tau = Q_p(|nabla log y|),  p = 0.987
#
# The symmetric gradient at time t over window h is:
#   nabla log y_t = 1/(2h) * (sum_{k=1}^{h} log y_{t+k} - sum_{k=1}^{h} log y_{t-k})
#
# h is fixed at 100 for all datasets except RH100 (i=10,23) where h=8.
# tau is set manually for i=10,23 due to their distinct signal structure.
#
# p = 0.987 was selected by minimising the total absolute difference in
# detected jump counts relative to expert-validated reference counts
# across 24 datasets (excluding RH100).
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# ── vectorised symmetric gradient (O(n) via cumsum) ──────────────────────────
mygrad <- function(x, maxh) {
  cs0 <- c(0, cumsum(x))
  idx <- (maxh + 1):(length(x) - maxh)
  fwd <- cs0[idx + maxh + 1] - cs0[idx + 1]
  bwd <- cs0[idx]             - cs0[idx - maxh]
  c(rep(NA, maxh), (fwd - bwd) / (2 * maxh), rep(NA, maxh))
}

# ── quantile-based threshold ──────────────────────────────────────────────────
quantile_tau <- function(log_y, h, jump_direction = "up", p = 0.987) {
  g <- mygrad(log_y, h)
  if (jump_direction == "down") g <- -g
  as.numeric(quantile(abs(g[!is.na(g)]), p))
}

# ── main setup function (drop-in replacement for auto_tune_v2 block) ──────────
#
# Usage in main_analysis.R:
#   source("./auto_tune_final.R")
#   tuning       <- setup_tuning(numdat, p = 0.987)
#   hs           <- tuning$hs
#   grad_cutoffs <- tuning$grad_cutoffs
#
setup_tuning <- function(numdat,
                         data_path = "./data/Dati_PERMANENT_reversible_all_ordered.xlsx",
                         p         = 0.987) {

  # ── h rules ────────────────────────────────────────────────────────────────
  # h = 100  default (sampling every 60*4 time units)
  # h = 8    RH100 datasets i=10,23 (sampling every 12*4 time units)
  hs           <- rep(100L, numdat)
  hs[c(10,23)] <- 8L

  # ── manual tau for RH100 (distinct signal structure) ───────────────────────
  tau_RH100 <- c("10" = 0.010, "23" = 0.020)

  grad_cutoffs <- numeric(numdat)

  cat(sprintf("Quantile-based auto-tuning (p = %.3f) for %d datasets...\n",
              p, numdat))

  for (i in seq_len(numdat)) {

    if (i %in% c(10L, 23L)) {
      grad_cutoffs[i] <- tau_RH100[as.character(i)]
      cat(sprintf("  Dataset %2d: MANUAL (RH100)  h = %3d  tau = %.4f\n",
                  i, hs[i], grad_cutoffs[i]))
    } else {
      dat <- readxl::read_excel(data_path, sheet = i)
      colnames(dat) <- c("time", "current")
      jd  <- "up"
      tau <- quantile_tau(log(dat$current), hs[i], jd, p)
      grad_cutoffs[i] <- tau
      cat(sprintf("  Dataset %2d: auto            h = %3d  tau = %.4f\n",
                  i, hs[i], grad_cutoffs[i]))
    }
  }
  cat("Auto-tuning complete.\n\n")

  list(hs = hs, grad_cutoffs = grad_cutoffs)
}
