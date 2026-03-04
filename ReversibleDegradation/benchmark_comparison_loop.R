# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # BENCHMARK COMPARISON LOOP – drop-in companion to the RFNR main script # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
# Prerequisites: run the RFNR main script first (or at least the setup section) so that
#   datanames, hs, grad_cutoffs, mygrad are already in the environment.
# This script re-uses the same ndx / dum / ind logic and produces a model_comp_bench
#   data frame with the same structure as model_comp from the RFNR loop.
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if ("pacman" %in% rownames(installed.packages()) == FALSE) install.packages("pacman")
pacman::p_load(foreach, parallel, dplyr, tidyr, ggplot2, readxl, optimx, laGP)

source("./RFNR_model_functions.R")
source("./benchmark_model_functions.R")   # <-- new file

# ── parallel cluster (reuse or recreate) ────────────────────────────────────────────────
unlink("./saved_results/log_bench.txt")
n.cores    <- parallel::detectCores()
my.cluster <- parallel::makeCluster(n.cores - 1, type = "PSOCK",
                                    outfile = "./saved_results/log_bench.txt")
doParallel::registerDoParallel(cl = my.cluster)

clusterExport(my.cluster,
              c("datanames", "models", "mygrad", "hs", "grad_cutoffs",
                "make_step_mat", "eval_fit_bench",
                "WienerDeg", "GammaDeg", "EKF_JD", "SparseGP",
                "predict.WienerDeg", "predict.GammaDeg",
                "predict.EKF_JD",    "predict.SparseGP",
                "coef.WienerDeg",    "coef.GammaDeg",
                "coef.EKF_JD",       "coef.SparseGP",
                "residuals.WienerDeg","residuals.GammaDeg",
                "residuals.EKF_JD",  "residuals.SparseGP",
                ".bench_plot", ".bench_summary"),
              envir = environment())

clusterEvalQ(my.cluster, {
  pacman::p_load(foreach, parallel, dplyr, tidyr, ggplot2, readxl, optimx, laGP)
  source("./RFNR_model_functions.R")
  source("./benchmark_model_functions.R")
})

# ── helper: shared data-prep for a single dataset index ─────────────────────────────────
.prep_dataset <- function(i) {
  data <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx", sheet = i)
  colnames(data) <- c("time", "current")
  data$time <- (data$time - 1) * 60
  
  jump_direction <- if (i == 10) "down" else "up"
  
  grad        <- mygrad(log(data$current), hs[i])
  grad_cutoff <- grad_cutoffs[i]
  if (i == 10) grad <- -grad
  ind <- grad > grad_cutoff
  for (k in seq_len(length(ind) - 1)) {
    if (!is.na(ind[k] * ind[k + 1]) && ind[k] && ind[k + 1]) ind[k] <- FALSE
  }
  if (i %in% c(10, 13, 26)) ind[1] <- TRUE
  
  x   <- data$time
  y   <- data$current
  ndx <- if (i %in% c(10, 23)) seq(1, length(x), 12 * 4) else seq(1, length(x), 60 * 4)
  xnew <- x[-ndx]
  
  dum <- numeric(length(ndx))
  for (k in which(ind)) dum[which.min(abs(k - ndx))] <- 1
  
  ind_train <- 1:round(0.7 * length(ndx))
  ind_test  <- (round(0.7 * length(ndx)) + 1):length(ndx)
  
  tmp_data_name <- if (i < 14) datanames[i] else datanames[i - 13]
  tmp_material  <- if (i < 14) "MEA1" else "MEA2"
  
  list(x = x, y = y, ndx = ndx, xnew = xnew, dum = dum,
       ind_train = ind_train, ind_test = ind_test,
       jump_direction = jump_direction,
       tmp_data_name = tmp_data_name, tmp_material = tmp_material)
}

# ── helper: fit one model on one dataset, return metric rows ─────────────────────────────
.fit_one <- function(model_fn, model_name, d) {
  out <- data.frame(data=NULL, material=NULL, model=NULL,
                    error=NULL, error_measure=NULL, value=NULL)
  
  # train/test
  modres_train <- try(model_fn(d$y, d$x, d$dum[d$ind_train],
                               d$ndx[d$ind_train], d$jump_direction), silent=TRUE)
  if (!inherits(modres_train, "try-error")) {
    pred_test <- try(predict(modres_train, xnew = d$x[d$ndx[d$ind_test]],
                             dumnew = d$dum[d$ind_test], type = "linear"), silent=TRUE)
    if (!inherits(pred_test, "try-error")) {
      out <- rbind(out, data.frame(
        data=d$tmp_data_name, material=d$tmp_material, model=model_name,
        error="train/test", error_measure=c("RMSE","Cor","MAE"),
        value=as.numeric(eval_fit_bench(log(d$y[d$ndx[d$ind_test]]), pred_test))))
    }
  }
  
  # full sample
  modres_full <- try(model_fn(d$y, d$x, d$dum, d$ndx, d$jump_direction), silent=TRUE)
  if (!inherits(modres_full, "try-error")) {
    out <- rbind(out, data.frame(
      data=d$tmp_data_name, material=d$tmp_material, model=model_name,
      error="sample fit", error_measure=c("RMSE","Cor","MAE"),
      value=as.numeric(eval_fit_bench(log(d$y[d$ndx]),
                                      predict(modres_full, type="linear")))))
    
    pred_withheld <- try(predict(modres_full, xnew=d$xnew, dumnew=NULL,
                                 type="linear"), silent=TRUE)
    if (!inherits(pred_withheld, "try-error")) {
      out <- rbind(out, data.frame(
        data=d$tmp_data_name, material=d$tmp_material, model=model_name,
        error="withheld fit", error_measure=c("RMSE","Cor","MAE"),
        value=as.numeric(eval_fit_bench(log(d$y[-d$ndx]), pred_withheld))))
    }
    
    # diagnostic plots
    for (plot_type in c("scatter","residuals","res-QQ","res-hist","res-acf")) {
      tmp_plot <- try(plot(modres_full, type=plot_type) + theme_classic(), silent=TRUE)
      if (!inherits(tmp_plot, "try-error")) {
        plot_dir  <- if (plot_type == "scatter") "./plots_rev_deg" else "./plots_residuals_rev_deg"
        tag       <- switch(plot_type,
                            "scatter"   = "Fitted",
                            "residuals" = "Residuals",
                            "res-QQ"   = "ResQQ",
                            "res-hist" = "ResHist",
                            "res-acf"  = "ResACF")
        ggsave(paste0(plot_dir, "/", tag, "_", model_name, "_",
                      d$tmp_data_name, ".pdf"),
               tmp_plot, width=8*0.6, height=5*0.6)
      }
    }
  }
  out
}

clusterExport(my.cluster, c(".prep_dataset", ".fit_one"), envir=environment())

# ── define parametric benchmark models (run in parallel) ────────────────────────────────
bench_models_par <- list(
  "WienerDeg" = function(y, x, dum, ndx, jump_direction)
    WienerDeg(y, x, dum, ndx, jump_direction=jump_direction,
              methods=c("nvm","lbfgsb3c")),
  
  "GammaDeg"  = function(y, x, dum, ndx, jump_direction)
    GammaDeg(y, x, dum, ndx, jump_direction=jump_direction,
             methods=c("nvm","lbfgsb3c")),
  
  "EKF_JD"    = function(y, x, dum, ndx, jump_direction)
    EKF_JD(y, x, dum, ndx, jump_direction=jump_direction)
)
clusterExport(my.cluster, "bench_models_par", envir=environment())

# ── PHASE 1: parallel loop for parametric models ─────────────────────────────────────────
cat("Phase 1: fitting WienerDeg, GammaDeg, EKF_JD in parallel...\n")
foreach(i = 1:length(datanames)) %dopar% {
  d <- .prep_dataset(i)
  model_comp_bench <- data.frame(data=NULL, material=NULL, model=NULL,
                                 error=NULL, error_measure=NULL, value=NULL)
  for (k in seq_along(bench_models_par)) {
    model_comp_bench <- rbind(model_comp_bench,
                              .fit_one(bench_models_par[[k]], names(bench_models_par)[k], d))
  }
  cat(sprintf("Finished parametric rep %d / %d at %s.\n", i, length(datanames), Sys.time()))
  saveRDS(model_comp_bench,
          paste0("./saved_results_loop/bench_par_results_i", i, ".rds"))
  model_comp_bench
}

parallel::stopCluster(my.cluster)
cat("Phase 1 complete.\n\n")

# ── PHASE 2: sequential loop for SparseGP (memory-heavy, not safe in PSOCK) ─────────────
cat("Phase 2: fitting SparseGP sequentially...\n")
for (i in 1:length(datanames)) {
  d <- .prep_dataset(i)
  gp_res <- .fit_one(
    function(y, x, dum, ndx, jump_direction)
      SparseGP(y, x, dum, ndx, jump_direction=jump_direction, m_local=50L),
    "SparseGP", d
  )
  cat(sprintf("Finished SparseGP rep %d / %d at %s.\n", i, length(datanames), Sys.time()))
  saveRDS(gp_res,
          paste0("./saved_results_loop/bench_gp_results_i", i, ".rds"))
}

# ── aggregate results (parametric + SparseGP) ───────────────────────────────────────────
cat("\nAggregating results...\n")
model_comp_bench_all <- data.frame()
for (i in seq_along(datanames)) {
  resi_par <- readRDS(paste0("./saved_results_loop/bench_par_results_i", i, ".rds"))
  resi_gp  <- readRDS(paste0("./saved_results_loop/bench_gp_results_i",  i, ".rds"))
  model_comp_bench_all <- rbind(model_comp_bench_all, resi_par, resi_gp)
}

# ── merge with RFNR results for joint comparison table ──────────────────────────────────
rfnr_results   <- readRDS("./saved_results/results_RFNR_final.rds")
model_comp_rfnr <- rfnr_results$model_comp %>%
  filter(model == "Full")   # compare against the best RFNR variant only

model_comp_all <- bind_rows(
  model_comp_rfnr %>% mutate(family = "RFNR"),
  model_comp_bench_all %>% mutate(family = "Benchmark")
)

saveRDS(model_comp_all, "./saved_results/results_benchmark_final.rds")

# ── quick summary table (RMSE by model × error type) ────────────────────────────────────
summary_table <- model_comp_all %>%
  filter(error_measure == "RMSE") %>%
  group_by(model, error) %>%
  summarise(mean_RMSE = mean(value, na.rm = TRUE),
            sd_RMSE   = sd(value,   na.rm = TRUE),
            .groups = "drop") %>%
  arrange(error, mean_RMSE)

print(summary_table)

# ── comparison plot ──────────────────────────────────────────────────────────────────────
p_comparison <- model_comp_all %>%
  filter(error_measure == "RMSE", error %in% c("train/test", "sample fit")) %>%
  mutate(error = factor(error, levels = c("sample fit", "train/test"))) %>%
  ggplot(aes(x = reorder(model, value), y = value, fill = family)) +
  geom_boxplot(outlier.size = 0.8) +
  facet_wrap(~error, scales = "free_y") +
  labs(x = NULL, y = "RMSE (log scale)", fill = "Model family",
       title = "Model comparison: RFNR Full vs. Benchmarks") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("./plots_rev_deg/ModelComparison_RMSE.pdf", p_comparison,
       width = 10, height = 5)

cat("\nDone. Results saved to ./saved_results/results_benchmark_final.rds\n")
cat("Comparison plot saved to ./plots_rev_deg/ModelComparison_RMSE.pdf\n")
