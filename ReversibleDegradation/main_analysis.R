# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # MAIN ANALYSIS SCRIPT — RFNR + Benchmark Models # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if ("pacman" %in% rownames(installed.packages()) == FALSE) install.packages("pacman")
pacman::p_load(foreach, parallel, doParallel,
               dplyr, tidyr, ggplot2, readxl,
               optimx, laGP, FKF,
               goftest, sandwich, Matrix, forecast)

source("./RFNR_help_functions.R")
source("./RFNR_model_functions.R")
source("./benchmark_model_functions.R")
source("./auto_tune_final.R")

# ── output folders ───────────────────────────────────────────────────────────────────────
dir.create("./saved_results",            showWarnings = FALSE)
dir.create("./saved_results_loop",       showWarnings = FALSE)
dir.create("./plots_rev_deg",            showWarnings = FALSE)
dir.create("./plots_residuals_rev_deg",  showWarnings = FALSE)

# ── dataset names ────────────────────────────────────────────────────────────────────────
datanames <- c("Pure air s1","Pure air s2",
               "Ambient air s1 t1","Ambient air s1 t2",
               "Ambient air s2 t1","Ambient air s2 t2",
               "Ambient air s2 t3",
               "PM filter s2","ALL filters s2",
               paste0("RH", c(100, 50, 30, 15)))
datanames <- c(datanames, paste0(datanames, " MEA2"))
numdat    <- length(datanames)

# ── RFNR model variants ──────────────────────────────────────────────────────────────────
models <- list(
  "Decay only" = function(y, x, dum, ndx, jump_direction = NULL)
    RFNR(y, x, dum, ndx, log_order=1, jump_direction="none", init=c(-1,-10)),
  
  "Decay curvature" = function(y, x, dum, ndx, jump_direction = NULL)
    RFNR(y, x, dum, ndx, log_order=2, jump_direction="none", init=c(-1,-10,-1,-10)),
  
  "Decay shocks" = function(y, x, dum, ndx, jump_direction)
    RFNR(y, x, dum, ndx, log_order=1, jump_direction=jump_direction, init=c(-1,-10,0,1)),
  
  "Full" = function(y, x, dum, ndx, jump_direction)
    RFNR(y, x, dum, ndx, log_order=2, jump_direction=jump_direction, init=c(-1,-10,-1,-10,0,1))
)

# ── load data (needed for automatic tuning) ──────────────────────────────────────────────
cat("Loading all datasets...\n")
all_log_y <- vector("list", numdat)
for (i in 1:numdat) {
  dat            <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx", sheet = i)
  colnames(dat)  <- c("time", "current")
  dat$time       <- (dat$time - 1) * 60
  all_log_y[[i]] <- log(dat$current)
  cat(sprintf("  sheet %2d OK\n", i))
}
cat("Done.\n\n")

# ── automatic tuning: quantile-based threshold ───────────────────────────────────────────
# tau = Q_p(|nabla log y|), p = 0.987
# h = 100 for all datasets; h = 8 for RH100 (i=10,23), tau set manually.
# p = 0.987 selected by minimising total |jump_diff| vs expert-validated counts (24 datasets).
tuning       <- setup_tuning(numdat, p = 0.987)
hs           <- tuning$hs
grad_cutoffs <- tuning$grad_cutoffs

# ═══════════════════════════════════════════════════════════════════════════════════════════
# PHASE 1 — RFNR (parallel)
# ═══════════════════════════════════════════════════════════════════════════════════════════
unlink("./saved_results/log.txt")
n.cores    <- parallel::detectCores()
my.cluster <- parallel::makeCluster(n.cores - 1, type = "PSOCK",
                                    outfile = "./saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)

clusterExport(my.cluster,
              c("datanames", "models", "mygrad", "hs", "grad_cutoffs"),
              envir = environment())
clusterEvalQ(my.cluster, {
  pacman::p_load(foreach, parallel, dplyr, tidyr, ggplot2, readxl,
                 optimx, goftest, sandwich, Matrix, forecast)
  source("./RFNR_help_functions.R")
  source("./RFNR_model_functions.R")
  source("./auto_tune_final.R")
})

i_done_rfnr <- as.integer(gsub(".*RFNR_results_i(\\d+)\\.rds", "\\1",
                               list.files("./saved_results_loop", pattern = "^RFNR_results_i")))
i_todo_rfnr <- setdiff(1:numdat, i_done_rfnr)
if (length(i_todo_rfnr) == 0) {
  cat("Phase 1: all RFNR results already present, skipping.\n")
} else {
  cat(sprintf("Phase 1: fitting RFNR for %d datasets: %s\n",
              length(i_todo_rfnr), paste(i_todo_rfnr, collapse = ", ")))
}

foreach(i = i_todo_rfnr) %dopar% {
  data <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx", sheet = i)
  colnames(data) <- c("time", "current")
  data$time <- (data$time - 1) * 60
  jump_direction <- if (i == 10) "down" else "up"
  
  grad        <- mygrad(log(data$current), hs[i])
  grad_cutoff <- grad_cutoffs[i]
  if (i == 10) grad <- -grad
  ind <- !is.na(grad) & grad > grad_cutoff
  for (k in seq_len(length(ind) - 1))
    if (!is.na(ind[k] * ind[k+1]) && ind[k] && ind[k+1]) ind[k] <- FALSE
  if (i %in% c(10, 13, 26)) ind[1] <- TRUE
  
  x   <- data$time
  y   <- data$current
  ndx <- if (i %in% c(10, 23)) seq(1, length(x), 12*4) else seq(1, length(x), 60*4)
  xnew     <- x[-ndx]
  time_ind <- 1 + (1:length(x))[-ndx] * (length(ndx)-1) / max(ndx)
  
  dum <- numeric(length(ndx))
  for (k in which(ind)) dum[which.min(abs(k - ndx))] <- 1
  
  ind_train <- 1:round(0.7 * length(ndx))
  ind_test  <- (round(0.7 * length(ndx)) + 1):length(ndx)
  
  tmp_data_name <- if (i < 14) datanames[i] else datanames[i - 13]
  tmp_material  <- if (i < 14) "MEA1" else "MEA2"
  
  model_comp <- data.frame(data=NULL, material=NULL, model=NULL,
                           error=NULL, error_measure=NULL, value=NULL)
  for (k in 1:length(models)) {
    modres <- try(models[[k]](y, x, dum[ind_train], ndx[ind_train], jump_direction))
    if (!inherits(modres, "try-error")) {
      model_comp <- rbind(model_comp,
                          data.frame(data=tmp_data_name, material=tmp_material, model=names(models)[k],
                                     error="train/test", error_measure=c("RMSE","Cor","MAE"),
                                     value=as.numeric(eval_fit(log(y[ndx[ind_test]]),
                                                               predict(modres, xnew=x[ndx[ind_test]],
                                                                       dumnew=dum[ind_test], type="linear")))))
    }
    modres <- try(models[[k]](y, x, dum, ndx, jump_direction))
    if (!inherits(modres, "try-error")) {
      model_comp <- rbind(model_comp,
                          data.frame(data=tmp_data_name, material=tmp_material, model=names(models)[k],
                                     error="sample fit", error_measure=c("RMSE","Cor","MAE"),
                                     value=as.numeric(eval_fit(log(y[ndx]), predict(modres, type="linear")))))
      model_comp <- rbind(model_comp,
                          data.frame(data=tmp_data_name, material=tmp_material, model=names(models)[k],
                                     error="withheld fit", error_measure=c("RMSE","Cor","MAE"),
                                     value=as.numeric(eval_fit(log(y[-ndx]),
                                                               predict(modres, xnew=xnew, time_ind=time_ind, type="linear")))))
    }
  }
  
  data_tmp <- data %>% mutate(id = 1:length(data$time)) %>% filter(id %in% ndx)
  
  tmp_plot <- try({
    plot(modres, "scatter") + aes(linetype=type) + labs(x="time (s)") + theme_classic()
  }, silent = TRUE)
  if (!inherits(tmp_plot, "try-error"))
    ggsave(paste0("./plots_rev_deg/Fitted_RFNR_", datanames[i], ".pdf"),
           tmp_plot, width=8*0.6, height=5*0.6)
  
  norm_tests <- list(
    "ad"      = goftest::ad.test(modres$log_residuals, "pnorm",
                                 mean=0, sd=modres$sigma, estimated=TRUE),
    "cvm"     = goftest::cvm.test(modres$log_residuals, null="pnorm",
                                  mean=0, sd=modres$sigma, estimated=TRUE),
    "shapiro" = shapiro.test(modres$log_residuals),
    "ks"      = ks.test(modres$log_residuals, "pnorm", mean=0, sd=modres$sigma)
  )
  
  model_fits <- numeric(21)
  names(model_fits) <- c("jumps","mean_time_jumps","rh",
                         "c1","lower","upper",
                         "lam1","lower","upper",
                         "c2","lower","upper",
                         "lam2","lower","upper",
                         "delta","lower","upper",
                         "rho","lower","upper")
  summod  <- summary(modres, corr=TRUE)
  alpha   <- 0.05
  myqnorm <- qnorm(alpha/2, lower.tail=FALSE)
  model_fits[1] <- length(which(ind))
  model_fits[2] <- mean(diff(data$time[which(ind)]))
  model_fits[3] <- c(rep(30,9),100,50,30,15,rep(30,9),100,50,30,15)[i]
  model_fits[4:6]   <- exp(summod$par_sig[1,1] + c(0,-1,1)*summod$par_sig[1,2]*myqnorm)
  model_fits[7:9]   <- exp(summod$par_sig[2,1] + c(0,-1,1)*summod$par_sig[2,2]*myqnorm)
  model_fits[10:12] <- exp(summod$par_sig[3,1] + c(0,-1,1)*summod$par_sig[3,2]*myqnorm)
  model_fits[13:15] <- exp(summod$par_sig[4,1] + c(0,-1,1)*summod$par_sig[4,2]*myqnorm)
  model_fits[16:18] <- exp(summod$par_sig[5,1] + c(0,-1,1)*summod$par_sig[5,2]*myqnorm)
  model_fits[19:21] <- modres$rho_upper /
    (1 + exp(-summod$par_sig[6,1] + c(0,1,-1)*summod$par_sig[6,2]*myqnorm))
  
  my_full_dat <- data.frame(data_tmp, data=tmp_data_name, material=tmp_material,
                            fit=as.numeric(modres$fitted_jumps))
  
  resi <- list(dataname    = datanames[i],
               model_fits  = model_fits,
               sum_list    = list(summary=summod, norm_tests=norm_tests),
               my_full_dat = my_full_dat,
               model_comp  = model_comp)
  cat(sprintf("Finished RFNR rep %d / %d at %s.\n", i, numdat, Sys.time()))
  saveRDS(resi, paste0("./saved_results_loop/RFNR_results_i", i, ".rds"))
  resi
}

parallel::stopCluster(my.cluster)
cat("Phase 1 (RFNR) complete.\n\n")

# ── aggregate RFNR results ────────────────────────────────────────────────────────────────
model_fits <- matrix(NA, numdat, 21)
row.names(model_fits) <- datanames
colnames(model_fits)  <- c("jumps","mean_time_jumps","rh",
                           "c1","lower","upper",
                           "lam1","lower","upper",
                           "c2","lower","upper",
                           "lam2","lower","upper",
                           "delta","lower","upper",
                           "rho","lower","upper")
sum_list    <- vector("list", numdat)
names(sum_list) <- datanames
model_comp  <- data.frame(data=NULL, model=NULL, error=NULL, error_measure=NULL, value=NULL)
my_full_dat <- data.frame(NULL)

for (i in 1:numdat) {
  resi           <- readRDS(paste0("./saved_results_loop/RFNR_results_i", i, ".rds"))
  model_fits[i,] <- resi$model_fits
  sum_list[[i]]  <- resi$sum_list
  my_full_dat    <- rbind(my_full_dat, resi$my_full_dat)
  model_comp     <- rbind(model_comp,  resi$model_comp)
}

saveRDS(list(model_fits  = model_fits,
             sum_list    = sum_list,
             my_full_dat = my_full_dat,
             datanames   = datanames,
             model_comp  = model_comp),
        "./saved_results/results_RFNR_final.rds")
cat("RFNR results saved.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════════════════
# PHASE 2 — BENCHMARK MODELS (parallel: Wiener/Gamma/EKF; sequential: SparseGP)
# ═══════════════════════════════════════════════════════════════════════════════════════════
unlink("./saved_results/log_bench.txt")
my.cluster <- parallel::makeCluster(n.cores - 1, type = "PSOCK",
                                    outfile = "./saved_results/log_bench.txt")
doParallel::registerDoParallel(cl = my.cluster)

clusterExport(my.cluster,
              c("datanames","mygrad","hs","grad_cutoffs",
                "make_step_mat","eval_fit_bench",
                "WienerDeg","GammaDeg","EKF_JD","SparseGP",
                "predict.WienerDeg","predict.GammaDeg",
                "predict.EKF_JD","predict.SparseGP",
                "coef.WienerDeg","coef.GammaDeg",
                "coef.EKF_JD","coef.SparseGP",
                "residuals.WienerDeg","residuals.GammaDeg",
                "residuals.EKF_JD","residuals.SparseGP",
                ".bench_plot",".bench_summary"),
              envir = environment())
clusterEvalQ(my.cluster, {
  pacman::p_load(foreach, parallel, dplyr, tidyr, ggplot2, readxl,
                 optimx, laGP, FKF, forecast)
  source("./RFNR_help_functions.R")
  source("./RFNR_model_functions.R")
  source("./benchmark_model_functions.R")
  source("./auto_tune_final.R")
})

.prep_dataset <- function(i) {
  data <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx", sheet = i)
  colnames(data) <- c("time", "current")
  data$time <- (data$time - 1) * 60
  jump_direction <- if (i == 10) "down" else "up"
  grad        <- mygrad(log(data$current), hs[i])
  grad_cutoff <- grad_cutoffs[i]
  if (i == 10) grad <- -grad
  ind <- !is.na(grad) & grad > grad_cutoff
  for (k in seq_len(length(ind) - 1))
    if (!is.na(ind[k]*ind[k+1]) && ind[k] && ind[k+1]) ind[k] <- FALSE
  if (i %in% c(10, 13, 26)) ind[1] <- TRUE
  x   <- data$time
  y   <- data$current
  ndx <- if (i %in% c(10, 23)) seq(1, length(x), 12*4) else seq(1, length(x), 60*4)
  xnew <- x[-ndx]
  dum  <- numeric(length(ndx))
  for (k in which(ind)) dum[which.min(abs(k - ndx))] <- 1
  ind_train     <- 1:round(0.7 * length(ndx))
  ind_test      <- (round(0.7 * length(ndx)) + 1):length(ndx)
  tmp_data_name <- if (i < 14) datanames[i] else datanames[i - 13]
  tmp_material  <- if (i < 14) "MEA1" else "MEA2"
  list(x=x, y=y, ndx=ndx, xnew=xnew, dum=dum,
       ind_train=ind_train, ind_test=ind_test,
       jump_direction=jump_direction,
       tmp_data_name=tmp_data_name, tmp_material=tmp_material)
}

.fit_one <- function(model_fn, model_name, d) {
  out <- data.frame(data=NULL, material=NULL, model=NULL,
                    error=NULL, error_measure=NULL, value=NULL)
  modres_train <- try(model_fn(d$y, d$x, d$dum[d$ind_train],
                               d$ndx[d$ind_train], d$jump_direction), silent=TRUE)
  if (!inherits(modres_train, "try-error")) {
    pred_test <- try(predict(modres_train, xnew=d$x[d$ndx[d$ind_test]],
                             dumnew=d$dum[d$ind_test], type="linear"), silent=TRUE)
    if (!inherits(pred_test, "try-error"))
      out <- rbind(out, data.frame(
        data=d$tmp_data_name, material=d$tmp_material, model=model_name,
        error="train/test", error_measure=c("RMSE","Cor","MAE"),
        value=as.numeric(eval_fit_bench(log(d$y[d$ndx[d$ind_test]]), pred_test))))
  }
  modres_full <- try(model_fn(d$y, d$x, d$dum, d$ndx, d$jump_direction), silent=TRUE)
  if (!inherits(modres_full, "try-error")) {
    out <- rbind(out, data.frame(
      data=d$tmp_data_name, material=d$tmp_material, model=model_name,
      error="sample fit", error_measure=c("RMSE","Cor","MAE"),
      value=as.numeric(eval_fit_bench(log(d$y[d$ndx]),
                                      predict(modres_full, type="linear")))))
    pred_withheld <- try(predict(modres_full, xnew=d$xnew, dumnew=NULL,
                                 type="linear"), silent=TRUE)
    if (!inherits(pred_withheld, "try-error"))
      out <- rbind(out, data.frame(
        data=d$tmp_data_name, material=d$tmp_material, model=model_name,
        error="withheld fit", error_measure=c("RMSE","Cor","MAE"),
        value=as.numeric(eval_fit_bench(log(d$y[-d$ndx]), pred_withheld))))
    tmp_plot <- try(plot(modres_full, type="scatter") + theme_classic(), silent=TRUE)
    if (!inherits(tmp_plot, "try-error"))
      ggsave(paste0("./plots_rev_deg/Fitted_", model_name, "_", d$tmp_data_name, ".pdf"),
             tmp_plot, width=8*0.6, height=5*0.6)
  }
  out
}

clusterExport(my.cluster, c(".prep_dataset", ".fit_one"), envir=environment())

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

i_done_par <- as.integer(gsub(".*bench_par_results_i(\\d+)\\.rds", "\\1",
                              list.files("./saved_results_loop", pattern="^bench_par_results_i")))
i_todo_par <- setdiff(1:numdat, i_done_par)
if (length(i_todo_par) == 0) {
  cat("Phase 2a: all parametric benchmark results already present, skipping.\n")
} else {
  cat(sprintf("Phase 2a: fitting WienerDeg, GammaDeg, EKF_JD for %d datasets: %s\n",
              length(i_todo_par), paste(i_todo_par, collapse=", ")))
}

foreach(i = i_todo_par) %dopar% {
  d <- .prep_dataset(i)
  model_comp_bench <- data.frame(data=NULL, material=NULL, model=NULL,
                                 error=NULL, error_measure=NULL, value=NULL)
  for (k in seq_along(bench_models_par))
    model_comp_bench <- rbind(model_comp_bench,
                              .fit_one(bench_models_par[[k]], names(bench_models_par)[k], d))
  cat(sprintf("Finished parametric rep %d / %d at %s.\n", i, numdat, Sys.time()))
  saveRDS(model_comp_bench,
          paste0("./saved_results_loop/bench_par_results_i", i, ".rds"))
  model_comp_bench
}

parallel::stopCluster(my.cluster)
cat("Phase 2a complete.\n\n")

# ── Phase 2b: SparseGP (sequential) ─────────────────────────────────────────────────────
i_done_gp <- as.integer(gsub(".*bench_gp_results_i(\\d+)\\.rds", "\\1",
                             list.files("./saved_results_loop", pattern="^bench_gp_results_i")))
i_todo_gp <- setdiff(1:numdat, i_done_gp)
if (length(i_todo_gp) == 0) {
  cat("Phase 2b: all SparseGP results already present, skipping.\n")
} else {
  cat(sprintf("Phase 2b: fitting SparseGP for %d datasets: %s\n",
              length(i_todo_gp), paste(i_todo_gp, collapse=", ")))
}

for (i in i_todo_gp) {
  d      <- .prep_dataset(i)
  gp_res <- .fit_one(
    function(y, x, dum, ndx, jump_direction)
      SparseGP(y, x, dum, ndx, jump_direction=jump_direction, m_local=50L),
    "SparseGP", d)
  cat(sprintf("Finished SparseGP rep %d / %d at %s.\n", i, numdat, Sys.time()))
  saveRDS(gp_res, paste0("./saved_results_loop/bench_gp_results_i", i, ".rds"))
}
cat("Phase 2b complete.\n\n")

# ═══════════════════════════════════════════════════════════════════════════════════════════
# PHASE 3 — AGGREGATE AND COMPARE
# ═══════════════════════════════════════════════════════════════════════════════════════════
cat("Aggregating all results...\n")
model_comp_bench_all <- data.frame()
for (i in seq_along(datanames)) {
  resi_par <- readRDS(paste0("./saved_results_loop/bench_par_results_i", i, ".rds"))
  resi_gp  <- readRDS(paste0("./saved_results_loop/bench_gp_results_i",  i, ".rds"))
  model_comp_bench_all <- rbind(model_comp_bench_all, resi_par, resi_gp)
}

rfnr_results    <- readRDS("./saved_results/results_RFNR_final.rds")
model_comp_rfnr <- rfnr_results$model_comp %>% filter(model == "Full")

model_comp_all <- bind_rows(
  model_comp_rfnr      %>% mutate(family = "RFNR"),
  model_comp_bench_all %>% mutate(family = "Benchmark")
)

saveRDS(model_comp_all, "./saved_results/results_benchmark_final.rds")

# ── summary table ────────────────────────────────────────────────────────────────────────
summary_table <- model_comp_all %>%
  filter(error_measure == "RMSE") %>%
  group_by(model, error) %>%
  summarise(mean_RMSE = mean(value, na.rm=TRUE),
            sd_RMSE   = sd(value,   na.rm=TRUE),
            .groups   = "drop") %>%
  arrange(error, mean_RMSE)
print(summary_table)

# ── RFNR variant comparison table ───────────────────────────────────────────────────────
model_comp_all <- model_comp_all %>%
  mutate(model = ifelse(model == "Full", "RFNR", model))

rfnr_tab <- rfnr_results$model_comp %>%
  filter(error_measure %in% c("RMSE","MAE"),
         error %in% c("train/test","sample fit","withheld fit"),
         is.finite(value)) %>%
  group_by(model, error, error_measure) %>%
  summarise(mean_val = round(mean(value, na.rm=TRUE), 3),
            sd_val   = round(sd(value,   na.rm=TRUE), 3),
            .groups  = "drop") %>%
  mutate(val = sprintf("%.3f (%.3f)", mean_val, sd_val)) %>%
  dplyr::select(model, error, error_measure, val) %>%
  pivot_wider(names_from=c("error","error_measure"), values_from=val) %>%
  arrange(`train/test_RMSE`)
print(rfnr_tab)

# ── comparison plots ─────────────────────────────────────────────────────────────────────
p_comparison <- model_comp_all %>%
  filter(error_measure == "RMSE", error %in% c("train/test","sample fit")) %>%
  mutate(error = factor(error, levels=c("sample fit","train/test"),
                        labels=c("In-sample","Out-of-sample"))) %>%
  ggplot(aes(x=reorder(model, value), y=value, fill=family)) +
  geom_boxplot(outlier.size=0.8) +
  facet_wrap(~error, scales="free_y") +
  labs(x=NULL, y="RMSE (log scale)", fill="Model family") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave("./plots_rev_deg/ModelComparison_RMSE.pdf", p_comparison, width=10, height=5)

p_comparison <- model_comp_all %>%
  filter(error_measure == "MAE", error %in% c("train/test","sample fit")) %>%
  mutate(error = factor(error, levels=c("sample fit","train/test"),
                        labels=c("In-sample","Out-of-sample"))) %>%
  ggplot(aes(x=reorder(model, value), y=value, fill=family)) +
  geom_boxplot(outlier.size=0.8) +
  facet_wrap(~error, scales="free_y") +
  labs(x=NULL, y="MAE (log scale)", fill="Model family") +
  theme_classic() +
  theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave("./plots_rev_deg/ModelComparison_MAE.pdf", p_comparison, width=10, height=5)

p_comparison <- model_comp_all %>%
  filter(error_measure == "RMSE",
         error %in% c("train/test","sample fit","withheld fit"),
         is.finite(value)) %>%
  mutate(error = factor(error,
                        levels=c("sample fit","train/test","withheld fit"),
                        labels=c("In-sample","Out-of-sample","Interpolation"))) %>%
  ggplot(aes(x=reorder(model, value), y=value, fill=family)) +
  geom_boxplot(outlier.size=0.8) +
  facet_wrap(~error, scales="free_y", nrow=1) +
  labs(x=NULL, y="RMSE (log scale)", fill="Model family") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=30, hjust=1))
ggsave("./plots_rev_deg/ModelComparison_RMSE2.pdf", p_comparison, width=12, height=5)

cat("\nDone.\n")
cat("RFNR results:      ./saved_results/results_RFNR_final.rds\n")
cat("Benchmark results: ./saved_results/results_benchmark_final.rds\n")
cat("Comparison plots:  ./plots_rev_deg/ModelComparison_RMSE.pdf\n")

sink("./saved_results/sessionInfo.txt")
print(sessionInfo())
sink()