# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # RESULTS — RFNR + Benchmark Models # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

pacman::p_load(dplyr, tidyr, ggplot2, knitr, ggrepel, kableExtra, MASS, Matrix)
source("./RFNR_help_functions.R")
source("./RFNR_model_functions.R")

# ── load results ─────────────────────────────────────────────────────────────────────────
res       <- readRDS("./saved_results/results_RFNR_final.rds")
res_bench <- readRDS("./saved_results/results_benchmark_final.rds")
datanames <- res$datanames

# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 1 — MODEL COMPARISON TABLES
# ═══════════════════════════════════════════════════════════════════════════════════════════

# ── full RFNR model comparison (all 4 variants, both materials) ───────────────────────────
comp_tab <- res$model_comp %>%
  filter(error_measure != "Cor") %>%
  pivot_wider(names_from=c("material","error","error_measure"), values_from="value")

colnames(comp_tab) <- c("Data","Model", rep(c("RMSE","MAE"), 6))
comp_kab <- kable(comp_tab[,-1], booktabs=TRUE, digits=3, format="latex") %>%
  add_header_above(header=c(" ","Test error"=2,"In-sample error"=2,"Interpol. error"=2,
                            "Test error"=2,"In-sample error"=2,"Interpol. error"=2)) %>%
  add_header_above(header=c(" ","MEA1 material"=6,"MEA2 material"=6))
for (i in 1:13)
  comp_kab <- comp_kab %>% group_rows(res$datanames[i], start_row=1+4*(i-1), 4+4*(i-1))
comp_kab

# ── single dataset table (Ambient air s1 t1, MEA1, no withheld error) ────────────────────
onedata_tab <- res$model_comp %>%
  filter(data=="Ambient air s1 t1", material=="MEA1",
         error != "withheld fit", error_measure != "Cor") %>%
  pivot_wider(names_from=c("error","error_measure"), values_from="value")

colnames(onedata_tab) <- c("Data","Material","Model", rep(c("RMSE","MAE"), 2))
kable(onedata_tab[-c(1,2)], booktabs=TRUE, format="latex", digits=3) %>%
  add_header_above(header=c(" ","Test error"=2,"In-sample error"=2))

# ── benchmark summary table (RFNR Full vs benchmarks) ────────────────────────────────────
bench_tab <- res_bench %>%
  filter(error_measure=="RMSE", error %in% c("train/test","sample fit")) %>%
  group_by(model, error) %>%
  summarise(mean_RMSE = round(mean(value, na.rm=TRUE), 4),
            sd_RMSE   = round(sd(value,   na.rm=TRUE), 4),
            .groups="drop") %>%
  arrange(error, mean_RMSE)

kable(bench_tab, booktabs=TRUE, format="latex", digits=4,
      col.names=c("Model","Error type","Mean RMSE","SD RMSE")) %>%
  collapse_rows(columns=2, latex_hline="major")

# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 2 — PARAMETER CI PLOTS
# ═══════════════════════════════════════════════════════════════════════════════════════════

parameters <- c("c1","lam1","c2","lam2","delta","rho")
df_parCIs  <- data.frame(NULL)

for (i in 1:6) {
  df_parCIs <- rbind(df_parCIs,
                     data.frame(data      = c(datanames[1:13], datanames[1:13]),
                                material  = c(rep("MEA1",13),  rep("MEA2",13)),
                                parameter = parameters[i],
                                estimate  = res$model_fits[, 4+3*(i-1)],
                                lower     = res$model_fits[, 5+3*(i-1)],
                                upper     = res$model_fits[, 6+3*(i-1)]))
}
df_parCIs$data      <- factor(df_parCIs$data, levels=datanames[13:1])
df_parCIs$parameter <- factor(df_parCIs$parameter,
                              levels=c("c1","lam1","c2","lam2","delta","rho"),
                              labels=c(expression(c[1]),  expression(lambda[1]),
                                       expression(c[2]),  expression(lambda[2]),
                                       expression(delta), expression(rho)))

# ── manual corrections for degenerate CIs ────────────────────────────────────────────────
# NOTE: re-check these after refitting with auto_tune_final.R — values may have changed.

# 1.1: RH15 MEA2 — CI essentially (0, Inf) for c1, lam1, c2, lam2
df_parCIs[26,    c(5,6)] <- NA
df_parCIs[26*2,  c(5,6)] <- NA
df_parCIs[26*3,  c(5,6)] <- NA
df_parCIs[26*4,  c(5,6)] <- NA

# 1.2: RH15 MEA1 — all CIs degenerate
df_parCIs[13,      c(5,6)] <- NA
df_parCIs[13+26,   c(5,6)] <- NA
df_parCIs[13+26*2, c(5,6)] <- NA
df_parCIs[13+26*3, c(5,6)] <- NA
df_parCIs[13+26*4, c(5,6)] <- NA
df_parCIs[13+26*5, c(5,6)] <- NA

# 1.3: RH15 MEA1 lam1 estimate numerically zero — remove
df_parCIs[13+26*3, 4] <- NA

# 2.1: PM filter MEA2 lam1
df_parCIs[8+13+26*1, c(5,6)] <- c(0, Inf)

# 2.2: PM filter MEA1 lam1
df_parCIs[8+26*1, c(5,6)] <- c(0, Inf)

# 3: Amb air s2 t1 MEA2 lam1
df_parCIs[5+13+26*1, c(5,6)] <- c(0, Inf)

# 4: Amb air s2 t2 MEA2 lam1
df_parCIs[6+13+26*1, c(5,6)] <- c(0, Inf)

# ── plot ──────────────────────────────────────────────────────────────────────────────────
df_parCIs %>%
  ggplot(aes(x=estimate, y=data, col=material)) +
  geom_point(position=position_nudge(y=ifelse(df_parCIs$material=="MEA1", 0.2, -0.2))) +
  geom_errorbar(aes(xmin=lower, xmax=upper, linetype=material),
                position=position_nudge(y=ifelse(df_parCIs$material=="MEA1", 0.2, -0.2))) +
  facet_wrap(.~parameter, nrow=3, scales="free_x", labeller=label_parsed) +
  scale_x_log10() +
  theme_bw() +
  labs(y=" ", x=" ") +
  scale_color_brewer(type="qual", direction=1, palette=1)
# ggsave("./plots_rev_deg/CIs_all_pars.pdf", width=10, height=12)

# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 3 — FITTED CURVES PLOTS
# ═══════════════════════════════════════════════════════════════════════════════════════════

my_full_dat <- res$my_full_dat

# ── non-RH, MEA1 ─────────────────────────────────────────────────────────────────────────
my_full_dat %>%
  filter(!startsWith(data,"RH"), material=="MEA1") %>%
  pivot_longer(c(fit,current), names_to="type", values_to="value") %>%
  ggplot(aes(x=time, y=value, color=type, linetype=type)) +
  geom_line() +
  facet_wrap(.~factor(data, levels=datanames)) +
  labs(y=expression(current~density~(A/cm^2)), x="time (s)") +
  theme_bw() +
  scale_color_brewer(type="qual", direction=1, palette=1)
# ggsave("./plots_rev_deg/ExpBaseline_noRH_original.pdf", width=10, height=6)

# ── non-RH, MEA2 ─────────────────────────────────────────────────────────────────────────
my_full_dat %>%
  filter(!startsWith(data,"RH"), material!="MEA1") %>%
  pivot_longer(c(fit,current), names_to="type", values_to="value") %>%
  ggplot(aes(x=time, y=value, color=type, linetype=type)) +
  geom_line() +
  facet_wrap(.~factor(data, levels=datanames)) +
  labs(y=expression(current~density~(A/cm^2)), x="time (s)") +
  theme_bw() +
  scale_color_brewer(type="qual", direction=1, palette=1)
# ggsave("./plots_rev_deg/ExpBaseline_noRH_MEA2.pdf", width=10, height=6)

# ── RH datasets, both materials ───────────────────────────────────────────────────────────
my_full_dat %>%
  filter(startsWith(data,"RH")) %>%
  pivot_longer(c(fit,current), names_to="type", values_to="value") %>%
  ggplot(aes(x=time, y=value, linetype=type,
             col=factor(data, levels=datanames[10:13]))) +
  geom_line() +
  facet_wrap(.~material) +
  labs(y=expression(current~density~(A/cm^2)), x="time (s)", col="rel. humidity") +
  theme_bw() +
  scale_color_brewer(type="seq", direction=-1, palette=1)
# ggsave("./plots_rev_deg/ExpBaseline_RH.pdf", width=10, height=4)

# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 4 — RELATIVE DEGRADATION CURVES
# ═══════════════════════════════════════════════════════════════════════════════════════════

deg_curves <- data.frame(NULL)
coefs      <- res$model_fits[, c(4,7,10,13)]
colnames(coefs) <- c("c1","lam1","c2","lam2")

for (i in 1:nrow(res$model_fits)) {
  tmp_name <- c(res$datanames[1:13], res$datanames[1:13])[i]
  tmp_mat  <- c(rep("MEA1",13), rep("MEA2",13))[i]
  tmp_t    <- res$my_full_dat$time[res$my_full_dat$data==tmp_name &
                                     res$my_full_dat$material==tmp_mat]
  type     <- c(rep("pure",2), rep("ambient",5), rep("filter",2), rep("rh",4),
                rep("pure",2), rep("ambient",5), rep("filter",2), rep("rh",4))[i]
  rcd      <- exp(-coefs[i,2]*tmp_t - coefs[i,3]*(1-exp(-coefs[i,4]*tmp_t)))
  deg_curves <- rbind(deg_curves,
                      data.frame(data=tmp_name, material=tmp_mat, time=tmp_t, type=type, rcd=rcd))
}

col_fac_labels        <- c("Ambient/filtered/pure air","Relative humidity")
names(col_fac_labels) <- c(FALSE, TRUE)

deg_curves <- deg_curves %>%
  mutate(isrh  = startsWith(data,"RH"),
         type2 = ifelse(isrh, data, type))
deg_curves$type2 <- factor(deg_curves$type2,
                           levels=c("ambient","filter","pure","RH100","RH50","RH30","RH15"))

deg_curves %>%
  filter(time > 1e3, data != "RH15") %>%
  ggplot(aes(x=1+time, y=rcd, col=type2, group=data)) +
  geom_line(alpha=0.6, linewidth=1) +
  scale_x_log10() +
  labs(y=expression(predicted~current~density/pcd[0]),
       col="type", linetype="type", x="time (s)") +
  facet_grid(material~isrh, scales="free_x",
             labeller=labeller(material=c("MEA1","MEA2"), isrh=col_fac_labels)) +
  geom_text(data=deg_curves %>%
              filter(time>1e5, data=="ALL filters s2", time>106001, time<106002),
            size=2, nudge_y=-0.05, nudge_x=0.05,
            aes(label=data), show.legend=FALSE) +
  theme_bw() +
  scale_linetype_manual(values=c(1:3,1:3)) +
  scale_color_manual(
    values=c(scale_color_brewer(type="seq",palette=1,direction= 1)$palette(3),
             scale_color_brewer(type="seq",palette=1,direction=-1)$palette(3)))
# ggsave("./plots_rev_deg/RelDegCurves_all.pdf", width=10, height=5)

# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 5 — BENCHMARK VS RFNR COMPARISON PLOT
# ═══════════════════════════════════════════════════════════════════════════════════════════

res_bench %>%
  filter(error_measure=="RMSE", error %in% c("train/test","sample fit")) %>%
  mutate(error = factor(error, levels=c("sample fit","train/test"))) %>%
  ggplot(aes(x=reorder(model, value), y=value, fill=family)) +
  geom_boxplot(outlier.size=0.8) +
  facet_wrap(~error, scales="free_y") +
  labs(x=NULL, y="RMSE (log scale)", fill="Model family",
       title="Model comparison: RFNR Full vs. Benchmarks") +
  theme_bw() +
  theme(axis.text.x=element_text(angle=30, hjust=1))
# ggsave("./plots_rev_deg/ModelComparison_RMSE.pdf", width=10, height=5)

# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 6 — CI ON DEGRADATION QUANTITIES (Monte Carlo, M = 100,000)
# ═══════════════════════════════════════════════════════════════════════════════════════════
# Computes via parametric bootstrap (MASS::mvrnorm from HAC asymptotic covariance):
#   - t_50: time until 50% of initial current density
#   - rd_5h: relative degradation after 5 hours

tmp_alpha <- 0.05
M <- 1e4 
# M         <- 1e5
t_rd      <- 0.5
rd_time   <- 5 * 3600
est_t     <- est_rd  <- numeric(26)
quant_t   <- quant_rd <- matrix(0, 26, 2)
set.seed(1234)

for (i in 1:nrow(res$model_fits)) {
  releg5h <- timesto70 <- numeric(M)
  se    <- res$sum_list[[i]]$summary$par_sig[, 2]
  sigma <- res$sum_list[[i]]$summary$par_cor * outer(se, se)
  
  if (any(is.nan(sigma))) {
    message(sprintf("Dataset %d: NaN in covariance, skipping.", i))
    next
  }
  sigma <- Matrix::nearPD(sigma)$mat
  coefj <- MASS::mvrnorm(M, mu=res$sum_list[[i]]$summary$par_sig[,1], Sigma=sigma)
  
  for (j in 1:M) {
    ur <- try(uniroot(
      function(tmp_t)
        (-exp(coefj[j,2])*tmp_t -
           exp(coefj[j,3])*(1-exp(-exp(coefj[j,4])*tmp_t))) - log(t_rd),
      interval=c(0,1e3), extendInt="downX", tol=1e-2), silent=TRUE)
    timesto70[j] <- if (inherits(ur,"try-error")) Inf else ur$root
    releg5h[j]   <- exp(-exp(coefj[j,2])*rd_time -
                          exp(coefj[j,3])*(1-exp(-exp(coefj[j,4])*rd_time)))
  }
  
  ur_est <- try(uniroot(
    function(tmp_t)
      (-exp(res$sum_list[[i]]$summary$par_sig[2,1])*tmp_t -
         exp(res$sum_list[[i]]$summary$par_sig[3,1]) *
         (1-exp(-exp(res$sum_list[[i]]$summary$par_sig[4,1])*tmp_t))) - log(t_rd),
    interval=c(0,1e3), extendInt="downX", tol=1e-2))
  
  est_t[i]    <- if (inherits(ur_est,"try-error")) Inf else ur_est$root
  quant_t[i,] <- quantile(timesto70, probs=c(tmp_alpha/2, 1-tmp_alpha/2))
  est_rd[i]   <- exp(
    -exp(res$sum_list[[i]]$summary$par_sig[2,1])*rd_time -
      exp(res$sum_list[[i]]$summary$par_sig[3,1]) *
      (1-exp(-exp(res$sum_list[[i]]$summary$par_sig[4,1])*rd_time)))
  quant_rd[i,] <- quantile(releg5h, probs=c(tmp_alpha/2, 1-tmp_alpha/2), na.rm=TRUE)
  
  message(sprintf("Finished rep %d / 26.", i))
}

df_degCIs <- bind_rows(
  data.frame(data=c(datanames[1:13],datanames[1:13]),
             material=c(rep("MEA1",13),rep("MEA2",13)),
             term="Time until 50% of pcd0",
             estimate=est_t, lower=quant_t[,1], upper=quant_t[,2]),
  data.frame(data=c(datanames[1:13],datanames[1:13]),
             material=c(rep("MEA1",13),rep("MEA2",13)),
             term="Rel. deg. after 5 hours",
             estimate=est_rd, lower=quant_rd[,1], upper=quant_rd[,2])
)
df_degCIs$data <- factor(df_degCIs$data, levels=datanames[13:1])
df_degCIs$term <- factor(df_degCIs$term,
                         levels=c("Time until 50% of pcd0","Rel. deg. after 5 hours"),
                         labels=c(expression(Time~until~"50%"~of~pcd[0]),
                                  expression(Rel.~deg.~after~5~hours)))

df_degCIs_plot <- df_degCIs %>% filter(data != "RH15")
df_degCIs_plot %>%
  ggplot(aes(x=estimate, y=data, col=material)) +
  geom_point(position=position_nudge(
    y=ifelse(df_degCIs_plot$material=="MEA1", 0.2, -0.2))) +
  geom_errorbar(aes(xmin=lower, xmax=upper, linetype=material), alpha=0.6,
                position=position_nudge(
                  y=ifelse(df_degCIs_plot$material=="MEA1", 0.2, -0.2))) +
  facet_wrap(.~term, nrow=1, scales="free_x", labeller=label_parsed) +
  labs(y=" ", x=" ", col="material") +
  theme_bw() +
  scale_color_brewer(type="qual", direction=1, palette=1)
# ggsave("./plots_rev_deg/CIs_rel_deg.pdf", width=10, height=4)



# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 7 — FITTED CURVES COMPARISON: RFNR vs BENCHMARKS (Ambient air s1 t1, MEA1)
# ═══════════════════════════════════════════════════════════════════════════════════════════

all_log_y <- readRDS("./data/all_log_y.rds")

i         <- 1
log_y     <- all_log_y[[i]]
y         <- exp(log_y)
x         <- seq(0, by=1/4, length.out=length(log_y))

grad        <- mygrad(log_y, hs[i])
grad_cutoff <- grad_cutoffs[i]
ind <- !is.na(grad) & grad > grad_cutoff
for (k in seq_len(length(ind)-1))
  if (!is.na(ind[k]*ind[k+1]) && ind[k] && ind[k+1]) ind[k] <- FALSE

ndx <- seq(1, length(x), 60*4)
dum <- numeric(length(ndx))
for (k in which(ind)) dum[which.min(abs(k - ndx))] <- 1

ind_train <- 1:round(0.7 * length(ndx))
ind_test  <- (round(0.7 * length(ndx)) + 1):length(ndx)
t_cutoff  <- x[ndx[max(ind_train)]]

# ── fit all models on training set ───────────────────────────────────────────────────────
fit_rfnr   <- RFNR(y, x, dum[ind_train], ndx[ind_train], log_order=2,
                   jump_direction="up", init=c(-1,-10,-1,-10,0,1))
fit_wiener <- WienerDeg(y, x, dum[ind_train], ndx[ind_train], jump_direction="up")
fit_gamma  <- GammaDeg(y, x, dum[ind_train], ndx[ind_train], jump_direction="up")
fit_ekf    <- EKF_JD(y, x, dum[ind_train], ndx[ind_train], jump_direction="up")
fit_gp     <- SparseGP(y, x, dum[ind_train], ndx[ind_train],
                       jump_direction="up", m_local=50L)

# ── collect in-sample predictions ────────────────────────────────────────────────────────
df_train <- data.frame(
  time    = x[ndx[ind_train]],
  current = y[ndx[ind_train]],
  RFNR      = exp(predict(fit_rfnr,   type="linear")),
  WienerDeg = exp(predict(fit_wiener, type="linear")),
  GammaDeg  = exp(predict(fit_gamma,  type="linear")),
  EKF_JD    = exp(predict(fit_ekf,    type="linear")),
  SparseGP  = exp(predict(fit_gp,     type="linear")),
  period    = "In-sample"
)

# ── collect out-of-sample predictions ────────────────────────────────────────────────────
df_test <- data.frame(
  time    = x[ndx[ind_test]],
  current = y[ndx[ind_test]],
  RFNR      = exp(predict(fit_rfnr,   xnew=x[ndx[ind_test]],
                          dumnew=dum[ind_test], type="linear")),
  WienerDeg = exp(predict(fit_wiener, xnew=x[ndx[ind_test]],
                          dumnew=dum[ind_test], type="linear")),
  GammaDeg  = exp(predict(fit_gamma,  xnew=x[ndx[ind_test]],
                          dumnew=dum[ind_test], type="linear")),
  EKF_JD    = exp(predict(fit_ekf,    xnew=x[ndx[ind_test]],
                          dumnew=dum[ind_test], type="linear")),
  SparseGP  = exp(predict(fit_gp,     xnew=x[ndx[ind_test]],
                          dumnew=dum[ind_test], type="linear")),
  period    = "Out-of-sample"
)

# ── combine and reshape ───────────────────────────────────────────────────────────────────
df_fits <- bind_rows(df_train, df_test) %>%
  pivot_longer(cols=c(RFNR, WienerDeg, GammaDeg, EKF_JD, SparseGP),
               names_to="model", values_to="fit") %>%
  mutate(model  = factor(model,
                         levels=c("RFNR","EKF_JD","SparseGP","WienerDeg","GammaDeg")),
         period = factor(period, levels=c("In-sample","Out-of-sample")))

# ── plot ──────────────────────────────────────────────────────────────────────────────────
p_fits <- df_fits %>%
  ggplot(aes(x=time)) +
  geom_line(aes(y=current), color="grey60", linewidth=0.4) +
  geom_line(aes(y=fit, linetype=period), color="steelblue", linewidth=0.8) +
  geom_vline(xintercept=t_cutoff, linetype="dotted",
             color="black", linewidth=0.5) +
  facet_wrap(~model, nrow=2) +
  scale_linetype_manual(values=c("In-sample"="solid",
                                 "Out-of-sample"="dashed")) +
  labs(y=expression(current~density~(A/cm^2)),
       x="time (s)",
       linetype=" ") +
  theme_bw() +
  theme(legend.position="bottom")
p_fits
ggsave("./plots_rev_deg/FittedComparison_i1.pdf", p_fits, width=10, height=6)


# ═══════════════════════════════════════════════════════════════════════════════════════════
# SECTION 8 — SUPPLEMENTARY TABLE 1: bandwidths h and thresholds tau
# ═══════════════════════════════════════════════════════════════════════════════════════════

# recompute tau for all datasets using auto_tune_final.R
all_log_y <- readRDS("./data/all_log_y.rds")

tau_table <- data.frame(
  dataset  = datanames,
  material = c(rep("MEA1", 13), rep("MEA2", 13)),
  h        = hs,
  tau_auto = NA_real_,
  tau_used = NA_real_,
  manual   = FALSE
)

for (i in 1:numdat) {
  log_y      <- all_log_y[[i]]
  jump_dir   <- if (i == 10) "down" else "up"
  tau_auto_i <- quantile_tau(log_y, hs[i], jump_dir, p=0.987)
  tau_table$tau_auto[i] <- round(tau_auto_i, 4)
  tau_table$tau_used[i] <- round(grad_cutoffs[i], 4)
  tau_table$manual[i]   <- i %in% c(10, 23)
}

# ── LaTeX table ───────────────────────────────────────────────────────────────────────────
tau_table_print <- tau_table %>%
  mutate(tau = ifelse(manual,
                      paste0(tau_used, "$^*$"),
                      as.character(tau_auto)),
         dataset = gsub(" MEA2$", "", dataset)) %>%
  dplyr::select(dataset, material, h, tau)

kable(tau_table_print,
      booktabs  = TRUE,
      format    = "latex",
      digits    = 4,
      escape    = FALSE,
      col.names = c("Dataset", "Material", "$h$", "$\\tau$"),
      caption   = "Bandwidths $h$ and jump detection thresholds $\\tau$ for all
                   26 datasets. $\\tau$ is the automatically computed 
                   quantile-based threshold ($p = 0.987$), except for the two 
                   RH100 datasets (marked with $^*$) where $\\tau$ was set 
                   manually due to their distinct signal structure.",
      label     = "") %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(general       = "* manually set threshold.",
           general_title = "",
           escape        = FALSE)


# ── Table 2 Supplementary: all 26 datasets, all 4 RFNR variants ──────────────────────────

comp_tab_rmse <- res$model_comp %>%
  filter(error_measure == "RMSE",
         error %in% c("train/test","sample fit","withheld fit")) %>%
  pivot_wider(names_from=c("material","error"), values_from="value") %>%
  rename(
    MEA1_test     = `MEA1_train/test`,
    MEA1_insample = `MEA1_sample fit`,
    MEA1_interp   = `MEA1_withheld fit`,
    MEA2_test     = `MEA2_train/test`,
    MEA2_insample = `MEA2_sample fit`,
    MEA2_interp   = `MEA2_withheld fit`
  ) %>%
  dplyr::select(data, model,
                MEA1_test, MEA1_insample, MEA1_interp,
                MEA2_test, MEA2_insample, MEA2_interp)

comp_tab_mae <- res$model_comp %>%
  filter(error_measure == "MAE",
         error %in% c("train/test","sample fit","withheld fit")) %>%
  pivot_wider(names_from=c("material","error"), values_from="value") %>%
  rename(
    MEA1_test     = `MEA1_train/test`,
    MEA1_insample = `MEA1_sample fit`,
    MEA1_interp   = `MEA1_withheld fit`,
    MEA2_test     = `MEA2_train/test`,
    MEA2_insample = `MEA2_sample fit`,
    MEA2_interp   = `MEA2_withheld fit`
  ) %>%
  dplyr::select(data, model,
                MEA1_test, MEA1_insample, MEA1_interp,
                MEA2_test, MEA2_insample, MEA2_interp)

comp_tab <- res$model_comp %>%
  filter(error_measure != "Cor") %>%
  pivot_wider(names_from=c("material","error","error_measure"), 
              values_from="value")

colnames(comp_tab) <- c("Data","Model", rep(c("RMSE","MAE"), 6))

comp_kab <- kable(comp_tab[,-1], booktabs=TRUE, digits=3, format="latex") %>%
  add_header_above(header=c(" ",
                            "Test error"=2,"In-sample error"=2,"Interpol. error"=2,
                            "Test error"=2,"In-sample error"=2,"Interpol. error"=2)) %>%
  add_header_above(header=c(" ","MEA1 material"=6,"MEA2 material"=6))

for (i in 1:13)
  comp_kab <- comp_kab %>%
  group_rows(res$datanames[i], start_row=1+4*(i-1), 4+4*(i-1))

comp_kab


# ── session info ──────────────────────────────────────────────────────────────────────────
sink("./saved_results/sessionInfo.txt")
print(sessionInfo())
sink()
cat("Session info saved to ./saved_results/sessionInfo.txt\n")