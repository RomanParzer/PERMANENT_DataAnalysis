pacman::p_load(R.matlab,dplyr,tidyr,ggplot2,mgcv,robustbase)

# read in data
data <- readMat("./PERMANENT_CCMG_REFERENCE2.mat")
exp_current <- ifelse(c(data$Voltage2) < 0.775,
                      mean(data$Current2[data$Voltage2 < 0.775]),
                      mean(data$Current2[data$Voltage2 >= 0.775]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = c(data$Current2),
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))


cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean)

# we model current as funciton of voltage and exponential degradation of time, time after cleaning and time in cycle
# current = C*f(voltage) * exp(beta_1*time_c + beta_2*time_clean + beta_3*time)
# betas will be negative; product form due to exponential decay, want to apply log transformation, log(C) intercept


# # model funciton of current; takes a long time
tstamp <- Sys.time()
gam_res <- mgcv::gam(log(current) ~ s(voltage) + time_c + time_clean + time,data=df_long)
time_gam <- as.numeric(Sys.time() - tstamp,units="secs")

# use indicator of voltage
tstamp <- Sys.time()
lm_res <- lm(log(current) ~ I(voltage<0.775) + time_c + time_clean + time,data=df_long)
time_lm <- as.numeric(Sys.time() - tstamp,units="secs")
# indicator might not be enough; R^2 = 0.9691 quite good; exact confidence bounds, all significant (too many data?)

# try to model relative residual curve (current - exp_current) / exp_current ;
# needs robustness, try LTS
# could we need actual expected value?
# df_long %>% filter(cycle==242) %>%
#   ggplot(aes(x=time_c,y=current/ exp_current - 1)) +
#   geom_line(alpha = 0.2) +
#   labs(x="time in cycle in s")

# compare non-robust result to robust LTS estimation
tstamp <- Sys.time()
reslm_res <- lm(I(log(current/ exp_current )) ~  time_c + time_clean + time,data=df_long)
time_reslm <- as.numeric(Sys.time() - tstamp,units="secs")

# also takes quite some time (too much)
# tstamp <- Sys.time()
# lts_res <- robustbase::ltsReg(I(log(current/ exp_current )) ~  time_c + time_clean + time,data=df_long,model=FALSE)
# time_lts <- as.numeric(Sys.time() - tstamp,units="secs")
time_lts <- NA
lts_res <- NULL

tstamp <- Sys.time()
lmrob_res <- robustbase::lmrob(I(log(current/ exp_current )) ~  time_c + time_clean + time,data=df_long,model=FALSE)
time_lmrob <- as.numeric(Sys.time() - tstamp,units="secs")

saveRDS(list(times = c(as.numeric(time_gam,units="secs"),time_lm,time_reslm,time_lts,time_lmrob),
             gam_res = gam_res,
             lm_res = lm_res,
             reslm_res = reslm_res,
             lts_res = lts_res,
             lmrob_res=lmrob_res),
        "./saved_results/model_res_all.rds")


# tstamp <- Sys.time()
# gamb_res <- mgcv::gam(current ~ s(voltage) + s(time_c) + s(time_clean) + s(time),data=df_long)
# time_gamb <- Sys.time() - tstamp
# 
# saveRDS(list(times = c(time_gamb),
#              gamb_res = gamb_res),
#         "./saved_results/model_gamb.rds")


# tstamp <- Sys.time()
# gamb_res <- mgcv::gam(log(current) ~ s(voltage) + s(time_c) + s(time_clean) + s(time),data=df_long)
# time_gamb <- Sys.time() - tstamp
# 
# saveRDS(list(times = c(time_gamb),
#              gamb_res = gamb_res),
#         "./saved_results/model_gamb_log.rds")

resgamb_log <- readRDS( "./saved_results/model_gamb_log.rds")

overview_list <- list(lm = list(time = time_lm,print = print(lm_res), 
                                summary = summary(lm_res), fitted = lm_res$fitted.values),
                      lm_res = list(time = time_reslm,print = print(reslm_res), 
                                    summary = summary(reslm_res), fitted = reslm_res$fitted.values),
                      lmRob_res = list(time = time_lmrob,print = print(lmrob_res), 
                                       summary = summary(lmrob_res), fitted = lmrob_res$fitted.values),
                      gam_volt = list(time = time_gam,print = print(gam_res), 
                                      summary = summary(gam_res), fitted = gam_res$fitted.values),
                      gam_all = list(time = resgamb_log$times[1],print = print(resgamb_log$gamb_res), 
                                     summary = summary(resgamb_log$gamb_res), fitted = resgamb_log$gamb_res$fitted.values))

saveRDS(overview_list,"./saved_results/all_res_overview.rds")

# add plot
# str(overview_list_small)

overview_list_small <- list(lm = list(time = time_lm, coefs = summary(lm_res)$coefficients,
                                     r2 = summary(lm_res)$adj.r.squared, fitted = lm_res$fitted.values),
                            lm_res = list(time = time_reslm,coefs = summary(reslm_res)$coefficients,
                                          r2 = summary(reslm_res)$adj.r.squared, fitted = reslm_res$fitted.values),
                            lmRob_res = list(time = time_lmrob,coefs = summary(lmrob_res)$coefficients,
                                             r2 = summary(lmrob_res)$adj.r.squared,fitted = lmrob_res$fitted.values),
                            gam_volt = list(time = time_gam,coefs = summary(gam_res)[1:4],
                                            r2 = summary(gam_res)$r.sq,
                                            fitted = gam_res$fitted.values,plot = plot(gam_res)),
                            gam_all = list(time = resgamb_log$times[1],
                                           coefs = summary(resgamb_log$gamb_res)$s.table,
                                           r2 = summary(resgamb_log$gamb_res)$r.sq,
                                           fitted = resgamb_log$gamb_res$fitted.values,
                                           plot = plot(resgamb_log$gamb_res)))

saveRDS(overview_list_small,"./saved_results/all_res_overview_sm.rds")

