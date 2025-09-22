
pacman::p_load(dplyr,tidyr,ggplot2,knitr,ggrepel,kableExtra)
source("./RFNR_model_functions.R")

# final version: 
res <- readRDS("./saved_results/results_RFNR_final.rds") # best at the moment
datanames <- res$datanames

# latest updates:
# - sigma <- nearestPD(sigma) in HAC cov summary
# - jump at start of RH15 cases
# - less downsampling for RH100 cases
# - also consider kkts and conv in selection of solver method
# - rho_upper changed to 1

# # # # model comparison
comp_tab <- res$model_comp  %>% filter(error_measure!="Cor") %>%
  pivot_wider(names_from = c("material","error","error_measure"),values_from = "value")

colnames(comp_tab) <- c("Data","Model",rep(c("RMSE","MAE"),6))
comp_kab <- kable(comp_tab[,-1],booktabs=TRUE,digits = 3,format = "latex") %>%
  add_header_above(header = c(" ","Test error"=2,"In-sample error"=2,"Interpol. error"=2,
                              "Test error"=2,"In-sample error"=2,"Interpol. error"=2)) %>%
  add_header_above(header = c(" ","MEA1 material"=6,"MEA2 material"=6))
for (i in 1:13) {
  comp_kab <- comp_kab %>% group_rows(res$datanames[i],start_row = 1 + 4*(i-1),4+4*(i-1))
}
comp_kab

# just one error type for one data set ; do not show interpolating error on withheld data
onedata_tab <- res$model_comp %>% filter(data=="Ambient air s1 t1",material=="MEA1",error!="withheld fit",error_measure!="Cor") %>% 
  pivot_wider(names_from = c("error","error_measure"),values_from = "value")

colnames(onedata_tab) <- c("Data","Material","Model",rep(c("RMSE","MAE"),2))
kable(onedata_tab[-c(1,2)],booktabs=TRUE,format="latex", digits=3) %>%
  add_header_above(header = c(" ","Test error"=2,"In-sample error"=2)) 


# # no norm tests in paper, since we use HAC cov
# str(res$sum_list$`Ambient air s1 t1`$norm_tests)
# norm_tests <- matrix(NA,4,3)
# colnames(norm_tests) <- c("test","statistic","p-value")
# for (i in 1:4) {
#   norm_tests[i,1] <- res$sum_list$`Ambient air s1 t1`$norm_tests[[i]]$method[1]
#   norm_tests[i,2] <- format(res$sum_list$`Ambient air s1 t1`$norm_tests[[i]]$statistic,digits=3)
#   norm_tests[i,3] <- res$sum_list$`Ambient air s1 t1`$norm_tests[[i]]$p.value
# }
# norm_tests[,3] <- format.pval(as.numeric(norm_tests[,3]) ,digits=3)
# kable(norm_tests,booktabs=TRUE,format="latex")


# CI plots of parameters
# first single plots,
# then one plot for all parameters and materials used in paper
# kable_res_tab

parameters <- c("c1","lam1","c2","lam2","delta","rho")
df_parCIs <- data.frame(NULL)
# i <- 3
for (i in 1:6) {
  df_parCIs <- rbind(df_parCIs,
                     data.frame(data=c(datanames[1:13],datanames[1:13]),material=c(rep("MEA1",13),rep("MEA2",13)),
                                parameter=parameters[i],
                                estimate=res$model_fits[,4+3*(i-1)],
                                lower=res$model_fits[,5+3*(i-1)],
                                upper=res$model_fits[,6+3*(i-1)]))
}
df_parCIs$data <- factor(df_parCIs$data,levels=datanames[13:1])
df_parCIs$parameter <- factor(df_parCIs$parameter,levels=c("c1","lam1","c2","lam2","delta","rho"),
                              labels = c(expression(c[1]),expression(lambda[1]),
                                         expression(c[2]),expression(lambda[2]),
                                         expression(delta),expression(rho)))


# new corrections pos pars:

# 1: RH15 high variance, practically Inf bounds
# 1.1: RH15  MEA2 for c1,lam1,c2,lam2:
df_parCIs[26,c(5,6)] <- NA
# data material parameter lower   upper
# RH15     MEA2 c[1]  2.467107e-14 1900556525
df_parCIs[26*2,c(5,6)] <- NA
# data material parameter   lower      upper
# RH15     MEA2 lambda[1] 3.688737e-99 5.664976e+86
df_parCIs[26*3,c(5,6)] <- NA
# data material parameter  lower  upper
# RH15     MEA2      c[2] 4.60199e-08 49040355
df_parCIs[26*4,c(5,6)] <- NA
# data material parameter   lower  upper
# RH15     MEA2 lambda[2] 8.44502e-11 12.25943

# 1.2: RH15  MEA1 for all (were all Inf/1 or 0)
df_parCIs[13,c(5,6)] <- NA
df_parCIs[13+26,c(5,6)] <- NA
df_parCIs[13+26*2,c(5,6)] <- NA
df_parCIs[13+26*3,c(5,6)] <- NA
df_parCIs[13+26*4,c(5,6)] <- NA
df_parCIs[13+26*5,c(5,6)] <- NA

# 1.3 RH15 MEA1 lam1 estimate too small, remove
df_parCIs[13+26*3,c(4)] <- NA # 8.367918e-12

# 2.1: PM filter MEA2 lam1 too big, set to 0, Inf
df_parCIs[8+13+26*1,c(5,6)] <- c(0,Inf) # 6.603397e-30 2.74845e+15

# 2.2: PM filter MEA1 lam1 too big, set to 0, Inf
df_parCIs[8+26*1,c(5,6)] <- c(0,Inf) # 2.018608e-11 0.0225249

# 3: Amb air s2 t1 MEA2 lam1 too big, set to Inf
df_parCIs[5+13+26*1,c(5,6)] <- c(0,Inf) # 1.36465e-18 5765.606

# 4: Amb air s2 t2 MEA2 lam1 too big, set to Inf
df_parCIs[6+13+26*1,c(5,6)] <- c(0,Inf) # 2.235798e-13 0.1700911


df_parCIs %>% 
  ggplot(aes(x=estimate,y=data,col=material)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper),alpha=0.6) +
  facet_wrap(.~parameter,nrow=3,scales="free_x",
             labeller = label_parsed) +
  scale_x_log10() +
  labs(y=" ",x="95% CIs")
# ggsave(paste0("./plots_rev_deg/CIs_all_pars.pdf"),width = 10,height = 12)


# # # # # # plots of fits 
my_full_dat <- res$my_full_dat
datanames <- res$datanames

tmp_plot <- my_full_dat %>%
  filter(!startsWith(data,"RH"),material=="MEA1") %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(data,levels=datanames)) +
  labs(y=expression(current~density~(A/cm^2)),x="time (s)")
tmp_plot
# ggsave(paste0("./plots_rev_deg/ExpBaseline_noRH_original.pdf"),tmp_plot,width = 10,height = 6)

tmp_plot <- my_full_dat %>%
  filter(!startsWith(data,"RH"),material!="MEA1") %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(data,levels=datanames)) +
  labs(y=expression(current~density~(A/cm^2)),x="time (s)")
tmp_plot
# ggsave(paste0("./plots_rev_deg/ExpBaseline_noRH_MEA2.pdf"),tmp_plot,width = 10,height = 6)


tmp_plot <- my_full_dat %>% 
  # mutate(material=ifelse(material=="original","MEA1","MEA2")) %>%
  filter(startsWith(data,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,linetype=type,col=factor(data,levels=c(datanames[10:13])))) +
  geom_line(alpha=0.8) +
  # facet_wrap(.~factor(material,levels=c("MEA1","MEA2"))) +
  facet_wrap(.~material) +
  labs(y=expression(current~density~(A/cm^2)),x="time (s)",col="rel. humidity")
tmp_plot
# ggsave(paste0("./plots_rev_deg/ExpBaseline_RH.pdf"),tmp_plot,width = 10,height = 4)


# # # # # plots of non-parametric degradation curves

# # # as function of time
deg_curves <- data.frame(NULL)
# log necessary? # no
coefs <- res$model_fits[,c(4,7,10,13)]
colnames(coefs) <- c("c1","lam1","c2","lam2")

for (i in 1:nrow(res$model_fits)) {
  tmp_name <- c(res$datanames[1:13],res$datanames[1:13])[i]
  tmp_mat <- c(rep("MEA1",13),rep("MEA2",13))[i]
  tmp_t <- res$my_full_dat$time[res$my_full_dat$data==tmp_name & res$my_full_dat$material==tmp_mat]
  type <- c(rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4),
    rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4))[i]
  
  rcd <- exp(-coefs[i,2]*tmp_t - coefs[i,3]*(1-exp(-coefs[i,4]*tmp_t)))
  deg_curves <- rbind(deg_curves,
                      data.frame(data=tmp_name,material = tmp_mat,
                                 time=tmp_t,type=type,
                                 rcd=rcd))
}

deg_curves <- deg_curves %>% mutate(isrh = startsWith(data,"RH"), type2 = ifelse(isrh,data,type)) 
deg_curves %>% filter(time>1e3, 
                      !(data=="RH15")
                      ) %>%
  ggplot(aes(x=1+time,y=rcd,col=factor(type2,levels=c("ambient","filter","pure","RH100","RH50","RH30","RH15")),
             group=data)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  labs(y=expression(predicted~current~density/pcd[0]),col="type",x="time (s)") +
  facet_grid(material~isrh,scales="free_x") +
  theme(strip.text.x = element_blank()) +
  geom_text(data=deg_curves %>% filter(time>1e5,data=="ALL filters s2",time>106001, time < 106002),
            size=2,nudge_y = -0.05,nudge_x = 0.05,
            aes(label=data),show.legend = FALSE)
# ggsave(paste0("./plots_rev_deg/RelDegCurves_all.pdf"),width = 10,height = 5)



# # # significance of 70% deg time and of rel deg after 5 hours (18000 sec)
# not used, too wide (up to Inf) Confidence intervals
tmp_alpha <- 0.05
M <- 1e5
# require(mvtnorm)

est_t <- est_rd <- numeric(26)
quant_t <- quant_rd <- matrix(c(0),26,2)
t_rd <- 0.5
rd_time <- 5*3600
set.seed(1234)
for (i in 1:nrow(res$model_fits)) {
  releg5h <- timesto70 <- numeric(M)
  se <- res$sum_list[[i]]$summary$par_sig[,2]
  sigma <- res$sum_list[[i]]$summary$par_cor * outer(se,se)
  # if (!isSymmetric(sigma)) {
  #   sigma <- (sigma + t(sigma)) / 2
  # }
  if (any(is.nan(sigma))) {
    next
  } else {
    sigma <- Matrix::nearPD(sigma)$mat
  }
  # sigma <- diag(6) * outer(se,se)
  coefj <- MASS::mvrnorm(M,mu = res$sum_list[[i]]$summary$par_sig[,1],Sigma=sigma)
  for (j in 1:M) {
    uniroot_res <- try(uniroot(function(tmp_t) {(-exp(coefj[j,2])*tmp_t - exp(coefj[j,3])*(1-exp(-exp(coefj[j,4])*tmp_t)))-log(t_rd)},
                               interval = c(0,1e3),extendInt = "downX",tol = 1e-2),
                       silent=TRUE)
    if (inherits(uniroot_res,what = "try-error")) {
      uniroot_res$root <- Inf
    }
    timesto70[j] <- uniroot_res$root
    releg5h[j] <- exp(-exp(coefj[j,2])*rd_time - exp(coefj[j,3])*(1-exp(-exp(coefj[j,4])*rd_time)))
  }

  uniroot_res <- try(uniroot(function(tmp_t) {(-exp(res$sum_list[[i]]$summary$par_sig[2,1])*tmp_t - exp(res$sum_list[[i]]$summary$par_sig[3,1])*(1-exp(-exp(res$sum_list[[i]]$summary$par_sig[4,1])*tmp_t)))-log(t_rd)},
                             interval = c(0,1e3),extendInt = "downX",tol = 1e-2))
  if (inherits(uniroot_res,what = "try-error")) {
    uniroot_res$root <- Inf
  }
  est_t[i] <- uniroot_res$root
  quant_t[i,] <- quantile(timesto70,probs = c(tmp_alpha/2,1-tmp_alpha/2))

  est_rd[i] <- exp(-exp(res$sum_list[[i]]$summary$par_sig[2,1])*rd_time - exp(res$sum_list[[i]]$summary$par_sig[3,1])*(1-exp(-exp(res$sum_list[[i]]$summary$par_sig[4,1])*rd_time)))
  quant_rd[i,] <- quantile(releg5h,probs = c(tmp_alpha/2,1-tmp_alpha/2),na.rm = TRUE)

  print(quantile(timesto70,probs=c(0.95,0.975,0.995)))
  print(quantile(releg5h,probs=c(0.95,0.975,0.995),na.rm=TRUE))
  message(sprintf("Finished rep %d / %d!",i,26))
}


df_degCIs <- data.frame(NULL)
df_degCIs <- rbind(df_degCIs,
                   data.frame(data=c(datanames[1:13],datanames[1:13]),material=c(rep("MEA1",13),rep("MEA2",13)),
                              term="Time until 50% of pcd0",
                              estimate=est_t,
                              lower=quant_t[,1],
                              upper=quant_t[,2]))

df_degCIs <- rbind(df_degCIs,
                   data.frame(data=c(datanames[1:13],datanames[1:13]),material=c(rep("MEA1",13),rep("MEA2",13)),
                              term="Rel. deg. after 5 hours",
                              estimate=est_rd,
                              lower=quant_rd[,1],
                              upper=quant_rd[,2]))


df_degCIs$data <- factor(df_degCIs$data,levels=datanames[13:1])
df_degCIs$term <- factor(df_degCIs$term,levels=c("Time until 50% of pcd0","Rel. deg. after 5 hours"),
                              labels = c(expression(Time~until~"50%"~of~pcd[0]),expression(Rel.~deg.~after~5~hours)))


df_degCIs %>% filter(data!="RH15") %>%
  ggplot(aes(x=estimate,y=data,col=material)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper),alpha=0.6) +
  facet_wrap(.~term,nrow=1,scales="free_x",
             labeller = label_parsed) +
  labs(y=" ",x="95% CIs",col="material")
# ggsave(paste0("./plots_rev_deg/CIs_rel_deg.pdf"),width = 10,height = 4)
