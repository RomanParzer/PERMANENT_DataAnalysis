
pacman::p_load(dplyr,tidyr,ggplot2,knitr,ggrepel,kableExtra)
source("./RFNR_model_functions.R")

res <- readRDS("../saved_results/results_RFNR_model_all_cor_CIs.rds")
datanames <- res$datanames

# # # # model comparison
comp_tab <- res$model_comp  %>% filter(error_measure!="Cor") %>%
  pivot_wider(names_from = c("material","error","error_measure"),values_from = "value")

colnames(comp_tab) <- c("Data","Model",rep(c("RMSE","MAE"),6))
comp_kab <- kable(comp_tab[,-1],booktabs=TRUE,digits = 3,format = "latex") %>%
  add_header_above(header = c(" ","Test error"=2,"In-sample error"=2,"Interpol. error"=2,
                              "Test error"=2,"In-sample error"=2,"Interpol. error"=2)) %>%
  add_header_above(header = c(" ","MEA1 material"=6,"MEA2 material"=6))
for (i in 1:13) {
  comp_kab <- comp_kab %>% group_rows(res$datanames[res$data_order[i]],start_row = 1 + 4*(i-1),4+4*(i-1))
}
comp_kab

# just one error type for one data set ; do not show interpolating error on withheld data
onedata_tab <- res$model_comp %>% filter(data=="Ambient air test 1",material=="original",error!="withheld fit",error_measure!="Cor") %>% 
  pivot_wider(names_from = c("error","error_measure"),values_from = "value")

colnames(onedata_tab) <- c("Data","Material","Model",rep(c("RMSE","MAE"),2))
kable(onedata_tab[-c(1,2)],booktabs=TRUE,format="latex", digits=3) %>%
  add_header_above(header = c(" ","Test error"=2,"In-sample error"=2)) 


str(res$sum_list$`Ambient air test 1`$norm_tests)

# only report KS and SW in text
# norm_tests <- matrix(NA,4,3)
# colnames(norm_tests) <- c("test","statistic","p-value")
# for (i in 1:4) {
#   norm_tests[i,1] <- res$sum_list$`Ambient air test 1`$norm_tests[[i]]$method[1]
#   norm_tests[i,2] <- format(res$sum_list$`Ambient air test 1`$norm_tests[[i]]$statistic,digits=3)
#   norm_tests[i,3] <- res$sum_list$`Ambient air test 1`$norm_tests[[i]]$p.value
# }
# norm_tests[,3] <- format.pval(as.numeric(norm_tests[,3]) ,digits=3)
# kable(norm_tests,booktabs=TRUE,format="latex")

# # # # table of parameter estimates and their significance

kable_res_tab <- matrix(" ",26,8)
row.names(kable_res_tab) <- c(row.names(res$model_fits)[1:13],row.names(res$model_fits)[1:13])
colnames(kable_res_tab) <- c("no. jumps","avg. time","c1 (99% CI)","lam1 (99% CI)","c2 (99% CI)","lam2 (99% CI)","delta (99% CI)","rho (99% CI)")

kable_res_tab[,c(1)] <- res$model_fits[,1]
kable_res_tab[,c(2)] <- round(res$model_fits[,c(2)],digits=1)

kable_res_tab[,c(3)] <- apply(format(res$model_fits[,c(4,5,6)],digits=3),1,function(tmp_row){
  paste0(tmp_row[1]," (",tmp_row[2],",",tmp_row[3],")")
})
kable_res_tab[,c(4)] <- apply(res$model_fits[,c(7:9)],1,function(tmp_row){
  paste0(format(tmp_row[1],digits=3)," (",format(tmp_row[2],digits=3),",",format(tmp_row[3],digits = 3),")")
})
kable_res_tab[,c(5)] <- apply(res$model_fits[,c(10:12)],1,function(tmp_row){
  paste0(format(tmp_row[1],digits=3)," (",format(tmp_row[2],digits=3),",",format(tmp_row[3],digits = 3),")")
})
kable_res_tab[,c(6)] <- apply(res$model_fits[,c(13:15)],1,function(tmp_row){
  paste0(format(tmp_row[1],digits=3)," (",format(tmp_row[2],digits=3),",",format(tmp_row[3],digits = 3),")")
})
kable_res_tab[,c(7)] <- apply(res$model_fits[,c(16:18)],1,function(tmp_row){
  paste0(format(tmp_row[1],digits=3)," (",format(tmp_row[2],digits=3),",",format(tmp_row[3],digits = 3),")")
})
kable_res_tab[,c(8)] <- apply(res$model_fits[,c(19:21)],1,function(tmp_row){
  paste0(format(tmp_row[1],digits=3)," (",format(tmp_row[2],digits=3),",",format(tmp_row[3],digits = 3),")")
})
# kable_res_tab
kable(kable_res_tab,booktabs=TRUE,format = "latex") %>%
  group_rows("MEA1 material",1,13) %>%
  group_rows("MEA2 material",14,26)

# CI plots of parameters

# kable_res_tab

parameters <- c("c1","lam1","c2","lam2","delta","rho")
df_parCIs <- data.frame(NULL)
name_order <- c(1,4,2,3,5,6,7,8,9,
                10,11,12,13)
for (i in 1:6) {
  df_parCIs <- rbind(df_parCIs,
                     data.frame(data=c(datanames[1:13],datanames[1:13]),material=c(rep("MEA1",13),rep("MEA2",13)),
                                parameter=parameters[i],
                                estimate=res$model_fits[,4+3*(i-1)],
                                lower=res$model_fits[,5+3*(i-1)],
                                upper=res$model_fits[,6+3*(i-1)]))
}
df_parCIs$data <- factor(df_parCIs$data,levels=datanames[name_order][13:1])

df_parCIs[10+13,c(5,6)] <- NA # was 0 Inf
df_parCIs %>% filter(parameter =="c1") %>%
  ggplot(aes(x=estimate,y=data)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper)) +
  facet_wrap(.~material,scales="free_x") +
  labs(y=" ",x=expression("99% CIs for c"[1]))
# ggsave(paste0("../plots_rev_deg/CIs_c1.pdf"),width = 8,height = 6)

res$sum_list$`RH100 MEA2`$summary$par_sig[1,1]
str(res$sum_list$`RH100 MEA2`$norm_tests)
res$sum_list$`RH100 MEA2`$summary$par_sig[1,1] - qt(0.005,719-6,lower.tail = FALSE)*res$sum_list$`RH100 MEA2`$summary$par_sig[1,2]
res$sum_list$`RH100 MEA2`$summary$par_sig[1,1] + qt(0.005,719-6,lower.tail = FALSE)*res$sum_list$`RH100 MEA2`$summary$par_sig[1,2]


df_parCIs[10+13+26,c(4,5,6)] <- NA # was 1e-4 -+=(1e-3)
df_parCIs %>% filter(parameter =="lam1") %>%
  ggplot(aes(x=estimate,y=data)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper)) +
  facet_wrap(.~material,scales="free_x") +
  labs(y=" ",x=expression("99% CIs for"~lambda[1]))
# ggsave(paste0("../plots_rev_deg/CIs_lam1.pdf"),width = 8,height = 6)

df_parCIs[10+13+26*2,c(4,5,6)] <- NA # was 150 -+=(1e3)
df_parCIs %>% filter(parameter =="c2") %>%
  ggplot(aes(x=estimate,y=data)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper)) +
  facet_wrap(.~material,scales="free_x") +
  labs(y=" ",x=expression("99% CIs for"~c[2]))
# ggsave(paste0("../plots_rev_deg/CIs_c2.pdf"),width = 8,height = 6)

df_parCIs %>% filter(parameter =="lam2") %>%
  ggplot(aes(x=estimate,y=data)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper)) +
  facet_wrap(.~material,scales="free_x") +
  labs(y=" ",x=expression("99% CIs for"~lambda[2]))
# ggsave(paste0("../plots_rev_deg/CIs_lam2.pdf"),width = 8,height = 6)

df_parCIs %>% filter(parameter =="delta") %>%
  ggplot(aes(x=estimate,y=data)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper)) +
  facet_wrap(.~material,scales="free_x") +
  labs(y=" ",x=expression("99% CIs for"~delta))
# ggsave(paste0("../plots_rev_deg/CIs_delta.pdf"),width = 8,height = 6)

df_parCIs %>% filter(parameter =="rho") %>%
  ggplot(aes(x=estimate,y=data)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper)) +
  facet_wrap(.~material,scales="free_x") +
  labs(y=" ",x=expression("99% CIs for"~rho))
# ggsave(paste0("../plots_rev_deg/CIs_rho.pdf"),width = 8,height = 6)

df_parCIs$parameter <- factor(df_parCIs$parameter,levels=c("c1","lam1","c2","lam2","delta","rho"),
                              labels = c(expression(c[1]),expression(lambda[1]),
                                         expression(c[2]),expression(lambda[2]),
                                         expression(delta),expression(rho)))


df_parCIs %>% 
  ggplot(aes(x=estimate,y=data,col=material)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper),alpha=0.6) +
  facet_wrap(.~parameter,nrow=3,scales="free_x",
             labeller = label_parsed) +
  labs(y=" ",x="99% CIs",col="Material")
ggsave(paste0("../plots_rev_deg/CIs_all_pars.pdf"),width = 10,height = 12)

# # # # # # plots of fits 
my_full_dat <- res$my_full_dat
datanames <- res$datanames
name_order <- c(1,4,2,3,5,6,7,8,9,
                10,11,12,13)
name_order <- c(name_order,13+name_order)

tmp_plot <- my_full_dat %>%
  filter(!startsWith(data,"RH"),material=="original") %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(data,levels=datanames[name_order])) +
  labs(y=expression(current~density~(A/cm^2)),x="time (s)")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_noRH_original.pdf"),tmp_plot,width = 10,height = 6)

tmp_plot <- my_full_dat %>%
  filter(!startsWith(data,"RH"),material!="original") %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(data,levels=datanames[name_order])) +
  labs(y=expression(current~density~(A/cm^2)),x="time (s)")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_noRH_MEA2.pdf"),tmp_plot,width = 10,height = 6)


tmp_plot <- my_full_dat %>% 
  mutate(material=ifelse(material=="original","MEA1","MEA2")) %>%
  filter(startsWith(data,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,linetype=type,col=factor(data,levels=c(datanames[10:13])))) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(material,levels=c("MEA1","MEA2"))) +
  labs(y=expression(current~density~(A/cm^2)),x="time (s)",col="rel. humidity")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_RH.pdf"),tmp_plot,width = 10,height = 4)


# # # # # plots of non-parametric degradation curves

# deg_times <- matrix(c(0),nrow(res$model_fits),9)
# row.names(deg_times) <- c(row.names(res$model_fits)[1:13],row.names(res$model_fits)[1:13])
# ratios <- 9:1/10
# colnames(deg_times) <- ratios
# 
# coefs <- res$model_fits[,c(4,7,10,13)]
# coefs[,1] <- log(coefs[,1])
# colnames(coefs) <- c("log_c1","lam1","c2","lam2")
# 
# for (i in 1:nrow(res$model_fits)) {
#   pred_start <- exp(lin_pred_RFNR(par = coefs[i,],x_use = 0,jump_direction = "none"))
#   for (k in 1:9) {
#     uniroot_res <- try(uniroot(function(tmp_t) {exp(lin_pred_RFNR(par = coefs[i,],x_use = tmp_t,jump_direction = "none"))-ratios[k]*pred_start},
#                            interval = c(0,1e3),extendInt = "downX",tol = 1e-2))
#     if (inherits(uniroot_res,what = "try-error")) {
#       uniroot_res$root <- Inf
#     }
#     deg_times[i,k] <- uniroot_res$root
#   }
# }
# 
# mydf <- data.frame(deg_times,type=c(rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4),
#                                     rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4)),
#                    material=c(rep("MEA1",13),rep("MEA2",13)))
# colnames(mydf)[1:9] <- ratios
# 
# mydf <- mydf %>% mutate(id=1:26,name = row.names(deg_times)) %>%
#   pivot_longer(c(1:9),names_to = "ratio",values_to = "time")
# mydf$ratio <- as.numeric(mydf$ratio)
# 
# 
# mydf <- mydf %>% mutate(isrh = startsWith(name,"RH"), type2 = ifelse(isrh,name,type)) 
# mydf %>%
#   ggplot(aes(x=time,y=ratio,col=factor(type2,levels=c("ambient","filter","pure","RH100","RH50","RH30","RH15")),
#              group=id)) +
#   geom_line(alpha=0.6,linewidth=1) +
#   scale_x_log10() +
#   # geom_text_repel(data = filter(mydf_all,ratio==0.5,!startsWith(name,"RH")),
#   #                 aes(x=time+1,y=ratio,label=name),
#   #                 show.legend = FALSE,
#   #                 size=2,
#   #                 max.overlaps = 20) +
#   labs(y=expression(current~density/pcd[0]),col="type",x="time (s)") +
#   facet_grid(material~isrh) +
#   theme(strip.text.x = element_blank())
# # ggsave(paste0("../plots_rev_deg/RelDegTimes_all.pdf"),width = 10,height = 5)


# # # same but as function of time
deg_curves <- data.frame(NULL)
coefs <- res$model_fits[,c(4,7,10,13)]
coefs[,1] <- log(coefs[,1])
colnames(coefs) <- c("log_c1","lam1","c2","lam2")

for (i in 1:nrow(res$model_fits)) {
  tmp_name <- c(res$datanames[1:13],res$datanames[1:13])[i]
  tmp_mat <- c(rep("original",13),rep("MEA 2",13))[i]
  tmp_t <- res$my_full_dat$time[res$my_full_dat$data==tmp_name & res$my_full_dat$material==tmp_mat]
  type <- c(rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4),
    rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4))[i]
  
  rcd <- exp(-coefs[i,2]*tmp_t - coefs[i,3]*(1-exp(-coefs[i,4]*tmp_t)))
  deg_curves <- rbind(deg_curves,
                      data.frame(data=tmp_name,material = ifelse(tmp_mat=="original","MEA1","MEA2"),
                                 time=tmp_t,type=type,
                                 rcd=rcd))
}

deg_curves <- deg_curves %>% mutate(isrh = startsWith(data,"RH"), type2 = ifelse(isrh,data,type)) 
deg_curves %>% filter(time>1e3) %>%
  ggplot(aes(x=1+time,y=rcd,col=factor(type2,levels=c("ambient","filter","pure","RH100","RH50","RH30","RH15")),
             group=data)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  labs(y=expression(current~density/pcd[0]),col="type",x="time (s)") +
  facet_grid(material~isrh,scales="free_x") +
  theme(strip.text.x = element_blank())
# ggsave(paste0("../plots_rev_deg/RelDegCurves_all.pdf"),width = 10,height = 5)



# # # significance of 70% deg time and of rel deg after 5 hours (18000 sec)



name_order <- c(1,4,2,3,5,6,7,8,9,
                10,11,12,13)
tmp_alpha <- 0.10
M <- 1e5
require(mvtnorm)

est_t70 <- est_rd5 <- numeric(26)
quant_t70 <- quant_rd5 <- matrix(c(0),26,2)
set.seed(1234)
for (i in 1:nrow(res$model_fits)) {
  releg5h <- timesto70 <- numeric(M)
  se <- res$sum_list[[i]]$summary$par_sig[,2]
  sigma <- res$sum_list[[i]]$summary$par_cor * outer(se,se)
  # sigma <- diag(6) * outer(se,se)
  coefj <- rmvt(M,sigma=sigma,delta = res$sum_list[[i]]$summary$par_sig[,1])
  for (j in 1:M) {
    uniroot_res <- try(uniroot(function(tmp_t) {(-coefj[j,2]*tmp_t - coefj[j,3]*(1-exp(-coefj[j,4]*tmp_t)))-log(0.7)},
                               interval = c(0,1e3),extendInt = "downX",tol = 1e-2),
                       silent=TRUE)
    if (inherits(uniroot_res,what = "try-error")) {
      uniroot_res$root <- Inf
    }
    timesto70[j] <- uniroot_res$root
    releg5h[j] <- exp(-coefj[j,2]*18000 - coefj[j,3]*(1-exp(-coefj[j,4]*18000)))
  }
  
  uniroot_res <- try(uniroot(function(tmp_t) {(-res$sum_list[[i]]$summary$par_sig[2,1]*tmp_t - res$sum_list[[i]]$summary$par_sig[3,1]*(1-exp(-res$sum_list[[i]]$summary$par_sig[4,1]*tmp_t)))-log(0.7)},
                             interval = c(0,1e3),extendInt = "downX",tol = 1e-2))
  if (inherits(uniroot_res,what = "try-error")) {
    uniroot_res$root <- Inf
  }
  est_t70[i] <- uniroot_res$root
  quant_t70[i,] <- quantile(timesto70,probs = c(tmp_alpha/2,1-tmp_alpha/2))
  
  est_rd5[i] <- exp(-res$sum_list[[i]]$summary$par_sig[2,1]*18000 - res$sum_list[[i]]$summary$par_sig[3,1]*(1-exp(-res$sum_list[[i]]$summary$par_sig[4,1]*18000)))
  quant_rd5[i,] <- quantile(releg5h,probs = c(tmp_alpha/2,1-tmp_alpha/2))
  
  print(quantile(timesto70,probs=c(0.95,0.975,0.995)))
  print(quantile(releg5h,probs=c(0.95,0.975,0.995)))
  message(sprintf("Finished rep %d / %d!",i,26))
}


# time to 70%
df_degCIs <- data.frame(NULL)
df_degCIs <- rbind(df_degCIs,
                   data.frame(data=c(datanames[1:13],datanames[1:13]),material=c(rep("MEA1",13),rep("MEA2",13)),
                              term="Time until 70% of pcd0",
                              estimate=est_t70,
                              lower=quant_t70[,1],
                              upper=quant_t70[,2]))

df_degCIs <- rbind(df_degCIs,
                   data.frame(data=c(datanames[1:13],datanames[1:13]),material=c(rep("MEA1",13),rep("MEA2",13)),
                              term="Rel. deg. after 5 hours",
                              estimate=est_rd5,
                              lower=quant_rd5[,1],
                              upper=quant_rd5[,2]))


df_degCIs$data <- factor(df_degCIs$data,levels=datanames[name_order][13:1])
df_degCIs$term <- factor(df_degCIs$term,levels=c("Time until 70% of pcd0","Rel. deg. after 5 hours"),
                              labels = c(expression(Time~until~"70%"~of~pcd[0]),expression(Rel.~deg.~after~5~hours)))

# RH100 MEA2
df_degCIs[10+13,c(5,6)] <- NA # both bounds higher than estimate
df_degCIs[10+13 + 26,c(5,6)] <- NA # both bounds higher than estimate

# RH15 MEA2 rel deg upper bound too high, set to Inf
df_degCIs[52,c(6)] <- Inf

# PM filter MEA2 upper bound too high, set to Inf
df_degCIs[c(21,47),c(6)] <- Inf

df_degCIs %>% 
  ggplot(aes(x=estimate,y=data,col=material)) +
  geom_point() +
  geom_errorbar(aes(xmin = lower,xmax=upper),alpha=0.6) +
  facet_wrap(.~term,nrow=1,scales="free_x",
             labeller = label_parsed) +
  labs(y=" ",x="90% CIs",col="Material")
ggsave(paste0("../plots_rev_deg/CIs_rel_deg.pdf"),width = 10,height = 4)
