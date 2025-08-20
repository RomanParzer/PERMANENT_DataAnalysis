
pacman::p_load(dplyr,tidyr,ggplot2,knitr,ggrepel)
source("./RFNR_model_functions.R")

res <- readRDS("../saved_results/results_RFNR_model_MEA2.rds")
res1 <- readRDS("../saved_results/results_RFNR_model.rds")

res_tab <- res$res_tab
kable_res_tab <- res_tab[,1:10]
kable_res_tab[,c(2)] <- round(res_tab[,c(2)],digits=1)
kable_res_tab[,c(4,7,9,10)] <- round(res_tab[,c(4,7,9,10)],digits=3)
kable_res_tab[,c(5,6,8)] <- format(res_tab[,c(5,6,8)],digits=3)
kable(kable_res_tab,booktabs=TRUE,format = "latex")

my_full_dat <- res$my_full_dat
datanames <- res$datanames

tmp_plot <- my_full_dat %>%
  filter(!startsWith(dataname,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(dataname,levels=datanames[name_order])) +
  labs(y="current",x="time in s")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_all_MEA2.pdf"),tmp_plot,width = 10,height = 6)


tmp_plot <- my_full_dat %>%
  filter(startsWith(dataname,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,linetype=type,col=factor(dataname,levels=c(datanames[9:12])))) +
  geom_line(alpha=0.8) +
  labs(y="current",x="time in s",col="rel. humidity")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_RH_MEA2.pdf"),tmp_plot,width = 8,height = 4)

tmp_plot <- rbind(data.frame(my_full_dat[,-3],material="MEA2"),
                  data.frame(res1$my_full_dat,material="orig")) %>%
  filter(startsWith(dataname,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,linetype=type,col=factor(dataname,levels=c(datanames[9:12])))) +
  geom_line(alpha=0.8) +
  facet_wrap(.~material) +
  labs(y="current",x="time in s",col="rel. humidity")
tmp_plot
ggsave(paste0("../plots_rev_deg/ExpBaseline_RH_all.pdf"),tmp_plot,width = 10,height = 4)



deg_times <- matrix(c(0),nrow(res_tab),9)
row.names(deg_times) <- row.names(res_tab)
ratios <- 9:1/10
colnames(deg_times) <- ratios

coefs <- res_tab[,5:8]
coefs[,1] <- log(coefs[,1])
colnames(coefs) <- c("log_c1","lam1","c2","lam2")

for (i in 1:nrow(res_tab)) {
  pred_start <- exp(lin_pred_RFNR(par = coefs[i,],x_use = 0,include_jumps = FALSE))
  for (k in 1:9) {
    uniroot_res <- try(uniroot(function(tmp_t) {exp(lin_pred_RFNR(par = coefs[i,],x_use = tmp_t,include_jumps = FALSE))-ratios[k]*pred_start},
                           interval = c(0,1e3),extendInt = "downX",tol = 1e-2))
    if (inherits(uniroot_res,what = "try-error")) {
      uniroot_res$root <- Inf
    }
    deg_times[i,k] <- uniroot_res$root
  }
}

mydf <- data.frame(deg_times,type=c(rep("pure",2),rep("ambient",4),rep("filter",2),rep("rh",4)))
colnames(mydf)[1:9] <- ratios

mydf <- mydf %>% mutate(id=1:12,name = row.names(res_tab)) %>%
  pivot_longer(c(1:9),names_to = "ratio",values_to = "time")
mydf$ratio <- as.numeric(mydf$ratio)

mydf %>% filter(!startsWith(name,"RH")) %>%
  ggplot(aes(x=time,y=ratio,col=type,group=id)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  geom_text_repel(data = filter(mydf,ratio==0.5,!startsWith(name,"RH")),
                  aes(x=time+1,y=ratio,label=name),
                  show.legend = FALSE,
                  size=2,
                  max.overlaps = 20) +
  labs(y="ratio of current at start")
# ggsave(paste0("../plots_rev_deg/RelDegTimes_all_MEA2.pdf"),width = 7,height = 4)

mydf %>% filter(startsWith(name,"RH")) %>%
  ggplot(aes(x=time,y=ratio,col=factor(name,levels=datanames[9:12]),group=id)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  labs(y="ratio of current at start",col="rel.hum.")
# ggsave(paste0("../plots_rev_deg/RelDegTimes_RH_MEA2.pdf"),width = 7,height = 4)


# can we say: under Cond A, the decay is significantly faster than under conditions B?


deg_times1 <- matrix(c(0),nrow(res1$res_tab),9)
row.names(deg_times1) <- row.names(res1$res_tab)
colnames(deg_times1) <- ratios

coefs1 <- res1$res_tab[,5:8]
coefs1[,1] <- log(coefs1[,1])
colnames(coefs1) <- c("log_c1","lam1","c2","lam2")

for (i in 1:nrow(res1$res_tab)) {
  pred_start <- exp(lin_pred_RFNR(par = coefs1[i,],x_use = 0,include_jumps = FALSE))
  for (k in 1:9) {
    uniroot_res <- try(uniroot(function(tmp_t) {exp(lin_pred_RFNR(par = coefs1[i,],x_use = tmp_t,include_jumps = FALSE))-ratios[k]*pred_start},
                               interval = c(0,1e3),extendInt = "downX",tol = 1e-2))
    if (inherits(uniroot_res,what = "try-error")) {
      uniroot_res$root <- Inf
    }
    deg_times1[i,k] <- uniroot_res$root
  }
}

mydf1 <- data.frame(deg_times1,type=c(rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4)))
colnames(mydf1)[1:9] <- ratios

mydf1 <- mydf1 %>% mutate(id=1:13,name = row.names(res1$res_tab)) %>%
  pivot_longer(c(1:9),names_to = "ratio",values_to = "time")
mydf1$ratio <- as.numeric(mydf1$ratio)


mydf_all <- rbind(data.frame(mydf,material="MEA2"),
                  data.frame(mydf1,material="orig"))



mydf_all <- mydf_all %>% mutate(isrh = startsWith(name,"RH"), type2 = ifelse(isrh,name,type)) 
mydf_all %>%
  ggplot(aes(x=time,y=ratio,col=factor(type2,levels=c("ambient","filter","pure","RH100","RH50","RH30","RH15")),
             group=id)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  # geom_text_repel(data = filter(mydf_all,ratio==0.5,!startsWith(name,"RH")),
  #                 aes(x=time+1,y=ratio,label=name),
  #                 show.legend = FALSE,
  #                 size=2,
  #                 max.overlaps = 20) +
  labs(y="ratio of current at start",col="type",) +
  facet_grid(material~isrh) +
  theme(strip.text.x = element_blank())
# ggsave(paste0("../plots_rev_deg/RelDegTimes_all.pdf"),width = 10,height = 5)

