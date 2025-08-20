
pacman::p_load(dplyr,tidyr,ggplot2,knitr,ggrepel)
source("./RFNR_model_functions.R")

res <- readRDS("../saved_results/results_RFNR_model.rds")

res_tab <- res$res_tab
kable_res_tab <- res_tab[,1:10]
kable_res_tab[,c(2)] <- round(res_tab[,c(2)],digits=1)
kable_res_tab[,c(4,5,7,9,10)] <- round(res_tab[,c(4,5,7,9,10)],digits=3)
kable_res_tab[,c(6,8)] <- format(res_tab[,c(6,8)],digits=3)
kable(kable_res_tab,booktabs=TRUE,format = "latex")

my_full_dat <- res$my_full_dat
datanames <- res$datanames
tmp_plot <- my_full_dat %>%
  filter(dataname %in% c("Pure Air","Ambient air test 1","PM filter sample 2")) %>%
  ggplot(aes(x=time,y=current,col=dataname)) +
  geom_line(alpha=0.8) +
  labs(y="current",x="time in s",col="data")
tmp_plot
# ggsave(paste0("../plots_rev_deg/Example_data.pdf"),tmp_plot,width = 7,height = 3)


tmp_plot <- my_full_dat %>%
  filter(!startsWith(dataname,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~factor(dataname,levels=datanames[name_order])) +
  labs(y="current",x="time in s")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_all.pdf"),tmp_plot,width = 10,height = 6)


tmp_plot <- my_full_dat %>%
  filter(startsWith(dataname,"RH")) %>%
  pivot_longer(c(fit,current),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,linetype=type,col=factor(dataname,levels=c(datanames[10:13])))) +
  geom_line(alpha=0.8) +
  labs(y="current",x="time in s",col="rel. humidity")
tmp_plot
# ggsave(paste0("../plots_rev_deg/ExpBaseline_RH.pdf"),tmp_plot,width = 8,height = 4)



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
    uniroot_res <- uniroot(function(tmp_t) {exp(lin_pred_RFNR(par = coefs[i,],x_use = tmp_t,include_jumps = FALSE))-ratios[k]*pred_start},
                           interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
    deg_times[i,k] <- uniroot_res$root
  }
}


mydf <- data.frame(deg_times,type=c(rep("pure",2),rep("ambient",5),rep("filter",2),rep("rh",4)))
colnames(mydf)[1:9] <- ratios

mydf <- mydf %>% mutate(id=1:13,name = row.names(res_tab)) %>%
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
# ggsave(paste0("../plots_rev_deg/RelDegTimes_all.pdf"),width = 7,height = 4)

mydf %>% filter(startsWith(name,"RH"),name!="RH15") %>%
  ggplot(aes(x=time,y=ratio,col=name,group=id)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  labs(y="ratio of current at start",col="rel.hum.")
# ggsave(paste0("../plots_rev_deg/RelDegTimes_RH.pdf"),width = 7,height = 4)


# TODOs: 
# - look at jumps towards the start: low start level for Pure Air, Amb Air 1 & sample 2, RH15
# - CI for rho and delta (delta or just monotone transf)
# - check why c2 for rh15 negative: also due to high jump towards begining
# - also lam1 negative for 2 datasets; force positive in optimization?
# - maybe log-order 3 for some?
# - maybe dummy not just 0 1, but 2 3 for higher jumps compared to the lower base jumps, 
# or jump height directly?
# - or try to let the algorithm find the jumps itself?

# can we say: under Cond A, the decay is significantly faster than under conditions B?



