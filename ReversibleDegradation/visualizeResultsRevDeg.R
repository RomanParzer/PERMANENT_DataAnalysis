pacman::p_load(dplyr,tidyr,ggplot2,knitr,ggrepel)
reslist <- readRDS("./saved_results/RevDegResults_adapted.rds")

datanames <- c("Pure Air","Ambient air test 1","Ambient air test 2",
               "Pure air sample 2","Ambient air sample 2","Ambient air sample 2 test2",
               "Ambient air sample 2 test3","PM filter sample 2","ALL Filters sample 2")
data_order <- c(1,4,2,3,5,6,7,8,9)

res_tab <- matrix(NA,9,6)
row.names(res_tab) <- datanames[data_order]
colnames(res_tab) <- c("jumps","cor","C1","lam1","C2","lam2")

# 1 for des
# 2 for les
mod <- 1
for (i in 1:9) {
  res_tab[data_order[i],1] <- reslist[[i]]$num_jumps
  res_tab[data_order[i],2] <- reslist[[i]]$model_res_list[[mod]]$cor
  res_tab[data_order[i],3] <- exp(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[1,1])
  res_tab[data_order[i],4] <- -(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[2,1])
  res_tab[data_order[i],5] <- (reslist[[i]]$model_res_list[[mod]]$summary$coefficients[3,1])
  res_tab[data_order[i],6] <- -(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[4,1])
}
res_tab
kable_res_tab <- res_tab
kable_res_tab[,c(2,3,5)] <- format(res_tab[,c(2,3,5)],digits=1)
kable_res_tab[,c(4,6)] <- format(res_tab[,c(4,6)],digits=3)
kable(kable_res_tab,booktabs=TRUE,format = "latex")


res_tab_jumps <- matrix(NA,9,9)
row.names(res_tab_jumps) <- datanames[data_order]
colnames(res_tab_jumps) <- c("jumps","cor","C1","lam1","C2","lam2","lam3","C3","lam4")

# 4 for DE
# 5 for LE
mod <- 4
for (i in 1:9) {
  if (i %in% c(6,7)) {
    res_tab_jumps[data_order[i],] <- NA
  } else {
    res_tab_jumps[data_order[i],1] <- reslist[[i]]$num_jumps
    res_tab_jumps[data_order[i],2] <- reslist[[i]]$model_res_list[[mod]]$cor
    res_tab_jumps[data_order[i],3] <- exp(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[3,1])
    res_tab_jumps[data_order[i],4] <- -(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[1,1])
    res_tab_jumps[data_order[i],5] <- reslist[[i]]$model_res_list[[mod]]$summary$coefficients[4,1] -
      reslist[[i]]$model_res_list[[mod]]$summary$coefficients[3,1]
    res_tab_jumps[data_order[i],6] <- exp(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[5,1])
    res_tab_jumps[data_order[i],7] <- -(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[2,1])
    res_tab_jumps[data_order[i],8] <- reslist[[i]]$model_res_list[[mod]]$summary$coefficients[6,1]
    res_tab_jumps[data_order[i],9] <- exp(reslist[[i]]$model_res_list[[mod]]$summary$coefficients[7,1])
  }
  
}
res_tab_jumps
kable_res_tab_jumps <- res_tab_jumps
kable_res_tab_jumps[,c(2,3,5,8)] <- format(res_tab_jumps[,c(2,3,5,8)],digits=1)
kable_res_tab_jumps[,c(4,6,7,9)] <- format(res_tab_jumps[,c(4,6,7,9)],digits=3)
kable(kable_res_tab_jumps,booktabs=TRUE,format = "latex")



deg_time_tab <- matrix(NA,9,10)
row.names(deg_time_tab) <- datanames[data_order]
colnames(deg_time_tab) <- c("start",paste0("tt",9:1*10,"%"))
# time to reach specific relative degradations in minutes

# 1 for des
# 2 for les
# 3 for gam (uncertainty quantification)
# 4 for DE (i6 and i7 NULL)
# 5 for LE (i6 NULL)
mod <- 4
for (i in 1:9) {
  if (i %in% c(6,7)) {
    deg_time_tab[data_order[i],1] <- NA
    deg_time_tab[data_order[i],2:10] <- NA
  } else {
    deg_time_tab[data_order[i],1] <- reslist[[i]]$model_res_list[[mod]]$pred_start
    deg_time_tab[data_order[i],2:10] <- reslist[[i]]$model_res_list[[mod]]$deg_times/60
  }
}
deg_time_tab
# kable(deg_time_tab,booktabs=TRUE,format = "latex")


# visualize
mydf <- data.frame(deg_time_tab,"1"=0,type=c(rep("pure",2),rep("ambient",5),rep("filter",2)))
colnames(mydf)[2:11] <- c(9:1/10,1)

mydf <- mydf %>% mutate(id=1:9,name = datanames[id]) %>%
  pivot_longer(c(2:11),names_to = "current_val",values_to = "time")
mydf$current_val <- as.numeric(mydf$current_val)

mydf %>% 
  ggplot(aes(x=time+1,y=current_val,col=type,group=id)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  geom_text_repel(data = filter(mydf,current_val==0.5),
                  aes(x=time+1,y=current_val,label=name),
                  show.legend = FALSE,
                  size=2,
                  max.overlaps = 20)
# ggsave("./plots_rev_deg/deg_time_des_relative.pdf",width = 6,height = 5)
# ggsave("./plots_rev_deg/deg_time_DE_relative.pdf",width = 6,height = 5)

mydf %>% 
  ggplot(aes(x=time+1,y=current_val*start,col=type,group=id)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() +
  geom_text_repel(data = filter(mydf,current_val==0.5),
                  aes(x=time+1,y=current_val*start,label=name),
                  show.legend = FALSE,
                  size=2,
                  max.overlaps = 20)
# ggsave("./plots_rev_deg/deg_time_des_absolute.pdf",width = 6,height = 5)
# ggsave("./plots_rev_deg/deg_time_DE_absolute.pdf",width = 6,height = 5)


# might depend on model, compare for i=1 Pure Air
i <- 1
data <- read_excel("./data/Dati_PERMANENT_reversible.xlsx",sheet=i)
colnames(data) <- c("time","current")
pred_start <- data$current[1]
deg_times_data <- numeric(9)
ratios <- 9:1/10
for (k in 1:9) {
  myt <- 0
  for (it in 1:length(data$time)) {
    if (data$current[it]/pred_start <= ratios[k]) {
      break
    }
  }
  if (it==length(data$time)) {
    deg_times_data[k] <- NA
  } else {
    deg_times_data[k] <- data$time[it]
  }
}

deg_time_mods_tab <- matrix(NA,6,10)
row.names(deg_time_mods_tab) <- c("des","les","gam","DE","LE","Data")
colnames(deg_time_mods_tab) <- c("start",paste0("tt",9:1*10,"%"))
# time to reach specific relative degradations in minutes

for (mod in 1:5) {
  deg_time_mods_tab[mod,1] <- reslist[[i]]$model_res_list[[mod]]$pred_start
  deg_time_mods_tab[mod,2:10] <- reslist[[i]]$model_res_list[[mod]]$deg_times/60
}
deg_time_mods_tab[6,1] <-  pred_start
deg_time_mods_tab[6,2:10] <- deg_times_data

deg_time_mods_tab

# visualize
mydf_mods <- data.frame(deg_time_mods_tab,"1"=0)
colnames(mydf_mods)[2:11] <- c(9:1/10,1)

mydf_mods <- mydf_mods %>% mutate(id=1:6,model = c("des","les","gam","DE","LE","Data")) %>%
  pivot_longer(c(2:11),names_to = "current_val",values_to = "time")
mydf_mods$current_val <- as.numeric(mydf_mods$current_val)

mydf_mods %>% 
  ggplot(aes(x=time+1,y=current_val,col=model)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10()
# ggsave("./plots_rev_deg/deg_time_mods_pureair_relative.pdf",width = 6,height = 5)

mydf_mods %>% 
  ggplot(aes(x=time+1,y=current_val*start,col=model)) +
  geom_line(alpha=0.6,linewidth=1) +
  scale_x_log10() 
# ggsave("./plots_rev_deg/deg_time_mods_pureair_absolute.pdf",width = 6,height = 5)
