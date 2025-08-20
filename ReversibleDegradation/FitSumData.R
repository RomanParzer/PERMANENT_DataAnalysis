if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(dplyr,tidyr,ggplot2,readxl,quantreg,robustbase)

datanames <- c("Pure Air","Ambient air test 1","Ambient air test 2",
               "Pure air sample 2","Ambient air sample 2","Ambient air sample 2 test2",
               "Ambient air sample 2 test3","PM filter sample 2","ALL Filters sample 2",
               paste0("RH",c(100,50,30,15))
               )

numdat <- length(datanames)
data_order <- c(1,3,4,2,5,6,7,8,9,
                10,11,12,13)
name_order <- c(1,4,2,3,5,6,7,8,9,
                10,11,12,13)

res_tab <- matrix(NA,numdat,10)
row.names(res_tab) <- datanames[name_order]
sum_list <- vector("list",numdat)
names(sum_list) <- datanames[name_order]
colnames(res_tab) <- c("jumps","mean_time_jumps","rh","cor_mod","C1","lam1","C2","lam2","h","grad_cutoff")


mygrad <- function(x,maxh) {
  myind <- (maxh+1):(length(x)-maxh)
  grad <- sapply(myind,function(j) sum(x[j+1:maxh] - x[j-1:maxh])/(2*maxh))
  return(c(rep(NA,maxh),
           grad,
           rep(NA,maxh)))
}

hs <- c(50,50,50,
        50,50,50,
        50,50,50,
        10,50,50,80)
grad_cutoffs <- c(0.0012,0.0008,0.0005,
                  0.0050,0.0030,0.0016,
                  0.0001,0.0032,0.0023,
                  0.005,0.0008,0.0008,0.0008)
i <- 6

my_full_sum_dat <- my_full_dat <- data.frame(NULL)

for (i in 1:numdat) {
  if (i<=9) {
    data <- read_excel("../data/Dati_PERMANENT_reversible.xlsx",sheet=i)
  } else {
    data <- read_excel("../data/Dati_PERMANENT_reversible_RH.xlsx",sheet=i-9)
  }
  colnames(data) <- c("time","current")
  # change to s
  data$time <- (data$time-1)*60
  
  # set correct h and cutoff 
  grad <- mygrad(data$current,hs[i])
  grad_cutoff <- grad_cutoffs[i]
  if (i==10) {
    ind <- grad <  -grad_cutoff
  } else {
    ind <- grad > grad_cutoff
  }
  
  # take first time point of positive gradient regions
  for (k in c(length(ind):2)) {
    if (!is.na(ind[k]*ind[k-1]) & ind[k] & ind[k-1]) {
      ind[k] <- FALSE
    }
  }
  # summary(ind)
  
  xlo <- 1000*0
  stepsize <- 1e3
  # check smoothness of grad and correct cutoff (only 1 per peak)
  tmp_plot <- data %>% mutate("grad"=grad) %>% 
    filter(time>xlo,
           time<xlo+stepsize) %>% 
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type)) +
    geom_line(alpha=0.8) + 
    geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1) +
    coord_cartesian(xlim=c(xlo,xlo+stepsize)) +
    geom_hline(yintercept = grad_cutoff,linetype=2)
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/JumpsClose_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  # check if all proper peaks are detected and not more or less
  tmp_plot <- data %>% mutate("grad"=grad+0.15,id=1:length(data$time)) %>% 
    filter(id%%100==0) %>% 
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type)) +
    geom_line(alpha=0.8) + 
    geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1) 
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  diffs <- diff(c(which(ind),length(data$time)+1))
  time_last_jump <- c(data$time[1:(which(ind)[1]-1)],rep(c(data$time)[c(which(ind))],diffs))
  jumps <- rep(c(0,1:sum(ind,na.rm = TRUE)),c(which(ind)[1]-1,diffs))
  data <- data %>% mutate(time_l_jump=time - time_last_jump,jumps=jumps)
  data$jumps <- factor(data$jumps)
  
  if (i==10) {
    data_jumps <- data %>% group_by(jumps) %>% summarize(max = max(current),
                                                         mean = mean(current),
                                                         med = median(current),
                                                         jump_start_time = min(time),
                                                         peak_time = time[which.max(current)],
                                                         min = min(current),
                                                         min_time = time[which.min(current)],
                                                         jump_duration = max(time) - min(time),
                                                         jump_height = max - min,
                                                         resp=max,
                                                         med_time= ifelse(any(time>=peak_time & current<=med),
                                                                          min(time[time>=peak_time & current<=med]),
                                                                          max(time[time<=peak_time])),
                                                         time=peak_time)
  } else {
    data_jumps <- data %>% group_by(jumps) %>% summarize(max = max(current),
                                                         mean = mean(current),
                                                         med = median(current),
                                                         jump_start_time = min(time),
                                                         peak_time = time[which.max(current)],
                                                         min = min(current[time>=peak_time]),
                                                         min_time = time[time>=peak_time][which.min(current[time>=peak_time])],
                                                         jump_duration = max(time) - min(time),
                                                         jump_height = max - head(current,1),
                                                         resp=min,
                                                         med_time= ifelse(any(time>=peak_time & current<=med),
                                                                          min(time[time>=peak_time & current<=med]),
                                                                          max(time[time<=peak_time])),
                                                         time=min_time)
  }
  
  
  if (nrow(data_jumps)<7) {
    data_jumps <- pivot_longer(data_jumps[,-11],cols = c(min,med),names_to = "resp_type",values_to = "resp") %>% 
      mutate("time"=ifelse(resp_type=="min",min_time,med_time))
  }
  
  # data_jumps
  # hist(data_jumps$jump_height)
  # hist(data_jumps$jump_duration)
  if (i==4) {
    b1start <- 0
  } else {
    b1start <- -2e-5
  }

  
  model_baseline <- try(robustbase::nlrob(log(resp) ~ b0 + b1*time + b2*exp(b3*time),
                                               data=data_jumps,
                                               start = c(b0=-1,b1=b1start,b2=0.4,b3=-2e-4)))
  if ("try-error" %in% class(model_baseline)) {
    model_baseline <- try(quantreg::nlrq(log(resp) ~ b0 + b1*time + b2*exp(b3*time),
                                     data=data_jumps,tau=0.25,
                                     start = c(b0=-1,b1=b1start,b2=0.4,b3=-2e-4)))
    model_baseline$rweights <- 1
  }
  
  sum_list[[data_order[i]]] <- summary(model_baseline)
  
  data_tmp <- data %>% 
    mutate(id=1:length(data$time),
           type="current") %>%
    filter(id%%500==0)
  
  
  my_full_dat <- rbind(my_full_dat,
                       data.frame(data_tmp,dataname=datanames[i]))
  
  my_full_sum_dat <- rbind(my_full_sum_dat,
                           data.frame(data_jumps %>% select(resp,time) %>%
                                        mutate("fit"=exp(predict(model_baseline)),
                                               "type"="fit"),dataname=datanames[i]))
  
  tmp_plot <- data_jumps %>% 
    mutate("fit"=exp(predict(model_baseline)),
           "type"="fit") %>%
    # pivot_longer(c(resp,fit),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=fit,color=type)) +
    geom_line(alpha=0.8) +
    geom_point(aes(x=time,y=resp),alpha=0.8,col="darkgrey") +
    geom_line(data=data_tmp,aes(x=time,y=current),alpha=0.8)
  # tmp_plot
  # ggsave(paste0("../plots_rev_deg/MinBaseline_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  tmp_plot <- data_jumps %>% 
    mutate("fit"=exp(predict(model_baseline)),
           "weight" = model_baseline$rweights) %>%
    ggplot(aes(x=time,y=resp-fit,col=weight)) +
    geom_point(alpha=0.8) +
    geom_hline(yintercept = 0,col=2,alpha=0.8)
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/MinBaselineResiduals_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  res_tab[data_order[i],1] <- tail(jumps,1)
  res_tab[data_order[i],2] <- mean(data_jumps$jump_duration)
  res_tab[data_order[i],3]<- c(rep(30,9),100,50,30,15)[i]
  res_tab[data_order[i],4] <- cor(exp(predict(model_baseline)),data_jumps$resp)
  res_tab[data_order[i],5] <- exp(summary(model_baseline)$coefficients[1,1])
  res_tab[data_order[i],6] <- -(summary(model_baseline)$coefficients[2,1])
  res_tab[data_order[i],7] <- (summary(model_baseline)$coefficients[3,1])
  res_tab[data_order[i],8] <- -(summary(model_baseline)$coefficients[4,1])
  res_tab[data_order[i],9] <- hs[i]
  res_tab[data_order[i],10] <- grad_cutoffs[i]
  
  cat(sprintf('Finished rep %d / %d at %s.\n',i,numdat,Sys.time()))
}

# sum_list

res_tab
kable_res_tab <- res_tab[,1:8]
kable_res_tab[,c(2)] <- round(res_tab[,c(2)],digits=1)
kable_res_tab[,c(4,5,7)] <- round(res_tab[,c(4,5,7)],digits=3)
kable_res_tab[,c(6,8)] <- format(res_tab[,c(6,8)],digits=3)
require(knitr)
kable(kable_res_tab,booktabs=TRUE,format = "latex")



tmp_plot <- my_full_sum_dat %>% 
  filter(!startsWith(dataname,"RH")) %>%
  ggplot(aes(x=time,y=fit,color=type)) +
  geom_line(alpha=0.8) +
  geom_point(aes(x=time,y=resp),alpha=0.8,col="darkgrey") +
  geom_line(data=filter(my_full_dat,!startsWith(dataname,"RH")),aes(x=time,y=current),alpha=0.8) +
  facet_wrap(.~factor(dataname,levels=datanames[name_order]))
tmp_plot
# ggsave(paste0("../plots_rev_deg/MinBaseline_all.pdf"),tmp_plot,width = 12,height = 8)

tmp_plot <- my_full_sum_dat %>% 
  filter(startsWith(dataname,"RH")) %>%
  ggplot(aes(x=time,y=fit,color=type)) +
  geom_line(alpha=0.8) +
  geom_point(aes(x=time,y=resp),alpha=0.8,col="darkgrey") +
  geom_line(data=filter(my_full_dat,startsWith(dataname,"RH")),aes(x=time,y=current),alpha=0.8) +
  facet_wrap(.~factor(dataname,levels=datanames[name_order]))
tmp_plot
# ggsave(paste0("../plots_rev_deg/MinBaseline_RH.pdf"),tmp_plot,width = 10,height = 7)
