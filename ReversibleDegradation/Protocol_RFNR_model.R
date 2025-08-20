if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(dplyr,tidyr,ggplot2,readxl)
source("./RFNR_model_functions.R")

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

res_tab <- matrix(NA,numdat,12)
row.names(res_tab) <- datanames[name_order]
sum_list <- vector("list",numdat)
names(sum_list) <- datanames[name_order]
colnames(res_tab) <- c("jumps","mean_time_jumps","rh","cor_mod","C1","lam1","C2","lam2","delta","rho","h","grad_cutoff")


mygrad <- function(x,maxh) {
  myind <- (maxh+1):(length(x)-maxh)
  grad <- sapply(myind,function(j) sum(x[j+1:maxh] - x[j-1:maxh])/(2*maxh))
  return(c(rep(NA,maxh),
           grad,
           rep(NA,maxh)))
}

hs <- c(200,200,200,
        200,200,200,
        200,200,200,
        8,200,200,200)
grad_cutoffs <- c(0.06,0.035,0.035,
                  0.035,0.035,0.035,
                  0.01,0.035,0.035,
                  0.01,0.035,0.035,0.035)
my_full_dat <- data.frame(NULL)

for (i in 1:numdat) {
  i <- 3
  # datanames[i]
  if (i<=9) {
    data <- read_excel("../data/Dati_PERMANENT_reversible.xlsx",sheet=i)
  } else {
    data <- read_excel("../data/Dati_PERMANENT_reversible_RH.xlsx",sheet=i-9)
  }
  if (i==10) {
    jump_direction <- "down"
  } else {
    jump_direction <- "up"
  }
  colnames(data) <- c("time","current")
  # change to s
  data$time <- (data$time-1)*60

  # set correct h and cutoff
  grad <- mygrad(log(data$current),hs[i])
  grad_cutoff <- grad_cutoffs[i]
  if (i==10) {
    grad <- -grad
  }
  ind <- grad > grad_cutoff

  # only take last time point of positive gradient regions
  for (k in c(1:(length(ind)-1))) {
    if (!is.na(ind[k]*ind[k+1]) &ind[k] & ind[k+1]) {
      ind[k] <- FALSE
    }
  }

  xlo <- 5000*0
  stepsize <- 5e3
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
  tmp_plot <- data %>% mutate("grad"=grad+0.1,id=1:length(data$time)) %>%
    filter(id%%100==0) %>%
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type)) +
    geom_line(alpha=0.6) +
    geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1)
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  x <- data$time
  y <- data$current
  ndx <- seq(1, length(x), 60*4) # one measurement every minute

  # for each jump, closest observation gets the dummy=1
  dum <- numeric(length(ndx))
  for (k in which(ind)) {
    dum[which.min(abs(k - ndx))] <- 1
  }
  log_order <- 2
  model_baseline <- try(RFNR(y,x,dum,ndx,log_order=log_order,jump_direction = jump_direction))
  if ("try-error" %in% class(model_baseline)) {
    model_baseline <- try(fit_exp_deg_filter(y,x,dum,ndx,log_order=1,jump_direction = jump_direction))
  }
  sum_list[[data_order[i]]] <- summary(model_baseline)

  data_tmp <- data %>%
    mutate(id=1:length(data$time)) %>%
    filter(id%in%ndx)

  my_full_dat <- rbind(my_full_dat,
                       data.frame(data_tmp,dataname=datanames[i],fit=as.numeric(model_baseline$fitted_jumps)))


  tmp_plot <- plot(model_baseline) +
    labs(x="time in s")
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/Fitted_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  tmp_plot <- plot(model_baseline,"residuals") +
    labs(x="time in s")
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/Residuals_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  res_tab[data_order[i],1] <- length(which(ind))
  res_tab[data_order[i],2] <- mean(diff(data$time[which(ind)]))
  res_tab[data_order[i],3]<- c(rep(30,9),100,50,30,15)[i]
  res_tab[data_order[i],4] <- cor(predict(model_baseline),data$current[ndx])
  res_tab[data_order[i],5] <- exp(model_baseline$pars$log_c1)
  res_tab[data_order[i],6] <- model_baseline$pars$lam1
  res_tab[data_order[i],7] <- model_baseline$pars$c2
  res_tab[data_order[i],8] <- model_baseline$pars$lam2
  res_tab[data_order[i],9] <- exp(model_baseline$pars$log_delta)
  res_tab[data_order[i],10] <- model_baseline$rho_upper / (1+exp(-model_baseline$pars$phi))
  res_tab[data_order[i],11] <- hs[i]
  res_tab[data_order[i],12] <- grad_cutoffs[i]

  cat(sprintf('Finished rep %d / %d at %s.\n',i,numdat,Sys.time()))
}

# sum_list

saveRDS(list(res_tab=res_tab,
             sum_list = sum_list,
             my_full_dat = my_full_dat,
             datanames = datanames,
             data_order = data_order),
        "../saved_results/results_RFNR_model.rds")
    
