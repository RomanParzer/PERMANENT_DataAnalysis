if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl)
source("./RFNR_model_functions.R")


# adapt data names
datanames <- c("Pure air s1","Pure air s2",
               "Ambient air s1 t1","Ambient air s1 t2",
               "Ambient air s2 t1","Ambient air s2 t2",
               "Ambient air s2 t3",
               "PM filter s2","ALL filters s2",
               paste0("RH",c(100,50,30,15))
               )
datanames <- c(datanames,paste0(datanames," MEA2"))

numdat <- length(datanames)

models <- list(
  "Decay only"=function(y,x,dum,ndx,jump_direction=NULL){
    RFNR(y,x,dum,ndx,log_order=1,jump_direction = "none",init=c(-1,-10))
  },
  "Decay curvature"=function(y,x,dum,ndx,jump_direction=NULL){
    RFNR(y,x,dum,ndx,log_order=2,jump_direction = "none",init=c(-1,-10,-1,-10))
  },
  "Decay shocks"=function(y,x,dum,ndx,jump_direction){
    RFNR(y,x,dum,ndx,log_order=1,jump_direction = jump_direction,init=c(-1,-10,0,1))
  },
  "Full"=function(y,x,dum,ndx,jump_direction){
    RFNR(y,x,dum,ndx,log_order=2,jump_direction = jump_direction,init=c(-1,-10,-1,-10,0,1))
  }
    )


mygrad <- function(x,maxh) {
  myind <- (maxh+1):(length(x)-maxh)
  grad <- sapply(myind,function(j) sum(x[j+1:maxh] - x[j-1:maxh])/(2*maxh))
  return(c(rep(NA,maxh),
           grad,
           rep(NA,maxh)))
}

hs <- c(100,100,100,
        100,100,100,
        100,100,100,
        8,100,100,100,
        100,100,100,
        100,100,100,
        100,100,100,
        8,80,100,100)
grad_cutoffs <- c(0.06,0.035,0.09,
                  0.035,0.035,0.07,
                  0.01,0.035,0.07,
                  0.01,0.035,0.035,0.04,
                  0.06,0.035,0.035,
                  0.12,0.07,0.05,
                  0.01,0.035,0.035,
                  0.02,0.028,0.035,0.04)
# just for supplementary table
# 
# mygraddf <- data.frame(datanames[1:13],hs[1:13],grad_cutoffs[1:13],hs[14:26],grad_cutoffs[14:26])
# colnames(mygraddf) <- c("data","h",expression(tau),"h",expression(tau))
# knitr::kable(mygraddf,format="latex",booktabs=TRUE) %>%
#   kableExtra::add_header_above(c(" "=1,"MEA1" = 2,"MEA2"=2))

unlink("./saved_results/log.txt")
n.cores <- parallel::detectCores() # adapt! 
my.cluster <- parallel::makeCluster(n.cores-1, type = "PSOCK", outfile = "./saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('datanames','models','mygrad','hs','grad_cutoffs'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl)
  source("./RFNR_model_functions.R")
})

# i <- 3
foreach(i=1:length(datanames)) %dopar% {
  data <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx",sheet=i)
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
  data <- data %>% mutate("grad"=grad)
  # add jump at start for following experiments
  if (i %in% c(10,13,26)) {
    ind[1] <- TRUE
  }
  
  xlo <- 5000*(5/5-1)
  stepsize <- 5e3
  # check smoothness of grad and correct cutoff (only 1 per peak)
  tmp_plot <- data %>%
    filter(time>xlo,
           time<xlo+stepsize) %>%
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type,linetype=type)) +
    geom_line(alpha=0.8) +
    geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1) +
    coord_cartesian(xlim=c(xlo,xlo+stepsize)) +
    geom_hline(yintercept = grad_cutoff,linetype=2) +
    labs(x="time (s)") +
    theme_classic() +
    scale_color_brewer(type="qual",direction = 1,palette=1)
  # tmp_plot
  ggsave(paste0("./plots_rev_deg/JumpsClose_",datanames[i],".pdf"),tmp_plot,width = 8*0.8,height = 5*0.7)

  # check if all proper peaks are detected and not more or less
  tmp_plot <- data %>% mutate(id=1:length(data$time)) %>%
    filter(id%%100==0) %>%
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type,linetype=type)) +
    geom_line(alpha=0.6) +
    geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1,alpha=0.4) +
    labs(x="time (s)") +
    theme_classic() +
    scale_color_brewer(type="qual",direction = 1,palette=1)
  # tmp_plot
  ggsave(paste0("./plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8*0.8,height = 5*0.7)

  # # application
  x <- data$time
  y <- data$current
  # less downsampling for RH100 cases
  if (i %in% c(10,23)) {
    ndx <- seq(1, length(x), 12*4) # one measurement every 12 s / 5 per minute
  } else {
    ndx <- seq(1, length(x), 60*4) # one measurement every minute
  }
  # str(y[ndx])
  
  xnew <- x[-ndx]
  time_ind <- 1 + (1:length(x))[-ndx] * (length(ndx)-1) / (max(ndx))
  
  # for each jump, closest observation gets the dummy=1
  dum <- numeric(length(ndx))
  for (k in which(ind)) {
    dum[which.min(abs(k - ndx))] <- 1
  }
  
  # compare 4 models on 70/30 train / test split and on withheld data on full
  ind_train <- 1:round(0.7*length(ndx))
  ind_test <- (round(0.7*length(ndx))+1):length(ndx)
  model_comp <-  data.frame(data=NULL,material=NULL,model=NULL,error=NULL,error_measure=NULL,value=NULL)
  if (i < 14) {
    tmp_data_name <- datanames[i]
    tmp_material <- "MEA1"
  } else {
    tmp_data_name <- datanames[i-13]
    tmp_material <- "MEA2"
  }
  
  for (k in 1:length(models)) {
    modres <- try(models[[k]](y,x,dum[ind_train],ndx[ind_train],jump_direction))
    if (!("try-error" %in% class(modres))) {
      model_comp <- rbind(model_comp,
                          data.frame(data=tmp_data_name,material=tmp_material,
                                     model=names(models)[k],
                                     error="train/test",error_measure=c("RMSE","Cor","MAE"),
                                     value=as.numeric(eval_fit(log(y[ndx[ind_test]]),predict(modres,xnew = x[ndx[ind_test]],dumnew = dum[ind_test],type = "linear")))))
    }
    modres <- try(models[[k]](y,x,dum,ndx,jump_direction))
    if (!("try-error" %in% class(modres))) {
      model_comp <- rbind(model_comp,
                          data.frame(data=tmp_data_name,material=tmp_material,model=names(models)[k],
                                     error="sample fit",error_measure=c("RMSE","Cor","MAE"),
                                     value=as.numeric(eval_fit(log(y[ndx]),predict(modres,type = "linear")))))
      model_comp <- rbind(model_comp,
                          data.frame(data=tmp_data_name,material=tmp_material,model=names(models)[k],
                                     error="withheld fit",error_measure=c("RMSE","Cor","MAE"),
                                     value=as.numeric(eval_fit(log(y[-ndx]),predict(modres,xnew = xnew,time_ind = time_ind,type = "linear")))))
    }
  }
  data_tmp <- data %>%
    mutate(id=1:length(data$time)) %>%
    filter(id%in%ndx)

  tmp_plot <- plot(modres) +
    aes(linetype=type) +
    labs(x="time (s)") +
    theme_classic()
  # tmp_plot
  ggsave(paste0("./plots_rev_deg/Fitted_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8*0.6,height = 5*0.6)

  tmp_plot <- plot(modres,"residuals") +
    labs(x="time (s)")+
    theme_classic()
  # tmp_plot
  ggsave(paste0("./plots_residuals_rev_deg/Residuals_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8*0.6,height = 5*0.6)

  tmp_plot <- plot(modres,"res-QQ")+
    theme_classic()
  # tmp_plot
  ggsave(paste0("./plots_residuals_rev_deg/ResQQ_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8*0.6,height = 5*0.6)
  
  tmp_plot <- plot(modres,"res-vs-fitted")+
    theme_classic()
  # tmp_plot
  ggsave(paste0("./plots_residuals_rev_deg/ResFitted_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8*0.6,height = 5*0.6)
  
  tmp_plot <- plot(modres,"res-hist")+
    theme_classic()
  # tmp_plot
  ggsave(paste0("./plots_residuals_rev_deg/ResHist_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8*0.6,height = 5*0.6)
  
  tmp_plot <- plot(modres,"res-acf")+
    theme_classic()
  # tmp_plot
  ggsave(paste0("./plots_residuals_rev_deg/ResACF_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8*0.6,height = 5*0.64)

  # normality tests
  norm_tests <- list(
    "ad" = goftest::ad.test(modres$log_residuals, "pnorm", mean=0, sd=modres$sigma,estimated = TRUE),
    "cvm" = goftest::cvm.test(modres$log_residuals, null = "pnorm", mean=0, sd=modres$sigma, estimated=TRUE),
    "shapiro" = shapiro.test(modres$log_residuals),
    "ks" = ks.test(modres$log_residuals,"pnorm",mean=0,sd=modres$sigma)
  )


  # estimated parameters and significance
  model_fits <- numeric(21)
  names(model_fits) <- c("jumps","mean_time_jumps","rh",
                            "c1","lower","upper",
                            "lam1","lower","upper",
                            "c2","lower","upper",
                            "lam2","lower","upper",
                            "delta","lower","upper",
                            "rho","lower","upper")
  summod <- summary(modres,corr=TRUE)
  alpha <- 0.05
  myqnorm <- qnorm(alpha/2,lower.tail = FALSE)
  model_fits[1] <- length(which(ind))
  model_fits[2] <- mean(diff(data$time[which(ind)]))
  model_fits[3] <- c(rep(30,9),100,50,30,15,rep(30,9),100,50,30,15)[i]
  
  model_fits[4:6] <- exp(summod$par_sig[1,1] + c(0,-1,1)*summod$par_sig[1,2]*myqnorm)
  model_fits[7:9] <- exp(summod$par_sig[2,c(1)] + c(0,-1,1)*summod$par_sig[2,c(2)]*myqnorm)
  model_fits[10:12] <- exp(summod$par_sig[3,c(1)] + c(0,-1,1)*summod$par_sig[3,c(2)]*myqnorm)
  model_fits[13:15] <- exp(summod$par_sig[4,c(1)] + c(0,-1,1)*summod$par_sig[4,c(2)]*myqnorm)
  model_fits[16:18] <- exp(summod$par_sig[5,c(1)] + c(0,-1,1)*summod$par_sig[5,c(2)]*myqnorm)
  
  model_fits[19:21] <- modres$rho_upper / (1+exp(-summod$par_sig[6,1] + c(0,1,-1)*summod$par_sig[6,2]*myqnorm))

  my_full_dat <- data.frame(data_tmp,data=tmp_data_name,material=tmp_material,fit=as.numeric(modres$fitted_jumps))
  
  resi <- list(dataname = datanames[i],
                model_fits=model_fits,
                sum_list = list(summary=summod,norm_tests = norm_tests),
                my_full_dat = my_full_dat,
                model_comp = model_comp)
  cat(sprintf('Finished rep %d / %d at %s.\n',i,numdat,Sys.time()))
  saveRDS(resi,paste0("./saved_results_loop/RFNR_results_i",i,".rds"))
  resi
}


model_fits <- matrix(NA,numdat,21)
row.names(model_fits) <- datanames
colnames(model_fits) <- c("jumps","mean_time_jumps","rh",
                          "c1","lower","upper",
                          "lam1","lower","upper",
                          "c2","lower","upper",
                          "lam2","lower","upper",
                          "delta","lower","upper",
                          "rho","lower","upper")

sum_list <- vector("list",numdat)
names(sum_list) <- datanames


model_comp <-  data.frame(data=NULL,model=NULL,error=NULL,error_measure=NULL,value=NULL)
my_full_dat <- data.frame(NULL)

for (i in 1:numdat) {
  resi <- readRDS(paste0("./saved_results_loop/RFNR_results_i",i,".rds"))
  model_fits[i,] <- resi$model_fits
  sum_list[[i]] <- resi$sum_list
  my_full_dat <- rbind(my_full_dat,
                       resi$my_full_dat)
  model_comp <- rbind(model_comp,
                       resi$model_comp)
}

saveRDS(list(model_fits=model_fits,
             sum_list = sum_list,
             my_full_dat = my_full_dat,
             datanames = datanames,
             model_comp = model_comp),
        "./saved_results/results_RFNR_final.rds")
    
