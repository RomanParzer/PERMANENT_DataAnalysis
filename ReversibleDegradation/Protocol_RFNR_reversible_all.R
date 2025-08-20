if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl)
source("./RFNR_model_functions.R")

datanames <- c("Pure air","Ambient air test 1","Ambient air test 2",
               "Pure air sample 2","Ambient air sample 2","Ambient air sample 2 test2",
               "Ambient air sample 2 test3",
               "PM filter sample 2","ALL Filters sample 2",
               paste0("RH",c(100,50,30,15))
               )
datanames <- c(datanames,paste0(datanames," MEA2"))

numdat <- length(datanames)

models <- list(
  "Decay only"=function(y,x,dum,ndx,jump_direction=NULL){
    RFNR(y,x,dum,ndx,log_order=1,jump_direction = "none")
  },
  "Decay curvature"=function(y,x,dum,ndx,jump_direction=NULL){
    RFNR(y,x,dum,ndx,log_order=2,jump_direction = "none")
  },
  "Decay shocks"=function(y,x,dum,ndx,jump_direction){
    RFNR(y,x,dum,ndx,log_order=1,jump_direction = jump_direction)
  },
  "Full"=function(y,x,dum,ndx,jump_direction){
    RFNR(y,x,dum,ndx,log_order=2,jump_direction = jump_direction)
  }
    )


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
        8,200,200,200,
        200,200,200,
        200,200,200,
        200,200,200,
        100,80,200,200)
grad_cutoffs <- c(0.06,0.035,0.035,
                  0.035,0.035,0.035,
                  0.01,0.035,0.035,
                  0.01,0.035,0.035,0.035,
                  0.06,0.035,0.032,
                  0.035,0.035,0.035,
                  0.01,0.035,0.035,
                  0.02,0.028,0.035,0.04)

unlink("../saved_results/log.txt")
my.cluster <- parallel::makeCluster(7, type = "PSOCK", outfile = "../saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('datanames','models','mygrad','hs','grad_cutoffs'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl)
  source("./RFNR_model_functions.R")
})
# i <- 13+10
foreach(i=1:numdat) %dopar% {
  data <- read_excel("../data/Dati_PERMANENT_reversible_all.xlsx",sheet=i)
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
  
  xlo <- 5000*(5/5-1)
  stepsize <- 5e3
  # check smoothness of grad and correct cutoff (only 1 per peak)
  tmp_plot <- data %>%
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
  tmp_plot <- data %>% mutate(id=1:length(data$time)) %>%
    filter(id%%100==0) %>%
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type)) +
    geom_line(alpha=0.6) +
    geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1)
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  
  # # application
  x <- data$time
  y <- data$current
  ndx <- seq(1, length(x), 60*4) # one measurement every minute

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
    tmp_material <- "original"
  } else {
    tmp_data_name <- datanames[i-13]
    tmp_material <- "MEA 2"
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
    labs(x="time in s")
  # tmp_plot
  ggsave(paste0("../plots_rev_deg/Fitted_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  tmp_plot <- plot(modres,"residuals") +
    labs(x="time in s")
  # tmp_plot
  ggsave(paste0("../plots_residuals_rev_deg/Residuals_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

  tmp_plot <- plot(modres,"res-QQ")
  # tmp_plot
  ggsave(paste0("../plots_residuals_rev_deg/ResQQ_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  tmp_plot <- plot(modres,"res-vs-fitted")
  # tmp_plot
  ggsave(paste0("../plots_residuals_rev_deg/ResFitted_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  tmp_plot <- plot(modres,"res-hist")
  # tmp_plot
  ggsave(paste0("../plots_residuals_rev_deg/ResHist_RFNR_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  pdf(paste0("../plots_residuals_rev_deg/ResACF_RFNR_",datanames[i],".pdf"),width = 8,height = 5)
  plot(res_obj,"res-acf")
  dev.off()
  
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
  alpha <- 0.01
  myqt <- qt(alpha,df = length(modres$y_use) - length(modres$pars),lower.tail = FALSE)
  model_fits[1] <- length(which(ind))
  model_fits[2] <- mean(diff(data$time[which(ind)]))
  model_fits[3] <- c(rep(30,9),100,50,30,15,rep(30,9),100,50,30,15)[i]
  
  model_fits[4] <- exp(summod$par_sig[1,1])
  model_fits[5] <- exp(summod$par_sig[1,1] - summod$par_sig[1,2]*myqt)
  model_fits[6] <- exp(summod$par_sig[1,1] + summod$par_sig[1,2]*qt(alpha,df = length(modres$y_use) - length(modres$pars),lower.tail = FALSE))
  
  model_fits[7:9] <- summod$par_sig[2,c(1)] + c(0,-1,1)*summod$par_sig[2,c(2)]
  model_fits[10:12] <- summod$par_sig[3,c(1)] + c(0,-1,1)*summod$par_sig[3,c(2)]
  model_fits[13:15] <- summod$par_sig[4,c(1)] + c(0,-1,1)*summod$par_sig[4,c(2)]
  
  model_fits[16] <- exp(summod$par_sig[5,1])
  model_fits[17] <- exp(summod$par_sig[5,1] - summod$par_sig[5,2]*myqt)
  model_fits[18] <- exp(summod$par_sig[5,1] + summod$par_sig[5,2]*myqt)
  
  model_fits[19] <- modres$rho_upper / (1+exp(-summod$par_sig[6,1]))
  model_fits[20] <- modres$rho_upper / (1+exp(-summod$par_sig[6,1]+ summod$par_sig[6,2]*myqt))
  model_fits[21] <- modres$rho_upper / (1+exp(-summod$par_sig[6,1]- summod$par_sig[6,2]*myqt))

  my_full_dat <- data.frame(data_tmp,data=tmp_data_name,material=tmp_material,fit=as.numeric(modres$fitted_jumps))
  
  resi <- list(dataname = datanames[i],
                model_fits=model_fits,
                sum_list = list(summary=summod,norm_tests = norm_tests),
                my_full_dat = my_full_dat,
                model_comp = model_comp)
  cat(sprintf('Finished rep %d / %d at %s.\n',i,numdat,Sys.time()))
  saveRDS(resi,paste0("../saved_results_loop/RFNR_results_i",i,".rds"))
  resi
}

data_order <- c(1,3,4,2,5,6,7,8,9,
                10,11,12,13)
data_order <- c(data_order,13+data_order)
name_order <- c(1,4,2,3,5,6,7,8,9,
                10,11,12,13)
name_order <- c(name_order,13+name_order)

model_fits <- matrix(NA,numdat,21)
row.names(model_fits) <- datanames[name_order]
colnames(model_fits) <- c("jumps","mean_time_jumps","rh",
                          "c1","lower","upper",
                          "lam1","lower","upper",
                          "c2","lower","upper",
                          "lam2","lower","upper",
                          "delta","lower","upper",
                          "rho","lower","upper")

sum_list <- vector("list",numdat)
names(sum_list) <- datanames[name_order]


model_comp <-  data.frame(data=NULL,model=NULL,error=NULL,error_measure=NULL,value=NULL)
my_full_dat <- data.frame(NULL)

for (i in 1:numdat) {
  resi <- readRDS(paste0("../saved_results_loop/RFNR_results_i",i,".rds"))
  model_fits[data_order[i],] <- resi$model_fits
  sum_list[[data_order[i]]] <- resi$sum_list
  my_full_dat <- rbind(my_full_dat,
                       resi$my_full_dat)
  model_comp <- rbind(model_comp,
                       resi$model_comp)
}

saveRDS(list(model_fits=model_fits,
             sum_list = sum_list,
             my_full_dat = my_full_dat,
             datanames = datanames,
             data_order = data_order,
             model_comp = model_comp),
        "../saved_results/results_RFNR_model_all_cor_CIs.rds")
    
