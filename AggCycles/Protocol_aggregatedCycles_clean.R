
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # start by reading in the .mat data file and explore + visualize the data # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(parallel,foreach,dplyr,tidyr,readxl,nlme)
source("./RFNR_model_functions_v2.R")


datanames <- c("REF","UPL","LPL","CYC","80","90","100",
               "Outlet","Inlet","Dry-Inlet","80","90","100","90v2")
material <- c(rep("MEA1",7),rep("MEA2",7))
numdat <- length(datanames)
stopifnot(length(material)==numdat)


unlink("./saved_results/log.txt")
n.cores <- parallel::detectCores() # adapt! 
my.cluster <- parallel::makeCluster(n.cores-1, type = "PSOCK", outfile = "./saved_results/log.txt")
doParallel::registerDoParallel(cl = my.cluster)
# foreach::getDoParRegistered()
clusterExport(my.cluster,c('datanames','material'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(dplyr,tidyr,readxl,nlme)
  source("./RFNR_model_functions_v2.R")
})

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # start of loop # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# problems with starting values: 3, 6, 7, 8, 9, 10, 13
# have been resolved 
# idat <- 6
foreach(idat=1:numdat) %dopar% {
  # read in data
  if (idat <= 7) {
    data <- read_xlsx("./data_avg_current/MEA_1_Currents_AST.xlsx",sheet= idat)
  } else {
    data <- read_xlsx("./data_avg_current/MEA_2_Currents_AST.xlsx",sheet= idat-7)
  }
  
  colnames(data) <- c("cycle","hours","current")
  
  # hist(diff(c(0,data$current)))
  # summary(diff(c(0,data$current)))
  # tail(sort(diff(c(0,data$current))),14)

  if (idat < 7) {
    cut_off <- 0.24
  } else if (idat<=8 || idat>=12) {
    cut_off <- 0.13
  } else {
    cut_off <- 0.15
  }
  
  # some exceptions that are not detected by the cut_off
  help_cl_cycles <- NULL
  if (idat==8) {
    help_cl_cycles <- c(200,399,800)
  }
  if (idat==9) {
    help_cl_cycles <- c(402,1003)
  }
  if (idat==10) {
    help_cl_cycles <- c(201,402)
  }
  if (idat ==13) {
    help_cl_cycles <- c(121)
  }
   
  # data %>% filter(cycle > 110,
  #                 cycle < 130) %>%
  #   ggplot(aes(x=cycle,y=current)) +
  #   geom_line()
  
  cleaned_cycle <- data %>% filter(diff(c(0,current)) > cut_off | cycle %in% help_cl_cycles) %>% select(c(cycle))
  # cleaned_cycle$cycle
  cl_cycle_indices <- cleaned_cycle$cycle
  num_clean <- length(cl_cycle_indices)
  # data$hours[cl_cycle_indices]
  time_clean <- data$hours - rep(c(data$hours[cl_cycle_indices]),times=c(cl_cycle_indices[-1],max(data$cycle)+1) - cl_cycle_indices)
  data <- data %>% mutate(time_clean=time_clean,
                          cleaning = rep(1:length(cl_cycle_indices),times=c(cl_cycle_indices[-1],max(data$cycle)+1) - cl_cycle_indices))
  

  # data %>% pivot_longer(c(3:5),names_to = "name",values_to = "value") %>%
  #   ggplot(aes(x=cycle,y=value)) +
  #   geom_line() +
  #   facet_wrap(name~.,scales="free_y") +
  #   geom_vline(xintercept = cl_cycle_indices,linetype=2,col="blue")
  
  data$cleaning <- factor(data$cleaning)
  
  # # # # # # # # # # # # # # # # # # # # # 
  # # # fit model for each cleaning
  # # # # # # # # # # # # # # # # # # # # # 
  
  myfitted <- NULL
  t_rd70 <- rd_t50 <- numeric(num_clean)
  quant_t_rd70 <- quant_rd_t50 <- matrix(c(0),num_clean,2)
  df_coefs <- data.frame(NULL)
  for (i in 1:num_clean) {
    tmp_data <- data %>% filter(cleaning==i)
    res <- RFNR(tmp_data$current,tmp_data$time_clean,jump_direction = "none",log_order = 2,
                init = c(0.5,-5,-3,-1))
    
    mysum <- summary(res,cov="standard",corr=TRUE)
    tmp_fitted <- predict(res)
    myfitted <- c(myfitted,tmp_fitted)

    df_coefs <- rbind(df_coefs,
                      data.frame(dataname = datanames[idat],
                                 material = material[idat],
                                 coef_name = c("lc1","llam1","lc2","llam2"),
                                 cleaning = i,
                                 coef = coef(res),
                                 se = mysum$par_sig[,2],
                                 tvalue = mysum$par_sig[,3],
                                 pvalue = mysum$par_sig[,4],
                                 sigma = res$sigma,
                                 df = length(tmp_data$current)-4,
                                 cor = cor(tmp_data$current,tmp_fitted)
                                 ))
    
    tmp_sigma <- mysum$par_cor*outer(mysum$par_sig[,2],mysum$par_sig[,2])
    M <- 10^5
    sim_values_t_rd70 <- sim_values_rd_t50 <- numeric(M)
    set.seed((1234+idat)^2)
    coefk <- MASS::mvrnorm(M,mu = coef(res),Sigma=tmp_sigma)
    
    # might need to change for some idat
    time_selection <- 50
    rd_selection <- 0.7
    # k <- 1
    for (k in 1:M) {
      uniroot_res <- try(uniroot(function(tmp_t) {(-exp(coefk[k,2])*tmp_t - exp(coefk[k,3])*(1-exp(-exp(coefk[k,4])*tmp_t)))-log(rd_selection)},
                                 interval = c(0,1e3),extendInt = "downX",tol = 1e-2),
                         silent=TRUE)
      if (inherits(uniroot_res,what = "try-error")) {
        uniroot_res$root <- Inf
      }
      sim_values_t_rd70[k] <- uniroot_res$root
      sim_values_rd_t50[k] <- exp(-exp(coefk[k,2])*time_selection - exp(coefk[k,3])*(1-exp(-exp(coefk[k,4])*time_selection)))
    }
    uniroot_res <- try(uniroot(function(tmp_t) {(-exp(coefk[k,2])*tmp_t - exp(coefk[k,3])*(1-exp(-exp(coefk[k,4])*tmp_t)))-log(rd_selection)},
                               interval = c(0,1e3),extendInt = "downX",tol = 1e-2))
    if (inherits(uniroot_res,what = "try-error")) {
      uniroot_res$root <- Inf
    }
    t_rd70[i] <- uniroot_res$root
    quant_t_rd70[i,] <- quantile(sim_values_t_rd70,probs = c(0.05/2,1-0.05/2))
    
    rd_t50[i] <- exp(-exp(coefk[k,2])*time_selection - exp(coefk[k,3])*(1-exp(-exp(coefk[k,4])*time_selection)))
    quant_rd_t50[i,] <- quantile(sim_values_rd_t50,probs = c(0.05/2,1-0.05/2))
    
  }
  df_data <- data %>% mutate(fitted=myfitted) 
  df_data <- df_data %>% group_by(cleaning) %>% mutate(relative_degradation = fitted/max(fitted),
                                                       relative_current = current/max(fitted),
                                                       dataname = datanames[idat],
                                                       material = material[idat]
                                                       )
  df_coefs$cleaning <- factor(df_coefs$cleaning,levels=c(1:num_clean))
  
  df_reldeg <- data.frame(dataname = datanames[idat],
                          material = material[idat],
                          cleaning = factor(c(1:num_clean,1:num_clean),levels=c(1:num_clean)),
                          name = rep(c("t_rd70","rd_t50"),each=num_clean),
                          value = c(t_rd70,rd_t50),
                          lower = c(quant_t_rd70[,1],quant_rd_t50[,1]),
                          upper = c(quant_t_rd70[,2],quant_rd_t50[,2]))
  
  resi <- list(dataname = datanames[idat],
               material = material[idat],
               df_data = df_data,
               df_coefs = df_coefs,
               df_reldeg = df_reldeg)
  cat(sprintf('Finished rep %d / %d at %s.\n',idat,numdat,Sys.time()))
  saveRDS(resi,paste0("./saved_results_loop/agg_cyc_results_i",idat,".rds"))
  resi
}

df_data <- df_coefs <- df_reldeg <- data.frame(NULL)
for (i in 1:numdat) {
  resi <- readRDS(paste0("./saved_results_loop/agg_cyc_results_i",i,".rds"))
  df_data <- rbind(df_data,
                   resi$df_data)
  df_coefs <- rbind(df_coefs,
                   resi$df_coefs)
  df_reldeg <- rbind(df_reldeg,
                   resi$df_reldeg)
}

saveRDS(list(df_data=df_data,
             df_coefs = df_coefs,
             df_reldeg = df_reldeg
             ),
        "./saved_results/results_agg_cycles.rds")
