if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl)


datanames <- c("Pure Air","Ambient air test 1","Ambient air test 2",
               "Pure air sample 2","Ambient air sample 2","Ambient air sample 2 test2",
               "Ambient air sample 2 test3","PM filter sample 2","ALL Filters sample 2")

numtest <- length(datanames)

reslist <- vector("list",numtest)
names(reslist) <- datanames

# report cor, significance
source("./models_reversible.R")

unlink("./saved_results/log2.txt")
# select number of cores eg detectCores()-1
my.cluster <- parallel::makeCluster(5, type = "PSOCK", outfile = "./saved_results/log2.txt")
doParallel::registerDoParallel(cl = my.cluster)
# # # # # # # # # # # # # # # # # # # # # # # # 
# # # start of loop # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # 
 
hs <- c(200,200,200,
        500,200,200,
        200,200,200)
grad_cutoffs <- c(0.0012,0.0008,0.0005,
                  0.0050,0.0030,0.0016,
                  0.0001,0.0032,0.0023)
clusterExport(my.cluster,c('hs','grad_cutoffs','datanames'), envir = environment())
clusterEvalQ(my.cluster, {  
  pacman::p_load(readxl)
  source("./models_reversible.R")
})

?fda::smooth.basis()
# second derivatve

# otherwise
# p-splines

reslist <- foreach(i=1:9) %dopar% {
  i <- 1
  k <- 1
  data <- read_excel("./data/Dati_PERMANENT_reversible.xlsx",sheet=i)
  colnames(data) <- c("time","current")
  # change to s
  data$time <- (data$time-1)*60
  
  # set correct h and cutoff 
  grad <- mygrad(data$current,hs[i])
  grad_cutoff <- grad_cutoffs[i]
  ind <- grad > grad_cutoff
  # only take last time point of positive gradient regions
  for (k in c(1:(length(ind)-1))) {
    if (!is.na(ind[k]*ind[k+1]) &ind[k] & ind[k+1]) {
      ind[k] <- FALSE
    }
  }
  # summary(ind)
  
  xlo <- 1000*0
  stepsize <- 1e4
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
  ggsave(paste0("./plots_rev_deg/JumpsClose_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  # check if all proper peaks are detected and not more or less
  tmp_plot <- data %>% mutate("grad"=grad+0.15,id=1:length(data$time)) %>% 
    filter(id%%100==0) %>% 
    pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
    ggplot(aes(x=time,y=value,color=type)) +
    geom_line(alpha=0.8) + 
    geom_vline(xintercept=data$time[ind],linetype=2,size=0.1) 
  # tmp_plot
  ggsave(paste0("./plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)
  
  diffs <- diff(c(which(ind),length(data$time)+1))
  time_last_jump <- c(data$time[1:(which(ind)[1]-1)],rep(c(data$time)[c(which(ind))],diffs))
  jumps <- rep(c(0,1:sum(ind,na.rm = TRUE)),c(which(ind)[1]-1,diffs))
  data <- data %>% mutate(time_l_jump=time - time_last_jump,jumps=jumps)
  data$jumps <- factor(data$jumps)
  
  # name <- datanames[i]
  loopres <- vector("list",length(model_list))
  for (k in 1:length(model_list)) {
    loopres[[k]] <- tryCatch(model_list[[k]](data,datanames[i]),error=function(error_message) {
      message(paste0("Error in model ",names(model_list)[k]," for ",datanames[i],": ",error_message))
      return(NULL)
    })
    saveRDS(loopres[[k]],paste0("./saved_results_loop/loopres_k",k,"_i",i,".rds"))
    cat(sprintf('Done model %d / %d in rep %d / %d at %s.\n',k,length(model_list),i,9,Sys.time()))
  }
  
  cat(sprintf('Finished rep %d / %d at %s.\n',i,9,Sys.time()))
  list(num_jumps = max(jumps),h=hs[i],grad_cutoff=grad_cutoffs[i],model_res_list=loopres)
}

saveRDS(reslist,"./saved_results/RevDegResults_adapted.rds")

parallel::stopCluster(cl = my.cluster)
