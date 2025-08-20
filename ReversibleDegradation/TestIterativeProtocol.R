if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl,fda,gnm)

datanames <- c("Pure Air","Ambient air test 1","Ambient air test 2",
               "Pure air sample 2","Ambient air sample 2","Ambient air sample 2 test2",
               "Ambient air sample 2 test3","PM filter sample 2","ALL Filters sample 2")

numtest <- length(datanames)

reslist <- vector("list",numtest)
names(reslist) <- datanames

# report cor, significance
source("./models_reversible.R")


hs <- c(50,200,200,
        500,200,200,
        200,200,200)
grad_cutoffs <- c(0.0012,0.0008,0.0005,
                  0.0050,0.0030,0.0016,
                  0.0001,0.0032,0.0023)


i <- 1

data <- read_excel("../data/Dati_PERMANENT_reversible.xlsx",sheet=i)
colnames(data) <- c("time","current")
# change to s
data$time <- (data$time-1)*60

# set correct h and cutoff 
grad <- mygrad(data$current,hs[i])
grad_cutoff <- grad_cutoffs[i]
ind <- grad > grad_cutoff
# take first time point of positive gradient regions
for (k in c(length(ind):2)) {
  if (!is.na(ind[k]*ind[k-1]) & ind[k] & ind[k-1]) {
    ind[k] <- FALSE
  }
}
# summary(ind)

xlo <- 1000*0 
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
tmp_plot
# ggsave(paste0("./plots_rev_deg/JumpsClose_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

# check if all proper peaks are detected and not more or less
tmp_plot <- data %>% mutate("grad"=grad+0.15,id=1:length(data$time)) %>% 
  filter(id%%100==0) %>% 
  pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) + 
  geom_vline(xintercept=data$time[ind],linetype=2,linewidth=0.1) 
tmp_plot
# ggsave(paste0("./plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

diffs <- diff(c(which(ind),length(data$time)+1))
time_last_jump <- c(data$time[1:(which(ind)[1]-1)],rep(c(data$time)[c(which(ind))],diffs))
jumps <- rep(c(0,1:sum(ind,na.rm = TRUE)),c(which(ind)[1]-1,diffs))
data <- data %>% mutate(time_l_jump=time - time_last_jump,jumps=jumps)
data$jumps <- factor(data$jumps)

data_jumps <- data %>% group_by(jumps) %>% summarize(max = max(current),
                                                     mean = mean(current),
                                                     med = median(current),
                                                     jump_start_time = min(time),
                                                     peak_time = time[which.max(current)],
                                                     min = min(current[time>peak_time]),
                                                     min_time = time[time>peak_time][which.min(current[time>peak_time])],
                                                     jump_duration = max(time) - min(time),
                                                     jump_height = max - head(current,1))

data_jumps

hist(data_jumps$jump_height)
hist(data_jumps$jump_duration)


model_baseline <- robustbase::nlrob(log(min) ~ b0 + b1*min_time + b3*exp(b4*min_time),
                      data=data_jumps,
                      start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4))

summary(model_baseline)


data_tmp <- data %>% 
  mutate(id=1:length(data$time),
         type="current") %>%
  filter(id%%100==0)
  

tmp_plot <- data_jumps %>% 
  mutate("fit"=exp(predict(model_baseline))) %>%
  pivot_longer(c(min,fit),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=min_time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  geom_line(data=data_tmp,aes(x=time,y=current),alpha=0.8)
tmp_plot
# ggsave(paste0("../plots_rev_deg/MinBaseline_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)


str(model_baseline)
data_jumps %>% 
  mutate("fit"=exp(predict(model_baseline)),
         "weight" = model_baseline$rweights) %>%
  ggplot(aes(x=min_time,y=min-fit,col=weight)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)
# ggsave(paste0("../plots_rev_deg/MinBaselineRsiduals_",datanames[i],".pdf"),width = 8,height = 5)

# quantile regression?

# find more complicated nls model without random effects,
# maybe iterative 

# # # one model for each jump period

data <- mutate(data,"res_current"=current/ exp(predict(model_baseline,newdata = data.frame(min_time=data$time))))



my_jump_function <- function(time,jump_info,b0,b1,b2,b3,a0,a1,a2,a3) {
  z <- (time - a2)/a3
  ifelse( time < jump_info$peak_time,
          a0 + (a1)/(1+exp(-z)),
          b0 + b1*(time-jump_info$peak_time) + b2*exp(b3*(time-jump_info$peak_time))
  )
}

njumps <- max(as.numeric(data_jumps$jumps)) + 1

summaries <- vector("list",njumps)
j <- 14
for (j in 0:(njumps-1)) {
  data_one_jump <- filter(data,jumps==as.character(j))
  jump_info <- filter(data_jumps,jumps==as.character(j))
  
  one_jump_model <- try(nls(log(res_current) ~ my_jump_function(time,jump_info,b0,b1,b2,b3,a0,a1,a2,a3),
                          data=data_one_jump,
                          start = c(b0=log(max(data_one_jump$res_current)*5.5/10),
                                    b1=-3e-5,
                                    b2=(max(data_one_jump$res_current)-tail(data_one_jump$res_current,1))*6.5/10,
                                    b3=-1.5e-2,
                                    a0=log(data_one_jump$res_current[1]),
                                    a1=jump_info$jump_height*8,
                                    a2=(jump_info$jump_start_time+jump_info$peak_time)/2+3,
                                    a3=(jump_info$peak_time-jump_info$jump_start_time)/18),
                        trace=FALSE,
                        algorithm = "port"
  ))
  
  
  # my_jump_function(data_one_jump$time,jump_info,
  #                  log(max(data_one_jump$res_current)*5.5/10),
  #                  -3e-5,
  #                  (max(data_one_jump$res_current)-tail(data_one_jump$res_current,1))*6.5/10,
  #                  -1.5e-2,
  #                  log(data_one_jump$res_current[1]),
  #                  jump_info$jump_height*8,
  #                  (jump_info$jump_start_time+jump_info$peak_time)/2+3,
  #                  (jump_info$peak_time-jump_info$jump_start_time)/18)

  if (!("try-error" == class(one_jump_model))) {
    data_one_jump %>% mutate("fit" = exp(predict(one_jump_model))) %>%
      # filter(time<24700) %>%
      pivot_longer(c(res_current,fit),names_to = "type",values_to = "value") %>%
      ggplot(aes(x=time,y=value,color=type)) +
      geom_line() +
      labs(title=paste0("Fit for jump j=",j))
    # ggsave(paste0("../plots_rev_deg/JumpModelClose_",datanames[i],".pdf"),width = 8,height = 5)
    summaries[[j+1]] <- summary(one_jump_model)
  }
  message(sprintf("Finished jump %d / %d!",j+1,njumps))
}

head(str(summaries))

jump_coefs <- matrix(NA,njumps,8)
for (j in 1:njumps) {
  if (length(summaries[[j]])!=0) {
    jump_coefs[j,] <- summaries[[j]]$coefficients[,1]
  }
}
j <- 3
colnames(jump_coefs) <- names(summaries[[j]]$coefficients[,1])
data_jumps <- cbind(data_jumps,jump_coefs)


data_jumps %>% 
  pivot_longer(c(jump_start_time,mean,jump_duration,jump_height),names_to = "covariates",values_to = "cov_values") %>%
  pivot_longer(c(b0,b1,b2,b3,a0,a1,a2,a3),names_to = "parameters",values_to = "par_values") %>%
ggplot(aes(x=cov_values,y=par_values)) +
  geom_point() +
  facet_wrap(covariates~parameters,scales="free")
# ggsave(paste0("../plots_rev_deg/JumpParametersCovariates_",datanames[i],".pdf"),width = 12,height = 10)



# # # one big model, no random effects; difficult with nls, try gnm


my_nonlin <- function(resp, conc){
  list(predictors = list(Vm = substitute(conc), K = 1),
       variables = list(substitute(resp), substitute(conc)),
       term = function(predictors, variables) {
         pred <- paste("(", predictors[1], "/(", predictors[2],
                       " + ", variables[2], "))", sep = "")
         pred <- paste("(", variables[1], " - ", pred, ")/sqrt(",
                       pred, ")", sep = "")
       })
}
class(my_nonlin) <- "nonlin"

## use to fitted weighted Michaelis-Menten model
Treated <- Puromycin[Puromycin$state == "treated", ]
Pur.wt.2 <- gnm( ~ -1 + my_nonlin(rate, conc), data = Treated,
                 start = c(Vm = 200, K = 0.1), verbose = FALSE)
Pur.wt.2



myExp <- function(arg){
  list(predictors = list(b0 = 1, b1 = 1),
       variables = list(substitute(arg)),
       term = function(predictors, variables) {
         pred <- paste(predictors[1], " * exp(", predictors[2],
                       " * ", variables[1], ")", sep = "")
         pred
       })
}
class(myExp) <- "nonlin"


x <- 1:100
y <- 2*exp( -x/20 +  3*exp(- x / 10))

mydata <- data.frame(current=y,time=x)
set.seed(4)
exp1 <- gnm(log(current) ~ time + myExp(time), 
            data=mydata,
            start = c(b0 = 1, b1 = 0.1),
            verbose = FALSE)

exp1
time <- mydata$time
tmp <- myExp(time)
tmp$predictors
tmp$variables
tmp$term(names(tmp$predictors),tmp$variables)

yhat <- exp(predict(exp1))
plot(x,y)
plot(x,yhat)
plot(x,y-yhat)



# # # old RE model

model_res <- nlmer(log(current) ~ mynonlin(time,time_l_jump,b1,b2,A,R1,lrc1,R2,lrc2)~R2|jumps,
                   data=data,
                   start = c(b1=-1e-5,b2=-1e-5,A=-2,R1=1,lrc1=-5,R2=1,lrc2=-5),
                   control=nlmerControl(optimizer = "Nelder_Mead",tolPwrss = 1e-10),
                   verbose=2)
name <- datanames[i]
preds <- exp(predict(model_res))
mycor <-  cor(preds,data$current)

pred_start <- preds[1]
deg_times <- numeric(9)
ratios <- 9:1/10

coefs <- summary(model_res)$coefficients[,1]
mypred <- function(tmp_t) {
  SS1 <- SSasymp(tmp_t,coefs[3],coefs[4],coefs[5])
  SS2 <- SSasymp(0,0,coefs[6],coefs[7])
  res <- coefs[1]*tmp_t +  SS1 + coefs[2]*0 + SS2
  exp(res)
}
for (k in 1:9) {
  uniroot_res <- uniroot(function(tmp_t) {mypred(tmp_t)-ratios[k]*pred_start},
                         interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
  deg_times[k] <- uniroot_res$root
}

plot(ranef(model_res)$jumps$R2)
tmp_df <- data %>% mutate(fitted=preds,
                          resp = current , 
                          fitted0 = mypred(data$time),
                          id = 1:length(data$time)) 
tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) 
tmp_plot
# ggsave(paste0("./plots_rev_deg/FitCurve_DER2_",name,".pdf"),tmp_plot,width = 8,height = 5)

show_time <- 1000*120
stepsize <- 1e4
tmp_plot <- tmp_df %>% filter(time>show_time,
                              time < show_time + stepsize) %>%
  pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
tmp_plot
# ggsave(paste0("./plots_rev_deg/FitCurveClose_DER2_",name,".pdf"),tmp_plot,width = 8,height = 5)

tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)
tmp_plot
# ggsave(paste0("./plots_rev_deg/Residuals_DER2_",name,".pdf"),tmp_plot,width = 8,height = 5)

list(summary = summary(model_res), 
            cor = mycor,
            deg_times = deg_times,
            pred_start = preds[1])

