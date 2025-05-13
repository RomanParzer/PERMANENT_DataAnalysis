if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(foreach,parallel,dplyr,tidyr,ggplot2,readxl,fda)

datanames <- c("Pure Air","Ambient air test 1","Ambient air test 2",
               "Pure air sample 2","Ambient air sample 2","Ambient air sample 2 test2",
               "Ambient air sample 2 test3","PM filter sample 2","ALL Filters sample 2")

numtest <- length(datanames)

reslist <- vector("list",numtest)
names(reslist) <- datanames

# report cor, significance
source("./models_reversible.R")


hs <- c(200,200,200,
        500,200,200,
        200,200,200)
grad_cutoffs <- c(0.0012,0.0008,0.0005,
                  0.0050,0.0030,0.0016,
                  0.0001,0.0032,0.0023)


i <- 1

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

xlo <- 1000*1
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
tmp_plot
# ggsave(paste0("./plots_rev_deg/JumpsClose_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

# check if all proper peaks are detected and not more or less
tmp_plot <- data %>% mutate("grad"=grad+0.15,id=1:length(data$time)) %>% 
  filter(id%%100==0) %>% 
  pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) + 
  geom_vline(xintercept=data$time[ind],linetype=2,size=0.1) 
tmp_plot
# ggsave(paste0("./plots_rev_deg/Jumps_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

diffs <- diff(c(which(ind),length(data$time)+1))
time_last_jump <- c(data$time[1:(which(ind)[1]-1)],rep(c(data$time)[c(which(ind))],diffs))
jumps <- rep(c(0,1:sum(ind,na.rm = TRUE)),c(which(ind)[1]-1,diffs))
data <- data %>% mutate(time_l_jump=time - time_last_jump,jumps=jumps)
data$jumps <- factor(data$jumps)


# monotone smooth

nbasis <- 20
basisobj <- create.bspline.basis(range(data$time),nbasis)
fdParobj <- fdPar(fdobj=basisobj, Lfdobj=2, lambda=1e10)
#  Smooth the data, outputting a list containing various quantities
smooth_current <- smooth.basis(data$time, data$current, fdParobj)
str(smooth_current)
plot(smooth_current)
# hist(diff(data$current))

smoothm_current <- smooth.monotone(data$time, data$current, fdParobj)
str(smoothm_current)
plot(smoothm_current)


data_tmp <- data %>% mutate("smooth"=eval.fd(data$time,smooth_current$fd),
                "smoothm"=eval.fd(data$time,smoothm_current$yhatfd),
                id=1:length(data$time))

model_de <- nls(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                 data=data_tmp,
                 start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4))

model_de_sm <- nls(log(smooth) ~ b0 + b1*time + b3*exp(b4*time),
                   data=data_tmp,
                   start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4))

model_de_smm <- nls(log(smoothm) ~ b0 + b1*time + b3*exp(b4*time),
                   data=data_tmp,
                   start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4))


tmp_plot <- data_tmp %>% 
  mutate("fit"=exp(predict(model_de)),
         "fit_sm"=exp(predict(model_de_sm)),
         "fit_smm"=exp(predict(model_de_smm))) %>%
  filter(id%%100==0) %>% 
  pivot_longer(c(current,smooth,smoothm,fit,fit_sm,fit_smm),names_to = "curves",values_to = "value") %>%
  mutate(type=ifelse(curves%in%c("current","smooth","smoothm"),"data","fitted"),
         data = ifelse(curves%in%c("current","fit"),"current",
                             ifelse(curves%in%c("smooth","fit_sm"),"smooth","smoothm"))) %>%
  ggplot(aes(x=time,y=value,color=data,linetype=type)) +
  geom_line(alpha=0.8)
tmp_plot
# ggsave(paste0("./plots_rev_deg/SmoothMonotone_",datanames[i],".pdf"),tmp_plot,width = 8,height = 5)

# # RE on R2 only

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


# # L E

model_res_LE <- nlmer(current ~ mynonlin(time,time_l_jump,b1,b2,A,R1,lrc1,R2,lrc2)~R2|jumps,
                   data=data,
                   start = c(b1=-1e-5,b2=-1e-7,A=0.6,R1=0.5,lrc1=-1,R2=0.5,lrc2=-6),
                   control=nlmerControl(optimizer = "Nelder_Mead",tolPwrss = 1e-10),
                   verbose=2)
preds <- predict(model_res_LE)
mycor <-  cor(preds,data$current)

coefs <- summary(model_res_LE)$coefficients[,1]
mypred <- function(tmp_t) {
  SS1 <- SSasymp(tmp_t,coefs[3],coefs[4],coefs[5])
  SS2 <- SSasymp(0,0,coefs[6],coefs[7])
  res <- coefs[1]*tmp_t +  SS1 + coefs[2]*0 + SS2
  res
}
pred_start <- preds[1]
deg_times <- numeric(9)
ratios <- 9:1/10
for (k in 1:9) {
  uniroot_res <- uniroot(function(tmp_t) {mypred(tmp_t)-ratios[k]*pred_start},
                         interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
  deg_times[k] <- uniroot_res$root
}

plot(ranef(model_res_LE)$jumps$R2)

tmp_df <- data %>% mutate(fitted=preds,
                          resp = current , 
                          fitted0 = mypred(data$time),
                          id = 1:length(data$time)) 
tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) 
tmp_plot
# ggsave(paste0("./plots_rev_deg/FitCurve_LER2_",name,".pdf"),tmp_plot,width = 8,height = 5)

show_time <- 1000*120
stepsize <- 1e4
tmp_plot <- tmp_df %>% filter(time>show_time,
                              time < show_time + stepsize) %>%
  pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
tmp_plot
# ggsave(paste0("./plots_rev_deg/FitCurveClose_LER2_",name,".pdf"),tmp_plot,width = 8,height = 5)

tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)
tmp_plot
# ggsave(paste0("./plots_rev_deg/Residuals_LER2_",name,".pdf"),tmp_plot,width = 8,height = 5)
list(summary = summary(model_res_LE), 
            cor = mycor,
            deg_times = deg_times,
            pred_start = preds[1])