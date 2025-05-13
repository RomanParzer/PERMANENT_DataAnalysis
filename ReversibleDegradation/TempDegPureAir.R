if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(dplyr,tidyr,ggplot2,readxl,robustbase,mgcv,lme4)


# # # # Purified Air 
data <- read_excel("./data/Dati_PERMANENT_reversible.xlsx",sheet=1)
# str(data)

colnames(data) <- c("time","current")
# change to s
data$time <- (data$time-1)*60

# # double exponential
# tmp_func <- function(lam) {
#   model_pure <- lm(log(current) ~ I(exp(-time/lam)) + time + I(time^2),data=data)
#   summary(model_pure)$adj.r.squared
# }
# opt_res_pure <- optimize(tmp_func,interval=c(1000,100000),maximum = TRUE)
# 
# opt_res_pure$maximum
# opt_res_pure$objective
# 
# model_pure <- lm(log(current) ~ I(exp(-time/opt_res_pure$maximum)) + time + I(time^2),data=data)
# 
# # try to optimize all parameters with optim
# # model_pure <- lm(log(current) ~  time + I(time^2),data=data)
# summary(model_pure)
# 
# tmp_df <- data %>% mutate(fitted=exp(model_pure$fitted.values),
#                               resp = current , 
#                           id = 1:length(data$time)) 
# cor(tmp_df$resp,tmp_df$fitted)
# tmp_gg <- tmp_df %>% filter(id%%100==0) %>%
#   pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
#   ggplot(aes(x=time,y=value,color=type)) +
#   geom_line(alpha=0.8)
# tmp_gg
# # ggsave("./plots/fit_temp_deg_pure.pdf",tmp_gg,height = 4,width = 7)

# same but one estimation
model_pure_all <- nls(log(current) ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                      data=data,
                      start = c(b0=-2,b1=-1e-4,b2=1e-8,b3=1,b4=-1e-5))
model_pure_all <- nlrob(log(current) ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                      data=data,
                      start = c(b0=-2,b1=-1e-4,b2=1e-8,b3=1,b4=-1e-5))
model_pure_all <- nlrob(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                        data=data,
                        start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-5))
model_pure_all
summary(model_pure_all)
# same results apart from doub exp rate
tmp_df <- data %>% mutate(fitted=exp(predict(model_pure_all)),
                          resp = current , 
                          id = 1:length(data$time)) 

cor(tmp_df$resp,tmp_df$fitted)
# same cor
tmp_gg <- tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)
tmp_gg
# ggsave("./plots/fit_temp_deg_pure.pdf",tmp_gg,height = 4,width = 7)


# lin + exp
model_pure_le <- nls(current ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                      data=data,
                      start = c(b0=0.3,b1=-1e-6,b2=1e-12,b3=0.2,b4=-1e-6),
                     algorithm="port")
model_pure_le <- nlrob(current ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                     data=data,
                     start = c(b0=0.3,b1=-1e-6,b2=1e-12,b3=0.2,b4=-1e-6),
                     algorithm="port")
model_pure_le <- nlrob(current ~ b0 + b1*time +  b3*exp(b4*time),
                       data=data,
                       start = c(b0=0.1,b1=-1e-7,b3=0.02,b4=-2e-6),
                       algorithm="port",trace=TRUE)
model_pure_le
summary(model_pure_le)
# same results apart from doub exp rate

tmp_df <- data %>% mutate(fitted=predict(model_pure_le),
                          resp = current , 
                          id = 1:length(data$time)) 

cor(tmp_df$resp,tmp_df$fitted)
# slightly higher cor
tmp_gg <- tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)
tmp_gg
# ggsave("./plots/fit_temp_deg_pure_le.pdf",tmp_gg,height = 4,width = 7)


# # # # Ambient air 1
data <- read_excel("./data/Dati_PERMANENT_reversible.xlsx",sheet=2)
# str(data)

colnames(data) <- c("time","current")
# change to s
data$time <- (data$time-1)*60


model_amb1 <- nls(log(current) ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                  data=data,
                  start = c(b0=-2,b1=-1e-4,b2=1e-8,b3=1,b4=-1e-5))
model_amb1 <- nlrob(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                  data=data,
                  start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-5))
summary(model_amb1)

tmp_df <- data %>% mutate(fitted=exp(predict(model_amb1)),
                          resp = current , 
                          id = 1:length(data$time)) 
cor(tmp_df$resp,tmp_df$fitted)
tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)

# lin + exp
model_amb1_le <- nls(current ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                     data=data,
                     start = c(b0=0.3,b1=-1e-6,b2=1e-12,b3=0.2,b4=-1e-6),
                     algorithm="port")
model_amb1_le <- nlrob(current ~ b0 + b1*time + b3*exp(b4*time),
                     data=data,
                     start = c(b0=0.3,b1=-1e-6,b3=0.2,b4=-1e-6),
                     algorithm="port")
summary(model_amb1_le)

tmp_df <- data %>% mutate(fitted=predict(model_amb1_le),
                          resp = current , 
                          id = 1:length(data$time)) 

cor(tmp_df$resp,tmp_df$fitted)
# slightly lower cor
tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)


# # # # Ambient air 2
data <- read_excel("./data/Dati_PERMANENT_reversible.xlsx",sheet=3)
# str(data)

colnames(data) <- c("time","current")
# change to s
data$time <- (data$time-1)*60

# try gam or manual residual modeling on amb2

mygrad <- function(x,maxh) {
  myind <- (maxh+1):(length(x)-maxh)
  grad <- sapply(myind,function(j) sum(x[j+1:maxh] - x[j-1:maxh])/(2*maxh))
  return(c(rep(NA,maxh),
           grad,
           rep(NA,maxh)))
}

# set correct h and cutoff 
grad <- mygrad(data$current,100)
grad_cutoff <- 0.0008

ind <- grad > grad_cutoff

# only take first time point of positive gradient regions
for (i in c(length(ind):2)) {
  if (!is.na(ind[i]*ind[i-1]) &ind[i] & ind[i-1]) {
    ind[i] <- FALSE
  }
}
summary(ind)

xlo <- 45000
stepsize <- 1e4
# check smoothness of grad and correct cutoff (only 1 per peak)
data %>% mutate("grad"=grad) %>% 
  filter(time>xlo,
         time<xlo+stepsize) %>% 
  pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) + 
  geom_vline(xintercept=data$time[ind],linetype=2,size=0.1) +
  coord_cartesian(xlim=c(xlo,xlo+stepsize)) +
  geom_hline(yintercept = grad_cutoff,linetype=2)


# check if all proper peaks are detected and not more or less
data %>% mutate("grad"=grad+0.15,id=1:length(data$time)) %>% 
  filter(id%%100==0) %>% 
  pivot_longer(c(current,grad),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) + 
  geom_vline(xintercept=data$time[ind],linetype=2,size=0.1) 

# add time since last jump as variable

# tmplog1 <- c(data$Voltage2)<cutoff
# tmplog2 <- c(tmplog1[1],head(tmplog1,-1))
# log_vec <- (tmplog1 | tmplog2) & !(tmplog1&tmplog2)
# ind <- c(1,which(log_vec))
diffs <- diff(c(which(ind),length(data$time)+1))
time_last_jump <- c(data$time[1:(which(ind)[1]-1)],rep(c(data$time)[c(which(ind))],diffs))
jumps <- rep(c(0,1:sum(ind,na.rm = TRUE)),c(which(ind)[1]-1,diffs))
data <- data %>% mutate(time_l_jump=time - time_last_jump,jumps=jumps)
data$jumps <- factor(data$jumps)


data %>% mutate(id=1:length(data$time)) %>% 
  filter(id%%100==0) %>% 
  ggplot(aes(x=time,y=time_l_jump)) +
  geom_line(alpha=0.8) + 
  geom_vline(xintercept=data$time[ind],linetype=2,size=0.1) 

data %>% mutate(id=1:length(data$time)) %>% 
  filter(id%%100==0) %>% 
  ggplot(aes(x=time,y=jumps)) +
  geom_line(alpha=0.8) + 
  geom_vline(xintercept=data$time[ind],linetype=2,size=0.1) 


# first try with nls, one height for all jumps doesnt work well

model_amb2 <- nls(log(current) ~ b0 + b1*time + b2*exp(b3*time) + b4*time_last_jump + b5*exp(b6*time_last_jump),
                  data=data,
                  start = c(b0=-3.4,b1=-7e-5,b2=1.4,b3=-1.6e-5,b4=6.3e-5,b5=0.19,b6=-2e-5),
                  trace=TRUE)


# therefore, add RF depending on jumps

# input1 <- data$time[1:100]
# input2 <- data$time[1:100]
# A <- -1
# R2 <- R1 <- 1
# lrc2 <- lrc1 <- -5
# b1 <- b2 <- 3
mynonlin <- function(input1,input2,b1,b2,A,R1,lrc1,R2,lrc2) {
  SS1 <- SSasymp(input1,A,R1,lrc1)
  A2 <- 0
  SS2 <- SSasymp(input2,A2,R2,lrc2)
  stopifnot(length(input1)==length(input2))
  res <- b1*input1 +  SS1 + b2*input2 + SS2
  grad <- array(c(0),dim=c(length(input1),7),dimnames=list(obs=NULL,vars=c("b1","b2","A","R1","lrc1","R2","lrc2")))
  grad[,1] <- input1
  grad[,2] <- input2
  grad[,3:5] <- attr(SS1,"gradient")
  grad[,6:7] <- attr(SS2,"gradient")[,2:3]
  attr(res,"gradient") <- grad
  return(res)
}


model_amb2 <- nlmer(log(current) ~ mynonlin(time,time_l_jump,b1,b2,A,R1,lrc1,R2,lrc2)~A|jumps,
                  data=data,
                  start = c(b1=-1e-5,b2=-1e-5,A=-1,R1=0.1,lrc1=-5,R2=0.1,lrc2=-5),
                  verbose=2)

# or try nls with manual matrix for interaction of intercept to jumps

# saveRDS(model_amb2,"./saved_results/model_amb2_nlmer.rds")
model_amb2 <- readRDS("./saved_results/model_amb2_nlmer.rds")

# undebug(nlmer)

re_amb2 <- ranef(model_amb2)
str(re_amb2)

summary(attr(re_amb2$jumps,"postVar")[1,1,])
plot(re_amb2$jumps$A)

summary(model_amb2)
tmp_df <- data %>% mutate(fitted=exp(predict(model_amb2)),
                          resp = current , 
                          id = 1:length(data$time)) 
cor(tmp_df$resp,tmp_df$fitted)

tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) 

show_time <- 1000*20
stepsize <- 1e4
tmp_df %>% filter(time>show_time,
                  time < show_time + stepsize) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8) +
  geom_vline(xintercept=data$time[ind&data$time>show_time & data$time<show_time+stepsize],linetype=2,size=0.1)

# maybe less smooth gradient?

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2)

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=log(fitted),y=log(resp))) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2)

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=time,y=log(resp)-log(fitted))) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)








# # # # #. exp model using nls / nlrob
model_amb2 <- nls(log(current) ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                  data=data,
                  start = c(b0=-2,b1=-1e-4,b2=1e-8,b3=1,b4=-1e-5))

model_amb2 <- nlrob(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                  data=data,
                  start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-5))
summary(model_amb2)

tmp_df <- data %>% mutate(fitted=exp(predict(model_amb2)),
                          resp = current , 
                          id = 1:length(data$time)) 
cor(tmp_df$resp,tmp_df$fitted)

tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)

# lin + exp
model_amb2_le <- nls(current ~ b0 + b1*time + b2*I(time^2)+ b3*exp(b4*time),
                     data=data,
                     start = c(b0=0.3,b1=-1e-6,b2=1e-12,b3=0.2,b4=-1e-6),
                     algorithm="port")

model_amb2_le <- nlrob(current ~ b0 + b1*time + b3*exp(b4*time),
                     data=data,
                     start = c(b0=0.3,b1=-1e-6,b3=0.2,b4=-1e-6),
                     algorithm="port")

summary(model_amb2_le)
tmp_df <- data %>% mutate(fitted=predict(model_amb2_le),
                          resp = current , 
                          id = 1:length(data$time)) 

cor(tmp_df$resp,tmp_df$fitted)
# lower cor
tmp_gg <- tmp_df %>% filter(id%%100==0) %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)
tmp_gg

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2)

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=log(fitted),y=log(resp))) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2)

tmp_df %>% filter(id%%100==0) %>%
  ggplot(aes(x=time,y=log(resp)-log(fitted))) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)

# ggsave("./plots/fit_temp_deg_amb2.pdf",tmp_gg,height = 4,width = 7)

model_amb2 <- nlrob(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                    data=data,
                    start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-5))
summary(model_amb2)


# common table TBD

coefs <- cbind(model_pure_all$coefficients,model_amb1$coefficients,model_amb2$coefficients)

require(knitr)
colnames(coefs) <- c("pure","ambient 1","ambient 2")
row.names(coefs) <- c("Intercept","time","Exp_time","Exp_rate")
coefs
kable(format(coefs,digits=3),format = "latex",booktabs=TRUE)
