
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # start by reading in the .mat data file and explore + visualize the data # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(R.matlab,dplyr,tidyr,ggplot2,mgcv,robustbase,lme4)


# # # # # # # # 1 Reference base case (upper voltage 0.85, lower voltage 0.7)
# read in data
data <- readMat("./data/PERMANENT_CCMG_REFERENCE2.mat")
cutoff <- (0.7+0.85)/2
exp_current <- ifelse(c(data$Voltage2) < cutoff,
                      mean(data$Current2[data$Voltage2 < cutoff]),
                      mean(data$Current2[data$Voltage2 >= cutoff]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = c(data$Current2),
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# # find time indices of cleanings
tmp_df <- df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage))
ggplot(tmp_df,aes(x=cycle,y=avg_current)) + 
  geom_line() +
  geom_hline(yintercept=0.95,col=2)
c(tmp_df %>% filter(avg_current > 0.95) %>% select(c(cycle)))
# 11 peaks
cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean,
                              cleaning = rep(1:length(cl_cycle_indices),times=11133*(c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

# stay with current as response variable


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of cycles approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

i <- 1
tmp_data <- df_long %>% filter(cycle==i)
tmp_data$cycle <- factor(tmp_data$cycle)
model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
model_cycles
summary(model_cycles)

tmp_df <- tmp_data %>% mutate(fitted=exp(model_cycles$fitted.values),
                              resp = current ) 
cor(tmp_df$resp,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2)

tmp_df %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)



# i <- 1
# tmp_data <- df_long %>% filter(cycle==i)
# X <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
# X <- Matrix(X,sparse=TRUE)
# X <- cbind(X,matrix(c(0),11133,1199))
# y <- log(tmp_data$current)
# for (i in 2:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   smallX <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   smallX <- cbind(smallX,matrix(c(0),11133,1199))
#   smallX[,4+i] <- 1
#   X <- rbind(X,smallX)
#   y <- c(y,log(tmp_data$current))
# }
# 
# colnames(X) <- c("(Intercept)","time_c","I(time_c^2)","I(voltage < cutoff)TRUE","time_c:I(voltage < cutoff)TRUE",paste0("cycle",2:1200))
# 
# XXinv <- solve(crossprod(X))
# betahat <- XXinv%*%crossprod(X,y)
# betahat
# yhat <- X%*%betahat
# k <- ncol(X)
# n <- nrow(X)
# se <- sqrt(sum((y-yhat)^2)/(n-k)*diag(XXinv))
# tval <- betahat/se
# pval <- 2*pt(q = as.numeric(abs(tval)),df=n-k,lower.tail = FALSE)
# 
# restab <- cbind(betahat,se,tval,pval)
# colnames(restab) <- c("Estimate","Std. Error","t-value", "p-value")
# restab
# saveRDS(restab,"./saved_results/cycle_tabmodel_summary.rds")
restab <-  readRDS("./saved_results/cycle_tabmodel_summary.rds")
restab[1:10,]


# coef_mat <- array(c(0),dim=c(1200,5,4))
# for (i in 1:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   coef_mat[i,,] <- summary(model_cycles)$coefficients
# }
# dimnames(coef_mat) <- c(list(cyles=1:1200),dimnames(summary(model_cycles)$coefficients))
# saveRDS(coef_mat,"./saved_results/cycle_models_summaries.rds")
coef_mat <-  readRDS("./saved_results/cycle_models_summaries.rds")
apply(coef_mat,c(2,3),mean)
apply(coef_mat,c(2,3),sd)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of aggregated cycles over time approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

df_cycles <- df_long %>% group_by(cycle) %>% 
  summarize(avg_hi=mean(current[voltage<cutoff]),
            avg_lo=mean(current[voltage>=cutoff]),
            min = min(current),
            max = max(current),
            Q5_hi = quantile(current[voltage<cutoff],probs = 0.5),
            Q25_hi = quantile(current[voltage<cutoff],probs = 0.25),
            med_hi = median(current[voltage<cutoff]),
            Q75_hi = quantile(current[voltage<cutoff],probs = 0.75),
            Q95_hi = quantile(current[voltage<cutoff],probs = 0.95),
            avg_volt = mean(voltage),
            min_volt = min(voltage),
            time=min(time)/60,
            time_clean=min(time_clean)/60,
            cleaning=mean(cleaning))
df_cycles
df_cycles$cleaning <- factor(df_cycles$cleaning)

# plot average current (and average voltage) for each cycle
df_cycles %>% pivot_longer(c(max,avg_hi,avg_lo,Q25_hi,Q75_hi,
                             Q95_hi,avg_volt,min_volt),
                           names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,col=measure)) +
  geom_line(alpha = 0.6)


# lme4, lmer() 
# long data format, one response, interactions, see https://stackoverflow.com/questions/78517026/multivariate-linear-mixed-model
# start with 3:

df_cycles_long <- df_cycles %>% pivot_longer(c(avg_hi,avg_lo,max),names_to = "variable",values_to = "value")
#   mv_model <-  lmer(log(value) ~ 0 + variable + variable:(time + time_clean + I(time_clean^2)) + (1|cleaning),
#                   data=df_cycles_long)
# 
# saveRDS(mv_model, "./saved_results/mv_model.rds")
mv_model <- readRDS("./saved_results/mv_model.rds")
summary(mv_model)

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model)))

cor(tmp_df$value,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,value),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~variable,scales = "free_y")

tmp_df %>%
  ggplot(aes(x=fitted,y=value)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable,scales = "free")

tmp_df %>%
  ggplot(aes(x=cycle,y=value-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable)


# # # # # # # # 2 UPL case (upper voltage shifted to 0.9V)
# read in data
data <- readMat("./data/PERMANENT_CCMG_UPL.mat")
cutoff <- 0.8
exp_current <- ifelse(c(data$Voltage2) < cutoff,
                      mean(data$Current2[data$Voltage2 < cutoff]),
                      mean(data$Current2[data$Voltage2 >= cutoff]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = c(data$Current2),
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# # find time indices of cleanings
tmp_df <- df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage))
ggplot(tmp_df,aes(x=cycle,y=avg_current)) + 
  geom_line() +
  geom_hline(yintercept=0.92,col=2)
c(tmp_df %>% filter(avg_current > 0.92) %>% select(c(cycle)))

# 11 peaks, same as before
cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean,
                              cleaning = rep(1:length(cl_cycle_indices),times=11133*(c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of cycles approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# i <- 1
# tmp_data <- df_long %>% filter(cycle==i)
# X <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
# X <- Matrix(X,sparse=TRUE)
# X <- cbind(X,matrix(c(0),11133,1199))
# y <- log(tmp_data$current)
# for (i in 2:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   smallX <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   smallX <- cbind(smallX,matrix(c(0),11133,1199))
#   smallX[,4+i] <- 1
#   X <- rbind(X,smallX)
#   y <- c(y,log(tmp_data$current))
# }
# 
# colnames(X) <- c("(Intercept)","time_c","I(time_c^2)","I(voltage < cutoff)TRUE","time_c:I(voltage < cutoff)TRUE",paste0("cycle",2:1200))
# 
# XXinv <- solve(crossprod(X))
# betahat <- XXinv%*%crossprod(X,y)
# betahat
# yhat <- X%*%betahat
# k <- ncol(X)
# n <- nrow(X)
# se <- sqrt(sum((y-yhat)^2)/(n-k)*diag(XXinv))
# tval <- betahat/se
# pval <- 2*pt(q = as.numeric(abs(tval)),df=n-k,lower.tail = FALSE)
# 
# restab <- cbind(betahat,se,tval,pval)
# colnames(restab) <- c("Estimate","Std. Error","t-value", "p-value")
# saveRDS(restab,"./saved_results/cycle_tabmodel_summary_upl.rds")
restab <-  readRDS("./saved_results/cycle_tabmodel_summary_upl.rds")
restab[1:10,]

# coef_mat <- array(c(0),dim=c(1200,5,4))
# for (i in 1:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   coef_mat[i,,] <- summary(model_cycles)$coefficients
# }
# dimnames(coef_mat) <- c(list(cyles=1:1200),dimnames(summary(model_cycles)$coefficients))
# saveRDS(coef_mat,"./saved_results/cycle_models_summaries_upl.rds")
coef_mat <-  readRDS("./saved_results/cycle_models_summaries_upl.rds")
apply(coef_mat,c(2,3),mean)
apply(coef_mat,c(2,3),sd)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of aggregated cycles over time approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

df_cycles <- df_long %>% group_by(cycle) %>% 
  summarize(avg_hi=mean(current[voltage<cutoff]),
            avg_lo=mean(current[voltage>=cutoff]),
            min = min(current),
            max = max(current),
            Q5_hi = quantile(current[voltage<cutoff],probs = 0.5),
            Q25_hi = quantile(current[voltage<cutoff],probs = 0.25),
            med_hi = median(current[voltage<cutoff]),
            Q75_hi = quantile(current[voltage<cutoff],probs = 0.75),
            Q95_hi = quantile(current[voltage<cutoff],probs = 0.95),
            avg_volt = mean(voltage),
            min_volt = min(voltage),
            time=min(time)/60,
            time_clean=min(time_clean)/60,
            cleaning=mean(cleaning))
df_cycles
df_cycles$cleaning <- factor(df_cycles$cleaning)

# plot average current (and average voltage) for each cycle
df_cycles %>% pivot_longer(c(max,avg_hi,avg_lo,Q25_hi,Q75_hi,
                             Q95_hi,avg_volt,min_volt),
                           names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,col=measure)) +
  geom_line(alpha = 0.6)


# lme4, lmer() 
# long data format, one response, interactions, see https://stackoverflow.com/questions/78517026/multivariate-linear-mixed-model
# start with 3:

df_cycles_long <- df_cycles %>% pivot_longer(c(avg_hi,avg_lo,max),names_to = "variable",values_to = "value")
# mv_model <-  lmer(log(value) ~ 0 + variable + variable:(time + time_clean + I(time_clean^2)) + (1|cleaning),
#                   data=df_cycles_long)
# 
# saveRDS(mv_model, "./saved_results/mv_model_upl.rds")
mv_model <- readRDS("./saved_results/mv_model_upl.rds")
summary(mv_model)

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model)))

cor(tmp_df$value,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,value),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~variable,scales = "free_y")

tmp_df %>%
  ggplot(aes(x=fitted,y=value)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable,scales = "free")

tmp_df %>%
  ggplot(aes(x=cycle,y=value-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable)


# # # # # # # # 3 LPL case (lower voltage shifted to 0.6V)
# read in data
data <- readMat("./data/PERMANENT_CCMG_LPL.mat")
cutoff <- (0.6+0.85)/2
exp_current <- ifelse(c(data$Voltage2) < cutoff,
                      mean(data$Current2[data$Voltage2 < cutoff]),
                      mean(data$Current2[data$Voltage2 >= cutoff]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = pmax(c(data$Current2),0.001), # avoid 0s
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# # find time indices of cleanings
tmp_df <- df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage))
ggplot(tmp_df,aes(x=cycle,y=avg_current)) + 
  geom_line() +
  geom_hline(yintercept=1.7,col=2)
c(tmp_df %>% filter(avg_current > 1.7) %>% select(c(cycle)))

# 11 peaks, same as before
cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean,
                              cleaning = rep(1:length(cl_cycle_indices),times=11133*(c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of cycles approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# i <- 1
# tmp_data <- df_long %>% filter(cycle==i)
# X <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
# X <- Matrix(X,sparse=TRUE)
# X <- cbind(X,matrix(c(0),11133,1199))
# y <- log(tmp_data$current)
# for (i in 2:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   smallX <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   smallX <- cbind(smallX,matrix(c(0),11133,1199))
#   smallX[,4+i] <- 1
#   X <- rbind(X,smallX)
#   y <- c(y,log(tmp_data$current))
# }
# 
# colnames(X) <- c("(Intercept)","time_c","I(time_c^2)","I(voltage < cutoff)TRUE","time_c:I(voltage < cutoff)TRUE",paste0("cycle",2:1200))
# 
# XXinv <- solve(crossprod(X))
# betahat <- XXinv%*%crossprod(X,y)
# yhat <- X%*%betahat
# k <- ncol(X)
# n <- nrow(X)
# se <- sqrt(sum((y-yhat)^2)/(n-k)*diag(XXinv))
# tval <- betahat/se
# pval <- 2*pt(q = as.numeric(abs(tval)),df=n-k,lower.tail = FALSE)
# 
# restab <- cbind(betahat,se,tval,pval)
# colnames(restab) <- c("Estimate","Std. Error","t-value", "p-value")
# saveRDS(restab,"./saved_results/cycle_tabmodel_summary_lpl.rds")
restab <-  readRDS("./saved_results/cycle_tabmodel_summary_lpl.rds")
restab[1:10,]

# coef_mat <- array(c(0),dim=c(1200,5,4))
# for (i in 1:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   coef_mat[i,,] <- summary(model_cycles)$coefficients
# }
# dimnames(coef_mat) <- c(list(cyles=1:1200),dimnames(summary(model_cycles)$coefficients))
# saveRDS(coef_mat,"./saved_results/cycle_models_summaries_lpl.rds")
coef_mat <-  readRDS("./saved_results/cycle_models_summaries_lpl.rds")
apply(coef_mat,c(2,3),mean)
apply(coef_mat,c(2,3),sd)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of aggregated cycles over time approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

df_cycles <- df_long %>% group_by(cycle) %>% 
  summarize(avg_hi=mean(current[voltage<cutoff]),
            avg_lo=mean(current[voltage>=cutoff]),
            min = min(current),
            max = max(current),
            Q5_hi = quantile(current[voltage<cutoff],probs = 0.5),
            Q25_hi = quantile(current[voltage<cutoff],probs = 0.25),
            med_hi = median(current[voltage<cutoff]),
            Q75_hi = quantile(current[voltage<cutoff],probs = 0.75),
            Q95_hi = quantile(current[voltage<cutoff],probs = 0.95),
            avg_volt = mean(voltage),
            min_volt = min(voltage),
            time=min(time)/60,
            time_clean=min(time_clean)/60,
            cleaning=mean(cleaning))
df_cycles
df_cycles$cleaning <- factor(df_cycles$cleaning)

# plot average current (and average voltage) for each cycle
df_cycles %>% pivot_longer(c(max,avg_hi,avg_lo,Q25_hi,Q75_hi,
                             Q95_hi,avg_volt,min_volt),
                           names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,col=measure)) +
  geom_line(alpha = 0.6)


# lme4, lmer() 
# long data format, one response, interactions, see https://stackoverflow.com/questions/78517026/multivariate-linear-mixed-model
# start with 3:

df_cycles_long <- df_cycles %>% pivot_longer(c(avg_hi,avg_lo,max),names_to = "variable",values_to = "value")
# mv_model <-  lmer(log(value) ~ 0 + variable + variable:(time + time_clean + I(time_clean^2)) + (1|cleaning),
#                   data=df_cycles_long)
# 
# saveRDS(mv_model, "./saved_results/mv_model_lpl.rds")
mv_model <- readRDS("./saved_results/mv_model_lpl.rds")
summary(mv_model)

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model)))

cor(tmp_df$value,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,value),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~variable,scales = "free_y")

tmp_df %>%
  ggplot(aes(x=fitted,y=value)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable,scales = "free")

tmp_df %>%
  ggplot(aes(x=cycle,y=value-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable)


# # # # # # # # 4 cycling case (constant cycling between 0.7 and 0.85)
# read in data
data <- readMat("./data/PERMANENT_CCMG_cycling.mat")
cutoff <- (0.7+0.85)/2
exp_current <- ifelse(c(data$Voltage2) < cutoff,
                      mean(data$Current2[data$Voltage2 < cutoff]),
                      mean(data$Current2[data$Voltage2 >= cutoff]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = pmax(c(data$Current2),0.001),
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# # find time indices of cleanings
tmp_df <- df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage))
ggplot(tmp_df,aes(x=cycle,y=avg_current)) + 
  geom_line() +
  geom_hline(yintercept=0.66,col=2)
c(tmp_df %>% filter(avg_current > 0.66) %>% select(c(cycle)))

# 11 peaks, same as before
cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean,
                              cleaning = rep(1:length(cl_cycle_indices),times=11133*(c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of cycles approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# i <- 1
# tmp_data <- df_long %>% filter(cycle==i)
# X <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
# X <- Matrix(X,sparse=TRUE)
# X <- cbind(X,matrix(c(0),11133,1199))
# y <- log(tmp_data$current)
# for (i in 2:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   smallX <- model.matrix(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   smallX <- cbind(smallX,matrix(c(0),11133,1199))
#   smallX[,4+i] <- 1
#   X <- rbind(X,smallX)
#   y <- c(y,log(tmp_data$current))
# }
# colnames(X) <- c("(Intercept)","time_c","I(time_c^2)","I(voltage < cutoff)TRUE","time_c:I(voltage < cutoff)TRUE",paste0("cycle",2:1200))
# 
# XXinv <- solve(crossprod(X))
# betahat <- XXinv%*%crossprod(X,y)
# yhat <- X%*%betahat
# k <- ncol(X)
# n <- nrow(X)
# se <- sqrt(sum((y-yhat)^2)/(n-k)*diag(XXinv))
# tval <- betahat/se
# pval <- 2*pt(q = as.numeric(abs(tval)),df=n-k,lower.tail = FALSE)
# restab <- cbind(betahat,se,tval,pval)
# colnames(restab) <- c("Estimate","Std. Error","t-value", "p-value")
# saveRDS(restab,"./saved_results/cycle_tabmodel_summary_cycling.rds")
restab <-  readRDS("./saved_results/cycle_tabmodel_summary_cycling.rds")
restab[1:10,]

# coef_mat <- array(c(0),dim=c(1200,5,4))
# for (i in 1:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
#   coef_mat[i,,] <- summary(model_cycles)$coefficients
# }
# dimnames(coef_mat) <- c(list(cyles=1:1200),dimnames(summary(model_cycles)$coefficients))
# saveRDS(coef_mat,"./saved_results/cycle_models_summaries_cycling.rds")
coef_mat <-  readRDS("./saved_results/cycle_models_summaries_cycling.rds")
apply(coef_mat,c(2,3),mean)
apply(coef_mat,c(2,3),sd)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of aggregated cycles over time approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

df_cycles <- df_long %>% group_by(cycle) %>% 
  summarize(avg_hi=mean(current[voltage<cutoff]),
            avg_lo=mean(current[voltage>=cutoff]),
            min = min(current),
            max = max(current),
            Q5_hi = quantile(current[voltage<cutoff],probs = 0.5),
            Q25_hi = quantile(current[voltage<cutoff],probs = 0.25),
            med_hi = median(current[voltage<cutoff]),
            Q75_hi = quantile(current[voltage<cutoff],probs = 0.75),
            Q95_hi = quantile(current[voltage<cutoff],probs = 0.95),
            avg_volt = mean(voltage),
            min_volt = min(voltage),
            time=min(time)/60,
            time_clean=min(time_clean)/60,
            cleaning=mean(cleaning))
df_cycles
df_cycles$cleaning <- factor(df_cycles$cleaning)

# plot average current (and average voltage) for each cycle
df_cycles %>% pivot_longer(c(max,avg_hi,avg_lo,Q25_hi,Q75_hi,
                             Q95_hi,avg_volt,min_volt),
                           names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,col=measure)) +
  geom_line(alpha = 0.6)


# lme4, lmer() 
# long data format, one response, interactions, see https://stackoverflow.com/questions/78517026/multivariate-linear-mixed-model
# start with 3:

df_cycles_long <- df_cycles %>% pivot_longer(c(avg_hi,avg_lo,max),names_to = "variable",values_to = "value")
# mv_model <-  lmer(log(value) ~ 0 + variable + variable:(time + time_clean + I(time_clean^2)) + (1|cleaning),
#                   data=df_cycles_long)
# 
# saveRDS(mv_model, "./saved_results/mv_model_cycling.rds")
mv_model <- readRDS("./saved_results/mv_model_cycling.rds")
summary(mv_model)

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model)))

cor(tmp_df$value,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,value),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~variable,scales = "free_y")

tmp_df %>%
  ggplot(aes(x=fitted,y=value)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable,scales = "free")

tmp_df %>%
  ggplot(aes(x=cycle,y=value-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)+
  facet_wrap(.~variable)


