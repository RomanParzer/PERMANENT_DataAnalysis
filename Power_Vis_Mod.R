
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # start by reading in the .mat data file and explore + visualize the data # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(R.matlab,dplyr,tidyr,ggplot2,mgcv,robustbase,lme4)

# read in data
data <- readMat("./PERMANENT_CCMG_REFERENCE2.mat")

exp_current <- ifelse(c(data$Voltage2) < 0.775,
                      mean(data$Current2[data$Voltage2 < 0.775]),
                      mean(data$Current2[data$Voltage2 >= 0.775]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = c(data$Current2),
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# # find time indices of cleanings
tmp_df <- df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage)) %>% filter(avg_current > 0.95) 
tmp_df$cycle

# 11 peaks
cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean,
                              cleaning = rep(1:length(cl_cycle_indices),times=11133*(c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


df_long %>% filter(cycle==242) %>% 
  ggplot(aes(x=time_c,y=current*voltage)) +
  geom_line(alpha = 0.8) +
  labs(x="time in cycle in s")

# still effect of voltage, so rather stay with current?



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of cycles approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# # adding indicator of cycles needs to much vector memory, try to manually adjust for each cycle
# avg_c <- df_long %>% group_by(cycle) %>% summarize(mean_c = mean(log(current))) 
# model_cycles <- lm(resp~time_c+ I(time_c^2)+I(voltage<0.775),data=df_long %>% mutate(resp = log(current)-rep(avg_c$mean_c,11133)))
# model_cycles
# summary(model_cycles)
# 
# tmp_df <- df_long %>% mutate(fitted=exp(model_cycles$fitted.values+rep(avg_c$mean_c,11133)),
#                      resp = current ) 
# tmp_df %>%
#   pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
#   filter(time>=239*11133/4,time<=241*11133/4) %>%
#   ggplot(aes(x=time,y=value,color=type)) +
#   geom_line(alpha=0.8)
# 
# tmp_df %>%
#   filter(time>=239*11133/4,time<=241*11133/4) %>%
#   ggplot(aes(x=fitted,y=resp)) +
#   geom_point(alpha=0.8)
# 
# tmp_df %>%
#   filter(time>=239*11133/4,time<=241*11133/4) %>%
#   ggplot(aes(x=time,y=resp-fitted)) +
#   geom_point(alpha=0.8)

# not a good fit
# try log differences

# instead, do for each cycle:

i <- 1
tmp_data <- df_long %>% filter(cycle==i)
model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<0.775)+time_c:I(voltage<0.775),data=tmp_data)
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

# coef_mat <- array(c(0),dim=c(1200,5,4))

# for (i in 1:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<0.775)+time_c:I(voltage<0.775),data=tmp_data)
#   coef_mat[i,,] <- 
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
  summarize(avg_hi=mean(current[voltage<0.775]),
            avg_lo=mean(current[voltage>=0.775]),
            min = min(current),
            max = max(current),
            Q5_hi = quantile(current[voltage<0.775],probs = 0.5),
            Q25_hi = quantile(current[voltage<0.775],probs = 0.25),
            med_hi = median(current[voltage<0.775]),
            Q75_hi = quantile(current[voltage<0.775],probs = 0.75),
            Q95_hi = quantile(current[voltage<0.775],probs = 0.95),
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

# find methods for multivariate regression, lm can do it
# for now, look at univariate lm:
model_avg_hi <- lm(log(avg_hi) ~ time + time_clean+I(time_clean^2) ,data=df_cycles)
model_avg_hi
# plot(model_avg_hi)
summary(model_avg_hi)
# R2 0.7168

tmp_df <- df_cycles %>% mutate(fitted=exp(model_avg_hi$fitted.values),
                     resp = avg_hi )
cor(tmp_df$resp,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8)

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)

tmp_df %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)

# include random effect for jump after cleaning

model_avg_hi_re <- gam(log(avg_hi) ~ s(cleaning,bs="re") + time+time_clean + I(time_clean^2) ,data=df_cycles)
model_avg_hi_re
# plot(model_avg_hi_re)
summary(model_avg_hi_re)
# way higher R2  0.96

tmp_df <- df_cycles %>% mutate(fitted=exp(model_avg_hi_re$fitted.values),
                     resp = avg_hi ) 
cor(tmp_df$resp,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8)
  # geom_abline(slope=-0.00359,intercept=2.05,col=1)

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)

tmp_df %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)


# try log differences, add higher orders

length(diff(df_cycles$avg_hi))
df_cycles <- df_cycles %>% mutate(log_avg_hi_diff = c(0,diff(log(df_cycles$avg_hi))))

model_avg_hi_diff <- gam(log_avg_hi_diff ~ s(cleaning,bs="re") + time+time_clean ,data=df_cycles)
model_avg_hi_diff
summary(model_avg_hi_diff)
# not very significant

tmp_df <- df_cycles %>% mutate(fitted=exp(model_avg_hi_diff$fitted.values + c(log(df_cycles$avg_hi[1]),log(df_cycles$avg_hi[-1200]))),
                               resp = avg_hi ) 

cor(tmp_df$resp,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8)
# geom_abline(slope=-0.00359,intercept=2.05,col=1)

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)

tmp_df %>%
  ggplot(aes(x=cycle,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)
# high resudials for first cycles after each cleaning, 

# combine first model with lagged model by using lagged variables as predictor

df_cycles <- df_cycles %>% mutate(log_avg_hi_lag = c(mean(log(df_cycles$avg_hi[-1200])),log(df_cycles$avg_hi[-1200])))

model_avg_hi_lag <- gam(log(avg_hi) ~ log_avg_hi_lag +  s(cleaning,bs="re") + time+time_clean+I(time_clean^2) ,data=df_cycles)
model_avg_hi_lag
summary(model_avg_hi_lag)
# way higher R2 0.96

tmp_df <- df_cycles %>% mutate(fitted=exp(model_avg_hi_lag$fitted.values),
                               resp = avg_hi ) 

cor(tmp_df$resp,tmp_df$fitted)
tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8)
# geom_abline(slope=-0.00359,intercept=2.05,col=1)

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(slope=1,intercept = 0,col=2,alpha=0.8)

tmp_df %>%
  ggplot(aes(x=cycle,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8)

# not really better than just RE for cleaning

# multivariate ? 

# lme4, lmer() 
# long data format, one response, interactions, see https://stackoverflow.com/questions/78517026/multivariate-linear-mixed-model
# start with 3:
df_cycles_long <- df_cycles %>% mutate(log_avg_lo_lag = c(mean(log(df_cycles$avg_lo[-1200])),log(df_cycles$avg_lo[-1200])),
                                       log_Q95_hi_lag = c(mean(log(df_cycles$Q95_hi[-1200])),log(df_cycles$Q95_hi[-1200]))) %>%
  pivot_longer(c(log_avg_hi_lag,log_avg_lo_lag,log_Q95_hi_lag),names_to = "variable_lag",values_to = "log_value_lag")
tmp_df <- pivot_longer(df_cycles,c(avg_hi,avg_lo,Q95_hi),names_to = "variable",values_to = "value")
df_cycles_long <- df_cycles_long %>% mutate(value = tmp_df$value,variable = tmp_df$variable)
  

mv_model_lag <-  lmer(log(value) ~ 0 + variable + variable:(log_value_lag + time + time_clean +I(time_clean^2) ) + (1|cleaning),
                  data=df_cycles_long)

mv_model_lag
summary(mv_model_lag)
# estimates fairly similar for the three different features, especially Q95 and avg hi
# way higher R2 

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model_lag)))

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


# same but without lagged variables:
mv_model <-  lmer(log(value) ~ 0 + variable + variable:(time + time_clean + I(time_clean^2)) + (1|cleaning),
                  data=df_cycles_long)

mv_model
summary(mv_model)
# estimates fairly similar for the three different features, especially Q95 and avg hi
# way higher R2 

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

# seems better than using lags

# or how else would you use a model for residuals?
