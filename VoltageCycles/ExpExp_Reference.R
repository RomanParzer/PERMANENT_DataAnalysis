
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # start by reading in the .mat data file and explore + visualize the data # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(R.matlab,dplyr,tidyr,ggplot2,mgcv,robustbase,lme4,knitr)


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
# str(df_long)

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


time_last_volt_change <- numeric(length(data$Time2))
tmplog1 <- c(data$Voltage2)<cutoff
tmplog2 <- c(tmplog1[1],head(tmplog1,-1))
log_vec <- (tmplog1 | tmplog2) & !(tmplog1&tmplog2)

ind <- c(1,which(log_vec))
diffs <- diff(c(ind,length(data$Time2)+1))

time_last_volt_change <- rep(c(data$Time2)[ind],diffs)
df_long <- df_long %>% mutate(time_v_ch=time - time_last_volt_change)


# check time since last voltage change
tmp_gg <- df_long %>% filter(cycle==242) %>% 
  ggplot(aes(x=time_c,y=time_v_ch)) +
  geom_line(alpha = 0.8) +
  labs(x="time in cycle in s")
tmp_gg
# ggsave("./plots/time_last_volt_change_cycle.pdf",tmp_gg,height = 4,width = 7)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # modeling of cycles approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

i <- 1
tmp_data <- df_long %>% filter(cycle==i)
tmp_data$cycle <- factor(tmp_data$cycle)
# model_cycles <- lm(log(current)~time_c+ I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff),data=tmp_data)
model_cycles <- lm(log(current)~time_c+I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff)+
                    time_v_ch + time_v_ch:I(voltage<cutoff),data=tmp_data)
model_cycles
summary(model_cycles)

tmp_df <- tmp_data %>% mutate(fitted=exp(model_cycles$fitted.values),
                              resp = current ) 


cor(tmp_df$resp,tmp_df$fitted)

tmp_gg <- tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)
tmp_gg

# tmp_gg <- tmp_df %>%
#   pivot_longer(c(voltage,current),names_to = "type",values_to = "value") %>%
#   ggplot(aes(x=time,y=value,color=type)) +
#   geom_line(alpha=0.8)
# tmp_gg
# ggsave("./plots/fit_one_cycle.pdf",tmp_gg,height = 4,width = 7)

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2)

tmp_df %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2)


# high correlation due to clustered data
cor(tmp_df$resp[tmp_data$voltage<cutoff],tmp_df$fitted[tmp_data$voltage<cutoff])
cor(tmp_df$resp[tmp_data$voltage>=cutoff],tmp_df$fitted[tmp_data$voltage>=cutoff])


# model fit teilen in hogh voltage, low voltage
tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2) +
  facet_wrap(.~voltage<cutoff,scale="free")

tmp_df %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2) +
  facet_wrap(.~voltage<cutoff,scale="free")

# try gam for more general shape

tmp_data <- tmp_data %>% mutate(indic=factor(voltage<cutoff))
model_cycles <- gam(log(current)~time_c+I(time_c^2)+indic+time_c:indic+ 
                      s(time_v_ch,by=indic),data=tmp_data)
model_cycles
summary(model_cycles)

plot(model_cycles)
tmp_df <- tmp_data %>% mutate(fitted=exp(model_cycles$fitted.values),
                              resp = current ) 


cor(tmp_df$resp,tmp_df$fitted)

tmp_gg <- tmp_df %>%
  pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=type)) +
  geom_line(alpha=0.8)
tmp_gg

cor(tmp_df$resp[tmp_data$voltage<cutoff],tmp_df$fitted[tmp_data$voltage<cutoff])
cor(tmp_df$resp[tmp_data$voltage>=cutoff],tmp_df$fitted[tmp_data$voltage>=cutoff])

# maybe two different models?

tmp_df %>%
  ggplot(aes(x=fitted,y=resp)) +
  geom_point(alpha=0.8) +
  geom_abline(intercept = 0,slope=1,alpha=0.8,col=2) +
  facet_wrap(.~voltage<cutoff,scale="free")

tmp_df %>%
  ggplot(aes(x=time,y=resp-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,alpha=0.8,col=2) +
  facet_wrap(.~voltage<cutoff,scale="free")

# one big model for all might not be possible anymore
# but can repeat this analysis for all cycles

# coef_mat <- array(c(0),dim=c(1200,8,4))
# for (i in 1:1200) {
#   tmp_data <- df_long %>% filter(cycle==i)
#   model_cycles <- lm(log(current)~time_c+I(time_c^2)+I(voltage<cutoff)+time_c:I(voltage<cutoff)+
#                        I(-exp(time_c/5e3)) + time_v_ch + time_v_ch:I(voltage<cutoff),data=tmp_data)
#   
#   coef_mat[i,,] <- summary(model_cycles)$coefficients
# }
# dimnames(coef_mat) <- c(list(cyles=1:1200),dimnames(summary(model_cycles)$coefficients))
# saveRDS(coef_mat,"./saved_results/cycle_models_summaries_expexp.rds")
coef_mat <-  readRDS("./saved_results/cycle_models_summaries_expexp.rds")
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
# df_cycles
df_cycles$cleaning <- factor(df_cycles$cleaning)

# plot average current (and average voltage) for each cycle
tmp_gg <- df_cycles %>% pivot_longer(c(max,avg_hi,avg_lo,Q25_hi,Q75_hi,
                             Q95_hi,avg_volt,min_volt),
                           names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,col=measure)) +
  geom_line(alpha = 0.6)
tmp_gg
# ggsave("./plots/cycle_summaries.pdf",tmp_gg,height = 4,width = 7)

# check source code for optimizer, to estimate parameters and double exp rate at the same time
#? try different rate parameters after each cleaning? for max?
# look at 95Q instead of max

# look at time effect of parameters for degradation 
# try to standardize curves after cleaning to look at rates (eg /volt(jump) or - volt(jump))

# lme4, lmer() 
# long data format, one response, interactions, see https://stackoverflow.com/questions/78517026/multivariate-linear-mixed-model
# start with 3:

df_cycles_long <- df_cycles %>% pivot_longer(c(avg_hi,avg_lo,Q95_hi),names_to = "variable",values_to = "value")
# could use nlmer for exponential term, allow rate to change after each cleaning
# delet time for now
mv_model <-  lmer(log(value) ~ 0 + variable + variable:( time_clean + I(time_clean^2) + time_clean:cleaning) + (1|cleaning:variable),
                  data=df_cycles_long)
# saveRDS(mv_model, "./saved_results/mv_model_expexp.rds")
# mv_model <- readRDS("./saved_results/mv_model_expexp.rds")
summary(mv_model)

# kable(format(summary(mv_model)$coefficients,digits=3),format = "latex",booktabs=TRUE)
# str(mv_model)
# print(mv_model, digits=7, ranef.comp="Var")
# 

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model))) 
cor(tmp_df$value,tmp_df$fitted)


tmp_gg <- tmp_df %>% 
  pivot_longer(c(fitted,value),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~variable,scales = "free_y")
tmp_gg
# ggsave("./plots/fit_mv_model.pdf",tmp_gg,height = 4,width = 7)

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

# these models look good

RE <- ranef(mv_model)
# RE
# sum(RE$cleaning)
# getME(mv_model, "b")
# # getME(mv_model, "Zt")
# 

randeffs <- data.frame(re=RE$`cleaning:variable`$`(Intercept)`,var=rep(c("avg_hi","avg_lo","Q95_hi"),11),ind=rep(1:11,each=3))

randeffs$re
randeffs %>% ggplot(aes(x=ind,y=re,col=var)) +
  geom_line()
# for avg_hi

tmp_jump <- df_cycles_long %>% group_by(variable,cleaning) %>% summarise(maxv=max(value),minv=min(value))
jump <- tmp_jump$maxv[-c(1,12,23)] - tmp_jump$minv[-c(11,22,33)]

coef_df <- data.frame(variable = rep(c("avg_hi","avg_lo","Q95_hi"),10),
                      coef=summary(mv_model)$coefficients[-(1:12),1],
                      cleaning=rep(2:11,each=3),
                      jump=jump)

coef_df %>% 
  ggplot(aes(x=cleaning,y=coef,col=variable)) +
  geom_line() +
  geom_hline(yintercept = 0,col=2,linetype=2) +
  labs(title="difference in rate to cleaning 1")

coef_df %>% 
  ggplot(aes(x=cleaning,y=jump,col=variable)) +
  geom_line() +
  # geom_hline(yintercept = 0,col=2,linetype=2) +
  labs(title="jump after cleaning")

tmp_df <- df_cycles_long %>% mutate(fitted=exp(predict(mv_model))) 

tmp_df %>% group_by(variable,cleaning) %>% 
  mutate(value=value/max(value),fitted=fitted/max(fitted)) %>% 
  pivot_longer(c(fitted,value),names_to = "type",values_to = "value") %>%
  ggplot(aes(x=time_clean,y=value,linetype=type,color=cleaning)) +
  geom_line(alpha=0.8) +
  facet_wrap(.~variable,scales = "free_y")

# could fit model to these scaled avg currents instead, there is a factor-offset
# how does the degradation change over time

# try without time?

