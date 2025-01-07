
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # start by reading in the .mat data file and explore + visualize the data # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(R.matlab,dplyr,tidyr,ggplot2,mgcv,robustbase)

# read in data
data <- readMat("./PERMANENT_CCMG_REFERENCE2.mat")
str(data)


exp_current <- ifelse(c(data$Voltage2) < 0.775,
                      mean(data$Current2[data$Voltage2 < 0.775]),
                      mean(data$Current2[data$Voltage2 >= 0.775]))
df_long <- data.frame(time=c(data$Time2), 
                      voltage = c(data$Voltage2),
                      current = c(data$Current2),
                      exp_current = exp_current,
                      cycle = rep(1:1200,each=11133),
                      time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# # find time indicesof cleanings
tmp_df <- df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage)) %>% filter(avg_current > 0.95) 
tmp_df$cycle

# 11 peaks
cl_cycle_indices <- c(1,61,121,241,361,481,601,721,841,961,1081)
time_clean <- c(scale(data$Time2,scale=FALSE,
                      center=rep(apply(data$Time2[,cl_cycle_indices],2,min),times=c(cl_cycle_indices[-1],1201) - cl_cycle_indices)))

df_long <- df_long %>% mutate(time_clean=time_clean)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# time kump between cycles
max(data$Time2[,1])
min(data$Time2[,2])

# time jump while cleaning
max(data$Time2[,240])
min(data$Time2[,241])
# always the same 240 s

# time in seconds, 
# roughly 4 measurement each second
# 1 cycle roughly 46 min (each column 11133 measurements), 
# 1200 cycles, time continues from column to column
# ca 1000 h in total
max(data$Time2[,1])/60
max(data$Time2)/60/60

# in each cycle, the following voltage is applied
# select cycle 1,...,1200
j <- 1
plot(data$Time2[,j]/60,data$Voltage2[,j],type="l",xlab="time in min",ylab="voltage")
# and the resulting current in the system is measured:
plot(data$Time2[,j]/60,data$Current2[,j],type="l",xlab="time in min",ylab="current")

# we can fix one timepoint and plot the current over the cycles
# on 1st peak
t <- 60
summary(c(data$Voltage[t+c(-50:50),]))
plot(data$Current2[t,],type="l",xlab="cycle",ylab="current")
# on 1st low
t <- 190
summary(c(data$Voltage[t+c(-50:50),]))
plot(data$Current2[t,],type="l",xlab="cycle",ylab="current")

# similar
# # on 12th peak
# t <- 4900
# summary(c(data$Voltage[t+c(-50:50),]))
# plot(data$Current2[t,],type="l",xlab="cycle",ylab="current")
# # on low afterwards
# t <- 5200
# summary(c(data$Voltage[t+c(-50:50),]))
# plot(data$Current2[t,],type="l",xlab="cycle",ylab="current")
# 
# # on last peak
# t <- 10475
# summary(c(data$Voltage[t+c(-50:50),]))
# plot(data$Current2[t,],type="l",xlab="cycle",ylab="current")
# # on low afterwards
# t <- 11000
# summary(c(data$Voltage[t+c(-50:50),]))
# plot(data$Current2[t,],type="l",xlab="cycle",ylab="current")


# plot voltage + current over one cycle
df_long %>% filter(cycle==242) %>% 
  pivot_longer(c(voltage,current),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time_c,y=value,color=measure)) +
  geom_line(alpha = 0.8) +
  labs(y=" ",x="time in cycle in s") +
  geom_segment(
    data=tibble(), 
    aes(
      x=154+0:7*360,
      y=c(1.57,1.47,1.44,1.43,1.4,1.39,1.39,1.37),
      xend=c(360+0:6*360,7*360+240),
      yend=c(1.57,1.47,1.44,1.43,1.4,1.39,1.39,1.37),
      colour = as.character(c(3:10))
    ), 
    size=1
  ) +
  scale_color_discrete(
    breaks = c("current","voltage")
  )
# ggsave("./plots/cycle.pdf",width = 8,height = 4)
# ggsave("./plots/cycle_avg.pdf",width = 8,height = 4)

df_long %>% filter(cycle==242) %>% 
  ggplot(aes(x=time_c,y=current / exp_current)) +
  geom_line(alpha = 0.2) +
  geom_smooth(alpha = 0.8,method="gam") +
  labs(x="time in cycle in s")
# ggsave("./plots/cycle_rescurr.pdf",width = 8,height = 4)
  
# plot average current (and average voltage) for each cycle
df_long %>% group_by(cycle) %>% summarize(avg_current=mean(current),avg_voltage=mean(voltage)) %>%
  ggplot(aes(x=cycle,y=avg_current)) +
  geom_line(alpha = 0.8)
# ggsave("./plots/avg_current_over_cycles.pdf",width = 8,height = 4)

df_long %>% filter(time<=min(data$Time2[,242]),
                   time>=max(data$Time2[,238])) %>%
  ggplot(aes(x=time/60,y=current)) +
  geom_line(alpha = 0.8) +
  labs(x = "time in min") +
  annotate("text", x = 11926, y = 1.6, label = "cycle 239") +
  annotate("text", x = 11976, y = 1.6, label = "cycle 240") +
  annotate("text", x = 12026, y = 1.6, label = "cycle 241") +
  annotate("text", x = 11955, y = 0.07, label = "restart?",col="red") +
  annotate("text", x = 12006, y = 0.07, label = "cleaning?",col="red") +
  geom_vline(xintercept = 11947,col="red") +
  geom_vline(xintercept = 11997,col="red")
# ggsave("./plots/current_protocol_avg.pdf",width = 8,height = 4)
# ggsave("./plots/current_protocol.pdf",width = 8,height = 4)

# Qs
# 1 are there results of experiments under different conditions eg materials, clean air
# 2 what is done after each cycle? cooldown?
# 3 what cleaning is done after some cylces before the jump?
# 4 how long do effects of cleaning processes last? into next cycle?
# 5 data not independent
# 6 feature generation, max, quantiles eg, portion over threshold, residuals to expected current
# 7 degradation effects of time overall, time after last cleaning, time in cycle
# 8 is current at peaks of voltage as relevant as at the lower timepoints

# research questions:
# 1 estimate temporal degradation between cleanings and overall degradation (trend + seasonality?)
# 2 significance tests for changing conditions (no data yet, eg purifiede air, different cleaning protocol)


# Qs from reading Thesis and paper
# 1 how is loss measured? degradation as a function of time? eg in Thesis Chapter 4



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # first modeling approaches # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# model current as function of s(voltage) + degradation of time overall/ after cleaning / in cycle
# three different time degradation C exp(-lambda t), estimate C and lambda from data, with confidance bands?
# maybe look at relative residual curves, might need robust (LTS) estimate, looks like fractale

overview_list <- readRDS("./saved_results/all_res_overview_sm.rds")
# look at results, summary, R2, plots of fitted values
# resobj <- readRDS( "./saved_results/model_res_all.rds")
# str(resobj)

# names(resobj$times) <- c("GAM s(voltage)","linear","linear_residual","lts_residual","lmRob_residual")
# resobj$times


# lm_res <- lm(log(current) ~ I(voltage<0.775) + time_c + time_clean + time,data=df_long)
# resobj$lm_res
# summary(resobj$lm_res)
tmp_mod <- overview_list$lm
tmp_mod$time
tmp_mod$coefs
tmp_mod$r2

# tmp_mod$plot

# R2 adj 0.9691, quite similar coef estimates
df_long %>% mutate(fitted = exp(tmp_mod$fitted)) %>%
  filter(cycle>=239,cycle<=241) %>%
  pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=measure)) +
  geom_line(alpha = 0.5) +
  labs(y=" ",x="time in cycle in s")
# not really satisfying?
# y vs yhat
df_long %>% mutate(fitted = exp(tmp_mod$fitted)) %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=current,y=fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 1, color=2)


# reslm_res <- lm(I(log(current/ exp_current )) ~  time_c + time_clean + time,data=df_long)
# resobj$reslm_res
# summary(resobj$reslm_res)
# low R2 0.5582. similar parameter estimates, plotting looks quite okay, similar to before

tmp_mod <- overview_list$lm_res
tmp_mod$time
tmp_mod$coefs
tmp_mod$r2

# tmp_mod$plot

tmp_df <- df_long %>% mutate(fitted = exp(tmp_mod$fitted)*exp_current)
cor(tmp_df$current,tmp_df$fitted) # actual "R2" 0.9649
tmp_df %>%
  filter(cycle>=239,cycle<=241) %>%
  pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=measure)) +
  geom_line(alpha = 0.5) +
  labs(y=" ",x="time in cycle in s")
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=time,y= log(current/exp_current) - log(fitted/exp_current))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=time,y= current - fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)


# # lts_res <- robustbase::ltsReg(I(log(current/ exp_current )) ~  time_c + time_clean + time,data=df_long,model=FALSE)
# resobj$lts_res
# summary(resobj$lts_res)
# # lower R2 adj 0.6565
# df_long %>% mutate(fitted = exp(resobj$lts_res$fitted.values)*exp_current) %>%
#   filter(cycle>=239,cycle<=241) %>%
#   pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
#   ggplot(aes(x=time,y=value,color=measure)) +
#   geom_line(alpha = 0.5) +
#   labs(y=" ",x="time in cycle in s")
# #  really unsatisfying, for now, NAs for log(cur/exp_curr -1), so now fit without -1, but try lmrob for faster runtime, this took almost 2h 
# df_long %>% mutate(fitted = exp(resobj$lts_res$fitted.values)*exp_current) %>%
#   filter(cycle<=241,cycle>=200) %>%
#   ggplot(aes(x=current,y=fitted)) +
#   geom_point(alpha = 0.5) +
#   geom_abline(intercept = 0, slope= 1, color=2)


# lmrob_res <- robustbase::lmrob(I(log(current/ exp_current )) ~  time_c + time_clean + time,data=df_long,model=FALSE)
# resobj$lmrob_res
# tmp_sum <- summary(resobj$lmrob_res)
# str(tmp_sum)
# again similar estimates, R2 = 0.728
tmp_mod <- overview_list$lmRob_res
tmp_mod$time
tmp_mod$coefs
tmp_mod$r2

# tmp_mod$plot

tmp_df <- df_long %>% mutate(fitted = exp(tmp_mod$fitted)*exp_current)
cor(tmp_df$current,tmp_df$fitted) # actual "R2" 0.9647
tmp_df %>%
  filter(cycle>=239,cycle<=241) %>%
  pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=measure)) +
  geom_line(alpha = 0.5) +
  labs(y=" ",x="time in cycle in s")
# decay in 3rd cycle is too strong, maybe not exponential?
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=time,y= log(current/exp_current) - log(fitted/exp_current))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=time,y= current - fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)


# add factor cycle? or effect of time in cycle in each cycle
# or instead of exp decay; fit s() smooth term


# # exponential decay moight be too high, try linear
# 
# tstamp <- Sys.time()
# lmrob2_res <- robustbase::lmrob(I(current/ exp_current ) ~  time_c + time_clean + time,data=df_long,model=FALSE)
# time_lmrob2 <- as.numeric(Sys.time() - tstamp,units="secs")
# lmrob2_res
# tmp_sum <- summary(lmrob2_res)
# str(tmp_sum)
# # R2 adj 0.713
# tmp_df <- df_long %>% mutate(fitted = lmrob2_res$fitted.values*exp_current)
# cor(tmp_df$current,tmp_df$fitted) # actual "R2" 0.96167, less than with log
# tmp_df %>%
#   filter(cycle>=239,cycle<=241) %>%
#   pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
#   ggplot(aes(x=time,y=value,color=measure)) +
#   geom_line(alpha = 0.5) +
#   labs(y=" ",x="time in cycle in s")
# # decay in 3rd cycle is too strong, maybe not exponential?
# tmp_df %>%
#   filter(cycle<=241,cycle>=220) %>%
#   ggplot(aes(x=time,y= (current/exp_current) - (fitted/exp_current))) +
#   geom_point(alpha = 0.5) +
#   geom_abline(intercept = 0, slope= 0, color=2)
# tmp_df %>%
#   filter(cycle<=241,cycle>=220) %>%
#   ggplot(aes(x=time,y= current - fitted)) +
#   geom_point(alpha = 0.5) +
#   geom_abline(intercept = 0, slope= 0, color=2)


# General Additive Models


# resobj$gam_res
# summary(resobj$gam_res)
# plot(resobj$gam_res)
# R2 adj 0.976, all time deg significant, cycle highest rate

tmp_mod <- overview_list$gam_volt
tmp_mod$time
tmp_mod$coefs
tmp_mod$r2

# tmp_mod$plot
tmp_df <- df_long %>% mutate(fitted = exp(tmp_mod$fitted))
cor(tmp_df$current,tmp_df$fitted) # actual "R2" 0.9871

tmp_df %>%
  filter(cycle>=239,cycle<=241) %>%
  pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=measure)) +
  geom_line(alpha = 0.5) +
  labs(y=" ",x="time in cycle in s")
# not really satisfying?
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=current,y=fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 1, color=2)
# residual plots
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=time,y= log(current) - log(fitted))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)
mean(resobj$gam_res$residuals)
tmp_df %>%
  filter(cycle<=241,cycle>=220) %>%
  ggplot(aes(x=time,y= current - fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)


# resgamb <- readRDS( "./saved_results/model_gamb.rds")
# str(resgamb)
# plot(resgamb$gamb_res)
# resgamb$gamb_res
# summary(resgamb$gamb_res)
# # R2 adj 0.94 with additive form

# gamb_res <- mgcv::gam(log(current) ~ s(voltage) + s(time_c) + s(time_clean) + s(time),data=df_long)
# resgamb_log <- readRDS( "./saved_results/model_gamb_log.rds")
# str(resgamb_log)
# plot(resgamb_log$gamb_res)
# resgamb_log$gamb_res
# summary(resgamb_log$gamb_res)
# R2 adj 0.99 with log additive form

tmp_mod <- overview_list$gam_all
tmp_mod$time
tmp_mod$coefs
tmp_mod$r2

# tmp_mod$plot


tmp_df <- df_long %>% mutate(fitted = exp(tmp_mod$fitted))
cor(tmp_df$current,tmp_df$fitted) # actual "R2" 0.9871
tmp_df %>%
  filter(cycle>=239,cycle<=241) %>%
  pivot_longer(c(current,fitted),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time,y=value,color=measure)) +
  geom_line(alpha = 0.5) +
  labs(y=" ",x="time in cycle in s")
tmp_df %>%
  filter(cycle<=241,cycle>=235) %>%
  ggplot(aes(x=time,y= log(current) - log(fitted))) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)
tmp_df %>%
  filter(cycle<=241,cycle>=235) %>%
  ggplot(aes(x=time,y= current - fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 0, color=2)
tmp_df %>%
  filter(cycle<=241,cycle>=235) %>%
  ggplot(aes(x=current,y=  fitted)) +
  geom_point(alpha = 0.5) +
  geom_abline(intercept = 0, slope= 1, color=2)



# if decay in each cycle is nnot relevant, we could compute a set of features of each cycle and employ
# multivariate regression methods on those features dependent of time overall and time after each cleaning

# too many data? is degradation within cycle of interest? 
# If not, we could summarize each cycle by set of features and only take cycles as observations
# and only estimate degradation over time after cleaning and time overall,
# in lm, we saw that degradation wrt time in cycle has the highest rate

# other modeling approach instead of exponential decay? all non parametric into gam?


# first differences of log current, heteroscedacity



# questions: overleaf project with plots and questions
# 0 compare their scatter plot with our plot? 
# 1 explain what happens between all cycles, and what is the different protocol to the cleaning process?
# 2 is degradation within cycle of interest? 
# 3 what is the main goal? which kind of degradation?
# 4 compare different materials/conditions (purified air)? need baseline? are there more datasets
# 5 is there a current paper, where stat analysis is needed?

# definition degradation? loss in current for same voltage
# formula expected voltage

# clear questions

