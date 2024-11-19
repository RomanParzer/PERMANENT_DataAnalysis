
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # # # # start by reading in the .mat data file and explore + visualize the data # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(R.matlab,dplyr,tidyr,ggplot2)

# read in data
data <- readMat("./PERMANENT_CCMG_REFERENCE2.mat")
str(data)

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

df_long <- data.frame(time=c(data$Time2), 
                        voltage = c(data$Voltage2),
                        current = c(data$Current2),
                        cycle = rep(1:1200,each=11133),
                        time_c = c(scale(data$Time2,scale=FALSE,center=apply(data$Time2,2,min))))

# plot voltage + current over one cycle
df_long %>% filter(cycle==400) %>% 
  pivot_longer(c(voltage,current),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=time_c,y=value,color=measure)) +
  geom_line(alpha = 0.8)
  
  
# plot average current (and average voltage) for each cycle
df_long %>% group_by(cycle) %>% summarize(current=mean(current),voltage=mean(voltage)) %>%
  pivot_longer(c(voltage,current),names_to = "measure",values_to = "value") %>%
  ggplot(aes(x=cycle,y=value,color=measure)) +
  geom_line(alpha = 0.8)

