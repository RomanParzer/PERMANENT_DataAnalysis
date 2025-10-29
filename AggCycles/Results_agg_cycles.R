
pacman::p_load(dplyr,tidyr,ggplot2,knitr,ggrepel,kableExtra)

# final version: 
res <- readRDS("./saved_results/results_agg_cycles.rds")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # Structure: 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 1) Reference + 3 other protocols
# # # # 2) outlet / inlet / dry inlet
# # # # 3) 80/90/100 degrees MEA1 + MEA2
# # # # 4) 90 v2 MEA2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # used model on each period between two cleaning protocols: 
# # # # log(avg current over cycle) = log(c1) - lam1*time_since_last_cleaning + c2*exp(-lam2*time_since_last_cleaning)
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# # Q: 
# # - Extend rel deg plots to same maximal length?
# # - different rd=0.7 or t=50h values of relative degradation for the different inspections? esp out/inlet
# # - what is the question for 90 vs 90v2 in MEA2?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 1) Reference + 3 other protocols
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data <- res$df_data %>% filter(dataname %in% c("REF","UPL","LPL","CYC"))
data$dataname <- factor(data$dataname,levels=c("REF","UPL","LPL","CYC"))

df_coefs <- res$df_coefs %>% filter(dataname %in% c("REF","UPL","LPL","CYC"))
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("REF","UPL","LPL","CYC"))



# # # plots of fit
tmp_data <- data %>% 
  pivot_longer(c(fitted,current),names_to = "type",values_to = "value")
tmp_data$type <- factor(tmp_data$type,levels=c("current","fitted"),
                        labels=c("data","fitted"))
tmp_gg <- tmp_data %>%
  ggplot(aes(x=hours,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.,scales="free_y") +
  labs(y=expression(avg.~current~density~(A/cm^2)),x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
tmp_gg
ggsave("./plots/REFs_fit.pdf",tmp_gg,height = 4.4,width = 8)


# # # plots of residuals
tmp_df_sigma <- df_coefs %>% group_by(cleaning,dataname) %>% summarize(sigma = mean(sigma))
mysigma <- function(cleaning,dataname) {
  n <- length(cleaning)
  res <- numeric(n)
  for (i in 1:n) {
    res[i] <- tmp_df_sigma$sigma[tmp_df_sigma$cleaning==cleaning[i] & tmp_df_sigma$dataname==dataname[i]]
  }
  return(res)
}

data <- data %>% group_by(cleaning,dataname) %>% mutate(sigma = mysigma(cleaning,dataname))

data %>%
  ggplot(aes(x=hours,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="stand. residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/REFs_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free_x") +
  labs(y="stand. residuals") +
  theme_bw()
ggsave("./plots/REFs_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# these models look good


# # # # plots of parameter

tab_coef <- df_coefs %>% pivot_wider(id_cols = c(dataname,cleaning),names_from = coef_name, values_from = c(coef,se),names_vary = "slowest")
kab_tab_coef <- matrix(c("0"),48,5)
colnames(kab_tab_coef) <- c("Clean. period","log(c1)","log(lambda1)","log(c2)","log(lambda2)")
kab_tab_coef[,1] <- as.numeric(tab_coef$cleaning)

tab_coef <- signif(as.matrix(tab_coef[,-c(1,2)]),digits=3)
for (i in 2:5) {
  kab_tab_coef[,i] <- paste0(tab_coef[,(i-2)*2 + 1]," (",tab_coef[,(i-2)*2 + 2],")")
}
# copy console output to latex
kab_tab_coef %>% kable(format = "latex",booktabs=TRUE) %>%
  group_rows( c("REF protocol"),1,12) %>%
  group_rows(c("UPL protocol"),13,24) %>%
  group_rows(c("LPL protocol"),25,36) %>%
  group_rows(c("CYC protocol"),37,48)

df_coefs$se[df_coefs$coef_name=="llam1"&df_coefs$dataname=="LPL"&df_coefs$cleaning==7] <- Inf
# was 12.73802
df_coefs$coef_name <- factor(df_coefs$coef_name,levels=c("lc1","llam1","lc2","llam2"),
                             labels = c(expression(log(c[1])),expression(log(lambda[1])),
                                        expression(log(c[2])),expression(log(lambda[2]))))

df_coefs %>% 
  ggplot(aes(x=cleaning,y=coef)) +
  geom_point() +
  geom_errorbar(aes(ymin = coef - qt(0.975,df)*se,
                    ymax = coef + qt(0.975,df)*se)) +
  facet_grid(coef_name~dataname,scales="free_y",
             labeller = label_parsed) +
  labs(y="coef. estimates and 95% CIs",x="cleaning period") +
  theme_bw()
ggsave("./plots/REFs_coefs.pdf",height = 5,width = 8)



# # # # plots of rel deg
data %>% filter(cleaning %in% c(1,3,5,7,9,11,12)) %>%
  ggplot(aes(x=time_clean,y=relative_degradation,color=cleaning)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="relative degradation",x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
ggsave("./plots/REFs_reldeg.pdf",height = 5*0.8^2,width = 8*0.8)

# # plot of rel deg values
df_reldeg <- res$df_reldeg %>% filter(dataname %in% c("REF","UPL","LPL","CYC"))
df_reldeg$dataname <- factor(df_reldeg$dataname,levels=c("REF","UPL","LPL","CYC"))

df_reldeg$name <- factor(df_reldeg$name,levels=c("t_rd70","rd_t50"),
                         labels = c(expression(Time~(h)~until~"70%"~of~start),expression(Rel.~deg.~after~50~hours)))


df_reldeg %>% 
  ggplot(aes(x=cleaning,y=value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower,ymax=upper)) +
  labs(x="cleaning period",y=" ") +
  facet_grid(name~dataname,scales="free_y",
             labeller = label_parsed) +
  scale_y_log10() +
  theme_bw()
ggsave("./plots/REFs_reldeg-values.pdf",height = 4,width = 8)

tab_rd <- df_reldeg %>% pivot_wider(names_from = 1,values_from = c(5,6,7),names_vary = "slowest") %>% select(-c(1,3))
kab_tab_rd <- matrix(c("0"),24,5)
colnames(kab_tab_rd) <- c("Clean. period","REF","UPL","LPL","CYC")
kab_tab_rd[,1] <- as.numeric(tab_rd[,1]$cleaning)
# tab_rd <- format(as.matrix(tab_rd[,-1]),digits = 3,trim=TRUE)
tab_rd <- signif(as.matrix(tab_rd[,-1]),digits=3)
tab_rd
for (i in 2:5) {
  kab_tab_rd[,i] <- paste0(tab_rd[,(i-2)*3 + 1]," (",tab_rd[,(i-2)*3 + 2],",",tab_rd[,(i-2)*3 + 3],")")
}
# copy console output to latex
kab_tab_rd %>% kable(format = "latex",booktabs=TRUE) %>%
  group_rows( c("Time (h) until 70% of start"),1,12) %>%
  group_rows(c("Relative degradation after 50 hours"),13,24)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 2) outlet / inlet / dry inlet
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

data <- res$df_data %>% filter(dataname %in% c("Outlet","Inlet","Dry-Inlet"))
data$dataname <- factor(data$dataname,levels=c("Outlet","Inlet","Dry-Inlet"))

df_coefs <- res$df_coefs %>% filter(dataname %in% c("Outlet","Inlet","Dry-Inlet"))
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("Outlet","Inlet","Dry-Inlet"))

# # # plots of fit
tmp_data <- data %>% 
  pivot_longer(c(fitted,current),names_to = "type",values_to = "value")
tmp_data$type <- factor(tmp_data$type,levels=c("current","fitted"),
                        labels=c("data","fitted"))
tmp_gg <- tmp_data %>%
  ggplot(aes(x=hours,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y=expression(avg.~current~density~(A/cm^2)),x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
tmp_gg
ggsave("./plots/OutIn_fit.pdf",tmp_gg,height = 3.4,width = 8)

# # # plots of residuals

tmp_df_sigma <- df_coefs %>% group_by(cleaning,dataname) %>% summarize(sigma = mean(sigma))
mysigma <- function(cleaning,dataname) {
  n <- length(cleaning)
  res <- numeric(n)
  for (i in 1:n) {
    res[i] <- tmp_df_sigma$sigma[tmp_df_sigma$cleaning==cleaning[i] & tmp_df_sigma$dataname==dataname[i]]
  }
  return(res)
}

data <- data %>% group_by(cleaning,dataname) %>% mutate(sigma = mysigma(cleaning,dataname))

data %>%
  ggplot(aes(x=hours,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free_y") +
  labs(y="stand. residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/OutIn_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free") +
  labs(y="stand. residuals") +
  theme_bw()
ggsave("./plots/OutIn_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# some unexplained structure


# # # # plots of parameter
tab_coef <- df_coefs %>% pivot_wider(id_cols = c(dataname,cleaning),names_from = coef_name, values_from = c(coef,se),names_vary = "slowest")
kab_tab_coef <- matrix(c("0"),15,5)
colnames(kab_tab_coef) <- c("Clean. period","log(c1)","log(lambda1)","log(c2)","log(lambda2)")
kab_tab_coef[,1] <- as.numeric(tab_coef$cleaning)
tab_coef <- signif(as.matrix(tab_coef[,-c(1,2)]),digits=3)
for (i in 2:5) {
  kab_tab_coef[,i] <- paste0(tab_coef[,(i-2)*2 + 1]," (",tab_coef[,(i-2)*2 + 2],")")
}
# copy console output to latex
kab_tab_coef %>% kable(format = "latex",booktabs=TRUE) %>%
  group_rows( c("Outlet protocol"),1,5) %>%
  group_rows(c("Inlet protocol"),6,10) %>%
  group_rows(c("Dry-Inlet protocol"),11,15)

df_coefs$coef_name <- factor(df_coefs$coef_name,levels=c("lc1","llam1","lc2","llam2"),
                             labels = c(expression(log(c[1])),expression(log(lambda[1])),
                                        expression(log(c[2])),expression(log(lambda[2]))))

df_coefs %>% 
  ggplot(aes(x=cleaning,y=coef)) +
  geom_point() +
  geom_errorbar(aes(ymin = coef - qt(0.975,df)*se,
                    ymax = coef + qt(0.975,df)*se)) +
  facet_grid(coef_name~dataname,scales="free_y",
             labeller = label_parsed) +
  labs(y="coef. estimates and 95% CIs",x="cleaning period") +
  theme_bw()
ggsave("./plots/OutIn_coefs.pdf",height = 4.6,width = 8)

# # # # plots of rel deg
data %>%
  ggplot(aes(x=time_clean,y=relative_degradation,color=cleaning)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="relative degradation",x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
ggsave("./plots/OutIn_reldeg.pdf",height = 5*0.8*0.8,width = 8*0.8)

# # plot of rel deg values
df_reldeg <- res$df_reldeg %>% filter(dataname %in% c("Outlet","Inlet","Dry-Inlet"))
df_reldeg$dataname <- factor(df_reldeg$dataname,levels=c("Outlet","Inlet","Dry-Inlet"))

df_reldeg$name <- factor(df_reldeg$name,levels=c("t_rd70","rd_t50"),
                         labels = c(expression(Time~(h)~until~"70%"~of~start),expression(Rel.~deg.~after~50~hours)))


df_reldeg %>% 
  ggplot(aes(x=cleaning,y=value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower,ymax=upper)) +
  labs(x="cleaning period",y=" ") +
  facet_grid(name~dataname,scales="free_y",
             labeller = label_parsed) +
  scale_y_log10() +
  theme_bw()
ggsave("./plots/OutIn_reldeg-values.pdf",height = 4,width = 8)

tab_rd <- df_reldeg %>% pivot_wider(names_from = 1,values_from = c(5,6,7),names_vary = "slowest") %>% select(-c(1,3))
kab_tab_rd <- matrix(c("0"),10,4)
colnames(kab_tab_rd) <- c("Clean. period","Outlet","Inlet","Dry-Inlet")
kab_tab_rd[,1] <- as.numeric(tab_rd[,1]$cleaning)
# tab_rd <- format(as.matrix(tab_rd[,-1]),digits = 3,trim=TRUE)
tab_rd <- signif(as.matrix(tab_rd[,-1]),digits=3)
for (i in 2:4) {
  kab_tab_rd[,i] <- paste0(tab_rd[,(i-2)*3 + 1]," (",tab_rd[,(i-2)*3 + 2],",",tab_rd[,(i-2)*3 + 3],")")
}
# copy console output to latex
kab_tab_rd %>% kable(format = "latex",booktabs=TRUE) %>%
  group_rows( c("Time (h) until 70% of start"),1,5) %>%
  group_rows(c("Relative degradation after 50 hours"),6,10)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 3) 80/90/100 degrees MEA1 + MEA2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


data <- res$df_data %>% filter(dataname %in% c("80","90","100"))
data$dataname <- factor(data$dataname,levels=c("80","90","100"))

df_coefs <- res$df_coefs %>% filter(dataname %in% c("80","90","100"))
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("80","90","100"))

# # # plots of fit
tmp_data <- data %>% 
  pivot_longer(c(fitted,current),names_to = "type",values_to = "value")
tmp_data$type <- factor(tmp_data$type,levels=c("current","fitted"),
                        labels=c("data","fitted"))
tmp_gg <- tmp_data %>%
  ggplot(aes(x=hours,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_grid(material~dataname,scale="free_y") +
  labs(y=expression(avg.~current~density~(A/cm^2)),x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
tmp_gg
ggsave("./plots/Temp_fit.pdf",tmp_gg,height = 4.6,width = 8)


# # # plots of residuals
tmp_df_sigma <- df_coefs %>% group_by(cleaning,dataname,material) %>% summarize(sigma = mean(sigma))
mysigma <- function(cleaning,dataname,material) {
  n <- length(cleaning)
  res <- numeric(n)
  for (i in 1:n) {
    res[i] <- tmp_df_sigma$sigma[tmp_df_sigma$cleaning==cleaning[i] & tmp_df_sigma$dataname==dataname[i] & tmp_df_sigma$material==material[i]]
  }
  return(res)
}

data <- data %>% group_by(cleaning,dataname,material) %>% mutate(sigma = mysigma(cleaning,dataname,material))

data %>%
  ggplot(aes(x=hours,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_grid(material~dataname,scales="free_y") +
  labs(y="stand. residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/Temp_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_grid(material~dataname,scales="free") +
  labs(y="stand. residuals") +
  theme_bw()
ggsave("./plots/Temp_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# these models look good


# # # # plots of parameter
tab_coef <- df_coefs %>% 
  pivot_wider(id_cols = c(dataname,cleaning),names_from = c(coef_name,material), values_from = c(coef,se),names_vary = "slowest")
kab_tab_coef <- matrix(c("0"),14,9)
colnames(kab_tab_coef) <- c("Clean. period",
                            "log(c1)","log(lambda1)","log(c2)","log(lambda2)",
                            "log(c1)","log(lambda1)","log(c2)","log(lambda2)")
kab_tab_coef[,1] <- as.numeric(tab_coef$cleaning)
tab_coef <- signif(as.matrix(tab_coef[,-c(1,2)]),digits=3)
for (i in 2:9) {
  kab_tab_coef[,i] <- paste0(tab_coef[,(i-2)*2 + 1]," (",tab_coef[,(i-2)*2 + 2],")")
}
# copy console output to latex
kab_tab_coef %>% kable(format = "latex",booktabs=TRUE) %>%
  add_header_above(header = c(" "=1, "MEA1"=4,"MEA2"=4)) %>%
  group_rows( c("80 degrees"),1,4) %>%
  group_rows(c("90 degrees"),5,9) %>%
  group_rows(c("100 degrees"),10,14)


df_coefs$coef_name <- factor(df_coefs$coef_name,levels=c("lc1","llam1","lc2","llam2"),
                             labels = c(expression(log(c[1])),expression(log(lambda[1])),
                                        expression(log(c[2])),expression(log(lambda[2]))))


df_coefs %>% filter(material=="MEA1") %>%
  ggplot(aes(x=cleaning,y=coef)) +
  geom_point() +
  geom_errorbar(aes(ymin = coef - qt(0.975,df)*se,
                    ymax = coef + qt(0.975,df)*se)) +
  facet_grid(coef_name~dataname,scales="free",
             labeller = label_parsed) +
  labs(y="coef. estimates and 95% CIs",x="cleaning period") +
  theme_bw()
ggsave("./plots/TempMEA1_coefs.pdf",height = 5,width = 8)

df_coefs %>% filter(material=="MEA2") %>%
  ggplot(aes(x=cleaning,y=coef)) +
  geom_point() +
  geom_errorbar(aes(ymin = coef - qt(0.975,df)*se,
                    ymax = coef + qt(0.975,df)*se)) +
  facet_grid(coef_name~dataname,scales="free",
             labeller = label_parsed) +
  labs(y="coef. estimates and 95% CIs",x="cleaning period") +
  theme_bw()
ggsave("./plots/TempMEA2_coefs.pdf",height = 5,width = 8)

# # # # plots of rel deg
data %>%
  ggplot(aes(x=time_clean,y=relative_degradation,color=cleaning)) +
  geom_line(alpha=0.8) +
  facet_grid(material~dataname) +
  labs(y="relative degradation",x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
ggsave("./plots/Temp_reldeg.pdf",height = 5*0.8*0.8,width = 8*0.8)

# # plot of rel deg values
df_reldeg <- res$df_reldeg %>% filter(dataname %in% c("80","90","100"))
df_reldeg$dataname <- factor(df_reldeg$dataname,levels=c("80","90","100"))

df_reldeg$name <- factor(df_reldeg$name,levels=c("t_rd70","rd_t50"),
                         labels = c(expression(Time~(h)~until~"70%"~of~start),expression(Rel.~deg.~after~50~hours)))


tab_rd <- df_reldeg %>% pivot_wider(names_from = c(1,2),values_from = c(5,6,7),names_vary = "slowest") %>% select(-c(2))
tab_rd <- tab_rd[c(1:4,9,5:8,10),]
kab_tab_rd <- matrix(c("0"),10,7)
colnames(kab_tab_rd) <- c("Clean. period","80 degrees","90 degrees","100 degrees","80 degrees","90 degrees","100 degrees")
kab_tab_rd[,1] <- as.numeric(tab_rd[,1]$cleaning)
tab_rd <- signif(as.matrix(tab_rd[,-1]),digits=3)
for (i in 2:7) {
  kab_tab_rd[,i] <- paste0(tab_rd[,(i-2)*3 + 1]," (",tab_rd[,(i-2)*3 + 2],",",tab_rd[,(i-2)*3 + 3],")")
}
# copy console output to latex
kab_tab_rd %>% kable(format = "latex",booktabs=TRUE) %>%
  add_header_above(header = c(" "=1, "MEA1"=3,"MEA2"=3)) %>%
  group_rows( c("Time (h) until 70% of start"),1,5) %>%
  group_rows(c("Relative degradation after 50 hours"),6,10)



df_reldeg$lower[df_reldeg$material=="MEA1"&df_reldeg$cleaning==4&df_reldeg$dataname=="100"] <- 0
# cleaning 4
# time lower was 0.7961848
# rd lower was 2.852945e-10

df_reldeg$lower[df_reldeg$material=="MEA1"&df_reldeg$cleaning==5&df_reldeg$dataname=="100"] <- 0
# cleaning 5
# time lower was 5.70101264
# rd lower was 0.09358461

df_reldeg %>% filter(material=="MEA1") %>%
  ggplot(aes(x=cleaning,y=value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower,ymax=upper)) +
  labs(x="cleaning period",y=" ") +
  facet_grid(name~dataname,scales="free_y",
             labeller = label_parsed) +
  scale_y_log10() +
  theme_bw()
ggsave("./plots/TempMEA1_reldeg-values.pdf",height = 4.4,width = 8)

df_reldeg %>% filter(material=="MEA2") %>%
  ggplot(aes(x=cleaning,y=value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower,ymax=upper)) +
  labs(x="cleaning period",y=" ") +
  facet_grid(name~dataname,scales="free_y",
             labeller = label_parsed) +
  scale_y_log10() +
  theme_bw()
ggsave("./plots/TempMEA2_reldeg-values.pdf",height = 4.4,width = 8)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 4) 90 v2 MEA2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
data <- res$df_data %>% filter(dataname %in% c("90","90v2"),material=="MEA2")
data$dataname <- factor(data$dataname,levels=c("90","90v2"),
                                             labels = c("90","`90v2`"))
# somehow the string "90v2" is acting weird (some escape?), solev for now by label "`90v2`"
df_coefs <- res$df_coefs %>% filter(dataname %in% c("90","90v2"),material=="MEA2")
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("90","90v2"),
                            labels = c("90","`90v2`"))

# # # plots of fit
tmp_data <- data %>% 
  pivot_longer(c(fitted,current),names_to = "type",values_to = "value")
tmp_data$type <- factor(tmp_data$type,levels=c("current","fitted"),
                        labels=c("data","fitted"))
tmp_gg <- tmp_data %>%
  ggplot(aes(x=hours,y=value,color=type)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y=expression(avg.~current~density~(A/cm^2)),x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
tmp_gg
ggsave("./plots/90_fit.pdf",tmp_gg,height = 3.4,width = 8)

# # # plots of residuals
tmp_df_sigma <- df_coefs %>% group_by(cleaning,dataname) %>% summarize(sigma = mean(sigma))
mysigma <- function(cleaning,dataname) {
  n <- length(cleaning)
  res <- numeric(n)
  for (i in 1:n) {
    res[i] <- tmp_df_sigma$sigma[tmp_df_sigma$cleaning==cleaning[i] & tmp_df_sigma$dataname==dataname[i]]
  }
  return(res)
}

data <- data %>% group_by(cleaning,dataname) %>% mutate(sigma = mysigma(cleaning,dataname))

data %>%
  ggplot(aes(x=hours,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="stand. residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/90_residuals.pdf",height = 5*0.8^2,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=(current-fitted)/sigma)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free_x") +
  labs(y="stand. residuals") +
  theme_bw()
ggsave("./plots/90_fitted-vs-residuals.pdf",height = 5*0.8^2,width = 8*0.8)
# these models look good


# # # # plots of parameter
tab_coef <- df_coefs %>% pivot_wider(id_cols = c(dataname,cleaning),names_from = coef_name, values_from = c(coef,se),names_vary = "slowest")
kab_tab_coef <- matrix(c("0"),5,5)
colnames(kab_tab_coef) <- c("Clean. period","log(c1)","log(lambda1)","log(c2)","log(lambda2)")
kab_tab_coef[,1] <- as.numeric(tab_coef$cleaning)
tab_coef <- signif(as.matrix(tab_coef[,-c(1,2)]),digits=3)
for (i in 2:5) {
  kab_tab_coef[,i] <- paste0(tab_coef[,(i-2)*2 + 1]," (",tab_coef[,(i-2)*2 + 2],")")
}
# copy console output to latex
kab_tab_coef %>% kable(format = "latex",booktabs=TRUE) %>%
  group_rows( c("90 degrees MEA2"),1,3) %>%
  group_rows(c("90 degrees v2 MEA2"),4,5)


df_coefs$coef_name <- factor(df_coefs$coef_name,levels=c("lc1","llam1","lc2","llam2"),
                             labels = c(expression(log(c[1])),expression(log(lambda[1])),
                                        expression(log(c[2])),expression(log(lambda[2]))))


df_coefs %>% 
  ggplot(aes(x=cleaning,y=coef)) +
  geom_point() +
  geom_errorbar(aes(ymin = coef - qt(0.975,df)*se,
                    ymax = coef + qt(0.975,df)*se)) +
  facet_grid(coef_name~dataname,scales="free",
             labeller = label_parsed) +
  labs(y="coef. estimates and 95% CIs",x="cleaning period") +
  theme_bw()
ggsave("./plots/90_coefs.pdf",height = 5,width = 8)

# # # # plots of rel deg
data %>% 
  ggplot(aes(x=time_clean,y=relative_degradation,color=cleaning)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="relative degradation",x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
ggsave("./plots/90_reldeg.pdf",height = 3.4,width = 8*0.8)

# # plot of rel deg values
df_reldeg <- res$df_reldeg %>% filter(dataname %in% c("90","90v2"),material=="MEA2")
df_reldeg$dataname <- factor(df_reldeg$dataname,levels=c("90","90v2"),
                            labels = c("90","`90v2`"))
df_reldeg$name <- factor(df_reldeg$name,levels=c("t_rd70","rd_t50"),
                         labels = c(expression(Time~(h)~until~"70%"~of~start),expression(Rel.~deg.~after~50~hours)))


df_reldeg %>% 
  ggplot(aes(x=cleaning,y=value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower,ymax=upper)) +
  labs(x="cleaning period",y=" ") +
  facet_grid(name~dataname,scales="free",
             labeller = label_parsed) +
  scale_y_log10() +
  theme_bw()
ggsave("./plots/90_reldeg-values.pdf",height = 4.4,width = 8*0.8)

tab_rd <- df_reldeg %>% pivot_wider(names_from = 1,values_from = c(5,6,7),names_vary = "slowest") %>% select(-c(1,3))
kab_tab_rd <- matrix(c("0"),6,3)
colnames(kab_tab_rd) <- c("Clean. period","90 degrees","90 degrees v2")
kab_tab_rd[,1] <- as.numeric(tab_rd[,1]$cleaning)
tab_rd <- signif(as.matrix(tab_rd[,-1]),digits=3)
for (i in 2:3) {
  kab_tab_rd[,i] <- paste0(tab_rd[,(i-2)*3 + 1]," (",tab_rd[,(i-2)*3 + 2],",",tab_rd[,(i-2)*3 + 3],")")
}
# copy console output to latex
kab_tab_rd %>% kable(format = "latex",booktabs=TRUE) %>%
  group_rows( c("Time (h) until 70% of start"),1,3) %>%
  group_rows(c("Relative degradation after 50 hours"),4,6)
