
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
# # . what is the question for 90 vs 90v2 in MEA2?

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 1) Reference + 3 other protocols
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

data <- res$df_data %>% filter(dataname %in% c("REF","UPL","LPL","CYC"))
data$dataname <- factor(data$dataname,levels=c("REF","UPL","LPL","CYC"))


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
ggsave("./plots/REFs_fit.pdf",tmp_gg,height = 5,width = 8)

# # # plots of residuals
data %>%
  ggplot(aes(x=hours,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/REFs_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free_x") +
  labs(y="residuals") +
  theme_bw()
ggsave("./plots/REFs_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# these models look good


# # # # plots of parameter

df_coefs <- res$df_coefs %>% filter(dataname %in% c("REF","UPL","LPL","CYC"))
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("REF","UPL","LPL","CYC"))
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
ggsave("./plots/REFs_reldeg.pdf",height = 5*0.8,width = 8*0.8)

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
ggsave("./plots/REFs_reldeg-values.pdf",height = 5,width = 8)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 2) outlet / inlet / dry inlet
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

data <- res$df_data %>% filter(dataname %in% c("Outlet","Inlet","Dry-Inlet"))
data$dataname <- factor(data$dataname,levels=c("Outlet","Inlet","Dry-Inlet"))


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
ggsave("./plots/OutIn_fit.pdf",tmp_gg,height = 5,width = 8)

# # # plots of residuals
data %>%
  ggplot(aes(x=hours,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free_y") +
  labs(y="residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/OutIn_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free") +
  labs(y="residuals") +
  theme_bw()
ggsave("./plots/OutIn_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# some unexplained structure


# # # # plots of parameter

df_coefs <- res$df_coefs %>% filter(dataname %in% c("Outlet","Inlet","Dry-Inlet"))
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("Outlet","Inlet","Dry-Inlet"))
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
ggsave("./plots/OutIn_coefs.pdf",height = 5,width = 8)

# # # # plots of rel deg
data %>%
  ggplot(aes(x=time_clean,y=relative_degradation,color=cleaning)) +
  geom_line(alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="relative degradation",x="time (h)") +
  theme_bw() +
  scale_color_brewer(type="qual",direction = 1,palette=1)
ggsave("./plots/OutIn_reldeg.pdf",height = 5*0.8,width = 8*0.8)

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
ggsave("./plots/OutIn_reldeg-values.pdf",height = 5,width = 8)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 3) 80/90/100 degrees MEA1 + MEA2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


data <- res$df_data %>% filter(dataname %in% c("80","90","100"))
data$dataname <- factor(data$dataname,levels=c("80","90","100"))


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
ggsave("./plots/Temp_fit.pdf",tmp_gg,height = 5,width = 8)

# # # plots of residuals
data %>%
  ggplot(aes(x=hours,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_grid(material~dataname,scales="free_y") +
  labs(y="residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/Temp_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_grid(material~dataname,scales="free") +
  labs(y="residuals") +
  theme_bw()
ggsave("./plots/Temp_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# these models look good


# # # # plots of parameter

df_coefs <- res$df_coefs %>% filter(dataname %in% c("80","90","100"))
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("80","90","100"))
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
ggsave("./plots/Temp_reldeg.pdf",height = 5*0.8,width = 8*0.8)

# # plot of rel deg values
df_reldeg <- res$df_reldeg %>% filter(dataname %in% c("80","90","100"))
df_reldeg$dataname <- factor(df_reldeg$dataname,levels=c("80","90","100"))

df_reldeg$name <- factor(df_reldeg$name,levels=c("t_rd70","rd_t50"),
                         labels = c(expression(Time~(h)~until~"70%"~of~start),expression(Rel.~deg.~after~50~hours)))

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
ggsave("./plots/TempMEA1_reldeg-values.pdf",height = 5,width = 8)

df_reldeg %>% filter(material=="MEA2") %>%
  ggplot(aes(x=cleaning,y=value)) +
  geom_point() + 
  geom_errorbar(aes(ymin = lower,ymax=upper)) +
  labs(x="cleaning period",y=" ") +
  facet_grid(name~dataname,scales="free_y",
             labeller = label_parsed) +
  scale_y_log10() +
  theme_bw()
ggsave("./plots/TempMEA2_reldeg-values.pdf",height = 5,width = 8)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # 4) 90 v2 MEA2
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


data <- res$df_data %>% filter(dataname %in% c("90","90v2"),material=="MEA2")

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
ggsave("./plots/90_fit.pdf",tmp_gg,height = 5,width = 8)

# # # plots of residuals
data %>%
  ggplot(aes(x=hours,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.) +
  labs(y="residuals",x="time (h)") +
  theme_bw()
ggsave("./plots/90_residuals.pdf",height = 5*0.8,width = 8*0.8)

data %>%
  ggplot(aes(x=fitted,y=current-fitted)) +
  geom_point(alpha=0.8) +
  geom_hline(yintercept = 0,col=2,alpha=0.8) +
  facet_wrap(dataname~.,scales="free_x") +
  labs(y="residuals") +
  theme_bw()
ggsave("./plots/90_fitted-vs-residuals.pdf",height = 5*0.8,width = 8*0.8)
# these models look good


# # # # plots of parameter

df_coefs <- res$df_coefs %>% filter(dataname %in% c("90","90v2"),material=="MEA2")
df_coefs$dataname <- factor(df_coefs$dataname,levels=c("90","90v2"),
                            labels = c("90","`90v2`"))
# df_coefs$se[df_coefs$coef_name=="llam1"&df_coefs$dataname=="LPL"&df_coefs$cleaning==7] <- Inf
# was 12.73802
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
ggsave("./plots/90_reldeg.pdf",height = 5*0.8,width = 8*0.8)

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
ggsave("./plots/90_reldeg-values.pdf",height = 5,width = 8)

