if("pacman" %in% rownames(installed.packages()) == FALSE) {install.packages("pacman")}
# install (if necessary) and load packages
pacman::p_load(dplyr,tidyr,ggplot2,readxl)

rhs <- c(paste0("RH",c(100,50,30,15)))
source("./models_reversible.R")

fits <- data.frame(time=NULL,current=NULL,fit=NULL,id=NULL,rh=NULL)
coefs <- matrix(NA,4,5)
row.names(coefs) <- rhs
colnames(coefs) <- c("cor","C1","lam1","C2","lam2")
for (i in 1:4) {
  data <- read_excel("./data/Dati_PERMANENT_reversible_RH.xlsx",sheet=i)
  colnames(data) <- c("time","current")
  data$time <- (data$time-1)*60
  
  model_res <- nls(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                  data=data,
                  start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4),
                  algorithm="port")
  
  preds <- exp(predict(model_res))
  coefs[i,1] <-  cor(preds,data$current)
  coefs[i,2] <- exp(summary(model_res)$coefficients[1,1])
  coefs[i,3] <- -(summary(model_res)$coefficients[2,1])
  coefs[i,4] <- (summary(model_res)$coefficients[3,1])
  coefs[i,5] <- -(summary(model_res)$coefficients[1,1])
  fits <- rbind(fits,
                data %>% mutate(fit = preds, id = 1:length(preds),rh=rhs[i]))
}


coefs
kable_coefs <- coefs
kable_coefs[,c(1,2,4,5)] <- format(coefs[,c(1,2,4,5)],digits=1)
kable_coefs[,c(3)] <- format(coefs[,c(3)],digits=3)

kable(kable_coefs,format = "latex",booktabs=TRUE)

fits %>% filter(id%%100==0) %>% 
  pivot_longer(c(2,3),names_to = "type",values_to = "current") %>%
  ggplot(aes(x=time,y=current,col=rh,linetype=type)) +
    geom_line()
# ggsave("./plots_rev_deg/relative_humidity_double_exponential.pdf",width = 8,height = 5)
