# # # # # # # # # # # # # # # # # # # # # # # # 
# # # some help functions # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # 
pacman::p_load(tidyr,dplyr,ggplot2,mgcv,lme4)

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

mygrad <- function(x,maxh) {
  myind <- (maxh+1):(length(x)-maxh)
  grad <- sapply(myind,function(j) sum(x[j+1:maxh] - x[j-1:maxh])/(2*maxh))
  return(c(rep(NA,maxh),
           grad,
           rep(NA,maxh)))
}

model_list <- list(
  "doub_exp_simple" = function(data,name,doplot=TRUE){
  model_res <- nls(log(current) ~ b0 + b1*time + b3*exp(b4*time),
                    data=data,
                    start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4),
                  algorithm="port")
  # model_res <- nls(log(current) ~ b0 + b1*time + b3*exp(b4*time),
  #                  data=data,
  #                  start = c(b0=-2,b1=-1e-4,b3=1,b4=-1e-4))
  preds <- exp(predict(model_res))
  mycor <-  cor(preds,data$current)
  pred_start <- preds[1]
  deg_times <- numeric(9)
  ratios <- 9:1/10
  for (k in 1:9) {
    uniroot_res <- uniroot(function(tmp_t) {exp(predict(model_res,newdata = data.frame(time=tmp_t,time_l_jump=0),))-ratios[k]*pred_start},
                           interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
    deg_times[k] <- uniroot_res$root
  }
  if (doplot) {
    tmp_df <- data %>% mutate(fitted=preds,
                              resp = current , 
                              id = 1:length(data$time)) 
    tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
      pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
      ggplot(aes(x=time,y=value,color=type)) +
      geom_line(alpha=0.8) 
    # tmp_plot
    ggsave(paste0("./plots_rev_deg/FitCurve_des_",name,".pdf"),tmp_plot,width = 8,height = 5)
    
    show_time <- 1000*120
    stepsize <- 1e4
    tmp_plot <- tmp_df %>% filter(time>show_time,
                                  time < show_time + stepsize) %>%
      pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
      ggplot(aes(x=time,y=value,color=type)) +
      geom_line(alpha=0.8) +
      coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
    # tmp_plot
    ggsave(paste0("./plots_rev_deg/FitCurveClose_des_",name,".pdf"),tmp_plot,width = 8,height = 5)
    
    tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
      ggplot(aes(x=time,y=resp-fitted)) +
      geom_point(alpha=0.8) +
      geom_hline(yintercept = 0,alpha=0.8,col=2)
    # tmp_plot
    ggsave(paste0("./plots_rev_deg/Residuals_des_",name,".pdf"),tmp_plot,width = 8,height = 5)
  }
  return(list(summary = summary(model_res), 
              cor = mycor,
              deg_times = deg_times,
              pred_start = pred_start))
  }, 
  "lin_exp_simple" = function(data,name,doplot=TRUE){
    model_res <- nls(current ~ b0 + b1*time + b3*exp(b4*time),
                     data=data,
                     start = c(b0=1,b1=-1e-4,b3=1,b4=-1e-4))
    preds <- predict(model_res)
    mycor <-  cor(preds,data$current)
    pred_start <- preds[1]
    deg_times <- numeric(9)
    ratios <- 9:1/10
    for (k in 1:9) {
      uniroot_res <- uniroot(function(tmp_t) {(predict(model_res,newdata = data.frame(time=tmp_t,time_l_jump=0),))-ratios[k]*pred_start},
                             interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
      deg_times[k] <- uniroot_res$root
    }
    
    if (doplot) {
      tmp_df <- data %>% mutate(fitted=preds,
                                resp = current , 
                                id = 1:length(data$time)) 
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurve_les_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      show_time <- 1000*120
      stepsize <- 1e4
      tmp_plot <- tmp_df %>% filter(time>show_time,
                                    time < show_time + stepsize) %>%
        pivot_longer(c(fitted,resp),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) +
        coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurveClose_les_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        ggplot(aes(x=time,y=resp-fitted)) +
        geom_point(alpha=0.8) +
        geom_hline(yintercept = 0,alpha=0.8,col=2)
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/Residuals_les_",name,".pdf"),tmp_plot,width = 8,height = 5)
    }
    return(list(summary = summary(model_res), 
                cor = mycor,
                deg_times = deg_times,
                pred_start = pred_start))
  },
  "gam" = function(data,name,doplot=TRUE){
    model_res <- gam(log(current) ~ s(time) + s(time_l_jump),data=data)
    preds <- exp(predict(model_res))
    
    mycor <-  cor(preds,data$current)
    pred_start <- preds[1]
    deg_times <- deg_times_lower <- deg_times_upper <- numeric(9)
    ratios <- 9:1/10
    for (k in 1:9) {
      uniroot_res <- uniroot(function(tmp_t) {exp(predict(model_res,newdata = data.frame(time=tmp_t,time_l_jump=0,jumps=0),))-ratios[k]*pred_start},
                             interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
      deg_times[k] <- uniroot_res$root
      uniroot_res_lower <- uniroot(function(tmp_t) {
        pred_res <- predict(model_res,newdata = data.frame(time=tmp_t,time_l_jump=0,jumps=0),se.fit = TRUE)
        exp(pred_res$fit-qnorm(0.975)*pred_res$se.fit)-ratios[k]*pred_start},
                             interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
      deg_times_lower[k] <- uniroot_res_lower$root
      uniroot_res_upper <- uniroot(function(tmp_t) {
        pred_res <- predict(model_res,newdata = data.frame(time=tmp_t,time_l_jump=0,jumps=0),se.fit = TRUE)
        exp(pred_res$fit+qnorm(0.975)*pred_res$se.fit)-ratios[k]*pred_start},
        interval = c(0,1e3),extendInt = "downX",tol = 1e-2)
      deg_times_upper[k] <- uniroot_res_upper$root
    }
    
    if (doplot) {
      pdf(paste0("./plots_rev_deg/Gam_",name,".pdf"),width = 8,height = 5)
      plot(model_res)
      dev.off()
      
      tmp_df <- data %>% mutate(fitted=preds,
                                resp = current, 
                                fitted0 = exp(predict(model_res,newdata = data.frame(time=data$time,time_l_jump=data$time,jumps=0),)),
                                id = 1:length(data$time)) 
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurve_gam_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      show_time <- 1000*120
      stepsize <- 1e4
      tmp_plot <- tmp_df %>% filter(time>show_time,
                                    time < show_time + stepsize) %>%
        pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) +
        coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurveClose_gam_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        ggplot(aes(x=time,y=resp-fitted)) +
        geom_point(alpha=0.8) +
        geom_hline(yintercept = 0,alpha=0.8,col=2)
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/Residuals_gam_",name,".pdf"),tmp_plot,width = 8,height = 5)
    }
    return(list(summary = summary(model_res), 
                cor = mycor,
                deg_times = deg_times,
                deg_times_upper = deg_times_upper,
                deg_times_lower = deg_times_lower,
                pred_start = pred_start))
  },
  "doub_exp" = function(data,name,doplot=TRUE){
    model_res <- nlmer(log(current) ~ mynonlin(time,time_l_jump,b1,b2,A,R1,lrc1,R2,lrc2)~A|jumps,
                                     data=data,
                                     start = c(b1=-1e-5,b2=-1e-5,A=-2,R1=1,lrc1=-5,R2=1,lrc2=-5),
                       control=nlmerControl(optimizer = "Nelder_Mead",tolPwrss = 1e-10),
                                     verbose=0)
    
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
    
    if (doplot) {
      pdf(paste0("./plots_rev_deg/RanEffsDE_",name,".pdf"),width = 8,height = 5)
      plot(ranef(model_res)$jumps$A)
      dev.off()
    
      tmp_df <- data %>% mutate(fitted=preds,
                                resp = current , 
                                fitted0 = mypred(data$time),
                                id = 1:length(data$time)) 
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurve_DE_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      show_time <- 1000*120
      stepsize <- 1e4
      tmp_plot <- tmp_df %>% filter(time>show_time,
                                    time < show_time + stepsize) %>%
        pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) +
        coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurveClose_DE_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        ggplot(aes(x=time,y=resp-fitted)) +
        geom_point(alpha=0.8) +
        geom_hline(yintercept = 0,alpha=0.8,col=2)
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/Residuals_DE_",name,".pdf"),tmp_plot,width = 8,height = 5)
    }
    return(list(summary = summary(model_res), 
                cor = mycor,
                deg_times = deg_times,
                pred_start = preds[1]))
  },
  "lin_exp" = function(data,name,doplot=TRUE){
    model_res <- nlmer(current ~ mynonlin(time,time_l_jump,b1,b2,A,R1,lrc1,R2,lrc2)~A|jumps,
                       data=data,
                       start = c(b1=-1e-5,b2=-1e-7,A=0.1,R1=0.5,lrc1=-1,R2=0.03,lrc2=-6),
                       control=nlmerControl(optimizer = "Nelder_Mead",tolPwrss = 1e-10),
                       verbose=0)
    preds <- predict(model_res)
    mycor <-  cor(preds,data$current)
    
    coefs <- summary(model_res)$coefficients[,1]
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
    
    if (doplot) {
      pdf(paste0("./plots_rev_deg/RanEffsLE_",name,".pdf"),width = 8,height = 5)
      plot(ranef(model_res)$jumps$A)
      dev.off()
      
      tmp_df <- data %>% mutate(fitted=preds,
                                resp = current , 
                                fitted0 = mypred(data$time),
                                id = 1:length(data$time)) 
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurve_LE_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      show_time <- 1000*120
      stepsize <- 1e4
      tmp_plot <- tmp_df %>% filter(time>show_time,
                                    time < show_time + stepsize) %>%
        pivot_longer(c(fitted,resp,fitted0),names_to = "type",values_to = "value") %>%
        ggplot(aes(x=time,y=value,color=type)) +
        geom_line(alpha=0.8) +
        coord_cartesian(xlim=c(show_time,show_time+stepsize)) 
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/FitCurveClose_LE_",name,".pdf"),tmp_plot,width = 8,height = 5)
      
      tmp_plot <- tmp_df %>% filter(id%%100==0) %>%
        ggplot(aes(x=time,y=resp-fitted)) +
        geom_point(alpha=0.8) +
        geom_hline(yintercept = 0,alpha=0.8,col=2)
      # tmp_plot
      ggsave(paste0("./plots_rev_deg/Residuals_LE_",name,".pdf"),tmp_plot,width = 8,height = 5)
    }
    return(list(summary = summary(model_res), 
                cor = mycor,
                deg_times = deg_times,
                pred_start = preds[1]))
  }
)
