# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # implementation of  Recursive Filter-Augmented Nonlinear Regression (RFNR) Model # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
pacman::p_load(optimx)

# help function to build linear predictor for each log level
lin_pred_RFNR <- function(par,x_use,log_order=2,dum_use=NULL,rho_upper=0.99,
                          jump_direction=c("up","down","none")) {
  jump_direction <- match.arg(jump_direction)
  if (jump_direction=="up") {
    stopifnot(!is.null(dum_use))
    jump_sign <- 1
    include_jumps <- TRUE
  } else if (jump_direction=="down") {
    stopifnot(!is.null(dum_use))
    jump_sign <- -1
    include_jumps <- TRUE
  } else {
    jump_sign <- 0
    include_jumps <- FALSE
  }
  if (log_order==0) {
    pred <- as.numeric(par["c1"]) - as.numeric(par["lam1"])*x_use
  } else if (log_order==1) {
    pred <- as.numeric(par["log_c1"]) - as.numeric(par["lam1"])*x_use
  } else if (log_order == 2) {
    pred <- as.numeric(par["log_c1"]) - as.numeric(par["lam1"])*x_use + as.numeric(par["c2"])*exp(-as.numeric(par["lam2"])*x_use)
  } else {
    if (log_order >3) {
      warning("Highest log order of 3 will be used!")
      log_order <- 3
    }
    pred <- as.numeric(par["log_c1"]) - as.numeric(par["lam1"])*x_use + 
      as.numeric(par["c2"])*exp(-as.numeric(par["lam2"])*x_use
                                + as.numeric(par["c3"])*exp(-as.numeric(par["lam3"])*x_use))
  }
  if (include_jumps) {
    stopifnot(!is.null(dum_use))
    pred <- pred + jump_sign*stats::filter(dum_use*exp(as.numeric(par["log_delta"])), rho_upper/(1+exp(-as.numeric(par["phi"]))), "recursive", init = 0)
  }
  return(pred)
}

# help function gradient for significance
grad_lin_pred_RFNR <- function(par,x_use,dum_use=NULL,log_order=2,rho_upper=0.99,jump_direction=c("up","down","none")) {
  grad <- matrix(c(0),length(x_use),length(par))
  colnames(grad) <- names(par)
  
  jump_direction <- match.arg(jump_direction)
  if (jump_direction=="up") {
    stopifnot(!is.null(dum_use))
    jump_sign <- 1
    include_jumps <- TRUE
  } else if (jump_direction=="down") {
    stopifnot(!is.null(dum_use))
    jump_sign <- -1
    include_jumps <- TRUE
  } else {
    jump_sign <- 0
    include_jumps <- FALSE
  }
  
  if (log_order==0) {
    grad[,"c1"] <- 1
    grad[,"lam1"] <- -x_use
  } else if (log_order==1) {
    grad[,"log_c1"] <- 1
    grad[,"lam1"] <- -x_use
  } else if (log_order == 2) {
    grad[,"log_c1"] <- 1
    grad[,"lam1"] <- -x_use
    grad[,"c2"] <- exp(-as.numeric(par["lam2"])*x_use)
    grad[,"lam2"] <- as.numeric(par["c2"])*exp(-as.numeric(par["lam2"])*x_use)*(-x_use)
  } else {
    if (log_order >3) {
      warning("Highest log order of 3 will be used!")
      log_order <- 3
    }
    grad[,"log_c1"] <- 1
    grad[,"lam1"] <- -x_use
    grad[,"c2"] <- exp(-as.numeric(par["lam2"])*x_use
                       + as.numeric(par["c3"])*exp(-as.numeric(par["lam3"])*x_use))
    grad[,"lam2"] <- as.numeric(par["c2"])*exp(-as.numeric(par["lam2"])*x_use
                                               + as.numeric(par["c3"])*exp(-as.numeric(par["lam3"])*x_use))*(-x_use)
    grad[,"c3"] <- as.numeric(par["c2"])*exp(-as.numeric(par["lam2"])*x_use + as.numeric(par["c3"])*exp(-as.numeric(par["lam3"])*x_use)) * 
      exp(-as.numeric(par["lam3"])*x_use)
    grad[,"lam3"] <- as.numeric(par["c2"])*exp(-as.numeric(par["lam2"])*x_use + as.numeric(par["c3"])*exp(-as.numeric(par["lam3"])*x_use)) * 
      as.numeric(par["c3"])*exp(-as.numeric(par["lam3"])*x_use) * (-x_use)
  }
  if (include_jumps) {
    rho <- rho_upper/(1+exp(-as.numeric(par["phi"])))
    rec_filt <- jump_sign*stats::filter(dum_use*exp(as.numeric(par["log_delta"])), rho, "recursive", init = 0)
    grad[,"log_delta"] <- rec_filt
    grad_rec_filt <- numeric(length(rec_filt))
    for (t in seq_along(grad_rec_filt)) {
      grad_rec_filt[t] <- sum(dum_use[1:t]*(t:1-1)*rho^(t:1-1))
    }
    grad[,"phi"] <- jump_sign*grad_rec_filt * ( exp(as.numeric(par["log_delta"])) / (1+exp(as.numeric(par["phi"]))) ) 
  }
  return(grad)
}

# # # test gradient
# res_obj <- res
# start <- res_obj$pars
# x_use <- res_obj$x_use
# dum_use <- res_obj$dum_use
# mygrad <- grad_lin_pred_RFNR(start,x_use,dum_use,log_order=2,jump_direction=res_obj$jump_direction)
# num_gr <- matrix(c(0),length(x_use),length(start))
# for (i in seq_along(x_use)) {
#   num_gr[i,] <- numDeriv::grad(function(par){
#     res <- lin_pred_RFNR(par,x_use=x_use,dum_use=dum_use,log_order=2,jump_direction=res_obj$jump_direction)
#     res[i]
#     },start,
#     method="simple",method.args = list(eps=1e-10))
# }
# absdiff <- numeric(length(x_use))
# for (i in seq_along(x_use)) {
#   absdiff[i] <- sum(abs(mygrad[i,]-num_gr[i,]))
# }
# summary(absdiff)
# hist(absdiff)
# which.max(absdiff)
# i <- which.max(absdiff)
# num_gr[i,]
# mygrad[i,]
# abs(mygrad[i,]-num_gr[i,])
# mean(abs(mygrad-num_gr))

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # main function to fit model  # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

RFNR <- function(y, x, 
                 dum = NULL,
                 ndx = 1:length(y),
                 log_order = 2,
                 rho_upper = 0.99,
                 jump_direction=c("up","down","none"),
                 init = NULL,
                 methods=c("Rvmmin","nvm"),
                 ...) {
  x_use <- x[ndx]
  jump_direction <- match.arg(jump_direction)
  if (jump_direction=="none") {
    include_jumps <- FALSE
    dum_use <- NULL
  } else {
    stopifnot(!is.null(dum))
    include_jumps <- TRUE
    if (length(dum)==length(y)) {
      dum_use <- dum[ndx]
      warning("Check that event indicators align with used ndx!")
    } else {
      stopifnot(length(ndx)==length(dum))
      dum_use <- dum
    }
    if (jump_direction=="up") {
      jump_sign <- 1
    } else {
      jump_sign <- -1
    }
  }
  
  # generate obj, gradient and start vec depending on log_order
  if (log_order == 0) {
    resp <- y[ndx]
    obj <- function(par) {
      mean(
        (
          resp - lin_pred_RFNR(par,x_use,log_order,dum_use,jump_direction=jump_direction)
        )^2
      )
    }
    gr_obj <- function(par) {
      grad <- numeric(length(par)) 
      names(grad) <- names(par)
      rho <- rho_upper/(1+exp(-par["phi"]))
      if (include_jumps) {
        rec_filt <- jump_sign*stats::filter(dum_use*exp(par["log_delta"]), rho, "recursive", init = 0)
      } else {
        rec_filt <- 0
      }
      
      res <- resp - par["c1"] + par["lam1"]*x_use - rec_filt
      grad["c1"] <- -2*mean(res)
      grad["lam1"] <- 2*mean(res*x_use)
      if (include_jumps) {
        grad["log_delta"] <- -2*mean(res*rec_filt)
        grad_rec_filt <- numeric(length(rec_filt))
        for (t in seq_along(grad_rec_filt)) {
          grad_rec_filt[t] <- sum(dum_use[1:t]*(t:1-1)*rho^(t:1-1))
        }
        grad_rec_filt <- jump_sign*grad_rec_filt * ( exp(par["log_delta"]) / (1+exp(par["phi"])) )
        grad["phi"] <- -2*mean(res*grad_rec_filt)
      }
      return(grad)
    }
    if (include_jumps) {
      start <- numeric(4)
      names(start) <- c("c1","lam1","log_delta","phi")
      if(!is.null(init)) {
        stopifnot(length(init) == length(start)) 
        start[1:4] <- init
      } else {
        start[1:4] <- c(0.5,-1e-4,log(2),-2) 
      }
    } else {
      start <- numeric(2)
      names(start) <- c("c1","lam1")
      if(!is.null(init)) {
        stopifnot(length(init) == length(start)) 
        start[1:2] <- init
      } else {
        start[1:2] <- c(0.5,-1e-4) 
      }
    }
    
  } else {
    resp <- log(y[ndx])
    if (log_order == 1) {
      obj <- function(par) {
        mean(
          (
            resp - lin_pred_RFNR(par,x_use,log_order,dum_use,jump_direction = jump_direction)
          )^2
        )
      }
      gr_obj <- function(par) {
        grad <- numeric(length(par)) 
        names(grad) <- names(par)
        rho <- rho_upper/(1+exp(-par["phi"]))
        if (include_jumps) {
          rec_filt <- jump_sign*stats::filter(dum_use*exp(par["log_delta"]), rho, "recursive", init = 0)
        } else {
          rec_filt <- 0
        }
        res <- resp - par["log_c1"] + par["lam1"]*x_use - rec_filt
        grad["log_c1"] <- -2*mean(res)
        grad["lam1"] <- 2*mean(res*x_use)
        if (include_jumps) {
          grad["log_delta"] <- -2*mean(res*rec_filt)
          grad_rec_filt <- numeric(length(rec_filt))
          for (t in seq_along(grad_rec_filt)) {
            grad_rec_filt[t] <- sum(dum_use[1:t]*(t:1-1)*rho^(t:1-1))
          }
          grad_rec_filt <- jump_sign*grad_rec_filt * ( exp(par["log_delta"]) / (1+exp(par["phi"])) )
          grad["phi"] <- -2*mean(res*grad_rec_filt)
        }
        return(grad)
      }
      if (include_jumps) {
        start <- numeric(4)
        names(start) <- c("log_c1","lam1","log_delta","phi")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:4] <- init
        } else {
          start[1:4] <- c(-2,-1e-4,log(2),-2) 
        }
      } else {
        start <- numeric(2)
        names(start) <- c("log_c1","lam1")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:2] <- init
        } else {
          start[1:2] <- c(-2,-1e-4) 
        }
      }
    } else if (log_order == 2) {
      obj <- function(par) {
        mean(
          (
            resp - lin_pred_RFNR(par,x_use,log_order,dum_use,jump_direction = jump_direction)
          )^2
        )
      }
      gr_obj <- function(par) {
        grad <- numeric(length(par)) 
        names(grad) <- names(par)
        rho <- rho_upper/(1+exp(-par["phi"]))
        if (include_jumps) {
          rec_filt <- jump_sign*stats::filter(dum_use*exp(par["log_delta"]), rho, "recursive", init = 0)
        } else {
          rec_filt <- 0
        }
        res <- resp - par["log_c1"] + par["lam1"]*x_use - par["c2"]*exp(-par["lam2"]*x_use) - rec_filt
        grad["log_c1"] <- -2*mean(res)
        grad["lam1"] <- 2*mean(res*x_use)
        grad["c2"] <- -2*mean(res*exp(-par["lam2"]*x_use))
        grad["lam2"] <- 2*mean(res*x_use*exp(-par["lam2"]*x_use)*par["c2"])
        if (include_jumps) {
          grad["log_delta"] <- -2*mean(res*rec_filt)
          grad_rec_filt <- numeric(length(rec_filt))
          for (t in seq_along(grad_rec_filt)) {
            grad_rec_filt[t] <- sum(dum_use[1:t]*(t:1-1)*rho^(t:1-1))
          }
          grad_rec_filt <- jump_sign*grad_rec_filt * ( exp(par["log_delta"]) / (1+exp(par["phi"])) )
          grad["phi"] <- -2*mean(res*grad_rec_filt)
        }
        return(grad)
      }
      if (include_jumps) {
        start <- numeric(6)
        names(start) <- c("log_c1","lam1","c2","lam2","log_delta","phi")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:6] <- init
        } else {
          start[1:6] <- c(-1,1e-4,0.5,1e-4,log(2),-2) 
        }
      } else {
        start <- numeric(4)
        names(start) <- c("log_c1","lam1","c2","lam2")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:4] <- init
        } else {
          start[1:4] <- c(-1,1e-4,0.5,1e-4) 
        }
      }
      
    } else {
      if (log_order > 3) {
        warning("Highest log order of 3 will be used!")
        log_order <- 3
      }
      obj <- function(par) {
        mean(
          (
            resp - lin_pred_RFNR(par,x_use,log_order,dum_use,jump_direction = jump_direction)
          )^2
        )
      }
      gr_obj <- function(par) {
        grad <- numeric(length(par)) 
        names(grad) <- names(par)
        rho <- rho_upper/(1+exp(-par["phi"]))
        if (include_jumps) {
          rec_filt <- jump_sign*stats::filter(dum_use*exp(par["log_delta"]), rho, "recursive", init = 0)
        } else {
          rec_filt <- 0
        }
        exp_term <- exp(-par["lam2"]*x_use + par["c3"]*exp(-par["lam3"]*x_use))
        res <- resp - par["log_c1"] + par["lam1"]*x_use - 
          par["c2"]*exp_term - rec_filt
        grad["log_c1"] <- -2*mean(res)
        grad["lam1"] <- 2*mean(res*x_use)
        grad["c2"] <- -2*mean(res*exp_term)
        grad["lam2"] <- 2*mean(res*x_use*exp_term*par["c2"])
        grad["c3"] <- -2*par["c2"]*mean(res*exp(-par["lam3"]*x_use)*exp_term)
        grad["lam3"] <- 2*par["c2"]*par["c3"]*mean(res*x_use*exp_term*exp(-par["lam3"]*x_use))
        if (include_jumps) {
          grad["log_delta"] <- -2*mean(res*rec_filt)
          grad_rec_filt <- numeric(length(rec_filt))
          for (t in seq_along(grad_rec_filt)) {
            grad_rec_filt[t] <- sum(dum_use[1:t]*(t:1-1)*rho^(t:1-1))
          }
          grad_rec_filt <- jump_sign*grad_rec_filt * ( exp(par["log_delta"]) / (1+exp(par["phi"])) )
          grad["phi"] <- -2*mean(res*grad_rec_filt)
        }
        return(grad)
      }
      if (include_jumps) {
        start <- numeric(8)
        names(start) <- c("log_c1","lam1","c2","lam2","c3","lam3","log_delta","phi")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:8] <- init
        } else {
          start[1:8] <- c(-1,1e-4,0.5,1e-4,0.5,1e-5,log(2),-2) 
        }
      } else {
        start <- numeric(6)
        names(start) <- c("log_c1","lam1","c2","lam2","c3","lam3")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:6] <- init
        } else {
          start[1:6] <- c(-1,1e-4,0.5,1e-4,0.5,1e-5) 
        }
      }
    }
  }
  
  
  # try all available and look what works best, grad for phi not correct
  opm_res <- opm(start,obj,gr=gr_obj,
                 method=methods,...) 
  
  # gr_obj(start)
  # num_gr <- numDeriv::grad(obj,start,method="Richardson",method.args = list(eps=1e-10))
  # num_gr
  # sum(abs(gr_obj(start)-num_gr))
  
  if (!any(opm_res$convergence==0)) {
    warning("No converged optimization method, try a different solver or adjust controls!")
  }
  # best methods for log_order 0: (nlminb) Rvmmin ucminf nvm
  # best methods for log_order 1: (nlminb) Rvmmin ucfminf (subplex) nvm
  # best methods for log_order 2: nlminb Rvmmin nvm 
  # best methods for log_order 3: BFGS Rvmmin nvm 
  
  ind_best <- which.min(opm_res$value)
  pars <- opm_res[ind_best,1:(ncol(opm_res)-8)]
  
  if (log_order > 0) {
    fitted_base <- exp(lin_pred_RFNR(pars,x_use,log_order,dum_use,jump_direction = "none"))
    fitted_jumps <- exp(lin_pred_RFNR(pars,x_use,log_order,dum_use,jump_direction = jump_direction))
  } else {
    fitted_base <- lin_pred_RFNR(pars,x_use,log_order,dum_use,jump_direction = "none")
    fitted_jumps <- lin_pred_RFNR(pars,x_use,log_order,dum_use,jump_direction = jump_direction)
  }
  
  if (log_order>0) {
    sigma <- sqrt(sum((log(y[ndx])-log(fitted_jumps))^2) / (length(x_use) - length(pars)))
  } else {
    sigma <- sqrt(sum((y[ndx]-fitted_jumps)^2) / (length(x_use) - length(pars)))
  }
  mypars <- as.numeric(pars)
  names(mypars) <- names(pars)
  res <- list(
    pars = mypars,
    fitted_base = fitted_base,
    fitted_jumps = fitted_jumps,
    x_use = x_use,
    y_use = y[ndx],
    dum_use = dum_use,
    opm_res = opm_res,
    rho_upper = rho_upper,
    log_order = log_order,
    jump_direction = jump_direction,
    sigma = sigma,
    residuals = y[ndx] - fitted_jumps,
    log_residuals = log(y[ndx]) - log(fitted_jumps)
    )
  class(res) <- "RFNR"
  return(res)
}

plot.RFNR <-  function(res_obj,type=c("scatter","residuals","res-vs-fitted","res-hist","res-QQ","res-acf")) {
  type <- match.arg(type)
  tmp_df <- data.frame(x=res_obj$x_use,
                       y=res_obj$y_use,
                       base=res_obj$fitted_base,
                       fitted=res_obj$fitted_jumps,
                       lres=as.numeric(log(res_obj$y_use) - log(res_obj$fitted_jumps)))
  
  if (type=="scatter") {
    tmp_df <- tmp_df %>% 
      pivot_longer(c(y,base,fitted),names_to = "type",values_to = "value")
    mygg <- ggplot(tmp_df,aes(x=x,y=value,col=type)) +
      geom_line() +
      labs(y=" ")
  } else if (type=="residuals") {
    mygg <- ggplot(tmp_df,aes(x=x,y=lres)) +
      geom_point() +
      labs(y="log residuals") +
      scale_y_continuous() +
      geom_hline(yintercept = 0,col=2)
  } else if (type=="res-vs-fitted") {
    mygg <- ggplot(tmp_df,aes(x=fitted,y=lres)) +
      geom_point() +
      labs(y="log residuals") +
      scale_y_continuous() +
      scale_x_log10() +
      geom_hline(yintercept = 0,col=2)
  } else if (type=="res-hist") {
    mysd <- res_obj$sigma
    mygg <- ggplot(tmp_df, aes(x=lres)) + 
      geom_histogram(aes(y=..density..), colour="black", fill="white")+
      stat_function(fun = dnorm, args = list(mean = 0, sd = mysd)) +
      labs(x="log residuals")
  } else if (type=="res-QQ") {
    mysd <- sqrt(sum(tmp_df$lres^2)/(length(tmp_df$lres)-length(res_obj$pars)))
    mygg <- ggplot(tmp_df, aes(sample = lres/mysd)) +
      stat_qq() +
      stat_qq_line() +
      labs(x="Theoretical quantiles",y="Standardised log-residuals")
  } else if (type=="res-acf") {
    mygg <- acf(res_obj$log_residuals,main=" ")
  }
  mygg
}
coef.RFNR <-  function(res_obj) {
  res_obj$pars
}
residuals.RFNR <- function(res_obj) {
  as.numeric(res_obj$log_residuals)
}
predict.RFNR <- function(res_obj,xnew=NULL,dumnew=NULL,time_ind=NULL,type=c("response","linear"),include_jumps=TRUE) {
  type <- match.arg(type)
  if (is.null(xnew)){ # return fitted values if no xnew provided
    if (include_jumps) {
      pred <- res_obj$fitted_jumps
    } else {
      pred <- res_obj$fitted_base
    }
    if (type=="linear" & res_obj$log_order>0) {
      pred <- log(pred)
    }
  } else { 
    if (!include_jumps) { # without jumps: just plug in xnew in explicit parametric formula
      pred <- lin_pred_RFNR(res_obj$pars,x_use = xnew,log_order = res_obj$log_order, jump_direction = "none")
      if (res_obj$log_order>0 & type=="response") {
        pred <- exp(pred)
      }
    } else { # with jumps: 
      if (is.null(time_ind)) { # xnew according to later times: need old dum and potential new dum
        if (is.null(dumnew)) {
          dumnew <- rep(0,length(xnew))
        }
        stopifnot(length(dumnew)==length(xnew))
        tmp_x <- c(res_obj$x_use,xnew) 
        tmp_dum <- c(res_obj$dum_use,dumnew)
        pred <- lin_pred_RFNR(res_obj$pars,x_use = tmp_x,log_order = res_obj$log_order,
                              dum_use = tmp_dum,rho_upper = res_obj$rho_upper, jump_direction = res_obj$jump_direction)
        if (res_obj$log_order>0 & type=="response") {
          pred <- exp(pred)
        }
        pred <- pred[-(1:length(res_obj$x_use))]
      } else { # xnew according to time in training range: no new jumps, interpolate rec filter values
        stopifnot(length(time_ind) == length(xnew))
        pred_base <- lin_pred_RFNR(res_obj$pars,x_use = xnew,log_order = res_obj$log_order, jump_direction = "none")
        
        if (res_obj$jump_direction != "none") {
          tmp_pars <- res_obj$pars
          tmp_pars[] <- 0
          tmp_pars["log_delta"] <- res_obj$pars["log_delta"]
          tmp_pars["phi"] <- res_obj$pars["phi"]
          jumps_use <- lin_pred_RFNR(tmp_pars,x_use = rep(0,length(res_obj$x_use)),dum_use = res_obj$dum_use,
                                     log_order = res_obj$log_order, jump_direction = res_obj$jump_direction,
                                     rho_upper = res_obj$rho_upper)
          
          pred_jumps <- approx(x=1:length(res_obj$x_use),y=jumps_use,xout=time_ind,
                               method = "linear",rule=2)
          pred <- pred_base + pred_jumps$y
        } else {
          pred <- pred_base
        }
        if (res_obj$log_order>0 & type=="response") {
          pred <- exp(pred)
        }
      }
    }
  }
  return(pred)
}

# theta_hat <- optimx(...)$par scores <- score_fn(theta_hat, data) hac_cov <- sandwich::vcovHAC(scores)
# optimx returns parameter estimates, but not score functions or covariance matrices. 
# To get a HAC matrix, you need to: Define and compute the score (first derivative of your objective function w.r.t. parameters) 
# at each observation. Use those scores in vcovHAC().
# see https://en.wikipedia.org/wiki/Newey–West_estimator

estfun.RFNR <- function(res_obj) {
  gradM_rhs <- grad_lin_pred_RFNR(res_obj$pars,res_obj$x_use,
                                  dum_use=res_obj$dum_use,
                                  log_order = res_obj$log_order,
                                  rho_upper = res_obj$rho_upper,
                                  jump_direction = res_obj$jump_direction)
  rval <- as.vector(res_obj$log_residuals) * gradM_rhs
  colnames(rval) <- names(coef(res_obj))
  rval
}

bread.RFNR <- function(res_obj) {
  if(!is.null(res_obj$na.action)) class(res_obj$na.action) <- "omit"
  gradM_rhs <- grad_lin_pred_RFNR(res_obj$pars,res_obj$x_use,
                                  dum_use=res_obj$dum_use,
                                  log_order = res_obj$log_order,
                                  rho_upper = res_obj$rho_upper,
                                  jump_direction = res_obj$jump_direction)
  QR <- qr(gradM_rhs,tol=1e-12)
  var_pars <- chol2inv(QR$qr[,1:QR$rank])
  return(var_pars * length(res_obj$y_use))
}


summary.RFNR <-  function(res_obj,corr=FALSE) {
  # # # old
  gradM_rhs <- grad_lin_pred_RFNR(res_obj$pars,res_obj$x_use,
                                  dum_use=res_obj$dum_use,
                                  log_order = res_obj$log_order,
                                  rho_upper = res_obj$rho_upper,
                                  jump_direction = res_obj$jump_direction)
  
  # not the same as grad of obj fcn (score), but = res*gradM
  
  tmp_res <- res_obj
  class(tmp_res) <- "nls"
  QR <- qr(gradM_rhs,tol=1e-12)
  tmp_res$m <- list(getAllPars = function() res_obj$pars,
                    resid = function() as.numeric(res_obj$log_residuals),
                    gradient = function() gradM_rhs,
                    deviance = function() sum(res_obj$log_residuals^2),
                    Rmat = function() qr.R(QR),
                    formula = function() "y~f(x)")
  str(tmp_res)
  
  undebug(sandwich::vcovHAC)
  HACcov <- sandwich::vcovHAC(tmp_res,diagnostics=TRUE)
  str(HACcov)
  HACcov
  HACcov2
  undebug(sandwich::vcovHAC)
  HACcov2 <- sandwich::vcovHAC(res_obj,diagnostics=TRUE)
  residuals(res_obj)
  
  getAnywhere(summary.nls)
  
  # use sqrt(variance x bias.correction), df t approx
  QR <- qr(gradM_rhs,tol=1e-12)
  var_pars <- chol2inv(QR$qr[,1:QR$rank])
  dimnames(var_pars) <- list(colnames(QR$qr)[1:QR$rank],colnames(QR$qr)[1:QR$rank])
  se <- sqrt(diag(var_pars))*res_obj$sigma
  param <- as.numeric(res_obj$pars[colnames(var_pars)])
  tval <- param/se
  param <- cbind(param, se, tval, 2 * pt(as.numeric(abs(tval)), length(res_obj$x_use) - length(res_obj$pars), lower.tail = FALSE))
  dimnames(param) <- list(colnames(QR$qr)[1:QR$rank], c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
  corr_mat <- NULL
  if (corr) {
    corr_mat <- var_pars * res_obj$sigma^2/outer(se, se)
  }
  if (QR$rank < length(res_obj$pars)) {
    warning("There were parameters with no variance leading to degenerate covariance!")
  }
  res <- list(par_sig = param,
              par_cor = corr_mat)
  class(res) <- "summary.RFNR"
  return(res)
}

print.summary.RFNR <- function(sum_obj) {
  print("Significance of parameters:")
  printCoefmat(sum_obj$par_sig) 
  if (!is.null(sum_obj$par_cor)) {
    print("Correlation of parameter estimators:")
    print(sum_obj$par_cor)
  }
}

print.RFNR <- function(res_obj) {
  print(summary(res_obj))
}

# # # # # my performance measures
eval_fit <- function(y,yhat) {
  RMSE <- sqrt(mean((y-yhat)^2))
  # NRMSE <- RMSE / (max(y)-min(y))
  Cor <- cor(y,yhat)
  MAE <- mean(abs(y-yhat))
  return(c(RMSE=RMSE,Cor=Cor,MAE=MAE))
}


# # # # # testing in development

pacman::p_load(tidyverse,readxl)
dt <- read_xlsx("../data/Dati_PERMANENT_reversible.xlsx", sheet = 3)
# dt <- read_excel("../data/Dati_PERMANENT_reversible_RH.xlsx",sheet=1)
names(dt) <- c("time", "current")
dt$time <- (dt$time - 1)*60

ndx <- seq(1, nrow(dt), 60*10)
x <- dt$time
y <- dt$current
dum <- c(0, diff(log(dt$current[ndx])) > 0.2)
# dum <- c(0, abs(diff(log(dt$current[ndx]))) > 0.1)

xnew <- x[-ndx]
time_ind <- 1 + (1:length(x))[-ndx] * (length(ndx)-1) / (max(ndx))

log_order <- 2
res <- RFNR(y,x,dum=dum,ndx=ndx,log_order=log_order,jump_direction = "up")
res_obj <- res
coef(res)
res$opm_res
summary(res,corr=TRUE)
res
# doesnt seem significant, but will be for higher n, s
# due to high noise / sigma, implies low t value

eval_fit(log(y[-ndx]),predict(res,xnew = xnew,time_ind = time_ind,type = "linear"))
eval_fit(log(y[ndx]),predict(res,xnew = x[ndx],time_ind = 1:length(ndx),type = "linear"))
res$sigma*sqrt((length(ndx)-6)/length(ndx))

plot(res) +
  labs(x="time in s")

plot(res,"residuals") +
  labs(x="time in s")

# looks fine, maybe estimate sd for significance of parameters

# # testing in dev of summary
# pred0 <- lin_pred_exp_dec_rec_filt(res_obj$pars,res_obj$x_use,log_order = log_order,dum_use=res_obj$dum_use)
# predup <- lin_pred_exp_dec_rec_filt(res_obj$pars + c(0,0,0,0,0,0,0.5,0),
#                                     res_obj$x_use,log_order = log_order,dum_use=res_obj$dum_use)
# predlo <- lin_pred_exp_dec_rec_filt(res_obj$pars + c(0,0,0,0,0,0,-0.5,0),
#                                     res_obj$x_use,log_order = log_order,dum_use=res_obj$dum_use)
# plot(pred0)
# lines(predup,col=2)
# lines(predlo,col=3)
# summary(predup - predlo)
# summary(gradM_rhs[,7])
# 
# pred0 <- lin_pred_exp_dec_rec_filt(res_obj$pars,res_obj$x_use,log_order = log_order,dum_use=res_obj$dum_use)
# predup <- lin_pred_exp_dec_rec_filt(res_obj$pars + c(0,0,0,0,0,0,0,0.5),
#                                     res_obj$x_use,log_order = log_order,dum_use=res_obj$dum_use)
# predlo <- lin_pred_exp_dec_rec_filt(res_obj$pars + c(0,0,0,0,0,0,0,-0.5),
#                                     res_obj$x_use,log_order = log_order,dum_use=res_obj$dum_use)
# plot(pred0)
# lines(predup,col=2)
# lines(predlo,col=3)
# summary(predup - predlo)
# summary(gradM_rhs[,8])
