# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # Structure:  # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # 1: build linear predictor
# # 2: gradient of linear predictor
# # 3: evaluation metrics
# # 4: build obj function and gradient

# # 1 help function to build linear predictor for each log level
lin_pred_RFNR <- function(par,x_use,log_order=2,dum_use=NULL,rho_upper=1,
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
    pred <- as.numeric(par["c1"]) - as.numeric(exp(par["log_lam1"]))*x_use
  } else if (log_order==1) {
    pred <- as.numeric(par["log_c1"]) - as.numeric(exp(par["log_lam1"]))*x_use
  } else if (log_order == 2) {
    pred <- as.numeric(par["log_c1"]) - as.numeric(exp(par["log_lam1"]))*x_use + as.numeric(exp(par["log_c2"]))*exp(-as.numeric(exp(par["log_lam2"]))*x_use)
  } else {
    if (log_order >3) {
      warning("Highest log order of 3 will be used!")
      log_order <- 3
    }
    pred <- as.numeric(par["log_c1"]) - as.numeric(exp(par["log_lam1"]))*x_use + 
      as.numeric(exp(par["log_c2"]))*exp(-as.numeric(exp(par["log_lam2"]))*x_use
                                         + as.numeric(exp(par["log_c3"]))*exp(-as.numeric(exp(par["log_lam3"]))*x_use))
  }
  if (include_jumps) {
    stopifnot(!is.null(dum_use))
    pred <- pred + jump_sign*stats::filter(dum_use*exp(as.numeric(par["log_delta"])), rho_upper/(1+exp(-as.numeric(par["phi"]))), "recursive", init = 0)
  }
  return(pred)
}

# # 2 help function gradient for significance
grad_lin_pred_RFNR <- function(par,x_use,dum_use=NULL,log_order=2,rho_upper=1,jump_direction=c("up","down","none")) {
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
    grad[,"log_lam1"] <- -x_use*as.numeric(exp(par["log_lam1"]))
  } else if (log_order==1) {
    grad[,"log_c1"] <- 1
    grad[,"log_lam1"] <- -x_use*as.numeric(exp(par["log_lam1"]))
  } else if (log_order == 2) {
    grad[,"log_c1"] <- 1
    grad[,"log_lam1"] <- -x_use*as.numeric(exp(par["log_lam1"]))
    grad[,"log_c2"] <- exp(-as.numeric(exp(par["log_lam2"]))*x_use)*as.numeric(exp(par["log_c2"]))
    grad[,"log_lam2"] <- as.numeric(exp(par["log_c2"]))*exp(-as.numeric(exp(par["log_lam2"]))*x_use)*(-x_use)*as.numeric(exp(par["log_lam2"]))
  } else { 
    if (log_order >3) {
      warning("Highest log order of 3 will be used!")
      log_order <- 3
    }
    grad[,"log_c1"] <- 1
    grad[,"log_lam1"] <- -x_use*as.numeric(exp(par["log_lam1"]))
    grad[,"log_c2"] <- (exp(-as.numeric(exp(par["log_lam2"]))*x_use
                            + as.numeric(exp(par["log_c3"]))*exp(-as.numeric(exp(par["log_lam3"]))*x_use)))*as.numeric(exp(par["log_c2"]))
    grad[,"log_lam2"] <- (as.numeric(exp(par["log_c2"]))*exp(-as.numeric(exp(par["log_lam2"]))*x_use
                                                             + as.numeric(exp(par["log_c3"]))*exp(-as.numeric(exp(par["log_lam3"]))*x_use))*(-x_use))*as.numeric(exp(par["log_lam2"]))
    
    grad[,"log_c3"] <- (as.numeric(exp(par["log_c2"]))*exp(-as.numeric(exp(par["log_lam2"]))*x_use + as.numeric(exp(par["log_c3"]))*exp(-as.numeric(exp(par["log_lam3"]))*x_use))) * 
      exp(-as.numeric(exp(par["log_c3"]))*x_use)
    grad[,"log_lam3"] <- (as.numeric(exp(par["log_c2"]))*exp(-as.numeric(exp(par["log_lam2"]))*x_use + as.numeric(exp(par["log_c3"]))*exp(-as.numeric(exp(par["log_lam3"]))*x_use)) * 
                            as.numeric(exp(par["log_c3"]))*exp(-as.numeric(exp(par["log_lam3"]))*x_use) * (-x_use))*exp(-as.numeric(exp(par["log_c3"]))*x_use)
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

# # # # test gradient
# res_obj <- res
# start <- res_obj$pars
# # start[c(2:4)] <- log(start[c(2:4)])
# # names(start)[2:4] <- c("log_lam1","log_c2","log_lam2")
# x_use <- res_obj$x_use
# dum_use <- res_obj$dum_use
# mygrad <- grad_lin_pred_RFNR(start,x_use,dum_use,log_order=2,jump_direction=res_obj$jump_direction)
# num_gr <- matrix(c(0),length(x_use),length(start))
# for (i in seq_along(x_use)) {
#   num_gr[i,] <- numDeriv::grad(function(par){
#     res <- lin_pred_RFNR(par,x_use=x_use,dum_use=dum_use,log_order=2,jump_direction=res_obj$jump_direction)
#     res[i]
#     },start,
#     method="Richardson",method.args = list(eps=1e-4))
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


# # 3 metrics
# # # # # my performance measures
eval_fit <- function(y,yhat) {
  RMSE <- sqrt(mean((y-yhat)^2))
  # NRMSE <- RMSE / (max(y)-min(y))
  Cor <- cor(y,yhat)
  MAE <- mean(abs(y-yhat))
  return(c(RMSE=RMSE,Cor=Cor,MAE=MAE))
}

# # 4 obj fcn and gradient

get_obj_fcn <- function(y_use,x_use,log_order,dum_use,jump_direction,init,rho_upper) {
  if (jump_direction =="none") {
    include_jumps <- FALSE
    jump_sign <- 0
  } else {
    include_jumps <- TRUE
    if (jump_direction=="up") {
      jump_sign <- 1
    } else {
      jump_sign <- -1
    }
  }
  
  if (log_order == 0) {
    resp <- y_use
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
      
      res <- resp - par["c1"] + exp(par["log_lam1"])*x_use - rec_filt
      grad["c1"] <- -2*mean(res)
      grad["log_lam1"] <- 2*mean(res*x_use)*exp(par["log_lam1"]) 
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
      names(start) <- c("c1","log_lam1","log_delta","phi")
      if(!is.null(init)) {
        stopifnot(length(init) == length(start)) 
        start[1:4] <- init
      } else {
        # see log_order=2 for details
        start_lm <- lm(resp~x_use)
        start[1:3] <- c(start_lm$coefficients[1],log(max(1e-10,-start_lm$coefficients[2])),log((max(start_lm$residuals) - min(start_lm$residuals))/2))
        start_rho <- optimize(function(phi) obj(c(start[1:3],phi=phi)),interval=c(-5,5))
        start[4] <- start_rho$minimum
      }
    } else {
      start <- numeric(2)
      names(start) <- c("c1","log_lam1")
      if(!is.null(init)) {
        stopifnot(length(init) == length(start)) 
        start[1:2] <- init
      } else {
        start_lm <- lm(resp~x_use)
        start[1:2] <- c(start_lm$coefficients[1],log(max(1e-10,-start_lm$coefficients[2])))
      }
    }
    
  } else {
    resp <- log(y_use)
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
        res <- resp - par["log_c1"] + exp(par["log_lam1"])*x_use - rec_filt
        grad["log_c1"] <- -2*mean(res)
        grad["log_lam1"] <- 2*mean(res*x_use)*exp(par["log_lam1"])
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
        names(start) <- c("log_c1","log_lam1","log_delta","phi")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:4] <- init
        } else {
          # see log_order=2 for details
          start_lm <- lm(resp~x_use)
          start[1:3] <- c(start_lm$coefficients[1],log(max(1e-10,-start_lm$coefficients[2])),log((max(start_lm$residuals) - min(start_lm$residuals))/2) )
          start_rho <- optimize(function(phi) obj(c(start[1:3],phi=phi)),interval=c(-5,5))
          start[4] <- start_rho$minimum
        }
      } else {
        start <- numeric(2)
        names(start) <- c("log_c1","log_lam1")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:2] <- init
        } else {
          start_lm <- lm(resp~x_use)
          start[1:2] <- c(start_lm$coefficients[1],log(max(1e-10,-start_lm$coefficients[2])))
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
        res <- resp - par["log_c1"] + exp(par["log_lam1"])*x_use - exp(par["log_c2"])*exp(-exp(par["log_lam2"])*x_use) - rec_filt
        grad["log_c1"] <- -2*mean(res)
        grad["log_lam1"] <- 2*mean(res*x_use)*exp(par["log_lam1"])
        grad["log_c2"] <- -2*mean(res*exp(-exp(par["log_lam2"])*x_use))*exp(par["log_c2"])
        grad["log_lam2"] <- 2*mean(res*x_use*exp(-exp(par["log_lam2"])*x_use)*exp(par["log_c2"]))*exp(par["log_lam2"])
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
        names(start) <- c("log_c1","log_lam1","log_c2","log_lam2","log_delta","phi")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:6] <- init
        } else {
          # do first one lm for c1 lam1,
          # then a second with x^2 for c2 lam2 
          # from 2nd order Taylor expansion of exp
          # beta_0 = log(c1)  + c2 -> c2
          # beta_2 = c2*lam2^2/2 -> lam2
          start_lm1 <- lm(resp~x_use)
          start[1] <- start_lm1$coefficients[1]
          start[2] <- log(-start_lm1$coefficients[2])
          start_lm2 <- lm(resp~x_use+I(x_use^2))
          start[3] <- log(max(1e-10,start_lm2$coefficients[1] - start_lm1$coefficients[1]))
          if (start_lm2$coefficients[3] / start[3] > 0) {
            start[4] <- log(sqrt(2*start_lm2$coefficients[3] / start[3]))
          } else {
            start[4] <- log(max(1e-10,-(start_lm2$coefficients[2] + start[2])/start[3]))
          }
          # delta from half range of residuals (right order of magnitude)
          start[5] <- log((max(start_lm2$residuals) - min(start_lm2$residuals))/2)
          # rho from optimal linesearch with current starting values
          start_rho <- optimize(function(phi) obj(c(start[1:5],phi=phi)),interval=c(-5,5))
          start[6] <- start_rho$minimum
        }
      } else {
        start <- numeric(4)
        names(start) <- c("log_c1","log_lam1","log_c2","log_lam2")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:4] <- init
        } else {
          start_lm1 <- lm(resp~x_use)
          start[1] <- start_lm1$coefficients[1]
          start[2] <- log(max(1e-10,-start_lm1$coefficients[2]))
          start_lm2 <- lm(resp~x_use+I(x_use^2))
          start[3] <- log(max(1e-10,start_lm2$coefficients[1] - start_lm1$coefficients[1]))
          if (start_lm2$coefficients[3] / start[3] > 0) {
            start[4] <- log(sqrt(2*start_lm2$coefficients[3] / start[3]))
          } else {
            start[4] <- log(max(1e-10,-(start_lm2$coefficients[2] + start[2])/start[3]))
          }
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
        exp_term <- exp(-exp(par["log_lam2"])*x_use + exp(par["log_c3"])*exp(-exp(par["log_lam3"])*x_use))
        res <- resp - par["log_c1"] + exp(par["log_lam1"])*x_use - 
          exp(par["log_c2"])*exp_term - rec_filt
        grad["log_c1"] <- -2*mean(res)
        grad["log_lam1"] <- 2*mean(res*x_use)*exp(par["log_lam1"])
        grad["log_c2"] <- -2*mean(res*exp_term)*exp(par["log_c2"])
        grad["log_lam2"] <- 2*mean(res*x_use*exp_term*exp(par["log_c2"]))*exp(par["log_lam2"])
        grad["log_c3"] <- -2*exp(par["log_c2"])*mean(res*exp(-exp(par["log_lam3"])*x_use)*exp_term)*exp(par["log_c3"])
        grad["log_lam3"] <- 2*exp(par["log_c2"])*exp(par["log_c3"])*mean(res*x_use*exp_term*exp(-exp(par["log_lam3"])*x_use))*exp(par["log_lam3"])
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
        names(start) <- c("log_c1","log_lam1","log_c2","log_lam2","log_c3","log_lam3","log_delta","phi")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:8] <- init
        } else {
          # do same as before, but set c3 = c2/2, lam3 = lam2/2, c2=c2/2, lam2 = lam2/2
          start_lm1 <- lm(resp~x_use)
          start[1] <- start_lm1$coefficients[1]
          start[2] <- log(max(1e-10,-start_lm1$coefficients[2]))
          start_lm2 <- lm(resp~x_use+I(x_use^2))
          start[3] <- log(max(1e-10,start_lm2$coefficients[1] - start_lm1$coefficients[1]))
          if (start_lm2$coefficients[3] / start[3] > 0) {
            start[4] <- log(sqrt(2*start_lm2$coefficients[3] / start[3]))
          } else {
            start[4] <- log(max(1e-10,-(start_lm2$coefficients[2] + start[2])/start[3]))
          }
          start[3:4] <- start[3:4] - log(2) # same as half of exp: log( exp(x)/2 = x - log(2)
          start[5:6] <- start[3:4]
          # delta from half range of residuals (right order of magnitude)
          start[7] <- log((max(start_lm2$residuals) - min(start_lm2$residuals))/2)
          # rho from optimal linesearch with current starting values
          start_rho <- optimize(function(phi) obj(c(start[1:7],phi=phi)),interval=c(-5,5))
          start[8] <- start_rho$minimum
        }
      } else {
        start <- numeric(6)
        names(start) <- c("log_c1","log_lam1","log_c2","log_lam2","log_c3","log_lam3")
        if(!is.null(init)) {
          stopifnot(length(init) == length(start)) 
          start[1:6] <- init
        } else {
          start_lm1 <- lm(resp~x_use)
          start[1] <- start_lm1$coefficients[1]
          start[2] <- log(max(1e-10,-start_lm1$coefficients[2]))
          start_lm2 <- lm(resp~x_use+I(x_use^2))
          start[3] <- log(max(1e-10,start_lm2$coefficients[1] - start_lm1$coefficients[1]))
          if (start_lm2$coefficients[3] / start[3] > 0) {
            start[4] <- log(sqrt(2*start_lm2$coefficients[3] / start[3]))
          } else {
            start[4] <- log(max(1e-10,-(start_lm2$coefficients[2] + start[2])/start[3]))
          }
          start[3:4] <- start[3:4] - log(2)
          start[5:6] <- start[3:4]
        }
      }
    }
  }
  return(list(obj=obj, grad=gr_obj, start=start))
}


