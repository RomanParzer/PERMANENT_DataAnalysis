# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # implementation of  Recursive Filter-Augmented Nonlinear Regression (RFNR) Model # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
pacman::p_load(optimx,ggplotify)

# to do: check all, adapt gradient by ... * exp(log thetaj)

# help function to build linear predictor for each log level
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

# help function gradient for significance
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # main function to fit model  # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# future: could add option to use standard COV matrix instead of HAC used here,
# but HAC is prob good for all time-series applications
RFNR <- function(y, x, 
                 dum = NULL,
                 ndx = 1:length(y),
                 log_order = 2,
                 rho_upper = 1,
                 jump_direction=c("up","down","none"),
                 init = NULL,
                 methods=c("nvm","subplex","mla","lbfgsb3c","spg"),
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
    if (sum(dum_use)==0) {
      warning("No jumps indicated; jump_direction is set to 'none'!")
      jump_direction <- "none"
      include_jumps <- FALSE
      dum_use <- NULL
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
  
  
  # try all available and look what works best
  opm_res <- opm(start,obj,gr=gr_obj,
                 method=methods,...) 
  
  # gr_obj(start)
  # num_gr <- numDeriv::grad(obj,start,method="Richardson",method.args = list(eps=1e-10))
  # num_gr
  # sum(abs(gr_obj(start)-num_gr))
  

  best_cand <- which(opm_res$kkt1 & opm_res$kkt2)
  if (length(best_cand)==0) {
    best_cand <- which(opm_res$kkt1 & opm_res$convergence==0)
  }
  if (length(best_cand)==0) {
    best_cand <- which(opm_res$convergence==0)
  }
  if (length(best_cand)==0) {
    warning("No converged optimization method, try a different solver or adjust controls!")
    best_cand <- 1:length(methods)
  }
  
  # # could also consider norm of found solutions, but just pick lowest obj value for now
  # values <- opm_res$value[best_cand]
  # norms <- apply(opm_res[best_cand,1:(ncol(opm_res)-8)],1,function(tmp_par) sqrt(mean(tmp_par^2)))
  # best_ind <- best_cand[which.min(values / mean(values^2) + norms/mean(norms^2))]
  best_ind <- best_cand[which.min(opm_res$value[best_cand])]
  pars <- opm_res[best_ind,1:(ncol(opm_res)-8)]
  
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
    start = start,
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
    # mygg <- acf(res_obj$log_residuals,main=" "))
    mygg <- forecast::ggAcf(res_obj$log_residuals) +
      labs(title=" ",x="lag",y="acf")
    mygg
  }
  mygg
}
coef.RFNR <-  function(res_obj) {
  res_obj$pars
}
residuals.RFNR <- function(res_obj,add_arg=NULL) {
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
          if (res_obj$jump_direction=="up") {
            jump_sign <- 1
          } else {
            jump_sign <- -1
          }
          jumps_use <- jump_sign*stats::filter(res_obj$dum_use*exp(as.numeric(res_obj$pars["log_delta"])), res_obj$rho_upper/(1+exp(-as.numeric(res_obj$pars["phi"]))), "recursive", init = 0)
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
  QR <- qr(gradM_rhs,tol=1e-16)
  var_pars <- chol2inv(QR$qr[,1:QR$rank])
  return(var_pars * length(res_obj$y_use))
}

summary.RFNR <-  function(res_obj,corr=FALSE) {
  # # # old
  # gradM_rhs <- grad_lin_pred_RFNR(res_obj$pars,res_obj$x_use,
  #                                 dum_use=res_obj$dum_use,
  #                                 log_order = res_obj$log_order,
  #                                 rho_upper = res_obj$rho_upper,
  #                                 jump_direction = res_obj$jump_direction)
  # # not the same as grad of obj fcn (score), but = res*gradM
  # QR <- qr(gradM_rhs,tol=1e-12)
  # var_pars <- chol2inv(QR$qr[,1:QR$rank])
  # dimnames(var_pars) <- list(colnames(QR$qr)[1:QR$rank],colnames(QR$qr)[1:QR$rank])
  # se <- sqrt(diag(var_pars))*res_obj$sigma
  # if (QR$rank < length(res_obj$pars)) {
  #   warning("There were parameters with no variance leading to degenerate covariance!")
  # }
  

  # HACcov <- sandwich::vcovHAC(res_obj,diagnostics=TRUE) 
  # uses by default cyclic weigths which can be negative
  
  # instead, use weave, weights based on residual acf
  HACcov <- sandwich::weave(res_obj,method="trunc",diagnostics=TRUE)
  HACcov <- as.matrix(Matrix::nearPD(HACcov)$mat)
  colnames(HACcov) <- row.names(HACcov) <- names(res_obj$pars)
  param <- as.numeric(res_obj$pars)
  se <- sqrt(diag(HACcov))
  zval <- param/se
  param <- cbind(param, se, zval, 2 * pnorm(as.numeric(abs(zval)), lower.tail = FALSE))
  dimnames(param) <- list(names(res_obj$pars), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
  corr_mat <- NULL
  if (corr) {
    corr_mat <- HACcov/outer(se, se)
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

# pacman::p_load(tidyverse,readxl)
# dt <- read_xlsx("../data/Dati_PERMANENT_reversible.xlsx", sheet = 1)
# # dt <- read_excel("../data/Dati_PERMANENT_reversible_RH.xlsx",sheet=1)
# names(dt) <- c("time", "current")
# dt$time <- (dt$time - 1)*60
# 
# ndx <- seq(1, nrow(dt), 60*10)
# x <- dt$time
# y <- dt$current
# summary(diff(log(dt$current[ndx])))
# dum <- c(0, diff(log(dt$current[ndx])) > 0.1)
# # dum <- c(0, diff(log(dt$current[ndx])) < -0.1)
# # dum <- c(0, abs(diff(log(dt$current[ndx]))) > 0.1)
# 
# 
# xnew <- x[-ndx]
# time_ind <- 1 + (1:length(x))[-ndx] * (length(ndx)-1) / (max(ndx))
# 
# log_order <- 2
# # for debugging:
# # jump_direction <- "up"
# # rho_upper <- 1
# # methods <- c("Rvmmin","nvm")
# # init <- NULL
# 
# res_obj <- RFNR(y,x,dum=dum,ndx=ndx,log_order=log_order,jump_direction = "up",
#             init=c(-1,-10,-1,-10,0,1),
#             methods=c("nvm","subplex","mla","lbfgsb3c","spg"))
# 
# res <- RFNR(y,x,dum=dum,ndx=ndx,log_order=log_order,jump_direction = "up",
#             init=c(-1,-10,-1,-10,0,1),
#             methods=c("nvm","Rvmmin",
#                       "nlminb","ucminf", "subplex", "ncg",
#                       "Rcgmin","mla",
#                       "lbfgsb3c", "Rtnmin","spg"))
# res$dum_use
# res_obj <- res
# coef(res)
# res$start
# res$opm_res
# which.min(res$opm_res$value)
# summary(res,corr=TRUE)
# res
# 
# # best candidates (lowest value and kkts): given init / automatic init
# # i=1: nvm, Rvmmin, nlminb, ucminf,  subplex, ncg, Rcgmin, mla, lbfgsb3c , Rtnmin, spg
# # i=2 subplex seems best
# # i=8: eg mla lbfgsb3c
# # RH100 spg, nvm
# 
# plot(res) +
#   labs(x="time in s")
# 
# plot(res,"residuals") +
#   labs(x="time in s")
# 

# eval_fit(log(y[-ndx]),predict(res,xnew = xnew,time_ind = time_ind,type = "linear"))
# eval_fit(log(y[ndx]),predict(res,xnew = x[ndx],time_ind = 1:length(ndx),type = "linear"))
# res$sigma*sqrt((length(ndx)-6)/length(ndx))



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
