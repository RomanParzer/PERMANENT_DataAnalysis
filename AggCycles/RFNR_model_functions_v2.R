# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # implementation of  Recursive Filter-Augmented Nonlinear Regression (RFNR) Model # # # # # # # # # # #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
pacman::p_load(optimx) 
source("./RFNR_help_functions.R")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # main function/class to fit model  # # # # # # # # # # # # 
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
  y_use <- y[ndx]
  jump_direction <- match.arg(jump_direction)
  if (jump_direction=="none") {
    dum_use <- NULL
  } else {
    stopifnot(!is.null(dum))
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
      dum_use <- NULL
    }
  }
  
  # generate obj, gradient and start vec depending on log_order
  myobjfcn <- get_obj_fcn(y_use,x_use,log_order,dum_use,jump_direction,init,rho_upper)
  obj <- myobjfcn$obj
  gr_obj <- myobjfcn$grad
  start <- myobjfcn$start
  
  
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
    y_use = y_use,
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

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # help functions for model class  # # # # # # # # # # # # #

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

summary.RFNR <-  function(res_obj,corr=FALSE,cov=c("HAC","standard")) {
  # new v2
  cov <- match.arg(cov)
  if (cov=="HAC") {
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
  } else {
    # # old
    gradM_rhs <- grad_lin_pred_RFNR(res_obj$pars,res_obj$x_use,
                                    dum_use=res_obj$dum_use,
                                    log_order = res_obj$log_order,
                                    rho_upper = res_obj$rho_upper,
                                    jump_direction = res_obj$jump_direction)
    # not the same as grad of obj fcn (score), but = res*gradM
    QR <- qr(gradM_rhs,tol=1e-12)
    var_pars <- chol2inv(QR$qr[,1:QR$rank])
    dimnames(var_pars) <- list(colnames(QR$qr)[1:QR$rank],colnames(QR$qr)[1:QR$rank])
    param <- as.numeric(res_obj$pars)
    se <- sqrt(diag(var_pars))*res_obj$sigma
    if (QR$rank < length(res_obj$pars)) {
      warning("There were parameters with no variance leading to degenerate covariance!")
    }
    zval <- param/se
    param <- cbind(param, se, zval, 2 * pnorm(as.numeric(abs(zval)), lower.tail = FALSE))
    dimnames(param) <- list(names(res_obj$pars), c("Estimate", "Std. Error", "z value", "Pr(>|z|)"))
    corr_mat <- NULL
    if (corr) {
      corr_mat <- var_pars*res_obj$sigma^2/outer(se, se)
    }
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
# res_obj$opm_res
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
