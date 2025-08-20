pacman::p_load(dplyr,tidyr,ggplot2,readxl)
source("./RFNR_model_functions.R")

hs <- c(200,200,200,
        200,200,200,
        200,200,
        100,80,200,200)
grad_cutoffs <- c(0.06,0.035,0.032,
                  0.035,0.035,0.035,
                  0.01,0.035,
                  0.02,0.028,0.035,0.04)

i <- 6
if (i<=8) {
  data <- read_excel("../data/Dati_PERMANENT_reversible_MEA_2.xlsx",sheet=i)
} else {
  data <- read_excel("../data/Dati_PERMANENT_reversible_RH_MEA_2.xlsx",sheet=i-8)
}

colnames(data) <- c("time","current")
# change to s
data$time <- (data$time-1)*60

if (i==1) { # remove higher level at end of experiment
  data <- filter(data, time < 200000)
}
mygrad <- function(x,maxh) {
  myind <- (maxh+1):(length(x)-maxh)
  grad <- sapply(myind,function(j) sum(x[j+1:maxh] - x[j-1:maxh])/(2*maxh))
  return(c(rep(NA,maxh),
           grad,
           rep(NA,maxh)))
}

# set correct h and cutoff
grad <- mygrad(log(data$current),hs[i])
grad_cutoff <- grad_cutoffs[i]
# if (i==9) {
#   grad <- -grad
# }
ind <- grad > grad_cutoff

# only take last time point of positive gradient regions
for (k in c(1:(length(ind)-1))) {
  if (!is.na(ind[k]*ind[k+1]) &ind[k] & ind[k+1]) {
    ind[k] <- FALSE
  }
}
data <- data %>% mutate("grad"=grad)

x <- data$time
y <- data$current
ndx <- seq(1, length(x), 60*40) 
myT <- length(ndx)
myT
# for each jump, closest observation gets the dummy=1
dum <- numeric(length(ndx))
for (k in which(ind)) {
  dum[which.min(abs(k - ndx))] <- 1
}
# compare log order 1 to lm
log_order <- 1
model_baseline <- RFNR(y,x,dum,ndx,log_order=log_order,jump_direction = "none")

model_lm <- lm(log(y)~I(-x),data=tibble(x=x[ndx],y=y[ndx]))
summary(model_lm)

model_bm <- nls(log(y)~ log_c1 - lam1*x,
                           data = tibble(x=x[ndx],y=y[ndx]),
                           start = list(log_c1=-2,lam1=1e-5))
mod_sum <- summary(model_baseline)$par_sig
model_baseline$sigma
summary(model_bm)
mod_sum
# very similar estimates, 


# compare log order 2 to nls
log_order <- 2
model_baseline <- RFNR(y,x,dum,ndx,log_order=log_order,jump_direction = "none")

model_bm <- nls(log(y)~ log_c1 - lam1*x + c2*exp(-lam2*x),
                data = tibble(x=x[ndx],y=y[ndx]),
                start = list(log_c1=-2,lam1=1e-5,c2=1,lam2=1e-5))

mod_sum <- summary(model_baseline)$par_sig
model_baseline$sigma
summary(model_bm)
mod_sum


# # # debugging test s in nls
# # at initialisation: all the same
# tmp_grad <- attr(rhs,"gradient")
# tmp_pars <- getPars()
# gradM_rhs <- grad_lin_pred_RFNR(tmp_pars,x[ndx],
#                                 dum_use=NULL,
#                                 log_order = 1,
#                                 rho_upper = 1,
#                                 jump_direction = "none")
# sum(abs(gradM_rhs-tmp_grad))
# cbind(gradM_rhs,tmp_grad)
# # same gradient
# 
# QR <- qr(gradM_rhs,tol=1e-12)
# var_pars <- chol2inv(QR$qr[,1:QR$rank])
# dimnames(var_pars) <- list(colnames(QR$qr)[1:QR$rank],colnames(QR$qr)[1:QR$rank])
# se <- sqrt(diag(var_pars))
# 
# XtXinv <- chol2inv(m$Rmat())
# se_nls <- sqrt(diag(XtXinv))
# 
# cbind(se,se_nls)
# sum(abs(se-se_nls))
# 
# sigma <- sqrt(sum((log(y[ndx])-lin_pred_RFNR(tmp_pars,x[ndx],
#                                              dum_use=NULL,
#                                              log_order = 1,
#                                              rho_upper = 1,
#                                              jump_direction = "none"))^2) / (length(ndx) - length(tmp_pars)))
# 
# sigma_nls <- sqrt(sum(resid^2)/(length(resid)-length(tmp_pars)))
# c(sigma,sigma_nls)
# abs(sigma-sigma_nls)
# 
# 
# # after convergence: all the same?
# tmp_grad <- nls.out$m$gradient()
# # tmp_pars <- nls.out$m$getPars()
# tmp_pars <- model_baseline$pars
# gradM_rhs <- grad_lin_pred_RFNR(tmp_pars,x[ndx],
#                                 dum_use=NULL,
#                                 log_order = 1,
#                                 rho_upper = 1,
#                                 jump_direction = "none")
# sum(abs(gradM_rhs-tmp_grad))
# cbind(gradM_rhs,tmp_grad)
# # similar gradient
# 
# QR <- qr(gradM_rhs,tol=1e-12)
# var_pars <- chol2inv(QR$qr[,1:QR$rank])
# dimnames(var_pars) <- list(colnames(QR$qr)[1:QR$rank],colnames(QR$qr)[1:QR$rank])
# se <- sqrt(diag(var_pars))
# 
# XtXinv <- chol2inv(nls.out$m$Rmat())
# se_nls <- sqrt(diag(XtXinv))
# 
# cbind(se,se_nls)
# sum(abs(se-se_nls))
# 
# sigma <- sqrt(sum((log(y[ndx])-lin_pred_RFNR(tmp_pars,x[ndx],
#                                              dum_use=NULL,
#                                              log_order = 1,
#                                              rho_upper = 1,
#                                              jump_direction = "none"))^2) / (length(ndx) - length(tmp_pars)))
# 
# sigma_nls <- sqrt(sum(nls.out$m$resid()^2)/(length(nls.out$m$resid())-length(tmp_pars)))
# c(sigma,sigma_nls)
# abs(sigma-sigma_nls)
# 
# se*sigma
# 
# model_baseline$sigma
# 
# undebug(nls)
# 
# getAnywhere(nlsModel)
# getAnywhere(summary.nls)
