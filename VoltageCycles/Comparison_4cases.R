
# compare data analysis of current given voltage and time (+cycles and cleanings)

# 1 temporal reversible degradation in each cycle

restab <-  readRDS("./saved_results/cycle_tabmodel_summary.rds")
restab_upl <-  readRDS("./saved_results/cycle_tabmodel_summary_upl.rds")
restab_lpl <-  readRDS("./saved_results/cycle_tabmodel_summary_lpl.rds")
restab_cycling <-  readRDS("./saved_results/cycle_tabmodel_summary_cycling.rds")

cycle_parameters <- data.frame(reference = restab[1:7,1],
                               upl = restab_upl[1:7,1],
                               lpl = restab_lpl[1:7,1],
                               cycling = restab_cycling[1:7,1])
cycle_parameters
require(kable)
kable(format(cycle_parameters,digits=3),format = "latex",booktabs=TRUE)

# test for difference
# e g , all differences are significant
z_val <- (restab[2,1]-restab_upl[2,1]) / sqrt(restab[2,2]^2+restab_upl[2,2]^2)
z_val
2*pnorm(abs(z_val),lower.tail = FALSE)


# 2 degradation after cleaning over time

mv_model <- readRDS("./saved_results/mv_model.rds")
mv_model_upl <-  readRDS("./saved_results/mv_model_upl.rds")
mv_model_lpl <-  readRDS("./saved_results/mv_model_lpl.rds")
mv_model_cycling <-  readRDS("./saved_results/mv_model_cycling.rds")

deg_parameters <- data.frame(reference = summary(mv_model)$coefficients[,1],
                               upl = summary(mv_model_upl)$coefficients[,1],
                               lpl = summary(mv_model_lpl)$coefficients[,1],
                               cycling = summary(mv_model_cycling)$coefficients[,1])
deg_parameters
kable(format(deg_parameters,digits=3),format = "latex",booktabs=TRUE)
# test for difference, on avg_hi
# e g , deg effect of time between ref and upl, not significant
z_val <- (summary(mv_model)$coefficients[4,1]-summary(mv_model_upl)$coefficients[4,1]) / sqrt(summary(mv_model)$coefficients[4,2]^2+summary(mv_model_upl)$coefficients[4,2]^2)
z_val
2*pnorm(abs(z_val),lower.tail = FALSE)

# e g , deg effect of time_clean between ref and cycl, not significant
z_val <- (summary(mv_model)$coefficients[7,1]-summary(mv_model_cycling)$coefficients[7,1]) / sqrt(summary(mv_model)$coefficients[7,2]^2+summary(mv_model_cycling)$coefficients[7,2]^2)
z_val
2*pnorm(abs(z_val),lower.tail = FALSE)

# e g , deg effect of time_clean between ref and upl, significant
z_val <- (summary(mv_model)$coefficients[7,1]-summary(mv_model_upl)$coefficients[7,1]) / sqrt(summary(mv_model)$coefficients[7,2]^2+summary(mv_model_upl)$coefficients[7,2]^2)
z_val
2*pnorm(abs(z_val),lower.tail = FALSE)


# reference case
# test for difference of coefficient for time_clean degradation between avg_hi and avg_lo

# use cov matrix of the one regression

# multiple comparison, adjust Bonferi alpha/k, k number of comparisons


