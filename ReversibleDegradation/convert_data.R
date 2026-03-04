library(readxl)

numdat    <- 26
all_log_y <- vector("list", numdat)

for (i in 1:numdat) {
  data           <- read_excel("./data/Dati_PERMANENT_reversible_all_ordered.xlsx", sheet = i)
  colnames(data) <- c("time", "current")
  data$time      <- (data$time - 1) * 60
  all_log_y[[i]] <- log(data$current)
}

saveRDS(all_log_y, "./data/all_log_y.rds")