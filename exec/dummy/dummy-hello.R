#!/usr/bin/env Rscript

##
# runner script example
#

devtools::load_all(".")


log <- function(...) {
  msg <- paste("#", ..., sep = " ")
  print(msg)
}

v <- function(...) cat(sprintf(...), "\n", sep = " ", file = stderr())
s <- function(...) do.call(paste, as.list(c(..., sep = ", ")))

scall <- function () {
  c(
    dmy_hello(),
    dmy_hello("Earth"),
    dmy_hello("Moon", "'Night")
  )
}

vcall <- function () {
  dmy_hello(c(
    "Mars",
    "Venus"
  ))
}

task <- function() {
  v("scall: %s", s(scall()))
  v("vcall: %s", s(vcall()))
  0
}

main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  log("> start:", paste(args, sep = " "))
  log("? args:", paste(commandArgs(), sep = ", "))
  rc <- 0
  print(elapsed <- system.time({
    rc <- task()
  }))
  log("< end:", rc, " -- ", summary(elapsed))
  rc
}

main()
