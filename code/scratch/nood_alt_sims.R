library(tidyverse)
library(segtest)
library(future)
library(doFuture)
library(foreach)
library(doRNG)
library(doParallel)

## Set up cluster
registerDoFuture()
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  nc <- 1
} else {
  eval(parse(text = args[[1]]))
}
if (nc == 1) {
  plan(sequential)
} else {
  plan(multisession, workers = nc)
}

## Simulations
nsamp <- 200
df <- expand_grid(
  n = c(20, 200),
  ploidy = c(4, 6, 8),
  rd = c(10, Inf),
  scenario = c("hard", "easy")
)
df$seed <- 1:nrow(df)

ret <- foreach(
  i = 1:nrow(df),
  .export = c("seg_lrt", "polymapr_test")
  ) %dorng% {
    set.seed(df$seed[[i]])
    ploidy <- df$ploidy[[i]]
    scenario <- df$scenario[[i]]
    n <- df$n[[i]]
    rd <- df$rd[[i]]

    if (scenario == "easy") {
      gf <- stats::runif(n = ploidy + 1)
      gf <- gf / sum(gf)
    } else if (scenario == "hard") {
      p1 <- stats::runif(n = ploidy / 2 + 1)
      p2 <- stats::runif(n = ploidy / 2 + 1)
      gf <- stats::convolve(x = p1, y = rev(p2), type = "open")
      gf[gf < 0] <- 0 ## for -1e-16
      gf <- gf / sum(gf)
    } else {
      stop("dub")
    }


    psims <- replicate(n = nsamp, expr = {
      x <- c(stats::rmultinom(n = 1, size = n, prob = gf))
      if (is.finite(rd)) {
        x <- simgl(nvec = x, rd = rd)
        pl <- t(apply(X = x, MARGIN = 1, FUN = \(x) exp(x - updog::log_sum_exp(x))))
      } else {
        pl <- x
      }
      sout <- seg_lrt(
        x = x,
        p1_ploidy = ploidy,
        p2_ploidy = ploidy,
        p1 = NULL,
        p2 = NULL,
        model = "seg",
        outlier = FALSE)
      pout <- polymapr_test(
        x = pl,
        g1 = NULL,
        g2 = NULL,
        type = "segtest")
      c(segtest = sout$p_value, polymapR = pout$p_value)
    })

    t(psims)
  }

## Unregister workers ----
if (nc > 1) {
  plan(sequential)
}

saveRDS(object = df, file = "./output/nood_altsims/nood_alt_paramdf.RDS")
saveRDS(object = ret, file = "./output/nood_altsims/nood_alt_pvalues.RDS")
