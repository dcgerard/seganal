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
  p1 = 0:8,
  p2 = 0:8,
  dr = c("low", "med", "high")
)
df |>
  filter(
    !(ploidy == 4 & (p1 > 4 | p2 > 4)),
    !(ploidy == 6 & (p1 > 6 | p2 > 6)),
    p2 > 1,
    p2 < (ploidy - 1)
  ) ->
  df
df$seed <- 1:nrow(df)

## Reorder rows for more even parallelization
set.seed(1)
df |>
  slice_sample(n = nrow(df), replace = FALSE) ->
  df

ret <- foreach(
  i = 1:nrow(df),
  .export = c("seg_lrt", "polymapr_test")
  ) %dorng% {
    set.seed(df$seed[[i]])
    ploidy <- df$ploidy[[i]]
    n <- df$n[[i]]
    rd <- df$rd[[i]]
    p1 <- df$p1[[i]]
    p2 <- df$p2[[i]]
    dr <- df$dr[[i]]

    if (dr == "low") {
      alpha <- rep(0, times = floor(ploidy / 4))
    } else if (dr == "med") {
      alpha <- drbounds(ploidy = ploidy) / 2
    } else if (dr == "high") {
      alpha <- drbounds(ploidy = ploidy)
    } else {
      stop("duh")
    }

    gf <- gf_freq(
      p1_g = p1,
      p1_ploidy = ploidy,
      p1_gamma = NULL,
      p1_alpha = alpha,
      p1_beta = NULL,
      p1_type = "polysomic",
      p1_add_dr = TRUE,
      p2_g = p2,
      p2_ploidy = ploidy,
      p2_gamma = NULL,
      p2_alpha = alpha,
      p2_beta = NULL,
      p2_type = "polysomic",
      p2_add_dr = TRUE,
      pi = 0,
      nudge = 0)


    psims <- replicate(n = nsamp, expr = {
      x <- c(stats::rmultinom(n = 1, size = n, prob = gf))
      if (is.finite(rd)) {
        x <- simgl(nvec = x, rd = rd)
        pl <- t(apply(X = x, MARGIN = 1, FUN = \(x) exp(x - updog::log_sum_exp(x))))
      } else {
        pl <- x
      }
      sout_nout <- seg_lrt(
        x = x,
        p1_ploidy = ploidy,
        p2_ploidy = ploidy,
        p1 = p1,
        p2 = p2,
        model = "seg",
        outlier = FALSE)
      sout <- seg_lrt(
        x = x,
        p1_ploidy = ploidy,
        p2_ploidy = ploidy,
        p1 = p1,
        p2 = p2,
        model = "seg",
        outlier = TRUE)
      pout <- polymapr_test(
        x = pl,
        g1 = p1,
        g2 = p2,
        type = "segtest")
      c(segtest = sout$p_value, segtest_nout = sout_nout$p_value, polymapR = pout$p_value)
    })

    t(psims)
  }

## Unregister workers ----
if (nc > 1) {
  plan(sequential)
}

saveRDS(object = df, file = "./output/dr_nullsims/dr_null_paramdf.RDS")
saveRDS(object = ret, file = "./output/dr_nullsims/dr_null_pvalues.RDS")
