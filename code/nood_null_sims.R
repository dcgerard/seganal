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

## set up simulation scenario
nsamp <- 200
df_4 <- expand_grid(
  ploidy = 4,
  p1 = 0:4,
  p2 = 0:4,
  n = c(20, 200),
  rd = c(10, Inf),
  pi = 0,
  gamma1 = list(NA, c(1, 0), c(0.5, 0.5), c(0, 1)),
  gamma2 = list(NA, c(1, 0), c(0.5, 0.5), c(0, 1)),
  beta1 = c(NA, 0, beta_bounds(4)/2, beta_bounds(4)),
  beta2 = c(NA, 0, beta_bounds(4)/2, beta_bounds(4))
)
df_4 |>
  filter(
    p2 >= p1,
    xor(p1 == 2, is.na(gamma1)),
    xor(p2 == 2, is.na(gamma2)),
    xor(p1 %in% c(1, 3), is.na(beta1)),
    xor(p2 %in% c(1, 3), is.na(beta2)),
    ) ->
  df_4

df_6 <- expand_grid(
  ploidy = 6,
  p1 = 0:6,
  p2 = 0:6,
  n = c(20, 200),
  rd = c(10, Inf),
  pi = 0,
  gamma1 = list(NA, c(1, 0), c(0.5, 0.5), c(0, 1)),
  gamma2 = list(NA, c(1, 0), c(0.5, 0.5), c(0, 1)),
  beta1 = c(NA, 0, beta_bounds(6)/2, beta_bounds(6)),
  beta2 = c(NA, 0, beta_bounds(6)/2, beta_bounds(6))
)
df_6 |>
  filter(
    p2 >= p1,
    xor(p1 %in% 2:4, is.na(gamma1)),
    xor(p2 %in% 2:4, is.na(gamma2)),
    xor(p1 %in% c(1, 5), is.na(beta1)),
    xor(p2 %in% c(1, 5), is.na(beta2)),
    ) ->
  df_6
df_8_1 <- expand_grid(
  ploidy = 8,
  p1 = c(0:3,5:8),
  p2 = c(0:3,5:8),
  n = c(20, 200),
  rd = c(10, Inf),
  pi = 0,
  gamma1 = list(NA, c(1, 0), c(0.5, 0.5), c(0, 1)),
  gamma2 = list(NA, c(1, 0), c(0.5, 0.5), c(0, 1)),
  beta1 = c(NA, 0, beta_bounds(8)/2, beta_bounds(8)),
  beta2 = c(NA, 0, beta_bounds(8)/2, beta_bounds(8))
)
df_8_1 |>
  filter(
    p2 >= p1,
    xor(p1 %in% c(2,3,5,6), is.na(gamma1)),
    xor(p2 %in% c(2,3,5,6), is.na(gamma2)),
    xor(p1 %in% c(1, 7), is.na(beta1)),
    xor(p2 %in% c(1, 7), is.na(beta2)),
    ) ->
  df_8_1
df_8_2 <- expand_grid(
  ploidy = 8,
  p1 = 4,
  p2 = 4,
  n = c(20, 200),
  rd = c(10, Inf),
  pi = 0,
  gamma1 = list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 0)/2, c(1, 0, 1)/2, c(0, 1, 1)/2, c(1, 1, 1)/3),
  gamma2 = list(c(1, 0, 0), c(0, 1, 0), c(0, 0, 1), c(1, 1, 0)/2, c(1, 0, 1)/2, c(0, 1, 1)/2, c(1, 1, 1)/3),
  beta1 = NA,
  beta2 = NA)

df <- bind_rows(df_4, df_6, df_8_1, df_8_2)

df |>
  mutate(
    gamma1 = if_else(!is.na(beta1), list(1), gamma1),
    gamma2 = if_else(!is.na(beta2), list(1), gamma2)) ->
  df

df$seed <- 1:nrow(df)

ret <- foreach(
  i = 1:nrow(df),
  .export = c("gf_freq", "seg_lrt", "polymapr_test")
  ) %dorng% {
    set.seed(df$seed[[i]])
    ploidy <- df$ploidy[[i]]
    p1 <- df$p1[[i]]
    p2 <- df$p2[[i]]
    pi <- df$pi[[i]]
    beta1 <- if (is.na(df$beta1[[i]])) NULL else df$beta1[[i]]
    beta2 <- if (is.na(df$beta1[[i]])) NULL else df$beta1[[i]]
    gamma1 <- if (is.na(df$gamma1[[i]][[1]])) NULL else df$gamma1[[i]]
    gamma2 <- if (is.na(df$gamma2[[i]][[1]])) NULL else df$gamma2[[i]]
    rd <- df$rd[[i]]
    n <- df$n[[i]]

    gf <- gf_freq(
      p1_g = p1,
      p1_ploidy = ploidy,
      p1_gamma = gamma1,
      p1_alpha = NULL,
      p1_beta = beta1,
      p1_type = "mix",
      p1_add_dr = TRUE,
      p2_g = p2,
      p2_ploidy = ploidy,
      p2_gamma = gamma2,
      p2_alpha = NULL,
      p2_beta = beta2,
      p2_type = "mix",
      p2_add_dr = TRUE,
      pi = pi,
      nudge = 0)


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
        p1 = p1,
        p2 = p2,
        model = "seg",
        outlier = FALSE)
      pout <- polymapr_test(
        x = pl,
        g1 = p1,
        g2 = p2,
        type = "segtest")
      c(segtest = sout$p_value, polymapR = pout$p_value)
    })

    t(psims)
  }

## Unregister workers ----
if (nc > 1) {
  plan(sequential)
}

saveRDS(object = df, file = "./output/nood_nullsims/nood_null_paramdf.RDS")
saveRDS(object = ret, file = "./output/nood_nullsims/nood_null_pvalues.RDS")
