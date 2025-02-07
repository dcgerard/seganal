## huge sample size under null to verify asymptotics
library(tidyverse)
library(segtest)
library(future)
library(foreach)
library(doRNG)
library(doParallel)

n <- 10000 ## sample size
nsamp <- 1000 ## number of replicates
df_4 <- expand_grid(
  ploidy = 4,
  p1 = 0:4,
  p2 = 0:4,
  model = c("seg", "auto", "auto_dr", "allo", "allo_pp", "auto_allo"),
  pi = c(0, 0.015, 0.03),
  alpha_1 = c(NA, 0, drbounds(4)/2, drbounds(4)),
  alpha_2 = c(NA, 0, drbounds(4)/2, drbounds(4)),
  gamma_1 = c(NA, 0, 1/2, 1),
  gamma_2 = c(NA, 0, 1/2, 1),
  beta_1 = c(NA, 0, beta_bounds(4)/2, beta_bounds(4)),
  beta_2 = c(NA, 0, beta_bounds(4)/2, beta_bounds(4))
)
df_4 |>
  filter(
    !(p1 %in% c(0, 4) & (!is.na(alpha_1) | !is.na(gamma_1) | !is.na(beta_1))),
    !(p2 %in% c(0, 4) & (!is.na(alpha_2) | !is.na(gamma_2) | !is.na(beta_2))),
    !(model == "seg" & (!is.na(alpha_1) | !is.na(alpha_2))),
    !(model == "seg" & !(p1 %in% c(1, 3)) & !is.na(beta_1)),
    !(model == "seg" & p1 %in% c(1, 3) & is.na(beta_1)),
    !(model == "seg" & !(p2 %in% c(1, 3)) & !is.na(beta_2)),
    !(model == "seg" & p2 %in% c(1, 3) & is.na(beta_2)),
    !(model == "seg" & p1 != 2 & !is.na(gamma_1)),
    !(model == "seg" & p1 == 2 & is.na(gamma_1)),
    !(model == "seg" & p2 != 2 & !is.na(gamma_2)),
    !(model == "seg" & p2 == 2 & is.na(gamma_2)),
    !(model == "auto" & (!is.na(gamma_1) | !is.na(gamma_2) | !is.na(beta_1) | !is.na(beta_2))),
    !(model == "auto" & p1 %in% 1:3 & alpha_1 > 0),
    !(model == "auto" & p2 %in% 1:3 & alpha_2 > 0),
    !(model == "auto" & p1 %in% 1:3 & is.na(alpha_1)),
    !(model == "auto" & p2 %in% 1:3 & is.na(alpha_2)),
    !(model %in% c("allo", "auto_allo") & (!is.na(alpha_1) | !is.na(alpha_2) | !is.na(beta_1) | !is.na(beta_2))),
    !(model %in% c("allo", "auto_allo") & p1 != 2 & !is.na(gamma_1)),
    !(model %in% c("allo", "auto_allo") & p2 != 2 & !is.na(gamma_2)),
    !(model %in% c("allo", "auto_allo") & p1 == 2 & !(gamma_1 %in% c(0, 1))),
    !(model %in% c("allo", "auto_allo") & p2 == 2 & !(gamma_2 %in% c(0, 1))),
    !(model == "auto_dr" & (!is.na(gamma_1) | !is.na(gamma_2) | !is.na(beta_1) | !is.na(beta_2))),
    !(model == "auto_dr" & p1 %in% c(1,2,3) & is.na(alpha_1)),
    !(model == "auto_dr" & p2 %in% c(1,2,3) & is.na(alpha_2)),
    !(model == "allo_pp" & (!is.na(beta_1) | !is.na(beta_2) | !is.na(alpha_1) | !is.na(alpha_2))),
    !(model == "allo_pp" & p1 != 2 & !is.na(gamma_1)),
    !(model == "allo_pp" & p1 == 2 & is.na(gamma_1)),
    !(model == "allo_pp" & p2 != 2 & !is.na(gamma_2)),
    !(model == "allo_pp" & p2 == 2 & is.na(gamma_2)),
  ) ->
  df_4

df_6 <- expand_grid(
  ploidy = 6,
  p1 = 0:6,
  p2 = 0:6,
  model = c("seg", "auto", "auto_dr", "allo", "allo_pp", "auto_allo"),
  pi = c(0, 0.015, 0.03),
  alpha_1 = c(NA, 0, drbounds(6)/2, drbounds(6)),
  alpha_2 = c(NA, 0, drbounds(6)/2, drbounds(6)),
  gamma_1 = c(NA, 0, 1/2, 1),
  gamma_2 = c(NA, 0, 1/2, 1),
  beta_1 = c(NA, 0, beta_bounds(6)/2, beta_bounds(6)),
  beta_2 = c(NA, 0, beta_bounds(6)/2, beta_bounds(6))
)
df_6 |>
  filter(
    !(p1 %in% c(0, 6) & (!is.na(alpha_1) | !is.na(gamma_1) | !is.na(beta_1))),
    !(p2 %in% c(0, 6) & (!is.na(alpha_2) | !is.na(gamma_2) | !is.na(beta_2))),
    !(model == "seg" & (!is.na(alpha_1) | !is.na(alpha_2))),
    !(model == "seg" & !(p1 %in% c(1, 5)) & !is.na(beta_1)),
    !(model == "seg" & p1 %in% c(1, 5) & is.na(beta_1)),
    !(model == "seg" & !(p2 %in% c(1, 5)) & !is.na(beta_2)),
    !(model == "seg" & p2 %in% c(1, 5) & is.na(beta_2)),
    !(model == "seg" & !(p1 %in% 2:4) & !is.na(gamma_1)),
    !(model == "seg" & p1 %in% 2:4 & is.na(gamma_1)),
    !(model == "seg" & !(p2 %in% 2:4) & !is.na(gamma_2)),
    !(model == "seg" & p2 %in% 2:4 & is.na(gamma_2)),
    !(model == "auto" & (!is.na(gamma_1) | !is.na(gamma_2) | !is.na(beta_1) | !is.na(beta_2))),
    !(model == "auto" & p1 %in% 1:5 & alpha_1 > 0),
    !(model == "auto" & p2 %in% 1:5 & alpha_2 > 0),
    !(model == "auto" & p1 %in% 1:5 & is.na(alpha_1)),
    !(model == "auto" & p2 %in% 1:5 & is.na(alpha_2)),
    !(model %in% c("allo", "auto_allo") & (!is.na(alpha_1) | !is.na(alpha_2) | !is.na(beta_1) | !is.na(beta_2))),
    !(model %in% c("allo", "auto_allo") & !(p1 %in% 2:4) & !is.na(gamma_1)),
    !(model %in% c("allo", "auto_allo") & !(p2 %in% 2:4) & !is.na(gamma_2)),
    !(model %in% c("allo", "auto_allo") & p1 %in% 2:4 & !(gamma_1 %in% c(0, 1))),
    !(model %in% c("allo", "auto_allo") & p2 %in% 2:4 & !(gamma_2 %in% c(0, 1))),
    !(model == "auto_dr" & (!is.na(gamma_1) | !is.na(gamma_2) | !is.na(beta_1) | !is.na(beta_2))),
    !(model == "auto_dr" & p1 %in% 1:5 & is.na(alpha_1)),
    !(model == "auto_dr" & p2 %in% 1:5 & is.na(alpha_2)),
    !(model == "allo_pp" & (!is.na(beta_1) | !is.na(beta_2) | !is.na(alpha_1) | !is.na(alpha_2))),
    !(model == "allo_pp" & !(p1 %in% 2:4) & !is.na(gamma_1)),
    !(model == "allo_pp" & p1 %in% 2:4 & is.na(gamma_1)),
    !(model == "allo_pp" & !(p2 %in% 2:4) & !is.na(gamma_2)),
    !(model == "allo_pp" & p2 %in% 2:4 & is.na(gamma_2)),
  ) ->
  df_6

df <- bind_rows(df_4, df_6)


registerDoParallel(cores = 12)
future::plan(future::multisession, workers = 12)

ret <- foreach(
  i = 1:nrow(df)
  ) %dorng% {
    ploidy <- df$ploidy[[i]]
    p1 <- df$p1[[i]]
    p2 <- df$p2[[i]]
    model <- df$model[[i]]
    pi <- df$pi[[i]]
    alpha_1 <- if (is.na(df$alpha_1[[i]])) NULL else df$alpha_1[[i]]
    alpha_2 <- if (is.na(df$alpha_2[[i]])) NULL else df$alpha_2[[i]]
    beta_1 <- if (is.na(df$beta_1[[i]])) NULL else df$beta_1[[i]]
    beta_2 <- if (is.na(df$beta_2[[i]])) NULL else df$beta_2[[i]]
    gamma_1 <- if (is.na(df$gamma_1[[i]])) NULL else df$gamma_1[[i]]
    gamma_2 <- if (is.na(df$gamma_2[[i]])) NULL else df$gamma_2[[i]]

    if (!is.null(gamma_1)) {
      gamma_1 <- c(gamma_1, 1 - gamma_1)
    }
    if (!is.null(gamma_2)) {
      gamma_2 <- c(gamma_2, 1 - gamma_2)
    }


    if (model == "seg") {
      type <- "mix"
      add_dr <- TRUE
      if (is.null(gamma_1)) {
        gamma_1 <- 1
      }
      if (is.null(gamma_2)) {
        gamma_2 <- 1
      }
    } else if (model == "allo_pp" || model == "allo" || model == "auto_allo") {
      type <- "mix"
      add_dr <- FALSE
      if (is.null(gamma_1)) {
        gamma_1 <- 1
      }
      if (is.null(gamma_2)) {
        gamma_2 <- 1
      }
    } else if (model == "auto_dr") {
      type <- "polysomic"
      add_dr <- FALSE
    } else if (model == "auto") {
      type <- "polysomic"
      add_dr <- FALSE
      alpha_1 <- 0
      alpha_2 <- 0
    }

    q <- gf_freq(
      p1_g = p1,
      p1_ploidy = ploidy,
      p1_gamma = gamma_1,
      p1_alpha = alpha_1,
      p1_beta = beta_1,
      p1_type = type,
      p1_add_dr = add_dr,
      p2_g = p2,
      p2_ploidy = ploidy,
      p2_gamma = gamma_2,
      p2_alpha = alpha_2,
      p2_beta = beta_2,
      p2_type = type,
      p2_add_dr = add_dr,
      pi = pi
    )


    pvec <- replicate(n = nsamp, expr = {
      x <- c(stats::rmultinom(n = 1, size = n, prob = q))
      seg_lrt(
        x = x,
        p1_ploidy = ploidy,
        p2_ploidy = ploidy,
        p1 = p1,
        p2 = p2,
        model = model,
        outlier = TRUE)$p_value
    })

    pvec
  }

future::plan(future::sequential)

attr(ret, "rng") <- NULL
attr(ret, "doRNG_version") <- NULL

saveRDS(object = ret, file = "../output/null_pvals.RDS")
