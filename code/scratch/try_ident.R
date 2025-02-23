## Get one example of non-identifiability
library(segtest)
ploidy <- 8
p1 <- 6
p2 <- 3
q1 <- gf_freq(
  p1_g = p1,
  p1_ploidy = ploidy,
  p1_gamma = c(1, 0),
  p1_alpha = NULL,
  p1_beta = NULL,
  p1_type = "mix",
  p1_add_dr = TRUE,
  p2_g = p2,
  p2_ploidy = ploidy,
  p2_gamma = c(1, 0),
  p2_alpha = NULL,
  p2_beta = NULL,
  p2_type = "mix",
  p2_add_dr = TRUE,
  pi = 0.03)

nsamp <- 1000
nvec <- c(stats::rmultinom(n = 1, size = nsamp, prob = q1))
sout <- seg_lrt(
  x = nvec,
  p1_ploidy = ploidy,
  p2_ploidy = ploidy,
  p1 = p1,
  p2 = p2,
  model = "seg",
  outlier = TRUE)
sout$null$df0
sout$null$gam[[1]]$gamma
sout$null$gam[[2]]$gamma
sout$null$gam[[3]]$pi

gamma1 <- c(0.976975533195293, 0.0230244668047073)
gamma2 <- c(0.977010812194103, 0.0229891878058969)
pi <- 0.0273334599991859

### New way
library(segtest)
ploidy <- 8
p1 <- 6
p2 <- 3
gamma1 <- c(0.97, 0.03)
gamma2 <- c(0.97, 0.03)
pi <- 0.027
q1 <- gf_freq(
  p1_g = p1,
  p1_ploidy = ploidy,
  p1_gamma = gamma1,
  p1_alpha = NULL,
  p1_beta = NULL,
  p1_type = "mix",
  p1_add_dr = TRUE,
  p2_g = p2,
  p2_ploidy = ploidy,
  p2_gamma = gamma2,
  p2_alpha = NULL,
  p2_beta = NULL,
  p2_type = "mix",
  p2_add_dr = TRUE,
  pi = pi)


obj <- function(pi, q1) {
q2 <- gf_freq(
  p1_g = p1,
  p1_ploidy = ploidy,
  p1_gamma = c(0.97, 0.03),
  p1_alpha = NULL,
  p1_beta = NULL,
  p1_type = "mix",
  p1_add_dr = TRUE,
  p2_g = p2,
  p2_ploidy = ploidy,
  p2_gamma = c(0.975, 0.025),
  p2_alpha = NULL,
  p2_beta = NULL,
  p2_type = "mix",
  p2_add_dr = TRUE,
  pi = pi)
  sum((q1 - q2)^2)
}

oout <- optim(
  par = 0.01,
  fn = obj,
  method = "Brent",
  lower = 0,
  upper = 0.03,
  q1 = q1)
pi <- oout$par
pi

