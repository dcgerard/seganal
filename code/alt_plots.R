library(tidyverse)
library(xtable)
df <- readRDS("./output/altsims/alt_paramdf.RDS")
pval <- readRDS("./output/altsims/alt_pvalues.RDS")

for (i in seq_along(pval)) {
  pval[[i]] <- as_tibble(pval[[i]])
  pval[[i]]$i <- i
}
df$i <- 1:nrow(df)

bind_rows(pval) |>
  left_join(df, by = join_by(i)) |>
  mutate(n = as.factor(n)) |>
  filter(rd == 10) |>
  select(-segtest_nout) |>
  pivot_longer(cols = c("segtest", "polymapR"), names_to = "method", values_to = "p") |>
  ggplot(aes(sample = p, color = method, lty = n)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(ploidy ~ scenario) +
  theme(strip.background = element_rect(fill = "white"))

bind_rows(pval) |>
  left_join(df, by = join_by(i)) |>
  mutate(n = as.factor(n)) |>
  filter(is.infinite(rd)) |>
  select(-segtest_nout) |>
  pivot_longer(cols = c("segtest", "polymapR"), names_to = "method", values_to = "p") |>
  ggplot(aes(sample = p, color = method, lty = n)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(ploidy ~ scenario) +
  theme(strip.background = element_rect(fill = "white"))

alpha <- 0.05
bind_rows(pval) |>
  left_join(df, by = join_by(i)) |>
  mutate(
    n = as.factor(n),
    ploidy = as.factor(ploidy),
    rd = as.factor(rd)) |>
  select(-segtest_nout) |>
  group_by(n, ploidy, rd, scenario) |>
  summarize(
    segtest = mean(segtest < alpha),
    polymapR = mean(polymapR < alpha)) |>
  xtable(
    label = "tab:alt.power",
    caption = "Results of the alternative simulations from Section \\ref{sec:alt.sims}. The table reports the power of \\texttt{segtest} and \\texttt{polymapR} under scenarios where segregation distortion is present (right two columns). Tests were performed at a nominal significance level of 0.05 across varying sample sizes ($n \\in \\{20, 200\\}$), ploidies ($K \\in \\{4, 6, 8\\}$), read depths (10 or infinite), and alternative scenarios (easy or hard). \\texttt{segtest} exhibits strong power for large $n$ and maintains reasonable power in most small $n$ scenarios. While \\texttt{segtest} has slightly lower power than \\texttt{polymapR} in most settings, this reduction is expected due to its ability to control Type I error, which \\texttt{polymapR} does not.") |>
  print(include.rownames = FALSE) |>
  cat(file = "./output/altsims/alt_tab.txt")
