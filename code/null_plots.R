library(tidyverse)
df <- readRDS("./output/nullsims/null_paramdf.RDS")
pval <- readRDS("./output/nullsims/null_pvalues.RDS")

# i <- sample(1:nrow(df), size = 1)
# pval[[i]] |>
#   as_tibble() |>
#   pivot_longer(cols = everything(), names_to = "method", values_to = "p") |>
#   ggplot(aes(sample = p, color = method)) +
#   geom_qq(distribution = qunif) +
#   geom_abline(slope = 1, intercept = 0, lty = 2) +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   theme_bw() +
#   ggtitle(label = i)
# paste(paste(names(df), "=",  df[i, ]), collapse = ", ")

for (i in seq_along(pval)) {
  pval[[i]] <- as_tibble(pval[[i]]) |>
    mutate(i = i)
}
pdf <- bind_rows(pval)
df$i <- 1:nrow(df)
pdf |>
  left_join(df, by = join_by(i)) |>
  pivot_longer(cols = c("segtest", "polymapR"), names_to = "method", values_to = "p") ->
  pdf

## All plots, for supplementary materials
for (rd_now in unique(pdf$rd)) {
  for (n_now in unique(pdf$n)) {
    for (ploidy_now in unique(pdf$ploidy)) {
      pdf |>
        filter(
          n == n_now,
          rd == rd_now,
          ploidy == ploidy_now,
          method == "segtest") |>
        mutate(i = as.factor(i)) |>
        ggplot(aes(sample = p, group = i)) +
        facet_wrap(~ p1 + p2) +
        geom_qq(distribution = qunif, geom = "line") +
        theme_bw() +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        theme(strip.background = element_rect(fill = "white")) ->
        pl_1
      ggsave(
        filename = paste0("./output/nullsims/qq_segtest_n", n_now, "_rd", rd_now, "_ploidy", ploidy_now, ".pdf"),
        plot = pl_1,
        height = 7,
        width = 6)
      pdf |>
        filter(
          n == n_now,
          rd == rd_now,
          ploidy == ploidy_now,
          method == "polymapR") |>
        mutate(i = as.factor(i)) |>
        ggplot(aes(sample = p, group = i)) +
        facet_wrap(~p1 + p2) +
        geom_qq(distribution = qunif, geom = "line") +
        theme_bw() +
        geom_abline(slope = 1, intercept = 0, color = "red") +
        theme(strip.background = element_rect(fill = "white")) ->
        pl_2
      ggsave(
        filename = paste0("./output/nullsims/qq_polymapr_n", n_now, "_rd", rd_now, "_ploidy", ploidy_now, ".pdf"),
        plot = pl_2,
        height = 7,
        width = 6)
    }
  }
}


## Example results from the null simulations. QQ plots
# (against the uniform distribution) of the p-values from segtest (blue)
# and polympaR (red). Plots are faceted by read-depth (columns) and
# the value of $\gamma_1$ (rows). This scenario had a ploidy of $K = 6$, with
# $\ell_1 = 0$ and $\ell_2 = 4$ and $n = 200$. Since the null is true,
# tests that control Type I error should lie at or above the y=x line (black).
# PolymapR's null assumption is only satisfied when $\gamma_1 = (1, 0)$
# (bottom row) and it can control type I error when genotypes are completely
# known (rd = Inf). However, polymapR cannot control type I error when
# there is not absolute preferential pairing (top row) because that scenario
# is not covered by the polymapR test. polymapR also does not control Type I
# error when genotypes are not known left column, even when their null
# scenario is satisfied ($\gamma_1 = (1,0))$. This is due to their ad hoc
# procedure of accounting for genotype uncertainty. The new segtest method
# controls Type I error in all scenarios.
pdf |>
  filter(ploidy == 6, p1 == 0, p2 == 4, n == 200, pi == 0) |>
  filter(
    (map_lgl(gamma2, \(x) all(x == c(0.5, 0.5)))) |
      (map_lgl(gamma2, \(x) all(x == c(1, 0))))
  ) |>
  mutate(gam = if_else(map_lgl(gamma2, \(x) all(x == c(0.5, 0.5))), "(0.5, 0.5)", "(1, 0)")) |>
  ggplot(aes(sample = p, color = method)) +
  geom_qq(distribution = qunif) +
  geom_abline(slope = 1, intercept = 0, color = "black") +
  facet_grid(gam ~ rd) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) ->
  pl
ggsave(filename = "./output/nullsims/qq_example.pdf", plot = pl, height = 4, width = 6)
