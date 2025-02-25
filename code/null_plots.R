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

## All plots, for supplementary materials ----
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
        mutate(pgeno = paste0("(", p1, ",", p2, ")")) |>
        ggplot(aes(sample = p, group = i)) +
        facet_wrap(~ pgeno) +
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
        mutate(i = as.factor(i))  |>
        mutate(pgeno = paste0("(", p1, ",", p2, ")")) |>
        ggplot(aes(sample = p, group = i)) +
        facet_wrap(~ pgeno) +
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

## Example plot ----
colvec <- palette.colors(n = 4, palette = "Okabe-Ito")[2:3]
pdf |>
  filter(ploidy == 8, p1 == 0, p2 == 6, n == 200, pi == 0) |>
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
  theme(strip.background = element_rect(fill = "white")) +
  scale_color_manual(values = colvec) ->
  pl
ggsave(filename = "./output/nullsims/qq_example.pdf", plot = pl, height = 4, width = 6)


## Look at just type I error ----
alpha <- 0.05
pdf |>
  group_by(i, method, n, rd) |>
  summarize(
    t1e = mean(p < alpha),
    nreject = sum(p < alpha)) |>
  ungroup() ->
  sumdf

# sumdf |>
#   filter(method == "segtest") |>
#   arrange(desc(t1e)) |>
#   view()

nsamp <- nrow(pval[[1]])
upper <- qbinom(p = 0.995, size = nsamp, prob = alpha) / nsamp

sumdf |>
  mutate(n = as.factor(n)) |>
  ggplot(aes(x = t1e)) +
  geom_histogram(fill = "black", color = "black", bins = 100) +
  facet_grid(n + rd ~ method, scales = "free") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_vline(xintercept = alpha, colour = "red", lty = 2) +
  geom_vline(xintercept = upper, colour = "blue", lty = 3) ->
  pl

ggsave(filename = "./output/nullsims/null_t1e.pdf", plot = pl, height = 6, width = 4)
