library(tidyverse)
library(xtable)
df <- readRDS("./output/dr_nullsims/dr_null_paramdf.RDS")
pval <- readRDS("./output/dr_nullsims/dr_null_pvalues.RDS")

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
  rename(segtestNO = segtest_nout) |>
  pivot_longer(cols = c("segtest", "segtestNO", "polymapR"), names_to = "method", values_to = "p") ->
  pdf

alpha <- 0.05
nsamp <- nrow(pval[[1]])
upper <- qbinom(p = 0.995, size = nsamp, prob = alpha) / nsamp

colvec <- palette.colors(n = 7, palette = "Okabe-Ito")[c(2, 3, 4)]

pdf |>
  group_by(n, ploidy, rd, p1, p2, dr, method) |>
  summarize(t1e = mean(p < alpha)) |>
  ungroup() |>
  mutate(dr = parse_factor(dr, levels = c("low", "med", "high")),
         method = parse_factor(method, levels = c("polymapR", "segtestNO", "segtest")),
         ploidy = as.factor(ploidy)) ->
  sumdf

sumdf |>
  mutate(dr = recode(dr, "low" = "none")) |>
  ggplot(aes(x = t1e, fill = ploidy)) +
  facet_grid(method ~ dr) +
  geom_histogram(bins = 30) +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  geom_vline(xintercept = alpha, col = "red", lty = 2) +
  geom_vline(xintercept = upper, colour = "blue", lty = 3) +
  scale_x_sqrt() +
  scale_fill_manual(values = colvec) +
  xlab("Type I Error Rate") ->
  pl
ggsave(filename = "./output/dr_nullsims/dr_null_hist.pdf", plot = pl, height = 5, width = 6)

sumdf |>
  filter(dr == "med", method == "segtest") |>
  arrange(desc(t1e)) |>
  filter(t1e > upper) |>
  mutate(n = factor(n)) |>
  select(-dr, -method) |>
  xtable(
    label = "tab:bad.dr.sims",
    caption = "Eight simulation scenarios for \\texttt{segtest} with poor type I error control under moderate double reduction (Section \\ref{sec:dr.sims}). All scenarios involve duplex-by-nullplex crosses, known genotypes, and large sample sizes."
    ) |>
  print(include.rownames = FALSE) |>
  cat(file = "./output/dr_nullsims/dr_bad_tab.txt")

sumdf |>
  filter(dr == "high", method == "segtest") |>
  arrange(desc(t1e)) |>
  filter(t1e > upper)

sumdf |>
  filter(method == "segtest") |>
  group_by(dr) |>
  summarize(nbad = sum(t1e > upper))
