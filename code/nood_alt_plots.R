library(tidyverse)
df <- readRDS("./output/nood_altsims/nood_alt_paramdf.RDS")
pval <- readRDS("./output/nood_altsims/nood_alt_pvalues.RDS")

for (i in seq_along(pval)) {
  pval[[i]] <- as_tibble(pval[[i]])
  pval[[i]]$i <- i
}
df$i <- 1:nrow(df)

bind_rows(pval) |>
  left_join(df, by = join_by(i)) |>
  mutate(n = as.factor(n)) |>
  filter(rd == 10) |>
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
  pivot_longer(cols = c("segtest", "polymapR"), names_to = "method", values_to = "p") |>
  ggplot(aes(sample = p, color = method, lty = n)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  facet_grid(ploidy ~ scenario) +
  theme(strip.background = element_rect(fill = "white"))

