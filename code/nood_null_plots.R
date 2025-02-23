library(tidyverse)
df <- readRDS("./output/nood_nullsims/nood_null_paramdf.RDS")
pval <- readRDS("./output/nood_nullsims/nood_null_pvalues.RDS")

i <- sample(1:nrow(df), size = 1)
pval[[i]] |>
  as_tibble() |>
  pivot_longer(cols = everything(), names_to = "method", values_to = "p") |>
  ggplot(aes(sample = p, color = method)) +
  geom_qq(distribution = qunif) +
  geom_abline(slope = 1, intercept = 0, lty = 2) +
  xlim(0, 1) +
  ylim(0, 1) +
  theme_bw() +
  ggtitle(label = i)
paste(paste(names(df), "=",  df[i, ]), collapse = ", ")

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
        theme(strip.background = element_rect(fill = "white"))
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
        theme(strip.background = element_rect(fill = "white"))
    }
  }
}
