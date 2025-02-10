library(tidyverse)
plist <- readRDS("./output/null_pvals.RDS")
df <- read_csv("./output/null_df.csv") |>
  rename(index = `...1`)

# i <- sample(1:length(plist), size = 1)
# suppressWarnings(ksout <- ks.test(x = plist[[i]], y = stats::punif)$p.value)
# tibble(p = plist[[i]]) |>
#   ggplot(aes(sample = p)) +
#   geom_qq(distribution = stats::qunif) +
#   geom_abline(slope = 1, intercept = 0, colour = 2, linetype = 2, linewidth = 1) +
#   xlim(0, 1) +
#   ylim(0, 1) +
#   ggtitle(paste0("i=", i, ", KS: ", round(ksout, digits = 6)))
# df[i, ]

names(plist) <- seq_along(plist)
as_tibble(plist) |>
  pivot_longer(cols = everything(), names_to = "index", values_to = "pvalue") |>
  mutate(index = parse_number(index)) |>
  left_join(df, by = join_by(index)) ->
  longdf

model_vec <- unique(longdf$model)
ploidy_vec <- unique(longdf$ploidy)

for (i in seq_along(model_vec)) {
  for (j in seq_along(ploidy_vec)) {
    ploidy_now <- ploidy_vec[[j]]
    model_now <- model_vec[[i]]
    longdf |>
      filter(model == model_now, ploidy == ploidy_now) |>
      ggplot(aes(sample = pvalue, group = index)) +
      geom_qq(distribution = qunif, geom = "line") +
      geom_abline(slope = 1, intercept = 0, colour = 2, linetype = 2) +
      facet_grid(p1 ~ p2) ->
      pl
    ggsave(filename = paste0("./output/qqplots/qq_", model_now, "_", ploidy_now, ".pdf"), plot = pl, height = 12, width = 12)
  }
}


longdf |>
  filter(ploidy == 4, model == "seg", p1 == 1, p2 == 3) |>
  mutate(pi = as.factor(round(pi, digits = 4)),
         beta_1 = as.factor(round(beta_1, digits = 4)),
         beta_2 = as.factor(round(beta_2, digits = 4))) |>
  ggplot(aes(sample = pvalue, group = index, color = pi)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, colour = 1, linetype = 2) +
  facet_grid(beta_1 ~ beta_2)

longdf |>
  filter(ploidy == 6, model == "seg", p1 == 2, p2 == 2) |>
  mutate(pi = as.factor(round(pi, digits = 4)),
         gamma_1 = as.factor(round(gamma_1, digits = 4)),
         gamma_2 = as.factor(round(gamma_2, digits = 4))) |>
  ggplot(aes(sample = pvalue, group = index, color = pi)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, colour = 1, linetype = 2) +
  facet_grid(gamma_1 ~ gamma_2)

longdf |>
  filter(ploidy == 6, model == "auto_dr", p1 == 2, p2 == 2) |>
  mutate(pi = as.factor(round(pi, digits = 4)),
         alpha_1 = as.factor(round(alpha_1, digits = 4)),
         alpha_2 = as.factor(round(alpha_2, digits = 4))) |>
  ggplot(aes(sample = pvalue, group = index, color = pi)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, colour = 1, linetype = 2) +
  facet_grid(alpha_1 ~ alpha_2)

longdf |>
  filter(ploidy == 4, model == "auto_dr", p1 == 1, p2 == 1) |>
  mutate(pi = as.factor(round(pi, digits = 4)),
         alpha_1 = as.factor(round(alpha_1, digits = 4)),
         alpha_2 = as.factor(round(alpha_2, digits = 4))) |>
  ggplot(aes(sample = pvalue, group = index, color = pi)) +
  geom_qq(distribution = qunif, geom = "line") +
  geom_abline(slope = 1, intercept = 0, colour = 1, linetype = 2) +
  facet_grid(alpha_1 ~ alpha_2)

