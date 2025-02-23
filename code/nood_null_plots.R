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

