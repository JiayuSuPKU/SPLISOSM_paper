library(tidyverse)
extrafont::loadfonts()

res_dir <- "~/Projects/SPLISOSM_paper/results/benchmark"

### Fig 2: Simulation benchmark results
## Fig 2A: Visualize simulation data and design
data_general <- read_csv(
  file.path(res_dir, "figures/source_data/sim_data_general_six.csv")
) %>%
  mutate(
    group_gene = factor(group_gene, levels = c("none", "donut"), labels = c("Random", "Regional")),
    group_iso = factor(group_iso, levels = c("none", "mvn", "donut"), labels = c("Random", "GP", "Regional")),
  )

# total gene expression pattern
m1.1 <- data_general %>%
  filter(group_iso == "Random") %>%
  ggplot(aes(x = x, y = y, fill = total_count)) +
  facet_wrap(~group_gene, nrow = 2, strip.position = "top") +
  labs(title = "Gene expression", fill = "Counts") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#628D56") +
  theme_void() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 0, 2, 0, "pt")),
    # strip.text = element_text(size = 12, family = "Arial", angle = 90, margin = margin(0,2,0,2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

# expected isoform ratio
m1.2 <- data_general %>%
  filter(group_gene == "Random") %>%
  select(c(x, y, group_gene, group_iso, starts_with("exp_ratio"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "exp_ratio_"), starts_with("exp_ratio")) %>%
  # capitalize the first letter of the isoform name
  pivot_longer(cols = starts_with("iso"), names_to = "isoform", values_to = "value") %>%
  mutate(isoform = str_to_title(isoform)) %>%
  # visualization
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_grid(isoform ~ group_iso, switch = "y") +
  labs(title = "Expected isoform usage ratio", fill = "Ratio") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#155289") +
  theme_void() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

# covariate design
m1.3 <- data_general %>%
  filter(
    (group_gene == "Regional") & (group_iso == "Regional")
  ) %>%
  select(c(x, y, group_gene, group_iso, starts_with("covar_"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "covar_"), starts_with("covar_")) %>%
  pivot_longer(cols = c(C1, C2, B1, B2), names_to = "covariate", values_to = "value") %>%
  mutate(covariate = factor(covariate, levels = c("B1", "B2", "C1", "C2"))) %>%
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_wrap(~covariate, nrow = 2, strip.position = "top") +
  labs(title = "Covariate design") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#94221F") +
  theme_void() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
    legend.position = "none"
  )

# observed isoform ratio
m1.4 <- data_general %>%
  select(c(x, y, group_gene, group_iso, starts_with("obs_ratio"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "obs_ratio_"), starts_with("obs_ratio")) %>%
  # capitalize the first letter of the isoform name
  pivot_longer(cols = starts_with("iso"), names_to = "isoform", values_to = "value") %>%
  mutate(
    isoform = str_to_title(isoform),
    G = group_gene,
    I = group_iso
  ) %>%
  # visualization
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_grid(isoform ~ G + I, switch = "y", labeller = labeller(.cols = label_both)) +
  labs(fill = "Observed ratio") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#155289") +
  theme_void() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

m1 <- cowplot::plot_grid(
  m1.1, NULL, m1.2, NULL, m1.3, NULL, m1.4,
  nrow = 1, align = "h", axis = "tbrl",
  rel_widths = c(0.55, 0.05, 1.1, 0.05, 0.75, 0.1, 2.2)
)
m1

ggsave(file.path(res_dir, "figures/sim_general_six.png"), m1, width = 16, height = 3, dpi = 300)
ggsave(file.path(res_dir, "figures/sim_general_six.pdf"), m1, width = 16, height = 3, dpi = 300)

## Fig 2B: SV test results of the general six scenarios
df_sv_res <- read_csv(file.path(res_dir, "sv_results/sv_general_six_1205.csv")) %>%
  mutate(
    group_gene = factor(group_gene, levels = c("none", "donut"), labels = c("Random", "Regional")),
    group_iso = factor(group_iso, levels = c("none", "mvn", "donut"), labels = c("Random", "GP", "Regional")),
  ) %>%
  pivot_longer(cols = starts_with("pvalue_"), names_to = "sv_test", values_to = "pval") %>%
  mutate(
    sv_test = str_remove(sv_test, "pvalue_"),
    sv_test = factor(
      sv_test,
      levels = c("spark-x", "hsic-gc", "hsic-ic", "hsic-ir"),
      labels = c("SPARK-X", "HSIC-GC", "HSIC-IC", "HSIC-IR")
    ),
  ) %>%
  # add expected pvalue under the null for each test by ranking pval
  group_by(group_gene, group_iso, sv_test) %>%
  mutate(
    pval_null = rank(pval),
    pval_null = pval_null / n(),
    Gene = group_gene,
    Isoform = group_iso
  )

m2 <- ggplot(
  df_sv_res,
  aes(x = pval_null, y = pval, group = sv_test, color = sv_test)
) +
  facet_wrap(~ Gene + Isoform, nrow = 1, scales = "free", labeller = labeller(.cols = label_both)) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  geom_line(linewidth = 1) +
  labs(x = "Expected p-value under the null", y = "Observed p-value", color = "Spatial variability test") +
  scale_x_continuous(transform = c("log10", "reverse")) +
  scale_y_continuous(transform = c("log10", "reverse")) +
  scale_color_manual(values = c(
    "HSIC-GC" = "#5AAA46", "HSIC-IR" = "#317EC2",
    "HSIC-IC" = "#825CA6", "SPARK-X" = "#F2B342"
  )) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
    legend.position = "right"
  )
m2

ggsave(file.path(res_dir, "figures/sv_general_six.png"), m2, width = 16, height = 3, dpi = 300)
ggsave(file.path(res_dir, "figures/sv_general_six.pdf"), m2, width = 16, height = 3, dpi = 300)

## Fig 2C: DU test results of the general six scenarios
# load and conbine du test results
du_np <- read_csv(file.path(res_dir, "du_results/du_np_general_six_1205.csv"))
du_p <- read_csv(file.path(res_dir, "du_results/du_p_general_six_1205.csv")) %>%
  mutate(gene = paste("gene", str_remove(gene, "gene_0+"), sep = "_"))
df_du_res <- full_join(
  du_np, du_p,
  by = join_by(gene, covariate, covar_type, group_gene, group_iso)
) %>%
  mutate(
    group_gene = factor(group_gene, levels = c("none", "donut"), labels = c("Random", "Regional")),
    group_iso = factor(group_iso, levels = c("none", "mvn", "donut"), labels = c("Random", "GP", "Regional")),
  ) %>%
  pivot_longer(cols = starts_with("pvalue_"), names_to = "du_test", values_to = "pval") %>%
  mutate(
    du_test = str_remove(du_test, "pvalue_"),
    du_test = factor(
      du_test,
      levels = c("t-fisher", "hsic", "hsic-gp", "glm", "glmm"),
      labels = c("T-test (Fisher's)", "Unconditional HSIC", "Conditional HSIC", "GLM", "GLMM")
    ),
  ) %>%
  # add expected pvalue under the null for each test by ranking pval
  group_by(group_gene, group_iso, du_test, covar_type) %>%
  mutate(
    pval_null = rank(pval, ties.method = "random"),
    pval_null = pval_null / n(),
    test_type = if_else(du_test %in% c("GLM", "GLMM"), "Parametric", "Non-parametric"),
    Gene = group_gene,
    Isoform = group_iso,
    Covariate = str_to_title(covar_type)
  )

# visualize data
m3 <- ggplot(
  df_du_res,
  aes(x = pval_null, y = pval + 1e-80, group = du_test, color = du_test)
) +
  facet_wrap(Covariate ~ Gene + Isoform, nrow = 2, scales = "free", labeller = labeller(.cols = label_both)) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  geom_line(linewidth = 1, aes(linetype = test_type)) +
  labs(x = "Expected p-value under the null", y = "Observed p-value", 
       color = "Differential usage test", linetype = '') +
  scale_x_continuous(transform = c("log10", "reverse")) +
  scale_y_continuous(transform = c("log10", "reverse")) +
  scale_color_manual(values = c(
    'GLMM' = '#00788C', 'GLM' = '#6ED3CF',
    'Conditional HSIC' = '#eb5c2f', 'Unconditional HSIC' = '#f9ca79',
    'T-test (Fisher\'s)' = '#4B0082')
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
    legend.position = "right"
  )
m3

ggsave(file.path(res_dir, "figures/du_general_six.png"), m3, width = 16, height = 6, dpi = 300)
ggsave(file.path(res_dir, "figures/du_general_six.pdf"), m3, width = 16, height = 6, dpi = 300)


### Fig S1: Power analysis of SV and DU tests
## Fig S1A: SV simulation design
data_sv_power <- read_csv(file.path(res_dir, "figures/source_data/sim_data_sv_power.csv")) %>%
  mutate(
    PctSpVar = paste0(prop_var_sp * 100, "%"),
  )

# expected isoform ratio
s1.1 <- data_sv_power %>%
  select(c(x, y, PctSpVar, starts_with("exp_ratio"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "exp_ratio_"), starts_with("exp_ratio")) %>%
  # capitalize the first letter of the isoform name
  pivot_longer(cols = starts_with("iso"), names_to = "isoform", values_to = "value") %>%
  mutate(
    isoform = str_to_title(isoform),
  ) %>%
  # visualization
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_grid(isoform ~ PctSpVar, switch = "y", labeller = labeller(.cols = label_both)) +
  labs(fill = "Expected ratio") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#155289") +
  # scale_fill_distiller(type = 'seq', palette = 'Blues', direction = 1) +
  theme_void() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    strip.text = element_text(size = 14, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 14, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

# observed isoform ratio
s1.2 <- data_sv_power %>%
  select(c(x, y, PctSpVar, starts_with("obs_ratio"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "obs_ratio_"), starts_with("obs_ratio")) %>%
  # capitalize the first letter of the isoform name
  pivot_longer(cols = starts_with("iso"), names_to = "isoform", values_to = "value") %>%
  mutate(
    isoform = str_to_title(isoform),
  ) %>%
  # visualization
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_grid(isoform ~ PctSpVar, switch = "y", labeller = labeller(.cols = label_both)) +
  labs(fill = "Observed ratio") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#155289") +
  # scale_fill_distiller(type = 'seq', palette = 'Blues', direction = 1) +
  theme_void() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    strip.text = element_text(size = 14, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 14, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

s1 <- cowplot::plot_grid(s1.1, s1.2, nrow = 1, align = "h", axis = "tbrl", rel_widths = c(1, 1))
s1

ggsave(file.path(res_dir, "figures/sim_sv_power.png"), s1, width = 20, height = 4.5, dpi = 300)

## ig S1B: Different HSIC-IR hyperparameters
df_sv_power <- read_csv(file.path(res_dir, "sv_results/sv_power_1205.csv")) %>%
  pivot_longer(cols = starts_with("hsic-ir_"), names_to = "config", values_to = "pval") %>%
  mutate(
    transformation = strsplit(config, "_") %>% map_chr(., 2),
    nan_filling = strsplit(config, "_") %>% map_chr(., 3),
  ) %>%
  mutate(
    transformation = factor(
      transformation,
      levels = c("none", "alr", "radial"),
      labels = c("None", "ALR", "Radial")
    ),
    nan_filling = factor(
      nan_filling,
      levels = c("none", "mean"),
      labels = c("None", "Mean")
    ),
  ) %>%
  # compute expected pvalue under the null for each test by ranking pval
  group_by(prop_var_sp, config) %>%
  mutate(
    pval_null = rank(pval),
    pval_null = pval_null / n(),
    PctSpVar = paste0(prop_var_sp * 100, "%"),
  )

# visualize the sv results
s2 <- ggplot(
  df_sv_power,
  aes(x = pval_null, y = pval, group = config, color = transformation, linetype = nan_filling)
) +
  facet_wrap(~PctSpVar, nrow = 1, scales = "free", labeller = labeller(.cols = label_both)) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  geom_line(linewidth = 1) +
  labs(
    x = "Expected p-value under the null", y = "Observed p-value",
    color = "Transformation", linetype = "NA filling"
  ) +
  scale_x_continuous(transform = c("log10", "reverse")) +
  scale_y_continuous(transform = c("log10", "reverse")) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
    legend.position = "right"
  )
s2

ggsave(file.path(res_dir, "figures/sv_power.png"), s2, width = 16, height = 3, dpi = 300)

## Fig S1C: DU simulation design
data_du_power <- read_csv(file.path(res_dir, "figures/source_data/sim_data_du_power.csv")) %>%
  mutate(
    EffectSize = beta_scale,
  )

# expected isoform ratio
s3.1 <- data_du_power %>%
  select(c(x, y, EffectSize, starts_with("exp_ratio"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "exp_ratio_"), starts_with("exp_ratio")) %>%
  # capitalize the first letter of the isoform name
  pivot_longer(cols = starts_with("iso"), names_to = "isoform", values_to = "value") %>%
  mutate(
    isoform = str_to_title(isoform),
  ) %>%
  # visualization
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_grid(isoform ~ EffectSize, switch = "y", labeller = labeller(.cols = label_both)) +
  labs(fill = "Expected ratio") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#155289") +
  # scale_fill_distiller(type = 'seq', palette = 'Blues', direction = 1) +
  theme_void() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    strip.text = element_text(size = 14, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 14, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

# observed isoform ratio
s3.2 <- data_du_power %>%
  select(c(x, y, EffectSize, starts_with("obs_ratio"))) %>%
  # remove the prefix from the column names
  rename_with(~ str_remove(., "obs_ratio_"), starts_with("obs_ratio")) %>%
  # capitalize the first letter of the isoform name
  pivot_longer(cols = starts_with("iso"), names_to = "isoform", values_to = "value") %>%
  mutate(
    isoform = str_to_title(isoform),
  ) %>%
  # visualization
  ggplot(aes(x = x, y = y, fill = value)) +
  facet_grid(isoform ~ EffectSize, switch = "y", labeller = labeller(.cols = label_both)) +
  labs(fill = "Observed ratio") +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "#155289") +
  # scale_fill_distiller(type = 'seq', palette = 'Blues', direction = 1) +
  theme_void() +
  theme(
    text = element_text(size = 14, family = "Arial"),
    strip.text = element_text(size = 14, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 14, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
  )

s3 <- cowplot::plot_grid(s3.1, s3.2, nrow = 1, align = "h", axis = "tbrl", rel_widths = c(1, 1))
s3

ggsave(file.path(res_dir, "figures/sim_du_power.png"), s3, width = 20, height = 4.5, dpi = 300)


## Fig S1D: Different non-parametric DU test configs
df_du_np_power <- read_csv(file.path(res_dir, "du_results/du_np_power_1205.csv")) %>%
  pivot_longer(cols = starts_with("pvalue_"), names_to = "du_test", values_to = "pval") %>%
  mutate(
    du_test = str_remove(du_test, "pvalue_"),
    du_test = factor(
      du_test,
      levels = c("t-fisher", "t-tippett", "hsic", "hsic-knn", "hsic-eps", "hsic-gp"),
      labels = c(
        "T-test (Fisher's)", "T-test (Tippett’s)", 
        "Uncond HSIC", "Cond HSIC (KNN residual)", 
        "Cond HSIC (fixed kernel)", "Cond HSIC (GP kernel)"
      )),
  ) %>%
  # compute expected pvalue under the null for each test by ranking pval
  group_by(beta_scale, covar_type, du_test) %>%
  mutate(
    pval_null = rank(pval, ties.method = "random"),
    pval_null = pval_null / n(),
    test_type = if_else(du_test %in% c("T-test (Fisher's)", "T-test (Tippett’s)"), "T-test", "HSIC-IR"),
    EffectSize = beta_scale,
    Covariate = str_to_title(covar_type)
  )

# visualize data
s4 <- ggplot(
  df_du_np_power,
  aes(x = pval_null, y = pval, group = du_test, color = du_test)
) +
  facet_wrap(Covariate ~ EffectSize, nrow = 2, scales = "free", labeller = labeller(.cols = label_both)) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  geom_line(linewidth = 1, aes(linetype = test_type)) +
  labs(x = "Expected p-value under the null", y = "Observed p-value", 
       color = "Differential usage test", linetype = "") +
  scale_x_continuous(transform = c("log10", "reverse")) +
  scale_y_continuous(transform = c("log10", "reverse")) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
    legend.position = "right"
  )
s4

ggsave(file.path(res_dir, "figures/du_np_power.png"), s4, width = 16, height = 5, dpi = 300)


## Fig S1E: Different parametric test configs
df_du_p_power <- read_csv(file.path(res_dir, "du_results/du_p_power_1205.csv")) %>%
  pivot_longer(cols = starts_with("pvalue_"), names_to = "du_test", values_to = "pval") %>%
  mutate(
    du_test = str_remove(du_test, "pvalue_"),
    model_type = strsplit(du_test, "-") %>% map_chr(., 1),
    test_type = strsplit(du_test, "-") %>% map_chr(., 2),
  ) %>%
  mutate(
    model_type = factor(model_type, levels = c("glm", "glmm"), labels = c("GLM", "GLMM")),
    test_type = factor(test_type, levels = c("score", "wald"), labels = c("Score", "Wald")),
  ) %>%
  # compute expected pvalue under the null for each test by ranking pval
  group_by(beta_scale, covar_type, du_test) %>%
  mutate(
    pval_null = rank(pval, ties.method = "random"),
    pval_null = pval_null / n(),
    EffectSize = beta_scale,
    Covariate = str_to_title(covar_type)
  )

# visualize data
s5 <- ggplot(
  df_du_p_power,
  aes(x = pval_null, y = pval + 1e-20, group = du_test)
) +
  facet_wrap(Covariate ~ EffectSize, nrow = 2, scales = "free", labeller = labeller(.cols = label_both)) +
  geom_abline(intercept = 0, slope = 1, color = "gray") +
  geom_line(linewidth = 1, aes(linetype = test_type, color = model_type)) +
  labs(x = "Expected p-value under the null", y = "Observed p-value", 
       color = "Parametric model", linetype = "Test statistic") +
  scale_x_continuous(transform = c("log10", "reverse")) +
  scale_y_continuous(transform = c("log10", "reverse")) +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = "Arial"),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = "Arial", margin = margin(2, 2, 2, 2, "pt")),
    plot.title = element_text(size = 12, family = "Arial", hjust = 0.5, margin = margin(0, 0, 5, 0, "pt")),
    legend.position = "right"
  )
s5

ggsave(file.path(res_dir, "figures/du_p_power.png"), s5, width = 16, height = 5, dpi = 300)

