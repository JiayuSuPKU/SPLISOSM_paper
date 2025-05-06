library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(rtracklayer)
library(ComplexHeatmap)

extrafont::loadfonts()

res_sc_dir <- "~/Projects/SPLISOSM_paper/results/sc_tasic_nature_18"
res_ont_dir <- "~/Projects/SPLISOSM_paper/results/sit_nar_23"
res_trend_dir <- "~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
ont_date <- "1107"
trend_date <- "1119"

rbp_with_motifs <- read_table("~/reference/cisbp-rna/mouse_pm/rbp_with_known_cleaned_motif.txt", col_names = 'rbp')
rbp_with_clip <- read_table("~/reference/POSTAR3/mouse.rbp_with_clip.txt", col_names = 'rbp')

### DU p-value correlation with RBFOX3
## Load DU results
# ONT
du_ont_rbp <- read_csv(sprintf("%s/figures/cbs/source_data/du_res_annot.csv", res_ont_dir)) %>%
  mutate(label = paste(gene, covariate, sep = ":")) %>%
  mutate(
    has_clip = covariate %in% rbp_with_clip$rbp,
    has_motif = covariate %in% rbp_with_motifs$rbp,
  ) %>%
  rowwise() %>%
  mutate(
    is_significant = max(
      `pvalue_hsic-gp_1`, `pvalue_hsic-gp_2`, 
      pvalue_glmm_1, pvalue_glmm_2
    ) < 0.01
  ) %>%
  ungroup()

# SR TREND
du_trend_rbp <- read_csv(sprintf("%s/figures/source_data/du_res_annot.csv", res_trend_dir)) %>%
  mutate(label = paste(gene, covariate, sep = ":")) %>%
  mutate(
    has_clip = covariate %in% rbp_with_clip$rbp,
    has_motif = covariate %in% rbp_with_motifs$rbp,
  ) %>%
  rowwise() %>%
  mutate(
    is_significant = max(
      `pvalue_hsic-gp`, pvalue_glmm
    ) < 0.01
  ) %>%
  ungroup()

## Calculate pvalue correlation with Rbfox3
# ONT
rbfox3_du_ont_corr <- data.frame(
  rbp = du_ont_rbp %>% pull(covariate) %>% unique() %>% sort(),
  spearman_rho = 0,
  spearman_pval = NA,
  n_both_sig = 0
) %>%
  mutate(
    has_clip = rbp %in% rbp_with_clip$rbp,
    has_motif = rbp %in% rbp_with_motifs$rbp
  )

for (rbp in rbfox3_du_ont_corr$rbp) {
  x1 <- du_ont_rbp %>% filter(covariate == 'Rbfox3') %>% pull(`pvalue_hsic-gp_2`)
  x2 <- du_ont_rbp %>% filter(covariate == rbp) %>% pull(`pvalue_hsic-gp_2`)
  res <- cor.test(x1, x2, method = 'spearman')
  rbfox3_du_ont_corr[
    rbfox3_du_ont_corr$rbp == rbp, 
    c('spearman_rho', 'spearman_pval', 'n_both_sig')
  ] <- c(res$estimate, res$p.value, sum(x1 < 0.01 & x2 < 0.01))
}

# SR TREND
rbfox3_du_trend_corr <- data.frame(
  rbp = du_trend_rbp %>% pull(covariate) %>% unique() %>% sort(),
  spearman_rho = 0,
  spearman_pval = NA,
  n_both_sig = 0
) %>%
  mutate(
    has_clip = rbp %in% rbp_with_clip$rbp,
    has_motif = rbp %in% rbp_with_motifs$rbp
  )

for (rbp in rbfox3_du_trend_corr$rbp) {
  x1 <- du_trend_rbp %>% filter(covariate == 'Rbfox3') %>% pull(`pvalue_hsic-gp`)
  x2 <- du_trend_rbp %>% filter(covariate == rbp) %>% pull(`pvalue_hsic-gp`)
  res <- cor.test(x1, x2, method = 'spearman')
  rbfox3_du_trend_corr[
    rbfox3_du_trend_corr$rbp == rbp, 
    c('spearman_rho', 'spearman_pval', 'n_both_sig')
  ] <- c(res$estimate, res$p.value, sum(x1 < 0.01 & x2 < 0.01))
}

## visualize top Rbfox3 co-regulators, ranked by number of overlaps
rbp_to_include <- c(
  union(rbp_with_clip$rbp, rbp_with_motifs$rbp),
  c('Rbfox3', 'Rbfox1', 'Rbfox2', 
    'Celf1', 'Celf2', 'Celf3', 'Celf4', 'Celf5', 'Celf6', 
    'Qk', 'Khdrbs3', 'Khdrbs2', 'Khdrbs1',
    'Arpp21', 'Adar', 'Adarb1', 'Adarb2')
)

# ONT
rbfox3_topk_ont <- rbfox3_du_ont_corr %>% 
  filter(rbp %in% rbp_to_include) %>%
  arrange(spearman_pval) %>% head(n=15) %>%
  mutate(rbp = factor(rbp, levels = rbp))

# barplot, rbp ranked by number of overlaps
p1.1 <- ggplot(rbfox3_topk_ont, aes(x = rbp, y = spearman_rho)) + 
  geom_bar(
    aes(fill = rbp %in% c(
      'Rbfox1', 'Rbfox2', 'Rbfox3', 
      'Celf1','Celf2', 'Celf3', 'Celf4', 'Celf5', 
      'Qk', 'Khdrbs3', 'Elavl2'
    )), 
    stat = 'identity') + 
  geom_text(aes(label = sprintf('%.0f', n_both_sig)), vjust = -0.05) + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
  labs(title = 'DU p-value correlation (ONT)', 
       x = '', y = 'Spearman rho') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = 'none'
  )

# SR TREND
rbfox3_topk_trend <- rbfox3_du_trend_corr %>% 
  filter(rbp %in% rbp_to_include) %>%
  arrange(spearman_pval) %>% head(n=15) %>%
  mutate(rbp = factor(rbp, levels = rbp))
p1.2 <- ggplot(rbfox3_topk_trend, aes(x = rbp, y = spearman_rho)) + 
  geom_bar(
    aes(fill = rbp %in% c(
      'Rbfox1', 'Rbfox2', 'Rbfox3', 
      'Celf1','Celf2', 'Celf3', 'Celf4', 'Celf5', 
      'Qk', 'Khdrbs3', 'Elavl2'
    )),
    stat = 'identity') + 
  geom_text(aes(label = sprintf('%.0f', n_both_sig)), vjust = -0.05) + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
  labs(title = 'DU p-value correlation (SR)', 
       x = '', y = 'Spearman rho') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = 'none'
  )

p1 <- cowplot::plot_grid(p1.1, p1.2, nrow = 1, align = 'hv')
p1

ggsave(sprintf("%s/figures/du_pval_corr_rbfox3.png", res_sc_dir), p1, width = 8, height = 4)

### Fig 5C: Pairwise coexpression vs DU pvalue correlation between top RBFOX3 co-regulators
## load expression data
rbp_sp_expr_ont <- read_csv(sprintf("%s/figures/cbs/source_data/rbp_expr.csv", res_ont_dir))
rbp_sp_expr_trend <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_trend_dir))

## ONT
rbp_to_plots <- rbfox3_topk_ont %>% pull(rbp)

# calculate Spearman R of DU pvalues
rbp_dupval_corr_ont <- data.frame(
  rbp1 = rep(rbp_to_plots, each = length(rbp_to_plots)),
  rbp2 = rep(rbp_to_plots, length(rbp_to_plots)),
  spearman_rho = 0,
  spearman_pval = NA
)

for (i in 1:nrow(rbp_dupval_corr_ont)) {
  x1 <- du_ont_rbp %>% filter(
    covariate == rbp_dupval_corr_ont$rbp1[i]
  ) %>% pull(`pvalue_hsic-gp_2`)
  x2 <- du_ont_rbp %>% filter(
    covariate == rbp_dupval_corr_ont$rbp2[i]
  ) %>% pull(`pvalue_hsic-gp_2`)
  res <- cor.test(x1, x2, method = 'spearman')
  rbp_dupval_corr_ont[i, c('spearman_rho', 'spearman_pval')] <- c(res$estimate, res$p.value)
}

# filter out insignificant correlations
rbp_dupval_corr_ont <- rbp_dupval_corr_ont %>%
  # filter out insignificant correlations
  mutate(spearman_rho = ifelse(spearman_pval > 0.01, NA, spearman_rho)) %>%
  mutate(rbp1 = factor(rbp1, levels = rbp_to_plots),
         rbp2 = factor(rbp2, levels = rbp_to_plots))

# ggplot(rbp_dupval_corr_ont, aes(x = rbp1, y = rbp2, fill = spearman_rho)) + 
#   geom_tile() + 
#   scale_fill_distiller(
#     type = 'div', palette = 'RdBu', na.value = 'gray',
#     limits = c(-1, 1), oob = scales::squish
#   ) +
#   labs(x = '', y = '', fill = 'Spearman rho', title = 'RBP co-regulation (DU ONT)') + 
#   theme_classic() + 
#   theme(
#     text = element_text(size = 12, family = 'Arial'),
#     axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
#     axis.text.y = element_text(size = 12, face = 'italic'),
#     legend.position = 'right'
#   )

# calculate Spearman R of expression log1p
rbp_expr_corr_ont <- data.frame(
  rbp1 = rep(rbp_to_plots, each = length(rbp_to_plots)),
  rbp2 = rep(rbp_to_plots, length(rbp_to_plots)),
  spearman_rho = 0,
  spearman_pval = NA
)

for (i in 1:nrow(rbp_expr_corr_ont)) {
  x1 <- rbp_sp_expr_ont %>% filter(
    gene == rbp_expr_corr_ont$rbp1[i],
    layer == 'log1p'
  ) %>% arrange(barcode) %>% pull(expression)
  x2 <- rbp_sp_expr_ont %>% filter(
    gene == rbp_expr_corr_ont$rbp2[i],
    layer == 'log1p'
  ) %>% arrange(barcode) %>% pull(expression)
  res <- cor.test(x1, x2, method = 'spearman')
  rbp_expr_corr_ont[i, c('spearman_rho', 'spearman_pval')] <- c(res$estimate, res$p.value)
}

# filter out insignificant correlations
rbp_expr_corr_ont <- rbp_expr_corr_ont %>%
  mutate(spearman_rho = ifelse(spearman_pval > 0.01, NA, spearman_rho)) %>%
  mutate(rbp1 = factor(rbp1, levels = rbp_to_plots),
         rbp2 = factor(rbp2, levels = rbp_to_plots))

# ggplot(rbp_expr_corr_ont, aes(x = rbp1, y = rbp2, fill = spearman_rho)) + 
#   geom_tile() + 
#   scale_fill_distiller(
#     type = 'div', palette = 'RdBu', na.value = 'gray',
#     limits = c(-1, 1), oob = scales::squish
#   ) +
#   labs(x = '', y = '', fill = 'Spearman rho', title = 'RBP co-expression (ONT-paired Visium)') + 
#   theme_classic() + 
#   theme(
#     text = element_text(size = 12, family = 'Arial'),
#     axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
#     axis.text.y = element_text(size = 12, face = 'italic'),
#     legend.position = 'right'
#   )

# combine the two heatmaps
rbp_pairwise_ont <- merge(
  rbp_dupval_corr_ont, rbp_expr_corr_ont, 
  by = c('rbp1', 'rbp2'), suffixes = c('.du', '.expr')
) %>%
  mutate(
    x = as.numeric(rbp1),
    y = as.numeric(rbp2),
    rho_display = ifelse(x < y, spearman_rho.du, spearman_rho.expr)
  )

p2.1 <- ggplot(rbp_pairwise_ont, aes(x = rbp1, y = rbp2, fill = rho_display)) + 
  geom_tile() + 
  scale_fill_distiller(
    type = 'div', palette = 'RdBu', na.value = 'gray',
    limits = c(-0.5, 0.5), oob = scales::squish
  ) +
  scale_y_discrete(limits = rev(levels(rbp_pairwise_ont$rbp2))) +
  labs(x = '', y = '', fill = 'Spearman rho', title = 'Rbfox3 potential co-regulators (ONT, AS)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
    axis.text.y = element_blank(),
    legend.position = 'bottom',
    plot.margin = margin(-1, 0, 0, 0),
  ) + 
  guides(fill = guide_colorbar(
    title.position = 'top', 
    title.hjust = 0.5, 
    title.vjust = 1, 
    barwidth = unit(2, 'inch')
  ))

p2.2 <- ggplot(rbfox3_topk_ont, aes(x = rbp, y = n_both_sig)) + 
  geom_bar(
    aes(fill = rbp %in% c(
      'Rbfox1', 'Rbfox2', 'Rbfox3', 
      'Celf1','Celf2', 'Celf3', 'Celf4', 'Celf5', 
      'Qk', 'Khdrbs3', 'Elavl2'
    )), 
    stat = 'identity') + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
  labs(y = '', x = '', title = '#shared targets') + 
  scale_x_discrete(limits = rev(levels(rbfox3_topk_ont$rbp)),
                   position = 'top') + 
  scale_y_continuous(
    transform = 'reverse', position = 'right'
  ) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = 'italic'),
    legend.position = 'none'
  )

p2 <- cowplot::plot_grid(
  p2.2, NULL, p2.1, nrow = 1, align = 'hv', 
  rel_widths = c(0.5, -0.1, 1)
)
p2
ggsave(sprintf("%s/figures/rbfox3_coreg_rbp_ont.pdf", res_sc_dir), p2, width = 5.8, height = 5.5)


## SR TREND
rbp_to_plots <- rbfox3_topk_trend %>% pull(rbp)

# calculate Spearman R of DU pvalues
rbp_dupval_corr_trend <- data.frame(
  rbp1 = rep(rbp_to_plots, each = length(rbp_to_plots)),
  rbp2 = rep(rbp_to_plots, length(rbp_to_plots)),
  spearman_rho = 0,
  spearman_pval = NA
)

for (i in 1:nrow(rbp_dupval_corr_trend)) {
  x1 <- du_trend_rbp %>% filter(
    covariate == rbp_dupval_corr_trend$rbp1[i]
  ) %>% pull(`pvalue_hsic-gp`)
  x2 <- du_trend_rbp %>% filter(
    covariate == rbp_dupval_corr_trend$rbp2[i]
  ) %>% pull(`pvalue_hsic-gp`)
  res <- cor.test(x1, x2, method = 'spearman')
  rbp_dupval_corr_trend[i, c('spearman_rho', 'spearman_pval')] <- c(res$estimate, res$p.value)
}

# filter out insignificant correlations
rbp_dupval_corr_trend <- rbp_dupval_corr_trend %>%
  # filter out insignificant correlations
  mutate(spearman_rho = ifelse(spearman_pval > 0.01, NA, spearman_rho)) %>%
  mutate(rbp1 = factor(rbp1, levels = rbp_to_plots),
         rbp2 = factor(rbp2, levels = rbp_to_plots))

# ggplot(rbp_dupval_corr_trend, aes(x = rbp1, y = rbp2, fill = spearman_rho)) +
#   geom_tile() +
#   scale_fill_distiller(
#     type = 'div', palette = 'RdBu', na.value = 'gray',
#     limits = c(-1, 1), oob = scales::squish
#   ) +
#   labs(x = '', y = '', fill = 'Spearman rho', title = 'RBP co-regulation (DU SR)') +
#   theme_classic() +
#   theme(
#     text = element_text(size = 12, family = 'Arial'),
#     axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
#     axis.text.y = element_text(size = 12, face = 'italic'),
#     legend.position = 'right'
#   )

# calculate Spearman R of expression log1p
rbp_expr_corr_trend <- data.frame(
  rbp1 = rep(rbp_to_plots, each = length(rbp_to_plots)),
  rbp2 = rep(rbp_to_plots, length(rbp_to_plots)),
  spearman_rho = 0,
  spearman_pval = NA
)

for (i in 1:nrow(rbp_expr_corr_trend)) {
  x1 <- rbp_sp_expr_trend %>% filter(
    gene == rbp_expr_corr_trend$rbp1[i],
    layer == 'log1p'
  ) %>% arrange(barcode) %>% pull(expression)
  x2 <- rbp_sp_expr_trend %>% filter(
    gene == rbp_expr_corr_trend$rbp2[i],
    layer == 'log1p'
  ) %>% arrange(barcode) %>% pull(expression)
  res <- cor.test(x1, x2, method = 'spearman')
  rbp_expr_corr_trend[i, c('spearman_rho', 'spearman_pval')] <- c(res$estimate, res$p.value)
}

# filter out insignificant correlations
rbp_expr_corr_trend <- rbp_expr_corr_trend %>%
  mutate(spearman_rho = ifelse(spearman_pval > 0.01, NA, spearman_rho)) %>%
  mutate(rbp1 = factor(rbp1, levels = rbp_to_plots),
         rbp2 = factor(rbp2, levels = rbp_to_plots))

# ggplot(rbp_expr_corr_trend, aes(x = rbp1, y = rbp2, fill = spearman_rho)) +
#   geom_tile() +
#   scale_fill_distiller(
#     type = 'div', palette = 'RdBu', na.value = 'gray',
#     limits = c(-1, 1), oob = scales::squish
#   ) +
#   labs(x = '', y = '', fill = 'Spearman rho', title = 'RBP co-expression (Visium)') +
#   theme_classic() +
#   theme(
#     text = element_text(size = 12, family = 'Arial'),
#     axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
#     axis.text.y = element_text(size = 12, face = 'italic'),
#     legend.position = 'right'
#   )

# combine the two heatmaps
rbp_pairwise_trend <- merge(
  rbp_dupval_corr_trend, rbp_expr_corr_trend, 
  by = c('rbp1', 'rbp2'), suffixes = c('.du', '.expr')
) %>%
  mutate(
    x = as.numeric(rbp1),
    y = as.numeric(rbp2),
    rho_display = ifelse(x < y, spearman_rho.du, spearman_rho.expr)
  )

p3.1 <- ggplot(rbp_pairwise_trend, aes(x = rbp1, y = rbp2, fill = rho_display)) + 
  geom_tile() + 
  scale_fill_distiller(
    type = 'div', palette = 'RdBu', na.value = 'gray',
    limits = c(-0.5, 0.5), oob = scales::squish
  ) +
  scale_y_discrete(limits = rev(levels(rbp_pairwise_trend$rbp2))) +
  labs(x = '', y = '', fill = 'Spearman rho', title = 'Rbfox3 potential co-regulators (SR, APA)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
    axis.text.y = element_blank(),
    legend.position = 'bottom',
    plot.margin = margin(-1, 0, 0, 0),
  ) + 
  guides(fill = guide_colorbar(
    title.position = 'top', 
    title.hjust = 0.5, 
    title.vjust = 1, 
    barwidth = unit(2, 'inch')
  ))

p3.2 <- ggplot(rbfox3_topk_trend, aes(x = rbp, y = n_both_sig)) + 
  geom_bar(
    aes(fill = rbp %in% c(
      'Rbfox1', 'Rbfox2', 'Rbfox3', 
      'Celf1','Celf2', 'Celf3', 'Celf4', 'Celf5', 
      'Qk', 'Khdrbs3', 'Elavl2'
    )), 
    stat = 'identity') + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
  labs(y = '', x = '', title = '#shared targets') + 
  scale_x_discrete(limits = rev(levels(rbfox3_topk_trend$rbp)),
                   position = 'top') + 
  scale_y_continuous(
    transform = 'reverse', position = 'right'
  ) + 
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12, face = 'italic'),
    legend.position = 'none'
  )

p3 <- cowplot::plot_grid(
  p3.2, NULL, p3.1, nrow = 1, align = 'hv', 
  rel_widths = c(0.5, -0.1, 1)
)
p3
ggsave(sprintf("%s/figures/rbfox3_coreg_rbp_sr.pdf", res_sc_dir), p3, width = 5.8, height = 5.5)


### Fig S4B: Heatmap of DU p-values
# rbp_to_plots <- c('Qk', 'Khdrbs3', paste0('Rbfox', 1:3), paste0('Celf', 1:5))
rbp_to_plots <- c('Qk', 'Khdrbs3', 'Elavl2', paste0('Celf', 1:5), paste0('Rbfox', 1:3))

# ONT, rank SVS by number of overlaps and DU p-values
df_ont_subset <- du_ont_rbp %>%
  rowwise() %>%
  mutate(
    # is_significant = max(`pvalue_hsic-gp_2`, pvalue_glmm_2) < 0.01,
    is_significant = max(
      `pvalue_hsic-gp_1`, `pvalue_hsic-gp_2`,
      pvalue_glmm_1, pvalue_glmm_2
    ) < 0.01
  ) %>%
  filter(
    is_significant,
    (covariate %in% rbp_to_plots)
  ) %>%
  # pad missing gene-covariate pairs with NA
  complete(gene, covariate, fill = list(is_significant = FALSE))
# order genes by the number of significant associations and p-values
gene_ont_order <- df_ont_subset %>% summarise(
  size = sum(is_significant), 
  pval = min(`pvalue_hsic-gp_2`, na.rm = TRUE),
  .by = gene
) %>% arrange(-size, pval) %>%
  filter(size >= 2)

# heatmap of DU p-values
x_color <- ifelse(
  gene_ont_order$gene %in% c('Clta', 'Cltb', 'Nnat', 'Myl6'), 'red', 'black'
)
p2 <- ggplot(
  df_ont_subset %>% 
    filter(gene %in% gene_ont_order$gene) %>%
    mutate(
      gene = factor(gene, levels = gene_ont_order$gene),
      covariate = factor(covariate, levels = rbp_to_plots)
    ), 
  aes(x = gene, y = covariate, fill = `pvalue_hsic-gp_2`)
) + 
  geom_tile(na.rm = FALSE) + 
  labs(x = '', y = '', fill = 'Conditional HSIC p-value (ONT)') +
  scale_fill_distiller(
    type = 'seq', palette = 'Reds', direction = -1, na.value = 'white',
    limits = c(0, 0.01), breaks = c(0, 0.005, 0.01)
  ) + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.y = element_text(size = 12, face = 'italic'),
    axis.text.x = element_text(
      size = 12, angle = 90, hjust = 1, vjust = 0.5, face = 'italic', color = x_color
    ),
    legend.position = 'top',
    legend.margin=margin(0,0,0,0),
    # legend.box.margin=margin(-15,0,0,0)
  ) + 
  guides(fill = guide_colorbar(
    title.position = 'top', 
    title.hjust = 0.5, 
    title.vjust = 1, 
    barwidth = unit(2, 'inch')
  ))
p2

ggsave(sprintf("%s/figures/du_pval_heatmap_ont.png", res_sc_dir), p2, width = 4, height = 5)

# SR, rank SVS by number of overlaps and DU p-values
df_trend_subset <- du_trend_rbp %>%
  filter(is_significant & (covariate %in% rbp_to_plots)) %>%
  # pad missing gene-covariate pairs with NA
  complete(gene, covariate, fill = list(is_significant = FALSE))
# order genes by the number of significant associations and p-values
gene_trend_order <- df_trend_subset %>% summarise(
  size = sum(is_significant), 
  pval = min(`pvalue_hsic-gp`, na.rm = TRUE),
  .by = gene
) %>% arrange(-size, pval) %>%
  filter(size >= 3)

## heatmap of DU p-values
x_color <- ifelse(
  gene_trend_order$gene %in% c('Myl6', 'Tln1', 'Klc1', 'Gnao1'), 'red', 'black'
)
p3 <- ggplot(
  df_trend_subset %>% 
    filter(gene %in% gene_trend_order$gene) %>%
    mutate(
      gene = factor(gene, levels = gene_trend_order$gene),
      covariate = factor(covariate, levels = rbp_to_plots)
    ), 
  aes(x = gene, y = covariate, fill = `pvalue_hsic-gp`)
) + 
  geom_tile(na.rm = FALSE) + 
  labs(x = '', y = '', fill = 'Conditional HSIC p-value (SR)') +
  scale_fill_distiller(
    type = 'seq', palette = 'Reds', direction = -1, na.value = 'white',
    limits = c(0, 0.01), breaks = c(0, 0.005, 0.01)
  ) + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.y = element_text(size = 12, face = 'italic'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5, 
                               face = 'italic', color = x_color),
    legend.position = 'top',
    legend.margin=margin(0,0,0,0),
    # legend.box.margin=margin(-15,0,0,0)
  ) + 
  guides(fill = guide_colorbar(
    title.position = 'top', 
    title.hjust = 0.5, 
    title.vjust = 1, 
    barwidth = unit(2, 'inch')
  ))
p3

ggsave(sprintf("%s/figures/du_pval_heatmap_trend.png", res_sc_dir), p3, width = 10, height = 5)


### Spatial co-expression patterns
rbp_sp_expr_ont <- read_csv(sprintf("%s/figures/cbs/source_data/rbp_expr.csv", res_ont_dir))
rbp_sp_expr_trend <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_trend_dir))

## ONT, spatial expression plot of selected example RBPs
p_list <- c()
for (rbp in c('Rbfox3', 'Celf5', 'Qk', 'Elavl2')){
  p <- rbp_sp_expr_ont %>% filter(
    gene == rbp, layer == 'log1p', array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = expression)) + 
    geom_point(size = 0.5) + 
    labs(color = '', title = rbp) + 
    scale_color_gradient2(low = 'white', high = ifelse(
      rbp %in% c('Qk'), '#0066CC', '#8B0000')
    ) +
    # scale_color_viridis_c(option = 'A') +
    theme_void() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      strip.text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      panel.background = element_rect(color = 'black', linewidth = 1),
      legend.position = 'right',
      legend.text = element_text(size = 12, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p4 <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
p4

ggsave(sprintf("%s/figures/sp_rbp_expr_ont.png", res_sc_dir), p4, width = 4.5, height = 4)

### Co-expression in the single-cell exon data from Tasic et al.
## load normalized expression of Rbfox, Celf, Qk and Khdrbs
df_rbp_sc_expr <- read_csv(sprintf("%s/rbp_expr.csv", res_sc_dir))

rbp_to_plots <- c('Qk', 'Khdrbs3', 'Elavl2', paste0('Rbfox', 1:3), paste0('Celf', 1:6))
df <- df_rbp_sc_expr %>% pivot_longer(
  cols = rbp_to_plots, names_to = 'rbp', values_to = 'expr'
)
ha <- rowAnnotation(
  Subtype = df_rbp_sc_expr$Subtype, 
  GE = df_rbp_sc_expr$GE,
  annotation_legend_param = list(Subtype = list(ncol = 2))
  # annotation_legend_param = list(
  #   Subtype = list(direction = "horizontal", nrow = 3),
  #   GE = list(direction = "horizontal", nrow = 1)
  # )
)
ht <- Heatmap(
  df_rbp_sc_expr %>% dplyr::select(rbp_to_plots) %>% as.matrix(),
  name = 'Log expression',
  column_title = 'RBP co-expression (Tasic et al.)',
  col = circlize::colorRamp2(c(-5, 0, 10), c('#0066CC', 'white', '#8B0000')),
  right_annotation = ha,
  cluster_rows = FALSE, cluster_columns = FALSE
)

png(sprintf("%s/figures/sup_sc_rbp_corr_heatmap.png", res_sc_dir), width = 6, height = 6, units = 'in', res = 300)
draw(ht, merge_legend = FALSE, heatmap_legend_side = "right", 
     annotation_legend_side = "right")
dev.off()


### Fig 5D Motif enrichment analysis
## Motif enrichment of ONT Rbfox SVS targets
dir.create(sprintf("%s/events/ont", res_sc_dir), recursive = TRUE, showWarnings = FALSE)

# select Rbfox3 targets to run motif enrichment analysis
ont_svs_gtf <- import(sprintf("%s/transcripts/cbs/cbs_svs.txid.gtf", res_ont_dir))
ont_rbfox_targets <- du_ont_rbp %>% filter(
  covariate %in% c('Rbfox3', 'Rbfox1', 'Rbfox2'),
  (`pvalue_hsic-gp_2` < 0.01) & (`pvalue_glmm_2` < 0.01)
) %>% pull(gene) %>% unique() %>% sort()
ont_rbfox_targets_gtf <- ont_svs_gtf[mcols(ont_svs_gtf)$gene_name %in% ont_rbfox_targets]
export(
  ont_rbfox_targets_gtf,
  sprintf("%s/events/ont/ont_rbfox_targets.gtf", res_sc_dir)
)

# load SEA and XSTREME results, exon + 150bp introns on each side
sea_res <- read_tsv(sprintf("%s/events/ont/sea.out/slop/sea.tsv", res_sc_dir), comment = '#')
ref_motifs <- read_meme("~/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme")

# visualize the top reference motifs
p5.1 <- view_motifs(ref_motifs %>% filter_motifs(name = sea_res$ID[1])) + 
  labs(title = sprintf("%s (Elavl family)\nSEA p-val=%s (%s/39 exons)", 
                       sea_res$ID[1], sea_res$PVALUE[1], sea_res$TP[1])) +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p5.2 <- view_motifs(ref_motifs %>% filter_motifs(name = sea_res$ID[2])) + 
  labs(title = sprintf("%s (Qk)\nSEA p-val=%s (%s/39 exons)", 
                       sea_res$ID[2], sea_res$PVALUE[2], sea_res$TP[2])) +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p5.3 <- view_motifs(ref_motifs %>% filter_motifs(name = sea_res$ID[4])) + 
  labs(title = sprintf("%s (Celf family)\nSEA p-val=%s (%s/39 exons)", 
                       sea_res$ID[4], sea_res$PVALUE[4], sea_res$TP[4])) +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

# add title
m_title <- ggplot() + 
  labs(title = 'Rbfox associated exons (±150bp, ONT)') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p5 <- cowplot::plot_grid(m_title, p5.1, p5.2, p5.3, nrow = 4, rel_heights = c(0.1, 0.5, 0.5, 0.5))
p5

ggsave(sprintf("%s/figures/motif_enrich_ont.png", res_sc_dir), p5, width = 3, height = 4)
ggsave(sprintf("%s/figures/motif_enrich_ont.pdf", res_sc_dir), p5, width = 2, height = 4)


## Motif enrichment of TREND Rbfox SVS targets
dir.create(sprintf("%s/events/trend", res_sc_dir), recursive = TRUE, showWarnings = FALSE)

# select Rbfox1/2/3 targets to run motif enrichment analysis
# the bed is in BED12 format
trend_svs_bed <- read.table(
  sprintf("%s/events/cbs_svs.exon.bed", res_trend_dir), 
  sep = '\t', header = FALSE,
  col.names = c(
    "chrom", "chromStart", "chromEnd", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "blockCount", "blockSizes", "blockStarts"
  )
)
trend_rbfox_targets <- du_trend_rbp %>% filter(
  covariate %in% c('Rbfox3', 'Rbfox1', 'Rbfox2'),
  (`pvalue_hsic-gp` < 0.01) & (`pvalue_glmm` < 0.01)
) %>% pull(gene) %>% unique() %>% sort()
trend_rbfox_targets_bed <- trend_svs_bed[
  trend_svs_bed$name %>% str_split_i(., ":", 1) %in% trend_rbfox_targets,
]
# # no duplicates
# trend_rbfox_targets_bed %>% 
#   duplicated(by = c('chrom', 'chromStart', 'chromEnd', 'strand')) %>% sum()
write.table(
  trend_rbfox_targets_bed,
  sprintf("%s/events/trend/trend_rbfox_targets.bed", res_sc_dir),
  sep = '\t', quote = FALSE, col.names = FALSE, row.names = FALSE
)

# load SEA and XSTREME results
sea_res <- read_tsv(sprintf("%s/events/trend/sea.out/exon/sea.tsv", res_sc_dir), comment = '#')
ref_motifs <- read_meme("~/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme")
xstreme_motifs <- read_meme(sprintf("%s/events/trend/xstreme.out/exon/combined.meme", res_sc_dir))

# visualize the top reference motifs
p6.1 <- view_motifs(ref_motifs %>% filter_motifs(name = sea_res$ID[1])) + 
  labs(title = sprintf("%s (Khdrbs3)\nSEA p-val=%s (%s/182 regions)", 
                       sea_res$ID[1], sea_res$PVALUE[1], sea_res$TP[1])) +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p6.2 <- view_motifs(ref_motifs %>% filter_motifs(name = sea_res$ID[3])) + 
  labs(title = sprintf("%s (Elavl family)\nSEA p-val=%s (%s/182 regions)", 
                       sea_res$ID[3], sea_res$PVALUE[3], sea_res$TP[3])) +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

# visualize the top de novo motif
p6.3 <- view_motifs(xstreme_motifs[[1]]) + 
  labs(title = "De novo motif 1\nSEA p-val=2.31e-7 (29/182 regions)") +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

m_title <- ggplot() + 
  labs(title = 'Rbfox associated TREND exons (SR)') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p6 <- cowplot::plot_grid(m_title, p6.1, p6.2, p6.3, nrow = 4, rel_heights = c(0.1, 0.5, 0.5, 0.5))
p6

ggsave(sprintf("%s/figures/motif_enrich_trend.png", res_sc_dir), p6, width = 4, height = 4)
ggsave(sprintf("%s/figures/motif_enrich_trend.pdf", res_sc_dir), p6, width = 2.2, height = 4)


## Motif enrichment of Rbfox targets in tKO data
# load and merge MATT results
motif_enrich_all <- c()
for (effect in c('strong', 'weak')){
  # load enrichment results
  res_exon <- read_tsv(
    sprintf("%s/matt_rbp_map/motif_enrichment.%s.exon.tab", res_sc_dir, effect),
    show_col_types = FALSE
  )
  res_up <- read_tsv(
    sprintf("%s/matt_rbp_map/motif_enrichment.%s.up150.tab", res_sc_dir, effect),
    show_col_types = FALSE
  )
  res_down <- read_tsv(
    sprintf("%s/matt_rbp_map/motif_enrichment.%s.down150.tab", res_sc_dir, effect),
    show_col_types = FALSE
  )
  
  # loop over RBPs
  for (rbp in c('Rbfox', 'Celf', 'Qk', 'Khdrbs', 'Elavl')){
    # loop over silenced and enhanced exons
    for (group in c('silenced', 'enhanced')){
      ## Exon pval and ratio
      df_rbp <- res_exon %>% filter(
        DATASETNAME_FG == group,
        grepl(tolower(rbp), GENES_IN_Mus_musculus)
      )
      min_all_index <- which.min(df_rbp$PVALUE_ALL)
      pval_exon <- df_rbp[min_all_index, c('PVALUE_ALL')]
      ratio_exon <- df_rbp[min_all_index, c('ENRICHMENT_ALL_FG_VS_BG')]
      
      ## Upstream
      df_rbp <- res_up %>% filter(
        DATASETNAME_FG == group,
        grepl(tolower(rbp), GENES_IN_Mus_musculus)
      )
      min_first_index <- which.min(df_rbp$PVALUE_FIRST)
      min_internal_index <- which.min(df_rbp$PVALUE_INTERNAL)
      min_end_index <- which.min(df_rbp$PVALUE_END)
      # keep the most significant p-value
      pval_up_first <- df_rbp[min_first_index, c('PVALUE_FIRST')]
      pval_up_internal <- df_rbp[min_internal_index, c('PVALUE_INTERNAL')]
      pval_up_end <- df_rbp[min_end_index, c('PVALUE_END')]
      # and the corresponding enrichment ratio
      ratio_up_first <- df_rbp[min_first_index, c('ENRICHMENT_FIRST_FG_VS_BG')]
      ratio_up_internal <- df_rbp[min_internal_index, c('ENRICHMENT_INTERNAL_FG_VS_BG')]
      ratio_up_end <- df_rbp[min_end_index, c('ENRICHMENT_END_FG_VS_BG')]
      
      ## Downstream
      df_rbp <- res_down %>% filter(
        DATASETNAME_FG == group,
        grepl(tolower(rbp), GENES_IN_Mus_musculus)
      )
      min_first_index <- which.min(df_rbp$PVALUE_FIRST)
      min_internal_index <- which.min(df_rbp$PVALUE_INTERNAL)
      min_end_index <- which.min(df_rbp$PVALUE_END)
      # keep the most significant p-value
      pval_down_first <- df_rbp[min_first_index, c('PVALUE_FIRST')]
      pval_down_internal <- df_rbp[min_internal_index, c('PVALUE_INTERNAL')]
      pval_down_end <- df_rbp[min_end_index, c('PVALUE_END')]
      # add the corresponding enrichment ratio
      ratio_down_first <- df_rbp[min_first_index, c('ENRICHMENT_FIRST_FG_VS_BG')]
      ratio_down_internal <- df_rbp[min_internal_index, c('ENRICHMENT_INTERNAL_FG_VS_BG')]
      ratio_down_end <- df_rbp[min_end_index, c('ENRICHMENT_END_FG_VS_BG')]
      
      # append to the dataframe
      motif_enrich_all <- rbind(
        motif_enrich_all,
        setNames(
          c(rbp, group, effect, 
            pval_exon, pval_up_first, pval_up_internal, pval_up_end, 
            pval_down_first, pval_down_internal, pval_down_end,
            ratio_exon, ratio_up_first, ratio_up_internal, ratio_up_end,
            ratio_down_first, ratio_down_internal, ratio_down_end
          ),
          c('NAME', 'group', 'effect', 
            'pval_exon', 'pval_up_first', 'pval_up_internal', 'pval_up_end', 
            'pval_down_first', 'pval_down_internal', 'pval_down_end',
            'ratio_exon', 'ratio_up_first', 'ratio_up_internal', 'ratio_up_end',
            'ratio_down_first', 'ratio_down_internal', 'ratio_down_end'
          )
        ) %>% data.frame()
      )
    }
  }
}

# reshape the dataframe
df <- motif_enrich_all %>%
  pivot_longer(
    cols = starts_with('pval_'), names_to = 'position1', values_to = 'pval'
  ) %>%
  pivot_longer(
    cols = starts_with('ratio_'), names_to = 'position2', values_to = 'ratio'
  ) %>%
  mutate(position1 = str_replace(position1, 'pval_', ''),
         position2 = str_replace(position2, 'ratio_', '')) %>%
  filter(position1 == position2) %>%
  mutate(
    ratio = ifelse(pval < 0.05, ratio, NA),
    position = factor(
      position1, levels = c(
        'up_first', 'up_internal', 'up_end', 
        'exon',
        'down_first', 'down_internal', 'down_end'),
      labels = c(
        '[-150, -100]', '[-100, -50]', '[-50, 0]',
        'Exon',
        '[0, 50]', '[50, 100]', '[100, 150]'
      )
    ),
    NAME = factor(NAME, levels = rev(c('Rbfox', 'Celf', 'Elavl', 'Khdrbs', 'Qk'))),
    effect = str_to_sentence(effect)
  )

# pvalue heatmap
ggplot(df, aes(x = position, y = NAME, fill = pval)) + 
  facet_wrap(effect~group, scales = 'free_y') +
  geom_tile() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  geom_vline(xintercept = 4.5, linetype = 'dashed') +
  scale_fill_distiller(
    type = 'seq', palette = 'Reds', direction = -1, na.value = 'white',
    limits = c(0, 0.05), breaks = c(0, 0.01, 0.05)
  ) +
  labs(x = '', y = '', fill = 'Motif enrichment p-value') +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12, face = 'italic'),
    legend.position = 'bottom',
  )

# keep only strong effect and visualize the enrichment ratios
p7 <- df %>%
  filter(effect == 'Strong') %>%
  mutate(
    group = factor(group, levels = c('enhanced', 'silenced'), 
                   labels = c('Enhanced exons (n=217)', 'Silenced exons (n=122)'))
  ) %>%
  ggplot(aes(x = position, y = NAME, fill = log(ratio))) + 
  facet_wrap(~group, scales = 'free_y', ncol = 1) +
  geom_tile() +
  geom_vline(xintercept = 3.5, linetype = 'dashed') +
  geom_vline(xintercept = 4.5, linetype = 'dashed') +
  scale_fill_gradient2(
    low = 'blue', mid = 'white', high = 'red', na.value = 'white',
    midpoint = 0, 
  ) +
  scale_x_discrete(
    labels = c('-125', '-75', '-25', 'Exon', '25', '75', '125')
  ) +
  labs(x = '', y = '', fill = 'LogRatio', title = 'Rbfox regulated exons (±150bp, triple KO)') +
  theme_classic() +
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    strip.text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, family = 'Arial'),
    # axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12, face = 'italic'),
    legend.position = 'bottom',
  )
p7
ggsave(sprintf("%s/figures/motif_ratio_tko.png", res_sc_dir), p7, width = 5, height = 4)
ggsave(sprintf("%s/figures/motif_ratio_tko.pdf", res_sc_dir), p7, width = 3, height = 4)

### RBFOX triple KO data
## Fig S4A: Expression change after Rbfox triple KO
# select RBPs to plot
rbp_to_include <- c(
  union(rbp_with_clip$rbp, rbp_with_motifs$rbp),
  c('Rbfox3', 'Rbfox1', 'Rbfox2', 
    'Celf1', 'Celf2', 'Celf3', 'Celf4', 'Celf5', 'Celf6', 
    'Qk', 'Khdrbs3', 'Khdrbs2', 'Khdrbs1',
    'Arpp21', 'Adar', 'Adarb1', 'Adarb2')
)
rbp_to_plots <- c(
  paste0('Rbfox', 1:3), paste0('Celf', 1:6),
  'Qk', 'Khdrbs3', 'Elavl2'
  # rbfox3_du_ont_corr %>% filter(rbp %in% rbp_to_include) %>% 
  #   arrange(desc(n_both_sig)) %>% head(n=15) %>% pull(rbp),
  # rbfox3_du_trend_corr %>% filter(rbp %in% rbp_to_include) %>% 
  #   arrange(desc(n_both_sig)) %>% head(n=15) %>% pull(rbp)
) %>% unique()

# rbp_to_plots <- c('Nova1', 'Nova2', 'Ptbp2', 'Septin8', 'Map4')

# load reference transcripts annotation
ensembl_gtf <- import.gff('~/reference/mm10/Mus_musculus.GRCm38.102.gtf')
rbp_iso_dict <- ensembl_gtf %>% as.data.frame() %>%
  filter(
    type == 'transcript', gene_name %in% rbp_to_plots
  ) %>%
  dplyr::select(gene_id, gene_name, transcript_id, transcript_name)

# load and prepare gene-level TPM
rbp_gene_tpm <- read_tsv(sprintf("%s/salmon_quant/all_samples.tpm.tsv", res_sc_dir)) %>%
  mutate(transcript_id = str_replace(Name, '\\..*', '')) %>%
  filter(transcript_id %in% rbp_iso_dict$transcript_id) %>%
  left_join(rbp_iso_dict, by = 'transcript_id') %>%
  pivot_longer(
    cols = starts_with('rbfox_'), names_to = 'sample', values_to = 'tpm'
  ) %>%
  mutate(
    group = ifelse(str_detect(sample, 'tko'), 'RBFOX tripleKO', 'WT')
    # gene_name = factor(gene_name, levels = rbp_to_plots)
  ) %>%
  # merge transcript level TPM to gene level
  group_by(gene_name, sample, group) %>%
  summarise(gene_tpm = sum(tpm))

# gene level differential expression
s2 <- ggplot(rbp_gene_tpm, aes(x = gene_name, y = gene_tpm, color = group)) + 
  geom_point(position = position_dodge(width = 0.5)) + 
  stat_compare_means(label = 'p.signif', method = 't.test', label.y.npc = 0.9) + 
  labs(
    x = '', y = 'TPM', color = 'Group', 
    title = 'Expression in day10 motor neurons'
  ) + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial'),
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5, face = 'italic'),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 14),
    legend.position = 'bottom'
  )
s2
ggsave(sprintf("%s/figures/sup_rbp_gene_expr_tko.png", res_sc_dir), s2, width = 4, height = 4)

# load salmon quantification results
rbp_iso_tpm <- read_tsv(sprintf("%s/salmon_quant/all_samples.tpm.tsv", res_sc_dir)) %>%
  mutate(transcript_id = str_replace(Name, '\\..*', '')) %>%
  filter(transcript_id %in% rbp_iso_dict$transcript_id) %>%
  left_join(rbp_iso_dict, by = 'transcript_id') %>%
  pivot_longer(
    cols = starts_with('rbfox_'), names_to = 'sample', values_to = 'tpm'
  ) %>%
  mutate(
    group = ifelse(str_detect(sample, 'tko'), 'tKO', 'WT')
    # gene_name = factor(gene_name, levels = rbp_to_plots)
  ) %>%
  # filter out transcripts that have mean WT tpm < 1
  group_by(gene_name, transcript_name) %>%
  filter(mean(tpm[group == 'WT']) > 5)

# isoform level differential expression
rbp_iso_tpm %>%
  ggplot(aes(x = transcript_name, y = tpm, color = group)) + 
  facet_wrap(~gene_name, scales = 'free', ncol = 4) +
  geom_point(position = position_dodge(width = 0.5)) + 
  stat_compare_means(label = 'p.signif', method = 't.test') + 
  labs(x = '', y = 'TPM', color = 'Group') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = 'bottom'
  )

s3 <- rbp_iso_tpm %>%
  group_by(gene_name, transcript_name, group) %>%
  summarise(
    n = n(),
    mean_tpm = mean(tpm),
    sd_tpm = sd(tpm),
    se_tpm = sd_tpm / sqrt(n)
  ) %>%
  ggplot(aes(x = transcript_name, y = mean_tpm, fill = group)) + 
  facet_wrap(~gene_name, scales = 'free', ncol = 5) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.8) +
  # geom_errorbar(aes(ymin = mean_tpm - se_tpm, ymax = mean_tpm + se_tpm), 
  #               width = 0.2, position = position_dodge(width = 0.8)) +
  geom_point(
    data = rbp_iso_tpm,
    mapping = aes(x = transcript_name, y = tpm),
    position = position_dodge(width = 0.8)
  ) + 
  stat_compare_means(
    data = rbp_iso_tpm, 
    mapping = aes(x = transcript_name, y = tpm),
    label = 'p.signif', method = 't.test',
    label.y.npc = 0.9
  ) +
  labs(x = '', y = 'TPM', fill = 'Day10 motor neuron condition') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial', face = 'italic'),
    strip.background = element_blank(),
    axis.text.x = element_text(size = 12, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 12),
    legend.position = 'bottom',
    legend.position.inside = c(0.9, 0.1)
  )

s3

ggsave(sprintf("%s/figures/sup_rbp_iso_expr_tko.png", res_sc_dir), s3, width = 15, height = 12)


