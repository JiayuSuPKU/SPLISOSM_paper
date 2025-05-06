library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(universalmotif)

extrafont::loadfonts()

res_dir <- "~/Projects/SPLISOSM_paper/results/sit_nar_23"
date <- "1107"

### MOB ONT figures
## make figure directory if doesn't exist
dir.create(sprintf("%s/figures/mob", res_dir), recursive = TRUE, showWarnings = FALSE)

## load SV data
df_sv_pval <- read_csv(sprintf("%s/sv_results/mob_sv_combined_%s.csv", res_dir, date)) %>%
  mutate(
    is_ont_svs = `pvalue_hsic-ir` < 0.05,
    is_ont_sve = `pvalue_hsic-gc` < 0.05,
    is_ont_expressed = count_avg > 0.5
  )

## Fig S2A: Data set stats
df_gene_stats <- read_csv(sprintf("%s/figures/mob/source_data/gene_stats_visium.csv", res_dir))
m1.1 <- ggplot(df_gene_stats, aes(x = ont_pass_qc, y = count_avg_visium)) + 
  geom_boxplot(width = 0.5) + 
  labs(x = 'ONT pass QC', y = 'Visium per-spot UMI') + 
  theme_classic() +
  scale_y_log10() + 
  labs(x = 'Pass ONT QC', y = 'Visium per-spot gene UMI') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    legend.position='none'
  )

df = df_gene_stats %>% filter(ont_pass_qc)
m1.2 <- ggplot(df, aes(x = count_avg, y = count_avg_visium)) + 
  stat_density_2d(color = 'black') + 
  scale_y_log10() + 
  scale_x_log10() + 
  geom_vline(xintercept = df['count_avg'] %>% pull() %>% median(na.rm = T), color = 'red') + 
  geom_hline(yintercept = df['count_avg_visium'] %>% pull() %>% median(na.rm = T), color = 'red') +
  geom_abline(slope = 1, intercept = log10(4), color = 'blue', linetype = 'dashed') + 
  geom_text(
    data = data.frame(x = 1, y = 10),
    mapping = aes(x = x, y = y), label = 'y = 4*x', color = 'blue'
  )+
  labs(y = 'Visium average per-spot UMI', x = 'ONT average per-spot UMI') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    legend.position='none'
  )
m1.3 <- ggplot(df, aes(y = pct_spot_on, x = count_avg)) + 
  geom_point(alpha = 0.5, size = 1) + 
  geom_vline(xintercept = df['count_avg'] %>% pull() %>% median(na.rm = T), color = 'red') + 
  scale_x_log10() + 
  labs(x = 'ONT average per-spot UMI', y = 'Pct spots expressed') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
  )

m1 <- cowplot::plot_grid(
  m1.1, m1.2, m1.3, nrow = 1, rel_widths = c(0.6, 1, 1)
)
m1
ggsave(sprintf("%s/figures/mob/mob_data_stats.png", res_dir), m1, width = 8, height = 3)

## Fig S2B: HSIC-GC vs HSIC-IR
df_ratio_comp_stats = data.frame(
  x = c(1e-54, 1e-54, 1e-62, 1e-62), 
  y = c(1e-3, 1e-5, 1e-3, 1e-5), 
  size = table(
    df_sv_pval$`pvalue_hsic-ir` < 0.05, 
    df_sv_pval$`pvalue_hsic-gc` < 0.05
  ) %>% as.vector() # FF, FT, TF, TT
)

m2 <- ggplot(df_sv_pval, aes(x = `pvalue_hsic-gc`, y = `pvalue_hsic-ir`)) +
  geom_point(color = 'gray', alpha = 0.2) +
  geom_hline(yintercept = 0.05, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.05, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_point(
    data = (df_sv_pval %>% filter(is_ont_svs & is_ont_expressed)),
    aes(color = perplexity)
  ) +
  geom_text_repel(
    data = (df_sv_pval %>% filter(is_ont_expressed) %>% top_n(n = -10, wt = `pvalue_hsic-ir`)),
    aes(label = gene),
    fontface = "italic"
  ) + 
  geom_raster(
    data = df_ratio_comp_stats,
    aes(x = x, y = y, fill = size)
  ) +
  geom_text(
    data = df_ratio_comp_stats,
    aes(x = x, y = y, label = size)
  ) + 
  annotate(
    geom = 'text', x = 1e-62, y = 1e-7, label = 'Number of genes', size = 12/.pt
  ) +
  labs(x = 'HSIC-GC p-value', y = 'HSIC-IR p-value', color = 'Perplexity') + 
  scale_color_viridis_b() + 
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = 1) + 
  scale_x_continuous(trans = c("log10", "reverse")) + 
  scale_y_continuous(trans = c("log10", "reverse")) + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    legend.position = c(0.8, 0.7)
  ) + 
  guides(fill = "none")
m2

ggsave(sprintf("%s/figures/mob/mob_sv_pval.png", res_dir), m2, width = 4.5, height = 3.5)

## Fig S2C: SVS gene expression and perplexity
m3.1 <- ggplot(df_sv_pval, aes(x = is_ont_svs, y = count_avg)) + 
  geom_boxplot(aes(fill = is_ont_svs), width = 0.5) + 
  scale_fill_manual(values = c('TRUE'='#FF9F1C', 'FALSE'='#718096')) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1.4, 
    label.y = 1
  ) +
  scale_y_log10() + 
  labs(x = 'ONT SVS', y = 'ONT per-spot UMI') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    legend.position='none'
  )
m3.2 <- ggplot(df_sv_pval, aes(x = is_ont_sve, y = count_avg)) + 
  geom_boxplot(aes(fill = is_ont_sve), width = 0.5) + 
  scale_fill_manual(values = c('TRUE'='#2EC4B6', 'FALSE'='#718096')) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1, 
    label.y = 0.5
  ) +
  scale_y_log10() + 
  labs(x = 'ONT SVE', y = '') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    legend.position='none'
  )
m3.3 <- ggplot(df_sv_pval, aes(x = is_ont_svs, y = perplexity)) + 
  geom_boxplot(aes(fill = is_ont_svs), width = 0.5) + 
  scale_fill_manual(values = c('TRUE'='#FF9F1C', 'FALSE'='#718096')) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1.3, 
    label.y = 2.5
  ) +
  labs(x = 'ONT SVS', y = 'Perplexity') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    legend.position='none'
  )

m3 <- cowplot::plot_grid(
  m3.1, m3.2, NULL, m3.3, 
  rel_widths = c(1, 1, 0.1, 1),
  align = c('h'), nrow = 1)
m3

ggsave(sprintf("%s/figures/mob/mob_svs_stats.png", res_dir), m3, width = 5.5, height = 3.5)

## Fig S2D: GO enrichment analysis
# load enrichment results
df_kegg <- read_csv(sprintf("%s/figures/mob/source_data/kegg_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    # geneset = factor(geneset, levels = c('SVE', 'SVS'), labels = c('SVENS', 'SVS'))
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced', 
                                'Spatially variably spliced'))
  ) %>%
  # keep the top 8 terms with the largest precision difference
  dplyr::slice((n()-15):n())

df_reac <- read_csv(sprintf("%s/figures/mob/source_data/reac_top20.csv", res_dir)) %>%
  # rename selected terms
  mutate(
    name = recode(
      name,
      'Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)' = 'NMD independent of the Exon Junction Complex (EJC)'
    )
  ) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced', 
                                'Spatially variably spliced'))
  ) %>%
  # keep the top 8 terms with the largest precision difference
  dplyr::slice((n()-15):n())

m4.1 <- ggplot(df_kegg, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = '', fill = 'Gene set', x = 'KEGG') + 
  coord_flip() +
  scale_fill_manual(values=c('Spatially variably spliced'= '#FF9F1C', 
                             'Spatially variably expressed but not spliced'= '#2EC4B6')) + 
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.position='inside',
    legend.position.inside = c(-1.6, -0.15),
    legend.background = element_blank()
  )
m4.2 <- ggplot(df_reac, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = 'Precision (proportion of term genes)', fill = 'Gene set',
       x = 'Reactome') +
  coord_flip() +
  scale_fill_manual(values=c('Spatially variably spliced'= '#FF9F1C', 
                             'Spatially variably expressed but not spliced'= '#2EC4B6')) + 
  theme_classic() +
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.text = element_text(size = 12, family = 'Arial'),
    legend.position='none',
    legend.background = element_blank()
  )
m4 <- cowplot::plot_grid(m4.1, m4.2, nrow = 2, align = 'hv', rel_heights = c(1,1))
m4

ggsave(sprintf("%s/figures/mob/mob_enrich_new.png", res_dir), m4, width = 7.5, height = 4.5)

## Fig S2H: Clustering
df_clu <- read_csv(sprintf("%s/figures/mob/source_data/clu_res.csv", res_dir))
clu_palette = c(
  'ONL'= '#6A0573',
  'EPL'= '#F2A03D',
  'GL'= '#009688',
  'GCL'= '#E63946',
  'MCL'= '#A8DADC'
)
m5.1 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = region_abbr)) + 
  scale_color_manual(values = clu_palette) + 
  labs(color = '', title = 'Short-reads-defined regions') + 
  theme_void() + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    legend.text = element_text(size = 14, family = 'Arial'),
    panel.background = element_rect(color = 'black', linewidth = 1),
    legend.key = element_rect(color = NA),
    plot.title = element_text(hjust = 0.5)
  )
m5.2 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_ont_sve)) + 
  scale_color_manual(values = clu_palette) + 
  labs(color = '', title = 'SVE (gene expression)') + 
  theme_void() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
m5.3 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_ont_svs)) + 
  scale_color_manual(values = clu_palette) + 
  labs(color = '', title = 'SVS (isoform ratio)') + 
  theme_void() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )

m5 <- cowplot::plot_grid(
  m5.1, m5.2, m5.3, rel_widths = c(1.25, 1, 1),
  align = c('h'), nrow = 1
)
m5

ggsave(sprintf("%s/figures/mob/mob_clusters_sv.png", res_dir), m5, width = 10, height = 3)

## Figure S2E: AS type
df_events_svs <- read_table(
  sprintf("%s/transcripts/mob/suppa.out/mob_svs/n_events.txt", res_dir),
  col_names = c('type', 'count')
)
df_events_svens <- read_table(
  sprintf("%s/transcripts/mob/suppa.out/mob_svens/n_events.txt", res_dir),
  col_names = c('type', 'count')
)
df_events <- rbind(
  df_events_svs %>% mutate(dataset = 'SVS'),
  df_events_svens %>% mutate(dataset = 'SVENS')
) %>%
  group_by(dataset) %>% 
  mutate(prop = count / sum(count))

pval <- chisq.test(
  matrix(c(df_events_svs$count, df_events_svens$count), ncol = 2)
)$p.value

m6 <- ggplot(df_events, aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity", width = 0.75) + 
  geom_text(aes(label = count, y = prop), position = position_stack(vjust = .5)) + 
  labs(y = 'Proportion', x = '', fill = 'AS event', 
       title = sprintf('Chi-square (contingency) p-value: %.3e', pval)) +
  scale_fill_brewer(type = 'qual', palette = 7) + 
  coord_flip() + 
  theme_classic() + 
  theme(
      legend.position = 'bottom',
      text = element_text(size = 14, family='arial'),
      axis.text = element_text(size = 14, family='arial'),
      plot.title = element_text(size = 14, hjust = 0.5)
  ) + 
  guides(fill = guide_legend(nrow = 1))
m6

ggsave(sprintf("%s/figures/mob/mob_as_type.png", res_dir), m6, width = 6, height = 2.2)

## Fig S2F: Matt exon features
df_matt <- read_tsv(sprintf("%s/transcripts/mob/matt.out/exon.with_efeatures.matt.tab", res_dir)) %>%
  mutate(GROUP = factor(GROUP, levels = c('mob_svens', 'mob_svs'), labels = c('SVENS', 'SVS')))

m7.1 <- ggplot(df_matt, aes(x = GROUP, y = EXON_MEDIANRELATIVERANK)) + 
  geom_boxplot(aes(fill = GROUP), width = 0.5) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1.4,
    label.y = 0.4
  ) +
  labs(y = 'Exon median rank', x = '') + 
  scale_fill_manual(values = c('SVS' = '#FF9F1C', 'SVENS' = '#2EC4B6')) +
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    legend.position = 'none'
  )
m7.2 <- ggplot(df_matt, aes(x = GROUP, y = PROP_EXON_IN_UTR)) + 
  geom_boxplot(aes(fill = GROUP), width = 0.5) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1.5,
    label.y = 0.4
  ) +
  labs(x = '', y = 'Proportion of UTR exons') +
  scale_fill_manual(values = c('SVS' = '#FF9F1C', 'SVENS' = '#2EC4B6')) +
  coord_flip() + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    legend.position = 'none'
  )

m7 <- cowplot::plot_grid(m7.1, m7.2, nrow = 1, rel_widths = c(1, 1))
m7

ggsave(sprintf("%s/figures/mob/mob_matt_ft.png", res_dir), m7, width = 6, height = 1.5)

## Extended Fig 2h: Xstreme de novo motif discovery
ref_motifs <- read_meme("~/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme")
exon_motifs <- read_meme(sprintf("%s/transcripts/mob/xstreme.out/exon/combined.meme", res_dir))
slop_motifs <- read_meme(sprintf("%s/transcripts/mob/xstreme.out/slop/combined.meme", res_dir))

# visualize the top de novo and reference motifs of exon sequences
motif_i <- exon_motifs[[1]]
motif_i@name <- "De novo motif 1\nSEA p-val=2.01e-7"
motif_i_ref <- (ref_motifs %>% filter_motifs(name = 'M298_0.6'))[[1]]
motif_i_ref@name <- "M298_0.6 (Rbfox family)"
m8.1 <- view_motifs(c(motif_i, motif_i_ref), min.overlap = 0.9) + 
  theme(plot.title = element_text(size = 14, family = 'Arial', hjust = 0.5),
        strip.text = element_text(size = 14, family = 'Arial', color = 'black'))

m8.2 <- view_motifs(exon_motifs[[2]]) + 
  labs(title = 'De novo motif 2\nSEA p-val=2.08e-4') +
  theme(plot.title = element_text(size = 14, family = 'Arial', hjust = 0.5))

m_title <- ggplot() + 
  labs(title = 'Variable exons (n=145)') + 
  theme(plot.title = element_text(size = 14, family = 'Arial', hjust = 0.5))
p1 <- cowplot::plot_grid(m_title, m8.1, m8.2, nrow = 3, rel_heights = c(0.1, 1, 0.5))
p1

# visualize the top de novo and reference motifs of exon plus 300bp intron sequences on each side
motif_i <- slop_motifs[[1]]
motif_i@name <- "De novo motif 1\nSEA p-val=1.25e-13"
motif_i_ref <- (ref_motifs %>% filter_motifs(name = 'M328_0.6'))[[1]]
motif_i_ref@name <- "M328_0.6 (Elavl family)"
m8.3 <- view_motifs(c(motif_i, motif_i_ref), min.overlap = 0.9) + 
  theme(plot.title = element_text(size = 14, family = 'Arial', hjust = 0.5),
        strip.text = element_text(size = 14, family = 'Arial', color = 'black'))

m8.4 <- view_motifs(slop_motifs[[2]]) + 
  labs(title = 'De novo motif 2\nSEA p-val=7.23e-12') +
  theme(plot.title = element_text(size = 14, family = 'Arial', hjust = 0.5))

m_title <- ggplot() + 
  labs(title = 'Variable exons ±300bp (n=145)') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p2 <- cowplot::plot_grid(m_title, m8.3, m8.4, nrow = 3, rel_heights = c(0.1, 1, 0.5))
p2

# combine and visualize
m8 <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1))
m8

ggsave(sprintf("%s/figures/mob/mob_xstreme_motif.png", res_dir), m8, width = 6, height = 4)


## Fig S2I: RBP expr clustering
m9 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) +
  geom_point(aes(color = clu_visium_rbp_sve)) +
  scale_color_manual(values = clu_palette) +
  labs(color = '', title = 'SVE RBP (gene expression)') +
  theme_void() +
  guides(color = guide_legend(override.aes = list(size = 2))) +
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.text = element_text(size = 12, family = 'Arial'),
    panel.background = element_rect(color = 'black', linewidth = 1),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
m9
ggsave(sprintf("%s/figures/mob/mob_clusters_rbp.png", res_dir), m9, width = 3, height = 3)

## Fig S2I: DU RBP results
rbp_with_motifs <- read_table("~/reference/cisbp-rna/mouse_pm/rbp_with_known_cleaned_motif.txt", col_names = 'rbp')
rbp_with_clip <- read_table("~/reference/POSTAR3/mouse.rbp_with_clip.txt", col_names = 'rbp')
df_du_rbp <- read_csv(sprintf("%s/figures/mob/source_data/du_res_annot.csv", res_dir)) %>%
  mutate(label = paste(gene, covariate, sep = ":")) %>%
  mutate(
    has_clip = covariate %in% rbp_with_clip$rbp,
    has_motif = covariate %in% rbp_with_motifs$rbp,
  ) %>%
  mutate(highlight = label %in% c(
    'Myl6:Rbfox3', 'Rps9:Cirbp'
  ))

# visualize p-values
m10.1 <- ggplot(df_du_rbp, aes(y = pvalue_glmm, x = `pvalue_hsic-gp`)) +
  geom_point(color = 'gray', alpha = 0.2) +
  geom_hline(yintercept = 0.05, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.05, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_text_repel(
    data = df_du_rbp %>% filter(has_clip & is_significant),
    aes(label = label, color = highlight),
    fontface = 'italic'
  ) + 
  geom_point(
    data = df_du_rbp %>% filter(has_clip & is_significant),
    aes(color = highlight)
  ) + 
  scale_color_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'black')) + 
  labs(y = 'GLMM-score p-value', x = 'Conditional HSIC p-value', 
       title = 'RBP-SVS validation with CLIP') + 
  scale_x_continuous(trans = c("log10", "reverse")) + 
  scale_y_continuous(trans = c("log10", "reverse")) + 
  theme_classic() + 
  theme(
    plot.title = element_text(size = 14, family = 'Arial'),
    text=element_text(size = 14, family = 'Arial'), 
    legend.position = 'none'
  )

# rank RBP by number of significant associations
df_rbp_rank <- df_du_rbp %>%
  summarise(count = sum(is_significant), .by = c(covariate, has_clip, has_motif)) %>%
  mutate(
    rank = rank(-count, ties.method = 'first'),
    covariate = factor(covariate, levels = covariate[order(count)])
  )

m10.2 <- ggplot(df_rbp_rank, aes(x = rank, y = count)) + 
  geom_line() + 
  geom_point(color = '#668586') +
  geom_label_repel(
    data = df_rbp_rank %>% filter((!has_clip) & (!has_motif)) %>% 
      slice_max(count, n = 10, with_ties = FALSE),
    aes(label = covariate), color = '#668586',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = -0.5, nudge_y = -2
  ) +
  geom_point(
    data = df_rbp_rank %>% filter(has_clip | has_motif) %>% slice_max(count, n = 7, with_ties = FALSE),
    color = '#8B0000',
  ) + 
  geom_label_repel(
    data = df_rbp_rank %>% filter(has_clip | has_motif) %>% slice_max(count, n = 7, with_ties = FALSE),
    aes(label = covariate), color = '#8B0000',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = 1, nudge_y = 2
  ) +
  scale_x_log10() + 
  labs(x = 'Rank', y = 'Number of significant associations') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
  )

m10 <- cowplot::plot_grid(m10.2, m10.1, nrow = 1, rel_widths = c(1, 1))
m10

ggsave(sprintf("%s/figures/mob/mob_du_selected_rbp.png", res_dir), m10, width = 7, height = 3.5)

## Fig S2K: RBP expression
df_rbp_expr <- read_csv(sprintf("%s/figures/mob/source_data/rbp_expr.csv", res_dir))
df_svs_ratio <- read_csv(sprintf("%s/figures/mob/source_data/svs_ratio.csv", res_dir))

# spatial expression of selected RBP
p_list <- c()
for (rbp in c('Rbfox3', 'Fmr1', 'Cirbp', 'Elavl4', 'Hnrnph2', 'Zfp36l1')){
  p <- df_rbp_expr %>% filter(
    gene == rbp, layer == 'log1p', array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = expression)) + 
    geom_point(size = 0.5) + 
    labs(color = '', title = rbp) + 
    scale_color_gradient2(low = 'white', high = ifelse(
      rbp %in% c('Rbfox3', 'Fmr1', 'Cirbp'),
      '#8B0000', '#668586')
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
      legend.key.width = unit(0.05, "in"),
      legend.key.height = unit(0.2, "in")
    )
  p_list <- c(p_list, list(p))
}
m11 <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
m11

ggsave(sprintf("%s/figures/mob/mob_selected_rbp_expr.png", res_dir), m11, width = 6.5, height = 3.5)

### Spatial isoform ratio of selected SVS genes
## Myl6
iso_to_plot <- c(
  'Myl6..ENSMUST00000164181.1' = 'Myl6-201',
  'Myl6..ENSMUST00000218127.1' = 'Myl6-206'
)
m_list <- c()
for (i in seq_along(iso_to_plot)){
  p <- df_svs_ratio %>% filter(
    isoform == names(iso_to_plot)[i], 
    layer == 'ratios_smoothed',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.8) +
    labs(color = '', title = iso_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
    # scale_color_gradient(low = 'white', high = '#008080') + 
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  m_list <- c(m_list, list(p))
}
p1 <- cowplot::plot_grid(plotlist = m_list, nrow = 1, align = 'hv')
p1

# expression vs ratio
df <- df_svs_ratio %>% filter(
  isoform %in% names(iso_to_plot),
  layer == 'ratios_smoothed'
) %>%
  # add transcript names according to iso_to_plot (names -> values)
  mutate(isoform = factor(isoform, levels = names(iso_to_plot), labels = iso_to_plot)) %>%
  left_join(
    df_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  facet_wrap(~ isoform, scales = 'free') +
  # geom_point(aes(color = ratio)) + 
  # scale_color_distiller(palette = 'Spectral') +
  geom_point(alpha = 0.5, color = 'gray') +
  stat_density_2d(
    data = df %>% filter(expression > 0), color = '#008080', bins = 5
  ) +
  labs(x = 'Rbfox3 log-normalized expression', y = 'Smoothed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 16, family = 'Arial'),
    strip.text = element_text(size = 16, family = 'Arial'),
    strip.background = element_blank()
  )

p <- cowplot::plot_grid(p1, NULL, p2, nrow = 1, rel_widths = c(0.6,0.05, 0.5))
p

ggsave(sprintf("%s/figures/mob/ratio_Myl6.png", res_dir), p, width = 11, height = 2.5)

# Cltb
iso_to_plot <- c(
  'Cltb..ENSMUST00000049575.7' = 'Cltb-201',
  'Cltb..ENSMUST00000091609.10' = 'Cltb-202'
)
m_list <- c()
for (i in seq_along(iso_to_plot)){
  p <- df_svs_ratio %>% filter(
    isoform == names(iso_to_plot)[i], 
    layer == 'ratios_smoothed',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.8) +
    labs(color = '', title = iso_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') + 
    # scale_color_gradient(low = 'white', high = '#0066CC') + 
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  m_list <- c(m_list, list(p))
}
p1 <- cowplot::plot_grid(plotlist = m_list, nrow = 1, align = 'hv')

# expression vs ratio
df <- df_svs_ratio %>% filter(
  isoform %in% names(iso_to_plot),
  layer == 'ratios_smoothed'
) %>%
  # add transcript names according to iso_to_plot (names -> values)
  mutate(isoform = factor(isoform, levels = names(iso_to_plot), labels = iso_to_plot)) %>%
  left_join(
    df_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  facet_wrap(~ isoform, scales = 'free') +
  geom_point(alpha = 0.5, color = 'gray') + 
  stat_density_2d(
    data = df %>% filter(expression > 0), color = '#0066CC', bins = 5
  ) + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'Smoothed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 16, family = 'Arial'),
    strip.text = element_text(size = 16, family = 'Arial'),
    strip.background = element_blank()
  )

p <- cowplot::plot_grid(p1, NULL, p2, nrow = 1, rel_widths = c(0.6,0.05, 0.5))
p

ggsave(sprintf("%s/figures/mob/ratio_Cltb.png", res_dir), p, width = 11, height = 2.5)


### MOB SVS-RBP pair genomic coordinates
source("~/Projects/SPLISOSM_paper/scripts/sit_data_analysis/visualization/plotgardener_isoform_structure.R")

## Load neccessary data for visualization
res_dir <- "~/Projects/SPLISOSM_paper/results/sit_nar_23"
fimo_dir <- sprintf('%s/transcripts/mob/fimo.out/', res_dir)

# Load MOB SVS transcript GTF
svs_gtf <- import.gff(sprintf('%s/transcripts/mob/mob_svs.txid.gtf', res_dir))
svs_gtf$gene_id <- gsub("\\..*", "", svs_gtf$gene_id)

# Find and load FIMO motif scanning results
fimo_rbp_dirs <- list.files(fimo_dir, pattern = 'mob_svs_', full.names = FALSE)
rbs_fimo_gtf <- lapply(fimo_rbp_dirs, function(rbp_dir){
  rbp_gtf <- import.gff(sprintf('%s/%s/fimo.gff', fimo_dir, rbp_dir))
  rbp_gtf$RBP_Name <- toupper(gsub('mob_svs_', '', rbp_dir))
  rbp_gtf$RBP_Name <- paste(rbp_gtf$RBP_Name, '(motif)')
  return(rbp_gtf)
}) %>% do.call(c, .)

# Load POSTAR RBP binding sites
rbs_clip_gtf <- import.gff('~/reference/POSTAR3/mouse.gtf')
rbs_clip_gtf$RBP_Name <- paste(rbs_clip_gtf$gene_id, '(CLIP)')

# Combine all RBP binding sites
rbs_all_gtf <- c(rbs_fimo_gtf, rbs_clip_gtf)


## Myl6 transcript structure
png(
  sprintf("%s/figures/mob/iso_struct_myl6.png", res_dir), 
  width = 7.25, height = 2, units = "in", res = 300
)
plot_transcripts(
  c('ENSMUST00000218127.1', 'ENSMUST00000164181.1'),
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'CELF2 (CLIP)', 'FMR1 (CLIP)', 'FMR1 (motif)'),
  rbs_gtf = rbs_all_gtf,
)
# Highlight exon
plotRect(
  x = 2.65, y = 0.2, width = 0.2, height = 0.4,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox2 KD – WT) = -0.26', fontsize = 8,
  x = 2.2 , y = 0.1,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# plot_transcripts(
#   c('ENSMUST00000218127.1', 'ENSMUST00000164181.1'),
#   ref_gtf = svs_gtf,
#   rbp_names = c('RBFOX2 (CLIP)', 'CELF2 (CLIP)', 'RBFOX3 (CLIP)', 'RBFOX3 (motif)'),
#   rbs_gtf = rbs_all_gtf,
# )

## Cltb transcripts structure
png(
  sprintf("%s/figures/mob/iso_struct_cltb.png", res_dir), 
  width = 7.25, height = 2, units = "in", res = 300
)
# rbp_to_plot <- rbs_all_gtf[seqnames(rbs_all_gtf) == "chr13" & 
#               start(rbs_all_gtf) >= 54592401 &
#               end(rbs_all_gtf) <= 54611344]$RBP_Name %>% unique()

plot_transcripts(
  c('ENSMUST00000049575.7', 'ENSMUST00000091609.10'),
  ref_gtf = svs_gtf,
  # rbp_names = rbp_to_plot,
  rbp_names = c('RBFOX3 (CLIP)', 'RBFOX3 (motif)', 'CELF4 (CLIP)', 'CELF4 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#0066CC'
)
# Highlight exon
plotRect(
  x = 2.4, y = 0.2, width = 0.2, height = 0.4,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKO – WT) = 0.15', fontsize = 8,
  x = 2 , y = 0.1,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# ## Rps9 transcripts and RBP binding sites
# png(
#   sprintf("%s/figures/mob/iso_struct_rps9.png", res_dir), 
#   width = 7.25, height = 2, units = "in", res = 300
# )
# plot_transcripts(
#   svs_gtf[svs_gtf$gene_name == 'Rps9']$transcript_id %>% unique(),
#   ref_gtf = svs_gtf,
#   rbp_names = c('CIRBP (CLIP)', 'RBM3 (CLIP)', 'RBFOX2 (CLIP)', 'LIN28A (CLIP)'),
#   rbs_gtf = rbs_all_gtf,
#   fill = '#8B0000'
# )
# pageGuideHide()
# dev.off()

### =======
# Anapc11
plot_transcripts(
  c('ENSMUST00000026128.9', 'ENSMUST00000093140.4'),
  ref_gtf = svs_gtf,
  rbp_names = c('ZFP36 (CLIP)', 'ZFP36L1 (motif)', 'CELF2 (CLIP)', 'CELF4 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Plp1
plot_transcripts(
  c('ENSMUST00000113085.1', 'ENSMUST00000033800.12'),
  ref_gtf = svs_gtf,
  rbp_names = c('ZFP36 (CLIP)', 'ZFP36L1 (motif)', 'FMR1 (CLIP)', 'FMR1 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Nnat
plot_transcripts(
  c('ENSMUST00000173839.1', 'ENSMUST00000153739.8',
    'ENSMUST00000134513.1', 'ENSMUST00000109526.1',
    'ENSMUST00000173595.1'),
  ref_gtf = svs_gtf,
  rbp_names = c('ZFP36 (CLIP)', 'ZFP36L1 (motif)', 'FMR1 (CLIP)', 'FMR1 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Clta
plot_transcripts(
  c('ENSMUST00000170241.7', 'ENSMUST00000107851.9',
    'ENSMUST00000107849.9', 'ENSMUST00000107847.9',
    'ENSMUST00000107845.3'),
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'RBFOX3 (motif)', 'CELF2 (CLIP)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Rps24
plot_transcripts(
  c('ENSMUST00000223999.1', 'ENSMUST00000224568.1',
    'ENSMUST00000169826.2', 'ENSMUST00000225023.1'),
  ref_gtf = svs_gtf,
  rbp_names = c('ZFP36 (CLIP)', 'ZFP36L1 (motif)', 'FMR1 (CLIP)', 'FMR1 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Rps9
plot_transcripts(
  svs_gtf[svs_gtf$gene_name == 'Rps9']$transcript_id %>% unique(),
  ref_gtf = svs_gtf,
  rbp_names = c('CIRBP (CLIP)', 'CIRBP (motif)', 'FMR1 (CLIP)', 'FMR1 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Pigp
plot_transcripts(
  svs_gtf[svs_gtf$gene_name == 'Pigp']$transcript_id %>% unique(),
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX2 (motif)', 'RBFOX3 (CLIP)', 'RBFOX3 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)

# Mrpl21
plot_transcripts(
  svs_gtf[svs_gtf$gene_name == 'Mrpl21']$transcript_id %>% unique(),
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX2 (motif)', 'RBFOX3 (CLIP)', 'RBFOX3 (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)


