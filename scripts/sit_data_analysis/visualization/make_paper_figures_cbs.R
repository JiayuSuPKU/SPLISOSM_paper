library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(universalmotif)

extrafont::loadfonts()

res_dir <- "~/Projects/SPLISOSM_paper/results/sit_nar_23"
date <- "1107"

### CBS main figures
## make figure directory if doesn't exist
dir.create(sprintf("%s/figures/cbs", res_dir), recursive = TRUE, showWarnings = FALSE)

## Fig 2a, 2d: Clustering
df_clu <- read_csv(sprintf("%s/figures/cbs/source_data/clu_res.csv", res_dir))

# match clusters to the reference annotation
library(clue)
match_labels <- function(ref_labels, alt_labels){
  ctable <- as.matrix(table(alt_labels, ref_labels))
  alignment <- solve_LSAP(ctable, maximum = TRUE)
  alt2ref_mapping <- cbind(
    rownames(ctable),
    colnames(ctable)[alignment]
  )
  return(alt2ref_mapping)
}

for (clu_col in c('clu_ont_sve', 'clu_ont_svs', 'clu_visium_rbp_sve')){
  alt2ref_mapping <- match_labels(df_clu$region, df_clu[[clu_col]])
  df_clu[[paste0(clu_col, '_anno')]] <- factor(
    df_clu[[clu_col]],
    levels = alt2ref_mapping[,1],
    labels = alt2ref_mapping[,2]
  )
}

# filter out outlier spots
df_clu <- df_clu %>% filter(array_row > 15)

# set color palette
clu_palette <- setNames(
  RColorBrewer::brewer.pal(12, 'Paired'),
  df_clu$region %>% unique() %>% sort()
)
clu_palette['Retrosplenial area'] <- '#2C3E50'

# reference palette
clu_palette2 <- setNames(
  RColorBrewer::brewer.pal(12, 'Set3'),
  df_clu$region %>% unique() %>% sort()
)
clu_palette2['CA3'] <- '#2C3E50'

# Short-reads-based reference annotation
m1 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = region), size = 1) + 
  scale_color_manual(values = clu_palette2) +
  # scale_color_brewer(palette = 'Set3') +
  labs(color = 'Regions based on gene expression') + 
  theme_void() + 
  guides(color = guide_legend(override.aes = list(size = 2), ncol = 2)) + 
  theme(
    # text = element_text(size = 12, family = 'Arial'),
    legend.title = element_text(size = 14, family = 'Arial'),
    legend.text = element_text(size = 14, family = 'Arial'),
    panel.background = element_rect(color = 'black', linewidth = 1),
    legend.key = element_rect(color = NA),
    plot.title = element_text(hjust = 0.5)
  )
m1

ggsave(sprintf("%s/figures/cbs/main_clu_ref_region.png", res_dir), m1, width = 7, height = 3)
ggsave(sprintf("%s/figures/cbs/main_clu_ref_region.pdf", res_dir), m1, width = 7, height = 3)

# ONT SVE and SVS clustering results
m2.1 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_ont_sve_anno), size = 1.5) + 
  scale_color_manual(values = clu_palette) +
  labs(color = '', title = 'SVE (gene expression)') + 
  theme_void() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
m2.2 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_ont_svs_anno), size = 1.5) + 
  scale_color_manual(values = clu_palette) +
  labs(color = '', title = 'SVS (isoform ratio)') + 
  theme_void() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
m2 <- cowplot::plot_grid(m2.1, m2.2, nrow = 1, align = 'h', rel_widths = c(1, 1))
m2

ggsave(sprintf("%s/figures/cbs/main_clu_sv.png", res_dir), m2, width = 7, height = 4)
ggsave(sprintf("%s/figures/cbs/main_clu_sv.pdf", res_dir), m2, width = 7, height = 4)

## load SV data
df_sv_cbs1 <- read_csv(sprintf("%s/sv_results/cbs1_sv_combined_%s.csv", res_dir, date)) %>%
  mutate(
    is_ont_svs = `padj_hsic-ir` < 0.05,
    is_ont_sve = `padj_hsic-gc` < 0.05,
    is_ont_expressed = count_avg > 0.5
  )
df_sv_cbs2 <- read_csv(sprintf("%s/sv_results/cbs2_sv_combined_%s.csv", res_dir, date)) %>%
  mutate(
    is_ont_svs = `padj_hsic-ir` < 0.05,
    is_ont_sve = `padj_hsic-gc` < 0.05,
    is_ont_expressed = count_avg > 0.5
  )

# combine two replicates
df_sv_pval <- inner_join(df_sv_cbs1, df_sv_cbs2, by = 'gene', suffix = c('.1', '.2')) %>%
  mutate(
    is_ont_svs = is_ont_svs.1 & is_ont_svs.2,
    is_ont_sve = is_ont_sve.1 & is_ont_sve.2,
    is_ont_expressed = is_ont_expressed.1 & is_ont_expressed.2
  )


## Gene clustering results
library(ggtree)

# load RV coefficients
cor_matrix <- read.csv(sprintf("%s/figures/cbs/source_data/rvcoeff.csv", res_dir), row.names = 1)
dist_matrix <- as.dist(1 - cor_matrix)

# hierarchical clustering
hc <- hclust(dist_matrix, method = "ward.D2")  # Use Ward's method
num_clusters <- 3  # Change this to the desired number of clusters
cluster_labels <- cutree(hc, k = num_clusters)

# group genes by clusters
g <- split(names(cluster_labels), cluster_labels)

# create ggtree object
p <- ggtree(hc, layout = "rectangular", linetype = "solid", size = 0.5)  # Root on the left by default

# identify and group clades for coloring
clades <- sapply(g, function(n) MRCA(p, n))  # Find the Most Recent Common Ancestor (MRCA)
p <- groupClade(p, clades, group_name = "subtree") + aes(color = subtree)

# create a data frame for additional tip annotations
d <- data.frame(
  label = names(cluster_labels),
  cluster = cluster_labels, 
  gene = rownames(cor_matrix)  # Add gene names from the correlation matrix
) %>% left_join(
  df_sv_pval %>% select(gene, `pvalue_hsic-ir.2`),
  by = 'gene'
) %>%
  mutate(logpval = log10(`pvalue_hsic-ir.2`))

# show only top5 genes per group
top_d <- d %>% 
  group_by(cluster) %>%
  slice_min(order_by = logpval, n = 5) %>%  # Select top 5 genes per cluster
  ungroup()
d <- d %>% mutate(label_to_display = ifelse(gene %in% top_d$gene, gene, "")) 

# Final dendrogram plot
m_tree <- p %<+% d + 
  # geom_tippoint(
  #   aes(fill = logpval),
  #   size = 1, shape = 21, color = "black"
  # ) + 
  labs(color = 'Programs') + 
  scale_y_continuous(transform = "reverse") +  # Reverse the y-axis
  geom_tiplab(mapping = aes(label = label_to_display), fontface = "italic") +  # Label tips with gene names
  scale_color_brewer(palette = "Set1") +  # Color clades with a color palette
  scale_fill_distiller(palette = 'Reds', direction = -1) +
  theme_tree2() +  # Add axes and grid lines
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = c(0.2, 0.4),
    legend.background = element_blank(),
  )  # Adjust plot margins
m_tree
ggsave(sprintf("%s/figures/cbs/gene_tree.pdf", res_dir), m_tree, width = 1.5, height = 4)

## Plot example genes per program
expr_ont <- read_csv(sprintf("%s/figures/cbs/source_data/gene_program_top5_ratio.csv", res_dir))

iso_to_plot <- c('Arpp19-201', 'Arpp19-203',
                 'Septin5-201', 'Septin5-205',
                 'Snap25-201', 'Snap25-202')
palette_list <- c(
  'Arpp19-201' = 'Greens',
  'Arpp19-203' = 'Greens',
  'Septin5-201' = 'Purples',
  'Septin5-205' = 'Purples',
  'Snap25-201' = 'Blues',
  'Snap25-202' = 'Blues'
)

# log1p
p_list <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    transcript_name == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    filter(array_row > 15) %>%
    ggplot(aes(x = array_col, y = - array_row, color = value)) + 
    geom_point(size = 0.2) +
    coord_fixed(2) +
    labs(color = 'LogExpr', title = sprintf("%s", iso_to_plot[i])) +
    scale_color_distiller(palette = palette_list[iso_to_plot[i]], direction = 1) +
    theme_void() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      strip.text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 12, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list <- c(p_list, list(p))
}

p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv', byrow = FALSE)
p
ggsave(sprintf("%s/figures/cbs/gene_program_top1.png", res_dir), p, width = 6, height = 4)
ggsave(sprintf("%s/figures/cbs/gene_program_top1.pdf", res_dir), p, width = 6, height = 4)


## Fig 2b: HSIC-GC vs HSIC-IR
df_ratio_comp_stats = data.frame(
  x = c(1e-20, 1e-20, 1e-45, 1e-45), 
  y = c(1e-165, 1e-185, 1e-165, 1e-185), 
  size = table(
    df_sv_pval$is_ont_svs, 
    df_sv_pval$is_ont_sve
  ) %>% as.vector() # FF, FT, TF, TT
)

m3 <- ggplot(df_sv_pval, aes(x = `pvalue_hsic-gc.2` + 1e-200, y = `pvalue_hsic-ir.2` + 1e-200)) +
  geom_point(color = 'gray', alpha = 0.2) +
  geom_hline(yintercept = 0.05, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.05, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_point(
    data = (df_sv_pval %>% filter(is_ont_svs & is_ont_expressed)),
    aes(color = perplexity.2)
  ) +
  geom_text_repel(
    data = (df_sv_pval %>% filter(is_ont_expressed) %>% top_n(n = -15, wt = `pvalue_hsic-ir.2`)),
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
    geom = 'text', x = 1e-45, y = 1e-205, label = 'Number of genes', size = 12/.pt
  ) +
  labs(x = 'HSIC-GC p-value (CBS2)', y = 'HSIC-IR p-value (CBS2)', color = 'Perplexity') + 
  scale_color_viridis_b() + 
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = 1) + 
  scale_x_continuous(trans = c("log10", "reverse")) + 
  scale_y_continuous(trans = c("log10", "reverse")) + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position.inside = c(0.8, 0.7)
  ) + 
  guides(fill = "none")
m3

ggsave(sprintf("%s/figures/cbs/main_gc_vs_ir.png", res_dir), m3, width = 5.25, height = 4)
ggsave(sprintf("%s/figures/cbs/main_gc_vs_ir.pdf", res_dir), m3, width = 5.25, height = 4)


## Fig 2c: functional enrichment analysis
# load enrichment results
df_go <- read_csv(sprintf("%s/figures/cbs/source_data/go_top20.csv", res_dir)) %>%
  mutate(
    name = str_to_sentence(name), # capitalize first letter,
    # name = stringr::str_wrap(name, 35), # auto line break
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = TRUE)])),
    geneset = factor(geneset, levels = c('SVS', 'SVE'), labels = c('SVS', 'SVENS'))
  ) %>%
  # keep the top 10 terms with the largest and smallest precision difference each
  dplyr::slice(1:10, (n()-9):n())
  
df_kegg <- read_csv(sprintf("%s/figures/cbs/source_data/kegg_top20.csv", res_dir)) %>%
  mutate(
    name = stringr::str_wrap(name, 35), # auto line break
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = TRUE)])),
    geneset = factor(geneset, levels = c('SVS', 'SVE'), labels = c('SVS', 'SVENS'))
  ) %>%
  # keep the top 10 terms with the largest and smallest precision difference each
  dplyr::slice(1:10, (n()-9):n())

df_reac <- read_csv(sprintf("%s/figures/cbs/source_data/reac_top20.csv", res_dir)) %>%
  # rename selected terms
  mutate(
    name = recode(
      name,
        'Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)' = 'NMD enhanced by the Exon Junction Complex (EJC)',
        'Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)' = 'NMD independent of the Exon Junction Complex (EJC)'
    )
  ) %>%
  mutate(
    name = stringr::str_wrap(name, 35), # auto line break
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = TRUE)])),
    geneset = factor(geneset, levels = c('SVS', 'SVE'), labels = c('SVS', 'SVENS'))
  ) %>%
  # keep the top 10 terms with the largest and smallest precision difference each
  dplyr::slice(1:10, (n()-9):n())

# # visualization
# m4.1 <- ggplot(df_go, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
#   geom_bar(stat='identity', position='dodge') + 
#   labs(y = 'Precision (proportion of term genes)', x = '', fill = 'Gene set', 
#        title = 'GO:BP') + 
#   geom_vline(xintercept = 5.5, linetype = 'dashed') + 
#   # coord_flip() + 
#   scale_fill_manual(values=c('SVS'= '#FF9F1C', 'SVENS'= '#2EC4B6')) + 
#   theme_classic() + 
#   theme(
#     text=element_text(size = 12, family = 'Arial'),
#     axis.text=element_text(size = 12, family = 'Arial'),
#     axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
#     plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
#     legend.position='inside',
#     legend.position.inside = c(0.2, 0.85)
#   )
# m4.2 <- ggplot(df_kegg, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
#   geom_bar(stat='identity', position='dodge') + 
#   labs(y = 'Precision (proportion of term genes)', x = '', fill = 'Gene set', 
#        title = 'KEGG') + 
#   geom_vline(xintercept = 5.5, linetype = 'dashed') + 
#   # coord_flip() + 
#   scale_fill_manual(values=c('SVS'= '#FF9F1C', 'SVENS'= '#2EC4B6')) + 
#   theme_classic() + 
#   theme(
#     text=element_text(size = 12, family = 'Arial'),
#     axis.text=element_text(size = 12, family = 'Arial'),
#     axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
#     plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
#     legend.position='inside',
#     legend.position.inside = c(0.3,0.85)
#   )
# m4.3 <- ggplot(df_reac, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
#   geom_bar(stat='identity', position='dodge') + 
#   labs(y = 'Precision (proportion of term genes)', x = '', fill = 'Gene set', 
#        title = 'Reactome') + 
#   geom_vline(xintercept = 5.5, linetype = 'dashed') + 
#   # coord_flip() + 
#   scale_fill_manual(values=c('SVS'= '#FF9F1C', 'SVENS'= '#2EC4B6')) + 
#   theme_classic() + 
#   theme(
#     text=element_text(size = 12, family = 'Arial'),
#     axis.text=element_text(size = 12, family = 'Arial'),
#     axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
#     plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
#     legend.position='inside',
#     legend.position.inside = c(0.3,0.85)
#   )
# 
# m4 <- cowplot::plot_grid(m4.1, m4.2, m4.3, nrow = 1, align = 'h', rel_widths = c(1,1.2,1.2))
# m4
# ggsave(sprintf("%s/figures/cbs/main_enrich.png", res_dir), m4, width = 12, height = 6)

df_kegg <- read_csv(sprintf("%s/figures/cbs/source_data/kegg_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced (SVENS)', 
                                'Spatially variably spliced (SVS)'))
  ) %>%
  # keep the top 5 terms with the largest precision difference
  dplyr::slice((n()-15):n())

df_reac <- read_csv(sprintf("%s/figures/cbs/source_data/reac_top20.csv", res_dir)) %>%
  # rename selected terms
  mutate(
    name = recode(
      name,
      'Nonsense Mediated Decay (NMD) enhanced by the Exon Junction Complex (EJC)' = 'NMD enhanced by the Exon Junction Complex (EJC)',
      'Nonsense Mediated Decay (NMD) independent of the Exon Junction Complex (EJC)' = 'NMD independent of the Exon Junction Complex (EJC)'
    )
  ) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced (SVENS)', 
                                'Spatially variably spliced (SVS)'))
  ) %>%
  # keep the top 5 terms with the largest precision difference
  dplyr::slice((n()-15):n())

m4.1 <- ggplot(df_kegg, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = '', fill = 'Gene set', x = 'KEGG') + 
  coord_flip() +
  scale_fill_manual(values=c('Spatially variably spliced (SVS)'= '#FF9F1C', 
                             'Spatially variably expressed but not spliced (SVENS)'= '#2EC4B6')) + 
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.position='none',
    legend.position.inside = c(0.3,0.85)
  )
m4.2 <- ggplot(df_reac, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = 'Precision (proportion of term genes)', fill = 'Gene set',
       x = 'Reactome') +
  coord_flip() +
  scale_fill_manual(values=c('Spatially variably spliced (SVS)'= '#FF9F1C', 
                             'Spatially variably expressed but not spliced (SVENS)'= '#2EC4B6')) + 
  theme_classic() +
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.position='inside',
    legend.position.inside = c(-1.35, 1.2),
    legend.background = element_blank()
  )
m4 <- cowplot::plot_grid(m4.1, m4.2, nrow = 2, align = 'hv', rel_heights = c(1,1))
m4

ggsave(sprintf("%s/figures/cbs/main_enrich_new.png", res_dir), m4, width = 7.5, height = 4)
ggsave(sprintf("%s/figures/cbs/main_enrich_new.pdf", res_dir), m4, width = 7, height = 4)

 ## Supfig: Visualize detailed KEGG enrichment terms
library(clusterProfiler)
library(org.Mm.eg.db)
eg <- bitr(
  df_sv_pval$gene[df_sv_pval$is_ont_svs] %>% unique(),
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = org.Mm.eg.db
)
kk <- enrichKEGG(
  gene = eg$ENTREZID,
  organism = 'mmu',
  keyType = 'kegg',
  pvalueCutoff = 0.01
)

# extract genes per KEGG term
kegg_glist <- list()
for (kegg_id in c('mmu04721', 'mmu05022', 'mmu05016', 'mmu05010', 'mmu03010')) {
  glist <- kk@result[kk@result$ID == kegg_id, 'geneID'] %>% str_split_1('/') %>% bitr(
    fromType = 'ENTREZID',
    toType = 'SYMBOL',
    OrgDb = org.Mm.eg.db
  )
  kegg_glist[[kegg_id]] <- df_sv_pval %>% filter(
    gene %in% glist$SYMBOL,
    is_ont_svs
  ) %>% pull(gene)
}

# visualize the gene list
df <- df_sv_pval %>% 
  mutate(rank = rank(`pvalue_hsic-ir.2`)) %>%
  filter(is_ont_svs)
p_list <- list()
for (kegg_id in c('mmu04721', 'mmu05022', 'mmu05016', 'mmu05010', 'mmu03010')) {
  p <- ggplot(df, aes(x = rank, y = `pvalue_hsic-ir.2`)) + 
    geom_point(alpha = 0.5, color = 'gray') + 
    scale_y_continuous(transform = c('log10', 'reverse')) +
    geom_point(
      data = df %>% filter(gene %in% kegg_glist[[kegg_id]]),
      color = 'black'
    ) +
    geom_text_repel(
      data = df %>% filter(gene %in% kegg_glist[[kegg_id]]),
      aes(label = gene, x = rank, y = `pvalue_hsic-ir.2`),
      force = 20, max.overlaps = 30, nudge_x = 0, nudge_y = 50,
      fontface = 'italic'
    ) +
    labs(x = 'Rank', y = 'HSIC-IR p-value (CBS2)', title = stringr::str_wrap(
      str_split_1(kk@result[kegg_id, 'Description'], ' - Mus')[1], 30
    )) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial')
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, ncol = 5, align = 'h')
p

ggsave(sprintf("%s/figures/cbs/sup_kegg_genes.png", res_dir), p, width = 15, height = 3)

## Fig 2e: AS type
df_events_svs <- read_table(
  sprintf("%s/transcripts/cbs/suppa.out/cbs_svs/n_events.txt", res_dir),
  col_names = c('type', 'count')
)
df_events_svens <- read_table(
  sprintf("%s/transcripts/cbs/suppa.out/cbs_svens/n_events.txt", res_dir),
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

m5 <- ggplot(df_events, aes(x = dataset, y = count, fill = type)) + 
  geom_bar(position="fill", stat="identity", width = 0.75) + 
  geom_text(aes(label = count, y = prop), position = position_stack(vjust = .5)) + 
  labs(y = 'Proportion', x = '', fill = 'AS event', 
       title = sprintf('Chi-square (contingency) p-value: %.3e', pval)) +
  scale_fill_brewer(type = 'qual', palette = 7) + 
  coord_flip() + 
  theme_classic() + 
  theme(
      legend.position = 'bottom',
      text = element_text(size = 12, family='Arial'),
      axis.text = element_text(size = 12, family='Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family='Arial')
  ) + 
  guides(fill = guide_legend(nrow = 1))
m5

ggsave(sprintf("%s/figures/cbs/main_as_type.png", res_dir), m5, width = 7, height = 2)
ggsave(sprintf("%s/figures/cbs/main_as_type.pdf", res_dir), m5, width = 7, height = 2.5)

## Fig 2f: Matt exon features
df_matt <- read_tsv(sprintf("%s/transcripts/cbs/matt.out/exon.with_efeatures.matt.tab", res_dir)) %>%
  mutate(GROUP = factor(GROUP, levels = c('cbs_svens', 'cbs_svs'), labels = c('SVENS', 'SVS')))

m6.1 <- ggplot(df_matt, aes(x = GROUP, y = EXON_MEDIANRELATIVERANK)) + 
  geom_boxplot(aes(fill = GROUP), width = 0.5) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1.5,
    label.y = 0.5
  ) +
  coord_flip() +
  labs(y = "Exon median relative rank (5'>3')", x = '') + 
  scale_fill_manual(values = c('SVS' = '#FF9F1C', 'SVENS' = '#2EC4B6')) +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'none'
  )
m6.2 <- ggplot(df_matt, aes(x = GROUP, y = MEDIAN_EXON_NUMBER)) + 
  geom_boxplot(aes(fill = GROUP), width = 0.5, outliers = FALSE) + 
  stat_compare_means(
    method = 't.test', 
    label = 'p.format', 
    label.x = 1.5,
    label.y = 7
  ) +
  coord_flip() +
  labs(x = '', y = 'Median exon number of transcripts\ncontaining the AS exon') +
  scale_fill_manual(values = c('SVS' = '#FF9F1C', 'SVENS' = '#2EC4B6')) +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'none'
  )

m6 <- cowplot::plot_grid(m6.1, m6.2, nrow = 1, align = 'h', rel_widths = c(1, 1))
m6

ggsave(sprintf("%s/figures/cbs/main_matt_ft.png", res_dir), m6, width = 7, height = 2)
ggsave(sprintf("%s/figures/cbs/main_matt_ft.pdf", res_dir), m6, width = 7, height = 2)

## Fig 2g: Xstreme de novo motif discovery
ref_motifs <- read_meme("~/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme")
exon_motifs <- read_meme(sprintf("%s/transcripts/cbs/xstreme.out/exon/combined.meme", res_dir))
slop_motifs <- read_meme(sprintf("%s/transcripts/cbs/xstreme.out/slop/combined.meme", res_dir))

# visualize the top de novo and reference motifs of exon sequences
motif_i <- exon_motifs[[1]]
motif_i@name <- "De novo motif 1\nSEA p-val=2.54e-10"
motif_i_ref <- (ref_motifs %>% filter_motifs(name = 'M102_0.6'))[[1]]
motif_i_ref@name <- "M102_0.6 (Srsf1/9)"
m7.1 <- view_motifs(c(motif_i, motif_i_ref), min.overlap = 0.9) + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
        strip.text = element_text(size = 12, family = 'Arial', color = 'black'))

m7.2 <- view_motifs(exon_motifs[[2]]) +
  labs(title = 'De novo motif 2\nSEA p-val=8.09e-6') +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

m_title <- ggplot() + 
  labs(title = 'SVS variable exons (n=214)') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p1 <- cowplot::plot_grid(m_title, m7.1, m7.2, nrow = 3, rel_heights = c(0.1, 1, 0.5))
p1

# visualize the top de novo and reference motifs of exon plus 300bp intron on both sides
motif_i <- slop_motifs[[1]]
motif_i@name <- "De novo motif 1\nSEA p-val=4.61e-11"
motif_i_ref <- (ref_motifs %>% filter_motifs(name = 'M065_0.6'))[[1]]
motif_i_ref@name <- "M065_0.6 (Srsf1/9)"
m7.3 <- view_motifs(c(motif_i, motif_i_ref), min.overlap = 0.9) + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
        strip.text = element_text(size = 12, family = 'Arial', color = 'black'))

m7.4 <- view_motifs(slop_motifs[[2]]) +
  labs(title = 'De novo motif 2\nSEA p-val=8.53e-8') +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

m_title <- ggplot() + 
  labs(title = 'SVS variable exons ±300bp (n=214)') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
p2 <- cowplot::plot_grid(m_title, m7.3, m7.4, nrow = 3, rel_heights = c(0.1, 1, 0.5))
p2

# combine and visualize
m7 <- cowplot::plot_grid(p1, p2, nrow = 1, rel_widths = c(1, 1))
m7

ggsave(sprintf("%s/figures/cbs/main_xstreme_motif.png", res_dir), m7, width = 6, height = 4)


## Fig 2h: RBP expr clustering
m8 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_visium_rbp_sve_anno)) + 
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
m8
ggsave(sprintf("%s/figures/cbs/main_clu_rbp.png", res_dir), m8, width = 3.5, height = 4)
ggsave(sprintf("%s/figures/cbs/main_clu_rbp.pdf", res_dir), m8, width = 3.2, height = 4)

## Fig 2i: DU RBP results
rbp_with_motifs <- read_table("~/reference/cisbp-rna/mouse_pm/rbp_with_known_cleaned_motif.txt", col_names = 'rbp')
rbp_with_clip <- read_table("~/reference/POSTAR3/mouse.rbp_with_clip.txt", col_names = 'rbp')
df_du_rbp <- read_csv(sprintf("%s/figures/cbs/source_data/du_res_annot.csv", res_dir)) %>%
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
df_du_rbp %>% filter(is_significant) %>% nrow()

# rank RBP by number of significant associations
df_rbp_rank <- df_du_rbp %>%
  summarise(count = sum(is_significant), .by = c(covariate, has_clip, has_motif)) %>%
  mutate(
    rank = rank(-count, ties.method = 'first'),
    covariate = factor(covariate, levels = covariate[order(count)])
  )

m9 <- ggplot(df_rbp_rank, aes(x = rank, y = count)) + 
  geom_line() + 
  geom_point(color = '#668586') +
  geom_label_repel(
    data = df_rbp_rank %>% filter((!has_clip) & (!has_motif)) %>% slice_max(count, n = 15, with_ties = FALSE),
    aes(label = covariate), color = '#668586',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = -0.5, nudge_y = -2
  ) +
  geom_point(
    data = df_rbp_rank %>% filter(has_clip | has_motif) %>% slice_max(count, n = 15, with_ties = FALSE),
    color = '#8B0000',
  ) + 
  geom_label_repel(
    data = df_rbp_rank %>% filter(has_clip | has_motif) %>% slice_max(count, n = 15, with_ties = FALSE),
    aes(label = covariate), color = '#8B0000',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = 1, nudge_y = 2
  ) +
  scale_x_log10() + 
  labs(x = 'Rank', y = 'Number of significant associations', color = 'In CLIPdb') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 12),
  )
m9

ggsave(sprintf("%s/figures/cbs/main_du_rbpr_rank.png", res_dir), m9, width = 4, height = 4)
ggsave(sprintf("%s/figures/cbs/main_du_rbpr_rank.pdf", res_dir), m9, width = 4, height = 4)

## visualize p-values for selected RBPs
# Rbfox3
m10.1 <- ggplot(
  df_du_rbp %>% filter(covariate == 'Rbfox3'), 
  aes(x = `pvalue_hsic-gp_2`, y = `pvalue_hsic-gp_1`)
  # aes(y = pvalue_glmm_2 + 1e-15, x = `pvalue_hsic-gp_2`)
) +
  geom_point(color = 'gray', alpha = 0.5) +
  geom_hline(yintercept = 0.01, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.01, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_point(
    data = df_du_rbp %>% 
      filter(covariate == 'Rbfox3', is_significant) %>%
      mutate(highlight = gene %in% c('Myl6', 'Gnas', 'Capzb', 'Pkm')),
    aes(color = highlight)
  ) + 
  geom_text_repel(
    data = df_du_rbp %>% 
      filter(covariate == 'Rbfox3', is_significant) %>%
      mutate(highlight = gene %in% c('Myl6', 'Gnas', 'Capzb', 'Pkm')),
    aes(label = label, color = highlight),
    fontface = 'italic',
    force = 5, force_pull = 0.1
  ) + 
  labs(y = 'Conditional HSIC p-value (CBS1)', x = 'Conditional HSIC p-value (CBS2)', 
       color = 'Has RBFOX binding', title = 'RBP-SVS validation with brain CLIP') + 
  scale_color_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'black')) +
  scale_x_continuous(trans = c("log10", "reverse")) + 
  scale_y_continuous(trans = c("log10", "reverse")) + 
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'), 
    legend.text = element_text(size = 12, family = 'Arial'),
    legend.position = 'inside',
    legend.background = element_blank(),
    legend.position.inside = c(0.78, 0.2)
  )
m10.1

# Celf5
m10.2 <- ggplot(
  df_du_rbp %>% filter(covariate == 'Celf5'), 
  aes(x = `pvalue_hsic-gp_2`, y = `pvalue_hsic-gp_1`)
  # aes(y = pvalue_glmm_2 + 1e-15, x = `pvalue_hsic-gp_2`)
) +
  # geom_text(aes(label = label), color = 'gray') + 
  geom_point(color = 'gray', alpha = 0.5) +
  geom_hline(yintercept = 0.01, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.01, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_point(
    data = df_du_rbp %>% 
      filter(covariate == 'Celf5', is_significant) %>%
      mutate(highlight = gene %in% c('Clta', 'Myl6', 'Nnat', 'Gnas', 'Pkm', 'S100a16', 'Cdc42', 'Cnih1', 'Caly')),
    aes(color = highlight)
  ) + 
  geom_text_repel(
    data = df_du_rbp %>% 
      filter(covariate == 'Celf5', is_significant) %>%
      mutate(highlight = gene %in% c('Clta', 'Myl6', 'Nnat', 'Gnas', 'Pkm', 'S100a16', 'Cdc42', 'Cnih1', 'Caly')),
    aes(label = label, color = highlight),
    fontface = 'italic',
    force = 5, force_pull = 0.1
  ) + 
  labs(y = 'Conditional HSIC p-value (CBS1)', x = 'Conditional HSIC p-value (CBS2)', 
       color = 'Has CELF binding', title = 'RBP-SVS validation with brain CLIP') + 
  scale_color_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'black')) +
  scale_x_continuous(trans = c("log10", "reverse")) + 
  scale_y_continuous(trans = c("log10", "reverse")) + 
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'), 
    legend.text = element_text(size = 12, family = 'Arial'),
    legend.position = 'inside',
    legend.background = element_blank(),
    legend.position.inside = c(0.25, 0.8)
  )
m10.2

m10 <- cowplot::plot_grid(m10.1, m10.2, nrow = 1, align = 'h', rel_widths = c(1, 1))
m10
ggsave(sprintf("%s/figures/cbs/main_du_rbfox3_celf5_pval.png", res_dir), m10, width = 8, height = 4)
ggsave(sprintf("%s/figures/cbs/main_du_rbfox3_celf5_pval.pdf", res_dir), m10, width = 7.5, height = 4)


## Fig 2j and 2i: Example expression and usage
df_rbp_expr <- read_csv(sprintf("%s/figures/cbs/source_data/rbp_expr.csv", res_dir))
df_svs_ratio <- read_csv(sprintf("%s/figures/cbs/source_data/svs_ratio.csv", res_dir))

# spatial expression of selected RBP
m12_list <- c()
for (rbp in c('Qk', 'Rbfox3', 'Celf5', 'Arpp21')){
  p <- df_rbp_expr %>% filter(
    gene == rbp, layer == 'log1p', array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = expression)) + 
    geom_point(size = 0.2) + 
    labs(color = '', title = rbp) + 
    scale_color_gradient2(low = 'white', high = ifelse(
      rbp %in% c('Qk', 'Rbfox3'),
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
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  m12_list <- c(m12_list, list(p))
}
m12 <- cowplot::plot_grid(plotlist = m12_list, nrow = 2, align = 'hv')
m12

ggsave(sprintf("%s/figures/cbs/main_rbp_expr.png", res_dir), m12, width = 4.5, height = 4)
ggsave(sprintf("%s/figures/cbs/main_rbp_expr.pdf", res_dir), m12, width = 4, height = 4)


### cbs SVS-RBP pair genomic coordinates
source("~/Projects/SPLISOSM_paper/scripts/sit_data_analysis/visualization/plotgardener_isoform_structure.R")

## Load necessary data for visualization
res_dir <- "~/Projects/SPLISOSM_paper/results/sit_nar_23"
fimo_dir <- sprintf('%s/transcripts/cbs/fimo.out/', res_dir)

# Load cbs SVS transcript GTF
svs_gtf <- import.gff(sprintf('%s/transcripts/cbs/cbs_svs.txid.gtf', res_dir))
svs_gtf$gene_id <- gsub("\\..*", "", svs_gtf$gene_id)

# Find and load FIMO motif scanning results
fimo_rbp_dirs <- list.files(fimo_dir, pattern = 'cbs_svs_', full.names = FALSE)
rbs_fimo_gtf <- lapply(fimo_rbp_dirs, function(rbp_dir){
  rbp_gtf <- import.gff(sprintf('%s/%s/fimo.gff', fimo_dir, rbp_dir))
  rbp_gtf$RBP_Name <- toupper(gsub('cbs_svs_', '', rbp_dir))
  rbp_gtf$RBP_Name <- paste(rbp_gtf$RBP_Name, '(motif)')
  return(rbp_gtf)
}) %>% do.call(c, .)

# Load POSTAR RBP binding sites
rbs_clip_gtf <- import.gff('~/reference/POSTAR3/mouse.gtf')
rbs_clip_gtf$RBP_Name <- paste(rbs_clip_gtf$gene_id, '(CLIP)')

# Combine all RBP binding sites
rbs_all_gtf <- c(rbs_fimo_gtf, rbs_clip_gtf)

## Fig 3l: Gnas 
gene_to_plot <- 'Gnas'
df <- df_svs_ratio %>% filter(
  gene == gene_to_plot, 
  transcript_name %in% c('Gnas-206', 'Gnas-208', 'Gnas-221'),
  layer == 'ratios_smoothed', array_row > 15
) %>% dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>%
  arrange(transcript_name)
iso_to_plot <- setNames(df$transcript_id, df$transcript_name)[c('Gnas-208', 'Gnas-206', 'Gnas-221')]


# transcripts and RBP binding sites
png(
  sprintf("%s/figures/cbs/iso_struct_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 2, units = "in", res = 300
)
# GNAS-L (~ Gnas-208) is a phenotypic MDS driver and encodes a hyperactive Gαs protein
plot_transcripts(
  iso_to_plot,
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF2 (CLIP)', 'CELF4 (CLIP)'),
  # rbp_names = c('RBFOX2 (CLIP)', 'CELF2 (CLIP)', 'QK (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)
plotRect(
  x = 3.75, y = 0.2, width = 0.15, height = 0.8,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'Exon3 inclusion drives the MDS phenotype in human HSPCs', fontsize = 8,
  x = 3.5 , y = 0.13,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial isoform ratio
p_list <- c()
for (i in seq_along(iso_to_plot)){
  p <- df_svs_ratio %>% filter(
    transcript_id == iso_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.8) +
    labs(color = '', title = names(iso_to_plot)[i]) + 
    # scale_color_viridis_c(option = 'H') +
    scale_color_distiller(palette = 'Reds', direction = 1) +
    # scale_color_gradient(
    #   low = 'white', high = '#8B0000'
    # ) +
    theme_void() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      strip.text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 12, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  # if (names(iso_to_plot)[i] == 'Gnas-208'){
  #   p <- p + scale_color_gradient2(
  #     low = '#8B0000', mid = 'white', high = '#5B3F8D', #'#0D98BA' 
  #     # limits = c(0.65, 1), oob = scales::squish,
  #     midpoint = 0.75
  #   )
  # } else if(names(iso_to_plot)[i] == 'Gnas-206'){
  #   p <- p + scale_color_gradient(low = 'white', high = '#8B0000')
  # } else {
  #   p <- p + scale_color_gradient(low = 'white', high = '#0066CC')
  # }
  
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p

ggsave(
  sprintf("%s/figures/cbs/ratio_%s.png", res_dir, gene_to_plot), 
  p, width = 9, height = 3
)
ggsave(
  sprintf("%s/figures/cbs/ratio_%s.pdf", res_dir, gene_to_plot), 
  p, width = 9, height = 3
)

## ========
plot_transcripts(
  svs_gtf[svs_gtf$gene_name == 'S100a16']$transcript_id %>% unique(),
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF2 (CLIP)', 'CELF4 (CLIP)', 
                'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  rbs_gtf = rbs_all_gtf,
  fill = '#8B0000'
)
pageGuideHide()