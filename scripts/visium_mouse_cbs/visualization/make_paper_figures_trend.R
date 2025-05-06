library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(universalmotif)
library(rtracklayer)

extrafont::loadfonts()

data_dir <- "~/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir <- "~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
date <- "1119"

### TREND main figures
## make figure directory if doesn't exist
dir.create(sprintf("%s/figures", res_dir), recursive = TRUE, showWarnings = FALSE)

## load SV data
# load short-reads-based TREND results
df_sv_trend <- read_csv(sprintf("%s/sv_results/cbs_peak_sv_combined_%s.csv", res_dir, date)) %>%
  mutate(
    is_trend_svs = `padj_hsic-ir` < 0.01,
    is_trend_sve = `padj_hsic-gc` < 0.01,
    is_trend_expressed = count_avg > 0.5
  )

# load and combine ONT results
df_sv_ont1 <- read_csv("~/Projects/SPLISOSM_paper/results/sit_nar_23/sv_results/cbs1_sv_combined_1107.csv") %>%
  mutate(
    is_ont_svs = `padj_hsic-ir` < 0.05,
    is_ont_sve = `padj_hsic-gc` < 0.05,
    is_ont_expressed = count_avg > 0.5
  )
df_sv_ont2 <- read_csv("~/Projects/SPLISOSM_paper/results/sit_nar_23/sv_results/cbs2_sv_combined_1107.csv") %>%
  mutate(
    is_ont_svs = `padj_hsic-ir` < 0.05,
    is_ont_sve = `padj_hsic-gc` < 0.05,
    is_ont_expressed = count_avg > 0.5
  )
df_sv_ont <- inner_join(df_sv_ont1, df_sv_ont2, by = 'gene', suffix = c('.1', '.2')) %>%
  mutate(
    is_ont_svs = is_ont_svs.1 & is_ont_svs.2,
    is_ont_sve = is_ont_sve.1 & is_ont_sve.2,
    is_ont_expressed = is_ont_expressed.1 & is_ont_expressed.2
  )

# load SlideseqV2 hippocampus results
df_sv_slideseq <- read_csv('~/Projects/SPLISOSM_paper/results/hippo_slideseqv2/slideseqv2_hippocampus_sv_combined_0307.csv') %>%
  mutate(
    is_slideseq_svs = `padj_hsic-ir` < 0.01,
    is_slideseq_sve = `padj_hsic-gc` < 0.01,
    is_slideseq_expressed = count_avg > 0.5
  )

## Fig 4B: HSIC-GC vs HSIC-IR
df_ratio_comp_stats = data.frame(
  x = c(1e-30, 1e-30, 1e-65, 1e-65), 
  y = c(1e-220, 1e-255, 1e-220, 1e-255), 
  size = table(
    df_sv_trend$is_trend_svs, 
    df_sv_trend$is_trend_sve
  ) %>% as.vector() # FF, FT, TF, TT
)

m1 <- ggplot(df_sv_trend, aes(x = `pvalue_hsic-gc` + 1e-300, y = `pvalue_hsic-ir` + 1e-300)) +
  geom_point(color = 'gray', alpha = 0.2) +
  geom_hline(yintercept = 0.05, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.05, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_point(
    data = (df_sv_trend %>% filter(is_trend_svs & is_trend_expressed)),
    aes(color = perplexity)
  ) +
  geom_text_repel(
    data = (df_sv_trend %>% filter(is_trend_expressed) %>% top_n(n = -20, wt = `pvalue_hsic-ir`)),
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
    geom = 'text', x = 1e-70, y = 1e-280, label = 'Number of genes', size = 12/.pt
  ) +
  labs(x = 'HSIC-GC p-value', y = 'HSIC-IR p-value', color = 'Perplexity') + 
  scale_color_viridis_b() + 
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = 1) + 
  scale_x_continuous(trans = c("log10", "reverse"), breaks = c(1e-05, 1e-50, 1e-100)) + 
  scale_y_continuous(trans = c("log10", "reverse"), breaks = c(1e-05, 1e-50, 1e-100)) + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position.inside = c(0.8, 0.7)
  ) + 
  guides(fill = "none")
m1

ggsave(sprintf("%s/figures/main_gc_vs_ir.png", res_dir), m1, width = 5.5, height = 4)
ggsave(sprintf("%s/figures/main_gc_vs_ir.pdf", res_dir), m1, width = 5.5, height = 4)

## Fig S3A: HSIC-IR Visium CBS vs SlideseqV2 hippocampus
df <- inner_join(df_sv_trend, df_sv_slideseq, by = 'gene', suffix = c('', '.ss')) 
table(df$is_slideseq_svs, df$is_trend_svs)
# FALSE TRUE
# FALSE   228  196
# TRUE     10   42

s1 <- ggplot(df, aes(x = `pvalue_hsic-ir.ss` + 1e-300, y = `pvalue_hsic-ir` + 1e-300)) +
  geom_point(color = 'gray', alpha = 0.2) +
  geom_hline(yintercept = 0.05, color = 'black', linetype = 2) + 
  geom_vline(xintercept = 0.05, color = 'black', linetype = 2) + 
  geom_abline(slope = 1, color = 'lightgray') + 
  geom_point(
    data = (df %>% filter(is_trend_svs & is_slideseq_svs)),
    aes(color = perplexity)
  ) +
  geom_text_repel(
    data = (df %>% filter(is_slideseq_svs & is_trend_svs) %>% top_n(n = -18, wt = `pvalue_hsic-ir.ss`)),
    aes(label = gene),
    fontface = "italic"
  ) + 
  labs(x = 'P-value (Hippocampus SlideseqV2)', y = 'P-value (CBS Visium)', 
       title = sprintf('HSIC-IR test agreement\n(Spearman R = %.2f)', cor(
         df$`pvalue_hsic-ir.ss`, df$`pvalue_hsic-ir`, method = 'spearman'
       )),
       color = 'Perplexity') + 
  scale_color_viridis_b() + 
  scale_fill_distiller(type = 'seq', palette = 'Reds', direction = 1) + 
  scale_x_continuous(trans = c("log10", "reverse"), breaks = c(1e-05, 1e-50, 1e-100)) + 
  scale_y_continuous(trans = c("log10", "reverse"), breaks = c(1e-05, 1e-50, 1e-100)) + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial', hjust = 0.5),
    legend.position = 'none',
    legend.position.inside = c(0.8, 0.7)
  ) + 
  guides(fill = "none")  

ggsave(sprintf("%s/figures/sup_hippo_slideseqv2.png", res_dir), s1, width = 4.5, height = 4)

## SupFig: Venn diagrams of SV genes
library(ggvenn)
venn_svs <- list(
  `Short-reads SVS` = df_sv_trend$gene[df_sv_trend$is_trend_svs] %>% unique(),
  `ONT SVS (rep1)` = df_sv_ont1$gene[df_sv_ont1$is_ont_svs] %>% unique(),
  `ONT SVS (rep2)` = df_sv_ont2$gene[df_sv_ont2$is_ont_svs] %>% unique()
)
venn_sve <- list(
  `Short-reads SVE` = df_sv_trend$gene[df_sv_trend$is_trend_sve] %>% unique(),
  `ONT SVE (rep1)` = df_sv_ont1$gene[df_sv_ont1$is_ont_sve] %>% unique(),
  `ONT SVE (rep2)` = df_sv_ont2$gene[df_sv_ont2$is_ont_sve] %>% unique()
)

m2.1 <- ggvenn(
  venn_svs, fill_color = c('#6B5B95', '#88D8C0', '#009999'), 
  show_percentage = TRUE, digits = 0,
  stroke_size = 0.5, text_size = 4, set_name_size = 4
)
m2.2 <- ggvenn(
  venn_sve, fill_color = c('#6B5B95', '#88D8C0', '#009999'), 
  show_percentage = TRUE, digits = 0,
  stroke_size = 0.5, text_size = 4, set_name_size = 4
)

m2 <- cowplot::plot_grid(m2.1, m2.2, nrow = 1, align = 'h', rel_widths = c(1, 1))
m2

ggsave(sprintf("%s/figures/main_venn.png", res_dir), m2, width = 6, height = 4)


## Fig 4C: Functional enrichment analysis
# load enrichment results
df_kegg <- read_csv(sprintf("%s/figures/source_data/kegg_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced (SVENS)', 
                                'Spatially variably spliced (SVS)'))
  ) %>%
  # keep the top 5 terms with the largest precision difference
  dplyr::slice((n()-15):n())

df_reac <- read_csv(sprintf("%s/figures/source_data/reac_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced (SVENS)', 
                                'Spatially variably spliced (SVS)'))
  ) %>%
  # keep the top 5 terms with the largest precision difference
  dplyr::slice((n()-15):n())

m3.1 <- ggplot(df_kegg, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
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
m3.2 <- ggplot(df_reac, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) +
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
    legend.position.inside = c(-1.2, 1.15),
    legend.background = element_blank()
  )
m3 <- cowplot::plot_grid(m3.1, m3.2, nrow = 2, align = 'hv', rel_heights = c(1,1))
m3

ggsave(sprintf("%s/figures/main_enrich_new.png", res_dir), m3, width = 7.5, height = 4)
ggsave(sprintf("%s/figures/main_enrich_new.pdf", res_dir), m3, width = 7, height = 4.2)


## SupFig: Visualize detailed KEGG enrichment terms
library(clusterProfiler)
eg <- bitr(
  df_sv_trend$gene[df_sv_trend$is_trend_svs] %>% unique(),
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
for (kegg_id in c('mmu04724', 'mmu04261', 'mmu04015', 'mmu04020', 
                  'mmu04024', 'mmu04010', 'mmu04510', 'mmu04810',
                  'mmu04970', 'mmu04611')) {
  glist <- kk@result[kk@result$ID == kegg_id, 'geneID'] %>% str_split_1('/') %>% bitr(
    fromType = 'ENTREZID',
    toType = 'SYMBOL',
    OrgDb = org.Mm.eg.db
  )
  kegg_glist[[kegg_id]] <- df_sv_trend %>% filter(
    gene %in% glist$SYMBOL,
    is_trend_svs
  ) %>% pull(gene)
}

# visualize the gene list
df <- df_sv_trend %>% 
  mutate(rank = rank(`pvalue_hsic-ir`)) %>%
  filter(is_trend_svs)
p_list <- list()
for (kegg_id in c('mmu04724', 'mmu04261', 'mmu04015', 'mmu04020', 
                  'mmu04024', 'mmu04010', 'mmu04510', 'mmu04810',
                  'mmu04970', 'mmu04611')) {
  p <- ggplot(df, aes(x = rank, y = `pvalue_hsic-ir`)) + 
    geom_point(alpha = 0.5, color = 'gray') + 
    scale_y_continuous(transform = c('log10', 'reverse')) +
    geom_point(
      data = df %>% filter(gene %in% kegg_glist[[kegg_id]]),
      color = 'black'
    ) +
    geom_text_repel(
      data = df %>% filter(gene %in% kegg_glist[[kegg_id]]),
      aes(label = gene, x = rank, y = `pvalue_hsic-ir`),
      force = 20, max.overlaps = 30, nudge_x = 0, nudge_y = 50,
      fontface = 'italic'
    ) +
    labs(x = 'Rank', y = 'HSIC-IR p-value', title = stringr::str_wrap(
      str_split_1(kk@result[kegg_id, 'Description'], ' - Mus')[1], 50
    )) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial')
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, ncol = 4, align = 'h')
p

ggsave(sprintf("%s/figures/sup_kegg_genes.png", res_dir), p, width = 15, height = 9)


## Fig 4F: Clustering results
df_clu <- read_csv(sprintf("%s/figures/source_data/clu_res.csv", res_dir)) %>%
  mutate(clu_trend_sve_anno = factor(clu_trend_sve))

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

for (clu_col in c('clu_trend_svs', 'clu_visium_rbp_sve')){
  alt2ref_mapping <- match_labels(df_clu$clu_trend_sve_anno, df_clu[[clu_col]])
  df_clu[[paste0(clu_col, '_anno')]] <- factor(
    df_clu[[clu_col]],
    levels = alt2ref_mapping[,1],
    labels = alt2ref_mapping[,2]
  )
}

# filter out outlier spots
df_clu <- df_clu %>% filter(array_col > 15)

# set color palette
clu_palette <- setNames(
  ggsci::pal_npg(palette = c("nrc"), alpha = 1)(10),
  df_clu$clu_trend_sve_anno %>% unique() %>% sort()
)
clu_palette <- setNames(
  c("#7E6148FF","#00A087FF","#DC0000FF","#4DBBD5FF", "#3C5488FF", "#F39B7FFF", "#8491B4FF", "#91D1C2FF"),
  df_clu$clu_trend_sve_anno %>% unique() %>% sort()
)

# "#E64B35FF" "#4DBBD5FF" "#00A087FF" "#3C5488FF" "#F39B7FFF" "#8491B4FF" "#91D1C2FF" "#DC0000FF" "#7E6148FF"

# TREND SVE and SVS clustering results
p1 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_trend_sve_anno), size = 0.3) + 
  coord_fixed(ratio = 1.5) +
  scale_color_manual(values = clu_palette) +
  labs(color = '', title = 'SVE (gene expression)') + 
  theme_void() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
p2 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_trend_svs_anno), size = 0.3) + 
  coord_fixed(ratio = 1.5) +
  scale_color_manual(values = clu_palette) +
  labs(color = '', title = 'SVS (isoform ratio)') + 
  theme_void() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )
p3 <- ggplot(df_clu, aes(x = array_col, y = - array_row)) + 
  geom_point(aes(color = clu_visium_rbp_sve_anno), size = 0.3) + 
  coord_fixed(ratio = 1.5) +
  scale_color_manual(values = clu_palette) +
  labs(color = '', title = 'SVE RBP (expression)') + 
  theme_void() + 
  guides(color = guide_legend(override.aes = list(size = 2))) + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    legend.text = element_text(size = 12, family = 'Arial'),
    panel.background = element_rect(color = 'black', linewidth = 1),
    legend.position = 'none',
    plot.title = element_text(hjust = 0.5)
  )

p <- cowplot::plot_grid(p1, p2, p3, nrow = 1, align = 'hv', rel_widths = c(1, 1, 1))
p

ggsave(sprintf("%s/figures/sup_clustering.png", res_dir), p, width = 6, height = 2.5)
ggsave(sprintf("%s/figures/sup_clustering.pdf", res_dir), p, width = 6, height = 2.4)

## Fig 4C: Xstreme de novo motif discovery
ref_motifs <- read_meme("~/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme")
xstreme_motifs <- read_meme(sprintf("%s/events/xstreme.out/combined.meme", res_dir))

# the poly-A motifs
m4.1 <- view_motifs((ref_motifs %>% filter_motifs(name = 'M042_0.6'))[[1]]) + 
  labs(title = 'M042_0.6 (Pabpc1/1l/2/4/6)\nSEA p-val=8.37e-7') +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
m4.2 <- view_motifs((ref_motifs %>% filter_motifs(name = 'M062_0.6'))[[1]]) +
  labs(title = 'M062_0.6 (Pabpc1/1l/2/4/6)\nSEA pval=5.71e-5') +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

# the de novo poly-A motif
mt1.denovo <- xstreme_motifs[[3]]
mt1.denovo@name <- "De novo motif 1\nSEA p-val=5.27e-9"
mt1.ref1 <- (ref_motifs %>% filter_motifs(name = 'M147_0.6'))[[1]]
mt1.ref1@name <- "M147_0.6 (Cnot4)"
m4.3 <- view_motifs(c(mt1.denovo, mt1.ref1), min.overlap = 5) + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
m4.3

# the de novo GAGGA motif
mt2.denovo <- xstreme_motifs[[7]]
mt2.denovo@name <- "De novo motif 2\nSEA p-val=6.26e-6"
mt2.ref1 <- (ref_motifs %>% filter_motifs(name = 'M154_0.6'))[[1]]
mt2.ref1@name <- "M154_0.6 (Srsf1/9)"
m4.4 <- view_motifs(c(mt2.denovo, mt2.ref1), min.overlap = 5) + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
m4.4

m4 <- cowplot::plot_grid(m4.1, m4.2, m4.3, m4.4, nrow = 2, rel_heights = c(0.6, 1))
m4

ggsave(sprintf("%s/figures/main_xstreme_motif.png", res_dir), m4, width = 6, height = 3.5)

# the de novo poly-A motif
m4.1 <- view_motifs(xstreme_motifs[[1]]) +
  labs(title = "De novo motif 1\nSEA p-val=5.84e-20") +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
# the de novo GAGGA motif
m4.2 <- view_motifs(xstreme_motifs[[7]]) + 
  labs(title = "De novo motif 2\nSEA p-val=6.26e-6") +
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
m4 <- cowplot::plot_grid(m4.1, m4.2, nrow = 1, align = 'hv', rel_widths = c(1, 1))
m4

ggsave(sprintf("%s/figures/main_xstreme_motif.pdf", res_dir), m4, width = 8, height = 1.5)


## Fig 4C: Peak annotation
# load sierra annotation outputs
# peak_annot <- read.table(sprintf("%s/peak_no_cutoff.annot.txt", data_dir), sep = '\t') %>%
#   rownames_to_column('peak_id')
peak_annot <- read.table(sprintf("%s/peak_no_cutoff.annot.overlap.txt", data_dir), sep = '\t') %>%
  rownames_to_column('peak_id')
# subset peaks for SVS and SVE
svs_peak <- import.bed(sprintf("%s/events/cbs_svs.peak.bed", res_dir))
svens_peak <- import.bed(sprintf("%s/events/cbs_svens.peak.bed", res_dir))
svs_peak_annot <- peak_annot %>% filter(peak_id %in% svs_peak$name) %>%
  # remove duplicated peaks by seqnames, start and end
  distinct(seqnames, start, end, .keep_all = TRUE)
svens_peak_annot <- peak_annot %>% filter(peak_id %in% svens_peak$name) %>%
  # remove duplicated peaks by seqnames, start and end
  distinct(seqnames, start, end, .keep_all = TRUE)

# count number of AS events in each group
n_svs <- nrow(svs_peak_annot)
n_sve <- nrow(svens_peak_annot)

df_event <- c()
for (event in c('UTR3','UTR5', 'intron', 'exon', 'CDS', 'Junctions', 'Alt-Exon')){
  if (event == 'Alt-Exon'){
    svs_counts <- sum(svs_peak_annot$exon == 'YES' & svs_peak_annot$intron == 'YES')
    sve_counts <- sum(svens_peak_annot$exon == 'YES' & svens_peak_annot$intron == 'YES')
  } else {
    svs_counts <- sum(svs_peak_annot[[event]] == ifelse(
      event == 'Junctions', 'across-junctions', 'YES'
    ))
    svens_counts <- sum(svens_peak_annot[[event]] == ifelse(
      event == 'Junctions', 'across-junctions', 'YES'
    ))
  }

  # test for differences between groups
  pval <- chisq.test(
    matrix(c(svs_counts, svens_counts, 
             n_svs - svs_counts, n_sve - svens_counts), nrow = 2)
  )$p.value
  
  # save results
  df_event <- rbind(
    df_event,
    data.frame(
      event = event,
      SVS = svs_counts,
      SVENS = svens_counts,
      pval = pval
    )
  )
}

df_event <- df_event %>%
  mutate(
    event = fct_recode(event, 'Intron' = 'intron', 'Exon' = 'exon'),
    event = factor(event, levels = c('Junctions', 'Exon', 'UTR3', 'Intron', 'Alt-Exon', 'CDS', 'UTR5')),
    label = ifelse(pval < 0.05, '*', 'ns'),
    label = ifelse(pval < 0.01, '**', label),
    label = ifelse(pval < 0.001, '***', label),
    label = ifelse(pval < 0.0001, '****', label)
  ) %>%
  pivot_longer(
    cols = c(SVS, SVENS), names_to = 'dataset', values_to = 'count'
  ) %>%
  mutate(
    prop = count / ifelse(dataset == 'SVS', n_svs, n_sve),
    dataset = factor(dataset, levels = c('SVS', 'SVENS'))  
  )


# plot event type frequency
m5.1 <- ggplot(df_event %>% filter(event == 'Junctions'), aes(x = event, y = prop)) + 
  geom_bar(aes(fill = dataset), position="dodge", stat="identity", width = 0.75) + 
  geom_text(aes(label = count, group = dataset, y = prop), position = position_dodge(width = 0.75)) +
  geom_text(data = df_event %>% filter(event == 'Junctions', dataset == 'SVS'),
            mapping = aes(x = event, y = 0.2, label = label)) +
  labs(y = 'Proportion of junction peaks', x = '', fill = 'Gene set') +
  scale_fill_manual(values = c('SVS' = '#FF9F1C', 'SVENS' = '#2EC4B6')) +
  theme_classic() + 
  theme(
    legend.position = 'none',
    text = element_text(size = 12, family='Arial'),
    axis.text = element_text(size = 12, family='Arial'),
    plot.title = element_text(hjust = 0.5, size = 12, family='Arial')
  )

m5.2 <- ggplot(df_event %>% filter(! event %in% c('Junctions', 'Intron')), 
               aes(x = event, y = prop)) + 
  geom_bar(aes(fill = dataset), position="dodge", stat="identity", width = 0.75) + 
  geom_text(aes(label = count, group = dataset, y = prop), position = position_dodge(width = 0.75)) +
  geom_text(data = df_event %>% filter(! event %in% c('Junctions', 'Intron'), dataset == 'SVS'),
            mapping = aes(x = event, y = 1, label = label)) +
  # coord_flip() + 
  labs(y = 'Proportion of peaks that overlap with', x = '', fill = 'Gene set') +
  scale_fill_manual(values = c('SVS' = '#FF9F1C', 'SVENS' = '#2EC4B6')) +
  theme_classic() + 
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.7, 0.7),
    text = element_text(size = 12, family='Arial'),
    axis.text = element_text(size = 12, family='Arial'),
    plot.title = element_text(hjust = 0.5, size = 12, family='Arial')
  )

m5 <- cowplot::plot_grid(m5.1, m5.2, nrow = 1, align = 'hv', rel_widths = c(1, 2.5))
m5

ggsave(sprintf("%s/figures/main_peak_annot.png", res_dir), m5, width = 6, height = 4)
ggsave(sprintf("%s/figures/main_peak_annot.pdf", res_dir), m5, width = 6, height = 3)

## Fig S3D: DU RBP results
rbp_with_motifs <- read_table("~/reference/cisbp-rna/mouse_pm/rbp_with_known_cleaned_motif.txt", col_names = 'rbp')
rbp_with_clip <- read_table("~/reference/POSTAR3/mouse.rbp_with_clip.txt", col_names = 'rbp')
df_du_rbp <- read_csv(sprintf("%s/figures/source_data/du_res_annot.csv", res_dir)) %>%
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
df_du_rbp %>% filter(is_significant) %>% nrow()

# rank RBP by number of significant associations
df_rbp_rank <- df_du_rbp %>%
  summarise(count = sum(is_significant), .by = c(covariate, has_clip, has_motif)) %>%
  mutate(
    rank = rank(-count, ties.method = 'first'),
    covariate = factor(covariate, levels = covariate[order(count)])
  )

m6.1 <- ggplot(df_rbp_rank, aes(x = rank, y = count)) + 
  geom_line() + 
  geom_point(color = '#668586') +
  geom_label_repel(
    data = df_rbp_rank %>% filter((!has_clip) & (!has_motif)) %>% slice_max(count, n = 15, with_ties = FALSE),
    aes(label = covariate), color = '#668586',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = -0.5, nudge_y = -20
  ) +
  geom_point(
    data = df_rbp_rank %>% filter(has_clip | has_motif) %>% slice_max(count, n = 15, with_ties = FALSE),
    color = '#8B0000',
  ) + 
  geom_label_repel(
    data = df_rbp_rank %>% filter(has_clip | has_motif) %>% slice_max(count, n = 15, with_ties = FALSE),
    aes(label = covariate), color = '#8B0000',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = 1, nudge_y = 20
  ) +
  scale_x_log10() + 
  labs(x = 'Rank', y = 'Number of significant associations', color = 'In CLIPdb',
       title = 'Less-studied vs CISBP-RNA RBPs') + 
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
  )

# Plot ranks for perturbed APA-related RBPs
df_apa_rbp <- read_csv(sprintf("%s/APA_related_RBP.csv", data_dir))
motif_rbp <- c(
  "Pabpc1", "Pabpc1l","Pabpc2", "Pabpc4", "Pabpc6", "Sart3", # poly-A
  "Zfp36", "Zfp36l1", "Zfp36l2", "Zfp36l3", # M269 AAAG
  "Cnot4", # M147 GACAGA
  "Srsf1", "Srsf9" # M154 GAGGA
)

m6.2 <- ggplot(df_rbp_rank, aes(x = rank, y = count)) + 
  geom_line() + 
  geom_point(color = '#668586') +
  geom_point(
    data = df_rbp_rank %>% filter(
      covariate %in% df_apa_rbp$ortholog_name,
      count > 0
    ),
    color = '#641E16',
  ) + 
  geom_label_repel(
    data = df_rbp_rank %>% filter(
      covariate %in% df_apa_rbp$ortholog_name,
      count > 0
    ),
    aes(label = covariate), color = '#641E16',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = 2, nudge_y = 30
  ) +
  geom_point(
    data = df_rbp_rank %>% filter(
      covariate %in% motif_rbp,
      count > 0
    ),
    color = '#0066CC',
  ) + 
  geom_label_repel(
    data = df_rbp_rank %>% filter(
      covariate %in% motif_rbp,
      count > 0
    ),
    aes(label = covariate), color = '#0066CC',
    fontface = 'italic',
    force = 20, force_pull = 0.1, max.overlaps = 20, nudge_x = -1, nudge_y = -10
  ) +
  scale_x_log10() + 
  labs(x = 'Rank', y = 'Number of significant associations', 
       title = 'Motif-based vs APA-related RBPs') +
  theme_classic() + 
  theme(
    text = element_text(family='Arial', size = 14),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
  )

m6 <- cowplot::plot_grid(m6.1, m6.2, nrow = 1, align = 'h', rel_widths = c(1, 1))
m6

ggsave(sprintf("%s/figures/main_du_rbpr_rank.png", res_dir), m6, width = 8, height = 4)


## Fig 4I: Example RBP expression
df_rbp_expr <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_dir)) %>%
  filter(array_col > 15)

# spatial expression of selected RBP
rbp_list <- c(
  'Arpp21', 'Celf5', 
  'Rbfox3', 'Qk', 'Celf2',
  'Zfp36l1', 'Pabpc4',
  'Cstf3', 'Tent4b', 'Pabpn1'
)
color_list <- c(
  'Arpp21' = '#668586', 'Celf5' = '#668586',
  'Rbfox3' = '#8B0000', 'Qk' = '#8B0000', 'Celf2' = '#8B0000', 
  'Zfp36l1' = '#0066CC', 'Pabpc4' = '#0066CC', 
  'Cstf3' = '#641E16', 'Tent4b' = '#641E16', 'Pabpn1' = '#641E16'
)
p_list <- c()
for (rbp in rbp_list){
  p <- df_rbp_expr %>% filter(
    gene == rbp, layer == 'log1p', array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = expression)) + 
    geom_point(size = 0.5) + 
    labs(color = '', title = rbp) + 
    scale_color_gradient2(low = 'white', high = color_list[rbp]) +
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial', face = 'italic'),
      panel.background = element_rect(color = 'black', linewidth = 1),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
m7 <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
m7

ggsave(sprintf("%s/figures/main_rbp_expr.png", res_dir), m7, width = 12, height = 4)


## Fig 4I: Plot selected SVS RBPs and their regulators
rbp_to_highlight <- union(
  union(rbp_with_clip$rbp, rbp_with_motifs$rbp),
  c('Celf5', 'Arpp21', 'Adar', 'Adarb1', 'Adarb2')
) %>% sort()

# Arpp21
m8.1 <- df_du_rbp %>% 
  filter(gene == 'Arpp21', covariate %in% rbp_to_highlight) %>%
  slice_min(pvalue_glmm, n = 5) %>%
  mutate(covariate = fct_reorder(covariate, pvalue_glmm, .desc = F)) %>%
  ggplot(aes(x = covariate, y = pvalue_glmm)) +
  geom_bar(stat = 'identity', aes(fill = covariate == 'Arpp21')) +
  geom_hline(yintercept = 0.01, linetype = 'dashed', color = 'black') +
  labs(title = 'Arpp21', x = '', y = 'DU test pvalue (GLMM)', fill = '') + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
  scale_y_continuous(trans = c('log10', 'reverse')) +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial',face = 'italic'),
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
    legend.position = 'none'
  )
# Celf2
m8.2 <- df_du_rbp %>% 
  filter(gene == 'Celf2', covariate %in% rbp_to_highlight) %>%
  slice_min(pvalue_glmm, n = 5) %>%
  mutate(covariate = fct_reorder(covariate, pvalue_glmm, .desc = F)) %>%
  ggplot(aes(x = covariate, y = pvalue_glmm)) +
  geom_bar(stat = 'identity', aes(fill = covariate == 'Celf5')) +
  geom_hline(yintercept = 0.01, linetype = 'dashed', color = 'black') +
  labs(title = 'Celf2', x = '', y = 'DU test pvalue (GLMM)', fill = '') + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
  scale_y_continuous(trans = c('log10', 'reverse')) +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial',face = 'italic'),
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
    legend.position = 'none'
  )
# Pcbp2
m8.3 <- df_du_rbp %>%
  filter(gene == 'Pcbp2') %>%
  slice_min(pvalue_glmm, n = 5) %>%
  mutate(covariate = fct_reorder(covariate, pvalue_glmm, .desc = F)) %>%
  ggplot(aes(x = covariate, y = pvalue_glmm)) +
  geom_bar(stat = 'identity', aes(fill = covariate == 'Pcbp2')) +
  geom_hline(yintercept = 0.01, linetype = 'dashed', color = 'black') +
  labs(title = 'Pcbp2', x = '', y = 'DU test pvalue (GLMM)', fill = '') + 
  scale_fill_manual(values = c('TRUE' = '#5B3F8D', 'FALSE' = 'gray')) +
  scale_y_continuous(trans = c('log10', 'reverse')) +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial',face = 'italic'),
    axis.text.x = element_text(angle = 45, hjust = 1, face = 'italic'),
    legend.position = 'none'
  )
m8 <- cowplot::plot_grid(m8.1, m8.2, m8.3, nrow = 1, align = 'h', rel_widths = c(1, 1, 1))
m8

ggsave(sprintf("%s/figures/du_arpp21_celf2_pcbp2.png", res_dir), m8, width = 6, height = 3)
ggsave(sprintf("%s/figures/du_arpp21_celf2_pcbp2.pdf", res_dir), m8, width = 6, height = 3)


## Fig S3B: validating neuronal APA genes
# Load cbs SVS peak bed
svs_bed <- import.bed(sprintf('%s/events/cbs_svs.peak.bed', res_dir)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1)
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order))  # Create the new name

# Load isoform ratio
df_svs_ratio <- read_csv(sprintf("%s/figures/source_data/svs_ratio.csv", res_dir)) %>%
  left_join(svs_bed %>% dplyr::select(name, peak_name), by = c('isoform' = 'name'))

## Fig 4E: Map4
gene_to_plot <- 'Map4'
peak_to_plot <- c('Map4-Event1', 'Map4-Event6')
peak_name_list <- c('Map4 - neuronal', 'Map4 - glial')
p_list <- c()
# spatial usage of selected peaks
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_name_list[i]) + 
    scale_color_distiller(palette = 'Spectral') +
    theme_void() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      strip.text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 12, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p

ggsave(sprintf("%s/figures/svs_cds/ratio_Map4.png", res_dir), p, width = 5, height = 2.5)

# spatial log1p expression of selected peaks
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.2) +
    labs(color = '', title = peak_name_list[i]) + 
    scale_color_distiller(palette = 'Purples', direction = 1) + 
    theme_void() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      strip.text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 12, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/log1p_Map4.pdf", res_dir), p, width = 4, height = 2)

# Klc1
gene_to_plot <- 'Klc1'
peak_to_plot <- c('Klc1-Event1', 'Klc1-Event2')
p_list <- c()
# spatial usage of selected peaks
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') + 
    scale_color_distiller(palette = 'Purples', direction = 1) +
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_Klc1.png", res_dir), p, width = 5, height = 2.5)

# Atp2a2
gene_to_plot <- 'Atp2a2'
peak_to_plot <- c('Atp2a2-Event1', 'Atp2a2-Event3')
p_list <- c()
# spatial usage of selected peaks
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') + 
    scale_color_distiller(palette = 'Purples', direction = 1) +
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_Atp2a2.png", res_dir), p, width = 5, height = 2.5)

# Itsn1
gene_to_plot <- 'Itsn1'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Itsn1-Event2', 'Itsn1-Event3', 'Itsn1-Event4')
p_list <- c()
# spatial usage of selected peaks
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') + 
    scale_color_distiller(palette = 'Blues', direction = 1) +
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_Itsn1.png", res_dir), p, width = 7.5, height = 2.5)

# Cdc42
gene_to_plot <- 'Cdc42'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Cdc42-Event1', 'Cdc42-Event4', 'Cdc42-Event5')
p_list <- c()
# spatial usage of selected peaks
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') + 
    scale_color_distiller(palette = 'Blues', direction = 1) +
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.25, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_Cdc42.png", res_dir), p, width = 7.5, height = 2.5)
