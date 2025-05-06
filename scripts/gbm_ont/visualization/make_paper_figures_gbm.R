library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(universalmotif)

extrafont::loadfonts()

res_dir_ont <- "~/Projects/SPLISOSM_paper/results/gbm_ont/"
res_dir_trend <- "~/Projects/SPLISOSM_paper/results/gbm_visium/"
res_dir_dlpfc <- "~/Projects/SPLISOSM_paper/results/human_dlpfc/"
date_ont <- "0104"
date_vis <- "0310"
date_dlpfc <- "0123"

# load GBM ONT SV results
df_ont_svs <- read.csv(paste0(res_dir_ont, "ont_svs_allsamples_", date_ont, ".csv"))
df_ont_sve <- read.csv(paste0(res_dir_ont, "ont_sve_allsamples_", date_ont, ".csv"))

# load GBM Visium SV results
df_trend_svs <- read.csv(paste0(res_dir_trend, "trend_svs_gbm_", date_vis, ".csv"))
df_trend_sve <- read.csv(paste0(res_dir_trend, "trend_sve_gbm_", date_vis, ".csv"))

# load human DLPFC SV results
df_dlpfc_svs <- read.csv(paste0(res_dir_dlpfc, "trend_svs_dlpfc_", date_dlpfc, ".csv"))
df_dlpfc_sve <- read.csv(paste0(res_dir_dlpfc, "trend_sve_dlpfc_", date_dlpfc, ".csv"))

## GBM main figures
## make figure directory if doesn't exist
dir.create(sprintf("%s/figures", res_dir_ont), recursive = TRUE, showWarnings = FALSE)
# dir.create(sprintf("%s/events", res_dir), recursive = TRUE, showWarnings = FALSE)

## Fig 7B: SVS results statistics
# ONT SVS gene ranking
df <- df_ont_svs %>% mutate(
  rank = row_number(),
  is_trend_svs = gene %in% df_trend_svs$gene,
  gbm_specific = !(gene %in% df_dlpfc_svs$gene)
) 
m1.1 <- ggplot(df, aes(x = rank, y = pvalue_hsic.ir)) +
  geom_point(aes(color = gbm_specific)) + 
  geom_text_repel(
    data = df %>% arrange(pvalue_hsic.ir) %>% head(n = 20), 
    force = 10, force_pull = 1,
    mapping = aes(label = gene, color = gbm_specific), nudge_y = 10, nudge_x = 0.1,
    fontface = "italic"
  ) + 
  labs(x = 'Rank (by occurrence)', y = 'Minimum HSIC-IR p-value', 
       title = 'SVP genes (ONT) in glioma',
       color = 'Disease specific') +
  scale_y_continuous(transform = c('log10', 'reverse')) + 
  scale_x_log10() +
  scale_color_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    legend.position = 'inside',
    legend.position.inside = c(0.8, 0.85),
    legend.background = element_blank()
  )

# TREND SVP gene ranking
df <- df_trend_svs %>% mutate(
  rank = row_number(),
  is_ont_svs = gene %in% df_ont_svs$gene,
  gbm_specific = !(gene %in% df_dlpfc_svs$gene)
)
m1.2 <- ggplot(df, aes(x = rank, y = pvalue_hsic.ir)) +
  geom_point(aes(color = gbm_specific)) + 
  geom_text_repel(
    data = df %>% arrange(pvalue_hsic.ir) %>% head(n = 20), 
    force = 10, force_pull = 1,
    mapping = aes(label = gene, color = gbm_specific), nudge_y = 15, nudge_x = 0.1,
    fontface = "italic"
  )+ 
  labs(x = 'Rank (by occurrence)', y = 'Minimum HSIC-IR p-value', 
       title = 'SVP genes (SR) in glioma',
       color = 'Disease specific') +
  scale_y_continuous(transform = c('log10', 'reverse')) + 
  scale_x_log10() + 
  scale_color_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    legend.position.inside = c(0.8, 0.85)
  )

m1 <- cowplot::plot_grid(m1.1, m1.2, ncol = 2, align = 'v')
m1

ggsave(sprintf("%s/figures/main_svs_pval.png", res_dir_ont), m1, width = 7.5, height = 3.5)
ggsave(sprintf("%s/figures/main_svs_pval.pdf", res_dir_ont), m1, width = 7.5, height = 3.5)


## Fig 7C: SVS gene recurrence
m2.1 <- rbind(
  (df_ont_svs %>% dplyr::select('gene', 'occurrence') %>% mutate(geneset = 'ONT')),
  (df_trend_svs %>% dplyr::select('gene', 'occurrence') %>% mutate(geneset = 'SR'))
) %>%
  mutate(
    group = factor(ifelse(occurrence > 1, 'More than one sample', 'Singleton'),
                   levels = c('Singleton', 'More than one sample')),
    geneset = factor(geneset, levels = c('SR', 'ONT'))
  ) %>%
  dplyr::count(geneset, group) %>%
  ggplot(aes(x = geneset, y = n, fill = group)) +
    geom_bar(stat = 'identity', position = 'fill') +
    scale_fill_manual(values = c('More than one sample' = '#FF9F1C', 'Singleton' = 'gray')) + 
    labs(x = '', y = 'Proportion', fill = '', title = 'SVP genes in glioma') + 
    geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
    coord_flip() + 
    theme_classic() + 
    theme(
      text=element_text(size = 12, family = 'Arial'),
      plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
      axis.title.x = element_blank(),
      legend.position = 'bottom',
      legend.margin = margin(-5, 0, 0, 0)
    )

m2.2 <- df_ont_svs %>%
  # filter(occurrence > 1) %>%
  mutate(
    gbm_specific = !(gene %in% df_dlpfc_svs$gene),
    occurrence = ifelse(occurrence >= 3, '3+', occurrence),
    occurrence = factor(occurrence, levels = c('3+', '2', '1'))
  ) %>%
  group_by(gbm_specific) %>%
  dplyr::count(occurrence) %>%
  ggplot(aes(x = occurrence, y = n, fill = gbm_specific)) +
  geom_bar(stat = 'identity', position = 'fill') +
  labs(x = 'N(ONT)', y = '', fill = 'Disease-specific') +
  coord_flip() +
  # scale_x_reverse() +
  scale_fill_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) + 
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
    legend.position = 'none',
    axis.title.x = element_blank(),
  )

m2.3 <- df_trend_svs %>% 
  # filter(occurrence > 1) %>%
  mutate(
    gbm_specific = !(gene %in% df_dlpfc_svs$gene),
    occurrence = ifelse(occurrence >= 7, '7+', occurrence),
    occurrence = factor(occurrence, levels = c('7+', '6', '5', '4', '3', '2', '1'))
  ) %>%
  group_by(gbm_specific) %>%
  dplyr::count(occurrence) %>%
  ggplot(aes(x = occurrence, y = n, fill = gbm_specific)) +
    geom_bar(stat = 'identity', position = 'fill') +
    labs(x = 'N(SR)', y = 'Proportion', fill = 'Disease-specific') +
    coord_flip() +
    # scale_x_reverse() +
    scale_fill_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) + 
    theme_classic() + 
    theme(
      text=element_text(size = 12, family = 'Arial'),
      plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
      legend.position = 'bottom',
      legend.margin = margin(-5, 0, 0, 0)
    )

# m2.2 <- rbind(
#   (df_ont_svs %>% dplyr::select('gene', 'occurrence') %>% mutate(geneset = 'ONT')),
#   (df_trend_svs %>% dplyr::select('gene', 'occurrence') %>% mutate(geneset = 'SR'))
# ) %>%
#   # filter(occurrence >= 2) %>%
#   mutate(
#     gbm_specific = !(gene %in% df_dlpfc_svs$gene)
#   ) %>%
#   dplyr::count(geneset, gbm_specific) %>%
#   ggplot(aes(x = geneset, y = n, fill = gbm_specific)) +
#   geom_bar(stat = 'identity', position = 'fill') +
#   scale_fill_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = '#f4d56e')) + 
#   labs(x = '', y = 'Proportion', fill = 'Disease-specific') + 
#   geom_text(aes(label = n), position = position_fill(vjust = 0.5)) +
#   coord_flip() + 
#   theme_classic() + 
#   theme(
#     text=element_text(size = 12, family = 'Arial'),
#     plot.title = element_text(size = 12, family = 'Arial'),
#     legend.position = 'top'
#   )

m2 <- cowplot::plot_grid(m2.1, m2.2, m2.3, ncol = 1, align = 'v', rel_heights = c(1.6, 0.9, 2))
m2

ggsave(sprintf("%s/figures/main_svs_recurrence.png", res_dir_ont), m2, width = 4, height = 4)
ggsave(sprintf("%s/figures/main_svs_recurrence.pdf", res_dir_ont), m2, width = 4, height = 4)

# the subset of overlapped GBM-specific SVS genes
intersect(
  df_trend_svs %>% filter(occurrence > 1, (! gene %in% df_dlpfc_svs$gene)) %>% pull(gene),
  df_ont_svs %>% filter(occurrence >= 1, (! gene %in% df_dlpfc_svs$gene)) %>% pull(gene)
)

intersect(
  df_ont_svs %>% filter(occurrence > 1, (! gene %in% df_dlpfc_svs$gene)) %>% pull(gene),
  df_trend_svs %>% filter(occurrence > 1, (! gene %in% df_dlpfc_svs$gene)) %>% pull(gene)
)


## SupFig: Venn diagrams of SV genes
library(ggvenn)
venn_svs_all <- list(
  `ONT SVS (GBM/DMG)` = df_ont_svs$gene %>% unique(),
  `SR SVS (GBM)` = df_trend_svs$gene %>% unique(),
  `SR SVS (DLPFC)` = df_dlpfc_svs$gene %>% unique()
)
venn_svs_rep <- list(
  `ONT SVS (GBM/DMG)` = df_ont_svs %>% filter(occurrence > 1) %>% pull(gene) %>% unique(),
  `SR SVS (GBM)` = df_trend_svs %>% filter(occurrence > 1) %>% pull(gene) %>% unique(),
  `SR SVS (DLPFC)` = df_dlpfc_svs %>% filter(occurrence > 1) %>% pull(gene) %>% unique()
)
venn_sve_all <- list(
  `ONT SVE (GBM/DMG)` = df_ont_sve$gene %>% unique(),
  `SR SVE (GBM)` = df_trend_sve$gene %>% unique(),
  `SR SVE (DLPFC)` = df_dlpfc_sve$gene %>% unique()
)
venn_sve_rep <- list(
  `ONT SVE (GBM/DMG)` = df_ont_sve %>% filter(occurrence > 1) %>% pull(gene) %>% unique(),
  `SR SVE (GBM)` = df_trend_sve %>% filter(occurrence > 1) %>% pull(gene) %>% unique(),
  `SR SVE (DLPFC)` = df_dlpfc_sve %>% filter(occurrence > 1) %>% pull(gene) %>% unique()
)

s1.1 <- ggvenn(
  venn_svs_all, 
  show_percentage = TRUE, digits = 0,
  stroke_size = 0.5, text_size = 4, set_name_size = 4
) + labs(title = 'Spatially variable in at least one sample') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
s1.2 <- ggvenn(
  venn_svs_rep,
  show_percentage = TRUE, digits = 0,
  stroke_size = 0.5, text_size = 4, set_name_size = 4
) + labs(title = 'Spatially variable in at least two samples') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
s1.3 <- ggvenn(
  venn_sve_all, 
  show_percentage = TRUE, digits = 0,
  stroke_size = 0.5, text_size = 4, set_name_size = 4
) + labs(title = 'Spatially variable in at least one sample') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))
s1.4 <- ggvenn(
  venn_sve_rep,
  show_percentage = TRUE, digits = 0,
  stroke_size = 0.5, text_size = 4, set_name_size = 4
) + labs(title = 'Spatially variable in at least two samples') + 
  theme(plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5))

s1 <- cowplot::plot_grid(s1.1, s1.2, s1.3, s1.4, nrow = 1, align = 'hv')
s1
ggsave(sprintf("%s/figures/sup_sv_venn.png", res_dir_ont), s1, width = 16, height = 4)

## Fig 7D: Functional enrichment results
df_enrich <- read_csv(sprintf("%s/figures/source_data/enrich_top20.csv", res_dir_ont)) %>%
  mutate(
    geneset = factor(geneset, levels = c('SVS', 'SVE'), labels = c('SVS', 'SVENS')),
    name = recode(
      name, 
      'SARS-CoV-2 activates/modulates innate and adaptive immune responses' = 'SARS-CoV-2 innate and adaptive immune responses'
    )
  )

# All recurrent SV genes, REAC
m3.1 <- df_enrich %>%
  filter(group == 'ont_all', source == 'REAC') %>%
  arrange(precision_diff, name) %>%
  mutate(
    # name = stringr::str_wrap(name, 40), # auto line break
    name_cat = factor(name, levels = unique(name)),
    geneset = factor(geneset, levels = c('SVENS', 'SVS'), 
                     labels = c('Recurrent SVENS', 'Recurrent SVS'))
  ) %>%
  # keep the top 10 terms of SVS genes
  dplyr::slice((n()-19):n()) %>%
  ggplot(aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = 'Precision (proportion of term genes)', x = 'Reactome (ONT)', fill = 'Gene set') + 
  # geom_vline(xintercept = 5.5, linetype = 'dashed') + 
  coord_flip() +
  scale_fill_manual(values=c('Recurrent SVENS' = '#2EC4B6', 'Recurrent SVS' = '#FF9F1C')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 12, family = 'Arial'),
    legend.position='none',
    legend.position.inside = c(0.7, 0.85),
    legend.background = element_blank()
  )
m3.2 <- df_enrich %>%
  filter(group == 'vis_all', source == 'REAC') %>%
  arrange(precision_diff, name) %>%
  mutate(
    # name = stringr::str_wrap(name, 40), # auto line break
    name_cat = factor(name, levels = unique(name)),
    geneset = factor(geneset, levels = c('SVENS', 'SVS'), 
                     labels = c('Recurrent SVENS', 'Recurrent SVS'))
  ) %>%
  # keep the top 10 terms of SVS genes
  dplyr::slice((n()-19):n()) %>%
  ggplot(aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = 'Precision (proportion of term genes)', fill = 'Gene set', 
       x = 'Reactome (SR)') + 
  # geom_vline(xintercept = 5.5, linetype = 'dashed') + 
  coord_flip() +
  scale_fill_manual(values=c('Recurrent SVENS' = '#2EC4B6', 'Recurrent SVS' = '#FF9F1C')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 12, family = 'Arial'),
    legend.position='inside',
    legend.position.inside = c(-0.8, 1.1),
    legend.background = element_blank()
  )

m3 <- cowplot::plot_grid(m3.1, m3.2, nrow = 2, align = 'hv')
m3
# ggsave(sprintf("%s/figures/enrich_reac_rec.png", res_dir_ont), m3, width = 8, height = 6)
ggsave(sprintf("%s/figures/enrich_reac_rec_new.png", res_dir_ont), m3, width = 8, height = 5)
ggsave(sprintf("%s/figures/enrich_reac_rec_new.pdf", res_dir_ont), m3, width = 8, height = 5)

## Fig S6E: All recurrent SV genes, KEGG and GO:BP
m4.1 <- df_enrich %>%
  filter(group == 'ont_all', source == 'KEGG') %>%
  arrange(desc(precision_diff), name) %>%
  mutate(
    # name = stringr::str_wrap(name, 30), # auto line break
    name_cat = factor(name, levels = unique(name)),
    geneset = factor(geneset, levels = c('SVENS', 'SVS'), 
                     labels = c('Recurrent SVENP', 'Recurrent SVP'))
  ) %>%
  # keep the top 5 terms of SVS genes
  dplyr::slice(1:10) %>%
  ggplot(aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = 'Precision', x = '', fill = 'Gene set', 
       title = 'KEGG (ONT)') + 
  coord_flip() +
  scale_fill_manual(values=c('Recurrent SVENP' = '#2EC4B6', 'Recurrent SVP' = '#FF9F1C')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial'),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 12, family = 'Arial'),
    legend.position='none',
    legend.position.inside = c(0.25, 0.85),
    legend.background = element_blank()
  )
m4.2 <- df_enrich %>%
  filter(group == 'vis_all', source == 'KEGG') %>%
  arrange(desc(precision_diff), name) %>%
  mutate(
    # name = stringr::str_wrap(name, 30), # auto line break
    name_cat = factor(name, levels = unique(name)),
    geneset = factor(geneset, levels = c('SVENS', 'SVS'), 
                     labels = c('Recurrent SVENP', 'Recurrent SVP'))
  ) %>%
  # keep the top 5 terms of SVS genes
  dplyr::slice(1:10) %>%
  ggplot(aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = 'Precision', x = '', fill = 'Gene set', 
       title = 'KEGG (SR)') + 
  coord_flip() +
  scale_fill_manual(values=c('Recurrent SVENP' = '#2EC4B6', 'Recurrent SVP' = '#FF9F1C')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial'),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 12, family = 'Arial'),
    legend.position='none',
    legend.position.inside = c(0.65, 0.85),
    legend.background = element_blank()
  )
m4.3 <- df_enrich %>%
  filter(group == 'ont_all', source == 'GO:BP') %>%
  arrange(desc(precision_diff), name) %>%
  mutate(
    # name = stringr::str_wrap(name, 30), # auto line break
    name_cat = factor(name, levels = unique(name)),
    geneset = factor(geneset, levels = c('SVENS', 'SVS'), 
                     labels = c('Recurrent SVENP', 'Recurrent SVP'))
  ) %>%
  # keep the top 5 terms of SVS genes
  dplyr::slice(1:10) %>%
  ggplot(aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = 'Precision', x = '', fill = 'Gene set', 
       title = 'GO:BP (ONT)') + 
  coord_flip() +
  scale_fill_manual(values=c('Recurrent SVENP' = '#2EC4B6', 'Recurrent SVP' = '#FF9F1C')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial'),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 12, family = 'Arial'),
    legend.position='none',
    legend.position.inside = c(0.25, 0.85),
    legend.background = element_blank()
  )
m4.4 <- df_enrich %>%
  filter(group == 'vis_all', source == 'GO:BP') %>%
  arrange(desc(precision_diff), name) %>%
  mutate(
    # name = stringr::str_wrap(name, 30), # auto line break
    name_cat = factor(name, levels = unique(name)),
    geneset = factor(geneset, levels = c('SVENS', 'SVS'), 
                     labels = c('Recurrent SVENP', 'Recurrent SVP'))
  ) %>%
  # keep the top 5 terms of SVS genes
  dplyr::slice(1:10) %>%
  ggplot(aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = 'Precision', x = '', fill = 'Gene set', 
       title = 'GO:BP (SR)') + 
  coord_flip() +
  scale_fill_manual(values=c('Recurrent SVENP' = '#2EC4B6', 'Recurrent SVP' = '#FF9F1C')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial'),
    axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 12, family = 'Arial'),
    legend.position='right',
    legend.position.inside = c(0.65, 0.85),
    legend.background = element_blank()
  )

m4 <- cowplot::plot_grid(m4.3, m4.4, m4.1, m4.2, nrow = 2, rel_widths = c(1, 1.3))
m4
ggsave(sprintf("%s/figures/enrich_go_kegg_rec.png", res_dir_ont), m4, width = 10, height = 4)

## Fig S6F-H: visualizing detailed enrichment terms
library(clusterProfiler)
library(org.Hs.eg.db)

# ONT recurrent SVS genes
eg_ont <- bitr(
  df_ont_svs$gene[df_ont_svs$occurrence >= 2] %>% unique(),
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db
)
# TREND recurrent SVS genes
eg_trend <- bitr(
  df_trend_svs$gene[df_trend_svs$occurrence >= 2] %>% unique(),
  fromType = 'SYMBOL',
  toType = 'ENTREZID',
  OrgDb = org.Hs.eg.db
)

# Reactome pathway enrichment results
library(ReactomePA)
reac_ont <- enrichPathway(
  gene = eg_ont$ENTREZID,
  pvalueCutoff = 0.01,
  readable = TRUE
)
reac_trend <- enrichPathway(
  gene = eg_trend$ENTREZID,
  pvalueCutoff = 0.01,
  readable = TRUE
)

# KEGG pathway enrichment result
kk_ont <- enrichKEGG(
  gene = eg_ont$ENTREZID,
  organism = 'hsa',
  keyType = 'kegg',
  pvalueCutoff = 0.01
)
kk_trend <- enrichKEGG(
  gene = eg_trend$ENTREZID,
  organism = 'hsa',
  keyType = 'kegg',
  pvalueCutoff = 0.01
)

# extract genes per KEGG term
kegg_ont_glist <- list()
for (kegg_id in c(kk_ont@result$ID)[1:20]) {
  glist <- kk_ont@result[kk_ont@result$ID == kegg_id, 'geneID'] %>% str_split_1('/') %>% bitr(
    fromType = 'ENTREZID',
    toType = 'SYMBOL',
    OrgDb = org.Hs.eg.db
  )
  kegg_ont_glist[[kegg_id]] <- df_ont_svs %>% filter(
    gene %in% glist$SYMBOL
  ) %>% pull(gene)
}

kegg_trend_glist <- list()
for (kegg_id in c(kk_trend@result$ID)[1:40]) {
  glist <- kk_trend@result[kk_trend@result$ID == kegg_id, 'geneID'] %>% str_split_1('/') %>% bitr(
    fromType = 'ENTREZID',
    toType = 'SYMBOL',
    OrgDb = org.Hs.eg.db
  )
  kegg_trend_glist[[kegg_id]] <- df_trend_svs %>% filter(
    gene %in% glist$SYMBOL
  ) %>% pull(gene)
}

# Fig S6F: ONT top KEGG pathways
df <- df_ont_svs %>% 
  mutate(
    rank = row_number(),
    gbm_specific = !(gene %in% df_dlpfc_svs$gene)
  )
  
p_list <- list()
for (kegg_id in c('hsa04612')) {
  p <- ggplot(df, aes(x = rank, y = `pvalue_hsic.ir`)) + 
    geom_point(alpha = 0.5, color = 'gray') + 
    scale_y_continuous(transform = c('log10', 'reverse')) +
    geom_point(
      data = df %>% filter(gene %in% kegg_ont_glist[[kegg_id]]),
      mapping = aes(color = gbm_specific)
    ) +
    geom_label_repel(
      data = df %>% filter(gene %in% kegg_ont_glist[[kegg_id]]),
      aes(label = gene, x = rank, y = `pvalue_hsic.ir`, color = gbm_specific),
      force = 20, max.overlaps = 30, nudge_x = 0, nudge_y = 30,
      fontface = 'italic'
    ) +
    scale_color_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) +
    labs(x = 'Rank', y = 'HSIC-IR p-value', color = 'Disease specific',
         title = stringr::str_wrap(
      paste(kk_ont@result[kegg_id, 'Description'], '-', kegg_id, '(ONT)'), 40
    )) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'inside',
      legend.position.inside = c(0.75, 0.8),
      legend.background = element_blank()
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, ncol = 1, align = 'h')
p

ggsave(sprintf("%s/figures/sup_kegg_ont_details.png", res_dir_ont), p, width = 3.5, height = 3.5)

# Fig S6G: TREND top KEGG pathways
df <- df_trend_svs %>% 
  mutate(
    rank = row_number(),
    gbm_specific = !(gene %in% df_dlpfc_svs$gene)
  )

# infection pathways
p_list <- list()
for (kegg_id in c('hsa05100', 'hsa05203', 'hsa04144'
                  # 'hsa05131', 'hsa05132'
                  # 'hsa04510', 'hsa04530',
)) {
  p <- ggplot(df, aes(x = rank, y = `pvalue_hsic.ir`)) + 
    geom_point(alpha = 0.5, color = 'gray') + 
    scale_y_continuous(transform = c('log10', 'reverse')) +
    geom_point(
      data = df %>% filter(gene %in% kegg_trend_glist[[kegg_id]]),
      color = 'black'
    ) +
    geom_text_repel(
      data = df %>% filter(gene %in% kegg_trend_glist[[kegg_id]]),
      aes(label = gene, x = rank, y = `pvalue_hsic.ir`, color = gbm_specific),
      force = 20, max.overlaps = 30, nudge_x = 20, nudge_y = 20,
      fontface = 'italic'
    ) +
    scale_color_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) +
    labs(x = 'Rank', y = 'HSIC-IR p-value', color = 'Disease specific',
         title = stringr::str_wrap(
           paste(kk_trend@result[kegg_id, 'Description'], '-', kegg_id, '(SR)'), 40
         )) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'inside',
      legend.position.inside = c(0.75, 0.85),
      legend.background = element_blank()
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, ncol = 3, align = 'h')
p

ggsave(sprintf("%s/figures/sup_kegg_trend_infection.png", res_dir_ont), p, width = 10.5, height = 3.5)

# Fig S6H: TREND Spliceosome 
p_list <- list()
for (kegg_id in c('hsa03040')) {
  p <- ggplot(df, aes(x = rank, y = `pvalue_hsic.ir`)) + 
    geom_point(alpha = 0.5, color = 'gray') + 
    scale_y_continuous(transform = c('log10', 'reverse')) +
    geom_point(
      data = df %>% filter(gene %in% kegg_trend_glist[[kegg_id]]),
      color = 'black'
    ) +
    geom_text_repel(
      data = df %>% filter(gene %in% kegg_trend_glist[[kegg_id]]),
      aes(label = gene, x = rank, y = `pvalue_hsic.ir`, color = gbm_specific),
      force = 20, max.overlaps = 30, nudge_x = 20, nudge_y = 20,
      fontface = 'italic'
    ) +
    scale_color_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) +
    labs(x = 'Rank', y = 'HSIC-IR p-value', color = 'Disease specific',
         title = stringr::str_wrap(
           paste(kk_trend@result[kegg_id, 'Description'], '-', kegg_id, '(SR)'), 40
         )) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'inside',
      legend.position.inside = c(0.75, 0.85),
      legend.background = element_blank()
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, ncol = 1, align = 'h')
p

ggsave(sprintf("%s/figures/sup_kegg_trend_spliceosome.png", res_dir_ont), p, width = 3.5, height = 3.5)

# Fig S6I: other important bio pathways
p_list <- list()
for (kegg_id in c('hsa04510', 'hsa04810', 'hsa04530', 'hsa05214',  'hsa04022'
                  )) {
  p <- ggplot(df, aes(x = rank, y = `pvalue_hsic.ir`)) + 
    geom_point(alpha = 0.5, color = 'gray') + 
    scale_y_continuous(transform = c('log10', 'reverse')) +
    geom_point(
      data = df %>% filter(gene %in% kegg_trend_glist[[kegg_id]]),
      color = 'black'
    ) +
    geom_text_repel(
      data = df %>% filter(gene %in% kegg_trend_glist[[kegg_id]]),
      aes(label = gene, x = rank, y = `pvalue_hsic.ir`, color = gbm_specific),
      force = 20, max.overlaps = 30, nudge_x = 20, nudge_y = 20,
      fontface = 'italic'
    ) +
    scale_color_manual(values = c('TRUE' = '#e65c1b', 'FALSE' = 'black')) +
    labs(x = 'Rank', y = 'HSIC-IR p-value', color = 'Disease specific',
         title = stringr::str_wrap(
           paste(kk_trend@result[kegg_id, 'Description'], '-', kegg_id, '(SR)'), 30
         )) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
      legend.position = 'inside',
      legend.position.inside = c(0.75, 0.85),
      legend.background = element_blank()
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, ncol = 5, align = 'h')
p

ggsave(sprintf("%s/figures/sup_kegg_trend_details.png", res_dir_ont), p, width = 17.5, height = 3.5)


## Fig S6D: AS type
# ONT AS type
df_as_types <- c()
for (sample_id in c('GBM_1', 'GBM_2', 'GBM_3', 'GBM_4', 'GBM_5', 'GBM_6',
                    'DMG_1', 'DMG_2', 'DMG_3', 'DMG_4', 'DMG_5')) {
  for (geneset in c('svs_iso', 'svens_iso')) {
    as_type_counts <- read_table(
      sprintf("%s/transcripts/%s/suppa.out/%s/n_events.txt", res_dir_ont, sample_id, geneset),
      col_names = c('type', 'count')
    ) %>%
      mutate(
        sample_id = sample_id, 
        geneset = ifelse(geneset == 'svs_iso', 'SVP', 'SVENP'),
        prop = count / sum(count)
      )
    df_as_types <- rbind(df_as_types, as_type_counts)
  }
}
df_as_types <- df_as_types %>%
  mutate(geneset = factor(geneset, levels = c('SVENP', 'SVP')))

s2 <- ggplot(df_as_types, aes(x = type, y = prop, fill = geneset)) + 
  geom_boxplot(position = 'dodge') + 
  stat_compare_means(label = 'p.signif', method = 'wilcox.test') +
  labs(y = 'Proportion', x = 'Alternative splicing event', fill = 'Gene set', 
       title = 'ONT AS event types per sample') + 
  scale_fill_manual(values = c('SVP'= '#FF9F1C', 'SVENP'= '#2EC4B6')) +
  theme_classic() + 
  theme(
    text=element_text(size = 14, family = 'Arial'),
    axis.text=element_text(size = 14, family = 'Arial'),
    plot.title=element_text(hjust = 0.5, size = 14, family = 'Arial'),
    legend.position='none'
  )
s2
ggsave(sprintf("%s/figures/sup_ont_as_types.png", res_dir_ont), s2, width = 3.5, height = 3.5)

## TREND type annotation
library(rtracklayer)

# helper function to count trend event types
count_trend_type <- function(sample_id) {
  # load sierra annotation outputs per sample
  data_dir_vis <- "/Users/jysumac/Projects/SPLISOSM_paper/data/gbm_visium_cell_24/"
  df_sample <-  read.table(
    sprintf("%s/GEO_sample_list.tsv", data_dir_vis), sep = '\t', 
    col.names = c('sample', 'GEO', 'group', 'description')
  )
  geo_id <- df_sample %>% filter(sample == sample_id) %>% pull(GEO)
  
  peak_annot <- read.table(sprintf(
    "%s/sierra_peaks_default_individual/%s/peak_default.annot.no_correction.txt", 
    data_dir_vis, geo_id
  ), sep = '\t') %>%
    rownames_to_column('peak_id')

  # load SVS and SVENS peaks
  svs_peak <- import.bed(sprintf("%s/events/%s/svs.exon.bed", res_dir_trend, sample_id))
  svens_peak <- import.bed(sprintf("%s/events/%s/svens.exon.bed", res_dir_trend, sample_id))
  
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
      label = ifelse(pval < 0.01, '**', label)
    ) %>%
    pivot_longer(
      cols = c(SVS, SVENS), names_to = 'dataset', values_to = 'count'
    ) %>%
    mutate(prop = count / ifelse(
      dataset == 'SVS', n_svs, n_sve
    )) %>%
    mutate(
      sample = sample_id,
      dataset = factor(dataset, levels = c('SVENS', 'SVS'), labels = c('SVENP', 'SVP')),
    )
  
  return(df_event)
}

# loop over samples and count event types
data_dir_vis <- "/Users/jysumac/Projects/SPLISOSM_paper/data/gbm_visium_cell_24/"
df_sample <- read.table(
  sprintf("%s/GEO_sample_list.tsv", data_dir_vis), sep = '\t', 
  col.names = c('sample', 'GEO', 'group', 'description')
) %>%
  filter(group == 'GBM')
df_events <- lapply(df_sample$sample, count_trend_type) %>%
  do.call(rbind, .)

## Fig S6D: TREND AS type
# plot event type frequency
s3.1 <- ggplot(
  df_events %>% filter(event == 'Junctions'), 
  aes(x = event, y = prop, group = interaction(event, dataset))
) + 
  geom_boxplot(aes(fill = dataset), position="dodge", width = 0.75) + 
  stat_compare_means(label = 'p.signif', method = 'wilcox.test') +
  labs(y = 'Proportion of junction peaks', x = '', fill = 'Gene set') +
  scale_fill_manual(values = c('SVP' = '#FF9F1C', 'SVENP' = '#2EC4B6')) +
  theme_classic() + 
  theme(
    legend.position = 'none',
    text = element_text(size = 14, family='Arial'),
    axis.text = element_text(size = 14, family='Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family='Arial')
  )

s3.2 <- ggplot(
  df_events %>% filter(! event %in% c('Junctions', 'Intron')), 
  aes(x = event, y = prop, group = interaction(event, dataset))
) + 
  geom_boxplot(aes(fill = dataset), position="dodge", width = 0.75) + 
  stat_compare_means(label = 'p.signif', method = 'wilcox.test') +
  # coord_flip() + 
  labs(y = 'Proportion of peaks that overlap with', x = '', fill = 'Gene set',
       title = 'TREND annotation per sample') +
  scale_fill_manual(values = c('SVP' = '#FF9F1C', 'SVENP' = '#2EC4B6')) +
  theme_classic() + 
  theme(
    legend.position = 'inside',
    legend.position.inside = c(0.7, 0.8),
    text = element_text(size = 14, family='Arial'),
    axis.text = element_text(size = 14, family='Arial'),
    plot.title = element_text(size = 14, family='Arial')
  )

s3 <- cowplot::plot_grid(s3.1, s3.2, nrow = 1, align = 'hv', rel_widths = c(1, 2.4))
s3

ggsave(sprintf("%s/figures/sup_peak_annot.png", res_dir_ont), s3, width = 5.5, height = 3.5)

## Fig 7E: HLA gene spatial expression
# load expression and annotation of GBM1
anno_gbm1 <- read_csv(sprintf("%s/figures/source_data/anno_GBM1.csv", res_dir_ont)) %>%
  mutate(region = ifelse(anno_hzy == 'Tumor', 'Tumor', 'Immune inf'))
hla_expr_gbm1 <- read_csv(sprintf("%s/figures/source_data/hla_expr_GBM1.csv", res_dir_ont)) %>%
  left_join(anno_gbm1 %>% dplyr::select('...1', 'region'), by = c( 'barcode' = '...1'))

# prepare separate dataframes for each HLA-DRB1 isoform with normalized alpha values
df <- hla_expr_gbm1 %>%
  filter(gene == 'HLA-DRB1', layer == 'counts') 

df_iso1 <- df %>%
  filter(isoform == 'HLA-DRB1-201', value > 0) %>%
  mutate(alpha_value = value / median(value), isoform_label = 'HLA-DRB1-201')

df_iso2 <- df %>%
  filter(isoform == 'HLA-DRB1_Iso_1', value > 0) %>%
  mutate(alpha_value = value / median(value), isoform_label = 'HLA-DRB1-Novel1')

df_iso3 <- df %>%
  filter(isoform == 'HLA-DRB1_Iso_2', value > 0) %>%
  mutate(alpha_value = value / median(value), isoform_label = 'HLA-DRB1-Novel2')

# Combine datasets for hexbin plotting
df_combined <- bind_rows(df_iso1, df_iso2, df_iso3)

df_combined_anno <- df_combined %>% 
  pivot_wider(
    id_cols = c(barcode, array_row, array_col),
    values_from = value,
    names_from = isoform_label
  ) %>%
  mutate(
    annotation = case_when(
      `HLA-DRB1-201` == pmax(`HLA-DRB1-201`, `HLA-DRB1-Novel1`, `HLA-DRB1-Novel2`, na.rm = TRUE) ~ 'HLA-DRB1-201',
      `HLA-DRB1-Novel1` == pmax(`HLA-DRB1-201`, `HLA-DRB1-Novel1`, `HLA-DRB1-Novel2`, na.rm = TRUE) ~ 'HLA-DRB1-Novel1',
      `HLA-DRB1-Novel2` == pmax(`HLA-DRB1-201`, `HLA-DRB1-Novel1`, `HLA-DRB1-Novel2`, na.rm = TRUE) ~ 'HLA-DRB1-Novel2',
      TRUE ~ NA_character_ # Handle cases where all values are NA
    )
  )

# merge annotation back to df_combined
df_combined <- left_join(
  df_combined, df_combined_anno %>% 
    dplyr::select(-c(`HLA-DRB1-201`, `HLA-DRB1-Novel1`, `HLA-DRB1-Novel2`)), 
  by = c('barcode', 'array_row', 'array_col'))

# Plot with geom_hex
m6.1 <- ggplot(df_combined, aes(x = -array_row, y = -array_col)) +
  geom_point(aes(color = annotation), 
             fill = NA, alpha = 1, shape = 21, size = 1,
             show.legend = FALSE) +
  geom_point(
    aes(fill = isoform_label, alpha = value),
    shape = 21, size = 1,
    show.legend = TRUE
  ) +
  scale_color_manual(
    values = c(
      'HLA-DRB1-201' = 'red',
      'HLA-DRB1-Novel1' = 'blue',
      'HLA-DRB1-Novel2' = '#4caf50'
    )
  ) +
  scale_fill_manual(
    values = c(
      'HLA-DRB1-201' = 'red',
      'HLA-DRB1-Novel1' = 'blue',
      'HLA-DRB1-Novel2' = '#4caf50'
    )
  ) +
  scale_alpha_continuous(range = c(0.1, 1)) +
  labs(
    title = 'HLA-DRB1 expression in GBM1 (ONT, IDH-mt)',
    color = '',
    fill = '',
    alpha = 'Raw counts'
  ) +
  theme_void() +
  theme(
    plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
    text = element_text(size = 12, family = 'Arial'),
    legend.position = 'bottom',
    legend.box="vertical", legend.margin=margin(),
    legend.text = element_text(size = 12, family = 'Arial')
  ) +
  guides(
    fill = guide_legend(override.aes = list(alpha = 1, size = 3)),
    alpha = guide_legend(override.aes = list(fill = 'gray', size = 3), ncol = 6)
  )
m6.1

ggsave(sprintf("%s/figures/main_HLA-DRB1_GBM1.png", res_dir_ont), m6.1, width = 3.5, height = 4.2)
ggsave(sprintf("%s/figures/main_HLA-DRB1_GBM1.pdf", res_dir_ont), m6.1, width = 3.5, height = 4.2)

# m6.1 <- anno_gbm1 %>%
#   ggplot(aes(x = -array_row, y = -array_col, color = region)) + 
#   geom_point(size = 1) +
#   labs(title = 'GBM1 (ONT, IDH-mt)', color = 'Region') +
#   scale_color_manual(values = c('Tumor' = '#f7c6a5', 'Immune inf' = '#6a4c93')) +
#   theme_void() + 
#   theme(
#     plot.title = element_text(size = 12, family = 'Arial', hjust = 0.5),
#     text = element_text(size = 12, family = 'Arial'),
#     legend.position = 'bottom',
#     legend.text = element_text(size = 12, family = 'Arial')
#   ) + 
#   guides(color = guide_legend(override.aes = list(size = 3)))


## Fig 7F: B2M and CD74 expression
# B2M
m6.2 <- hla_expr_gbm1 %>% filter(
  transcript_id %in% c('B2M-212', 'B2M-214'),
) %>% mutate(
  transcript_id = factor(transcript_id, levels = c('B2M-214', 'B2M-212'), 
                         labels = c('B2M-211\n(coding)', 'B2M-210\n(retained intron)'))
) %>% filter(
  layer == 'counts',
) %>% ggplot(aes(x = transcript_id, y = value + 0.01, group = interaction(region, transcript_id))) + 
  geom_boxplot(aes(fill = region)) +
  scale_y_log10() +
  scale_fill_manual(values = c('Tumor' = '#f7c6a5', 'Immune inf' = '#6a4c93')) +
  labs(x = '', y = 'Raw transcript count', title = 'B2M (ONT)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    title = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
  )

m6.3 <- hla_expr_gbm1 %>% filter(
  transcript_id %in% c('B2M-212', 'B2M-214'),
) %>% mutate(
  transcript_id = factor(transcript_id, levels = c('B2M-214', 'B2M-212'), 
                         labels = c('B2M-211\n(coding)', 'B2M-210\n(retained intron)'))
) %>% filter(
  layer == 'ratios_obs',
) %>% ggplot(aes(x = transcript_id, y = value, group = interaction(region, transcript_id))) + 
  geom_boxplot(aes(fill = region)) +
  scale_fill_manual(values = c('Tumor' = '#f7c6a5', 'Immune inf' = '#6a4c93')) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format') +
  labs(x = '', y = 'Observed transcript ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    title = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
  )

# CD74
m6.4 <- hla_expr_gbm1 %>% filter(
  gene == 'CD74',
  isoform %in% c('CD74-202', 'CD74_Iso_2', 'CD74_Iso_3', 'CD74_Iso_5') # intron retention isoforms
) %>% mutate(
  t_name = ifelse(transcript_id == 'CD74-202', 'CD74-202\n(coding)', 'CD74-novel\n(retained intron)')
) %>%
  group_by(t_name, layer, region, barcode) %>%
  summarize(value = sum(value)) %>%
  filter(
    layer == 'counts'
  ) %>% 
  ggplot(aes(x = t_name, y = value + 0.01, group = interaction(region, t_name))) + 
  geom_boxplot(aes(fill = region)) +
  scale_fill_manual(values = c('Tumor' = '#f7c6a5', 'Immune inf' = '#6a4c93')) +
  scale_y_log10() + 
  labs(x = '', y = 'Raw transcript count', title = 'CD74 (ONT)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    title = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
  )

m6.5 <- hla_expr_gbm1 %>% filter(
  gene == 'CD74',
  isoform %in% c('CD74-202', 'CD74_Iso_2', 'CD74_Iso_3', 'CD74_Iso_5') # intron retention isoforms
) %>% mutate(
  t_name = ifelse(transcript_id == 'CD74-202', 'CD74-202\n(coding)', 'CD74-novel\n(retained intron)')
) %>%
  group_by(t_name, layer, region, barcode) %>%
  summarize(value = sum(value)) %>%
  filter(
    layer == 'ratios_smoothed'
  ) %>% 
  ggplot(aes(x = t_name, y = value, group = interaction(region, t_name))) + 
  geom_boxplot(aes(fill = region)) +
  scale_fill_manual(values = c('Tumor' = '#f7c6a5', 'Immune inf' = '#6a4c93')) +
  stat_compare_means(method = 'wilcox.test', label = 'p.format') +
  labs(x = '', y = 'Smoothed transcript ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    title = element_text(size = 12, family = 'Arial'),
    axis.text.x = element_text(size = 12, family = 'Arial'),
    legend.position = 'none',
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
  )

m6 <- cowplot::plot_grid(m6.2, m6.3, m6.4, m6.5, align = 'hv', nrow = 1, rel_widths = c(1, 1, 1, 1))
m6
ggsave(sprintf("%s/figures/main_B2M_CD74_GBM1.png", res_dir_ont), m6, width = 12, height = 4)
ggsave(sprintf("%s/figures/main_B2M_CD74_GBM1.pdf", res_dir_ont), m6, width = 8, height = 4)

