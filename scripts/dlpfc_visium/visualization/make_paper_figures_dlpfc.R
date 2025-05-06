library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(universalmotif)

extrafont::loadfonts()

data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/human_dlpfc/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/human_dlpfc/"
date = '0123'

### TREND main figures
## make figure directory if doesn't exist
dir.create(sprintf("%s/figures", res_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(sprintf("%s/events", res_dir), recursive = TRUE, showWarnings = FALSE)

## load mouse SV results
mm_sv_cbs <- read.csv(
  "/Users/jysumac/Projects/SPLISOSM_paper/results/visium_mouse_cbs/sv_results/cbs_peak_sv_combined_1119.csv", 
) %>%
  mutate(is_trend_svs = padj_hsic.ir < 0.01,
         is_trend_sve = padj_hsic.gc < 0.01,
         is_trend_expressed = count_avg > 0.5)

## convert mouse genes to human homologs
library("biomaRt")
# use archived ensembl 98 for renamed genes (Sept -> Septin)
mouse_98 <- useEnsembl(
  biomart = "ensembl", dataset = "mmusculus_gene_ensembl", version = 98
)
m2h_98 <- getBM(
  attributes=c('ensembl_gene_id',
               'external_gene_name',
               'hsapiens_homolog_ensembl_gene',
               'hsapiens_homolog_associated_gene_name'),
  filters = "mgi_symbol",
  values = mm_sv_cbs$gene, mart = mouse_98, uniqueRows = T
) %>%
  mutate(across(where(is.character), ~ na_if(.,"")))
write.csv(m2h_98, paste0(res_dir, "events/m2h_sv_ensembl98.csv"), 
          row.names = FALSE, quote = FALSE)

# load and merge with ensembl 113
m2h_113 <- read_csv(paste0(res_dir, "events/m2h_sv_ensembl113.csv"))
m2h <- rbind(m2h_98, m2h_113) %>%
  drop_na() %>%
  distinct(external_gene_name, hsapiens_homolog_associated_gene_name, .keep_all = TRUE)
write.csv(m2h, paste0(res_dir, "events/m2h_sv_combined.csv"), 
          row.names = FALSE, quote = FALSE)

## load human SV results
m2h_sv_cbs <- mm_sv_cbs %>%
  left_join(m2h, by = c('gene' = 'external_gene_name'))

df_svs_pval <- read.csv(paste0(res_dir, "trend_svs_dlpfc_", date, ".csv")) %>%
  mutate(is_mouse_svs = gene %in% (m2h_sv_cbs %>% filter(is_trend_svs) %>% 
                                     pull(hsapiens_homolog_associated_gene_name)))
df_sve_pval <- read.csv(paste0(res_dir, "trend_sve_dlpfc_", date, ".csv")) %>%
  mutate(is_mouse_sve = gene %in% (m2h_sv_cbs %>% filter(is_trend_sve) %>% 
                                     pull(hsapiens_homolog_associated_gene_name)))
table(df_svs_pval$is_mouse_svs)
table(df_sve_pval$is_mouse_sve)

## Fig 6D: SVS results statistics
# SVS gene ranking
df = df_svs_pval %>%
  mutate(rank = row_number())

m1 <- ggplot(df, aes(x = rank, y = pvalue_hsic.ir)) +
  geom_point(aes(color = is_mouse_svs)) + 
  geom_text_repel(
    data = df %>% arrange(pvalue_hsic.ir) %>% head(n = 15), 
    force = 10, force_pull = 1,
    mapping = aes(label = gene, color = is_mouse_svs), nudge_y = 10, nudge_x = 0.1,
    fontface = "italic"
  )+ 
  labs(x = 'Rank (by occurrence)', y = 'Minimum HSIC-IR p-value', 
       title = 'SVS genes in human DLPFC',
       color = 'SVS in mouse') +
  scale_y_continuous(transform = c('log10', 'reverse')) + 
  scale_x_log10() + 
  scale_color_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'darkgray')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    legend.position = 'inside',
    legend.position.inside = c(0.85, 0.85)
  )
m1

ggsave(sprintf("%s/figures/main_svs_pval.png", res_dir), m1, width = 4, height = 3.5)
ggsave(sprintf("%s/figures/main_svs_pval.pdf", res_dir), m1, width = 4, height = 3.5)

## Fig 6B-C and Fig S5A-C: Expression and SV test power
# prepare the data
df_seq_depth <- data.frame(
  sample = "mouse",
  n_svs = sum(mm_sv_cbs$padj_hsic.ir < 0.01), # Replace '-' with '.'
  n_sve = sum(mm_sv_cbs$padj_hsic.gc < 0.01), # Replace '-' with '.'
  n_gene_tested = nrow(mm_sv_cbs),
  sum_count_avg = sum(mm_sv_cbs$count_avg, na.rm = TRUE),
  median_count_avg = median(mm_sv_cbs$count_avg, na.rm = TRUE),
  median_pct_spot_on = median(mm_sv_cbs$pct_spot_on, na.rm = TRUE)
)
sample_list <- c(
  "151507", "151508", "151509", "151510", 
  "151669", "151670", "151671", "151672", 
  "151673", "151674", "151675", "151676"
)
# Loop through each sample
for (sample_id in sample_list) {
  # Read the CSV file
  df <- read.csv(file.path(res_dir, "sv_results", paste0(sample_id, ".sv_combined_", date, ".csv")))
  
  # Create a temporary data frame
  temp_df <- data.frame(
    sample = sample_id,
    n_svs = sum(df$padj_hsic.ir < 0.01, na.rm = TRUE), # Replace '-' with '.'
    n_sve = sum(df$padj_hsic.gc < 0.01), # Replace '-' with '.'
    n_gene_tested = nrow(df),
    sum_count_avg = sum(df$count_avg, na.rm = TRUE),
    median_count_avg = median(df$count_avg, na.rm = TRUE),
    median_pct_spot_on = median(df$pct_spot_on, na.rm = TRUE)
  )
  
  # Append the temporary data frame to the main data frame
  df_seq_depth <- rbind(df_seq_depth, temp_df)
}

# main, depth vs n_svs
m2.1 <- ggplot(df_seq_depth %>% filter(sample != 'mouse'), 
       aes(x = sum_count_avg, y = n_svs)) +
  geom_smooth(method = 'lm', color = 'darkgray', fullrange=TRUE) +
  geom_point() +
  geom_text(data = data.frame(x = 3500, y = 400), 
            aes(x = x, y = y), label = 'Human DLPFC', color = 'black') +
  geom_point(data = df_seq_depth %>% filter(sample == 'mouse'), 
             aes(x = sum_count_avg, y = n_svs), color = '#8B0000') +
  geom_text(data = data.frame(x = 7000, y = 760), 
            aes(x = x, y = y), label = 'Mouse CBS', color = '#8B0000') +
  labs(x = 'Avg total UMI per spot', title = 'Number of SVS genes', y = '') +
  # scale_x_log10() +
  theme_classic() +
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial')
  )

# main, depth vs n_sve
m2.2 <- ggplot(df_seq_depth %>% filter(sample != 'mouse'), 
               aes(x = sum_count_avg, y = n_sve)) +
  geom_smooth(method = 'lm', color = 'darkgray', fullrange=TRUE) +
  geom_point() +
  geom_text(data = data.frame(x = 3500, y = 7500), 
            aes(x = x, y = y), label = 'Human DLPFC', color = 'black') +
  geom_point(data = df_seq_depth %>% filter(sample == 'mouse'), 
             aes(x = sum_count_avg, y = n_sve), color = '#8B0000') +
  geom_text(data = data.frame(x = 7000, y = 4000), 
            aes(x = x, y = y), label = 'Mouse CBS', color = '#8B0000') +
  labs(x = 'Avg total UMI per spot', title = 'Number of SVE genes', y = '') +
  # scale_x_log10() +
  theme_classic() +
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial')
  )

m2 <- cowplot::plot_grid(m2.1, m2.2, ncol = 1, align = 'hv')
m2

ggsave(sprintf("%s/figures/main_nsig_dataset.png", res_dir), m2, width = 2.5, height = 4.5)
ggsave(sprintf("%s/figures/main_nsig_dataset.pdf", res_dir), m2, width = 2.5, height = 4.5)

# Fig S5A: depth vs proportion of SVS
s1.1 <- ggplot(df_seq_depth %>% filter(sample != 'mouse'), 
               aes(x = sum_count_avg, y = n_svs / n_gene_tested)) +
  geom_smooth(method = 'lm', color = 'darkgray', fullrange=TRUE) +
  geom_point() +
  geom_text(data = data.frame(x = 3500, y = 0.1), 
            aes(x = x, y = y), label = 'Human DLPFC', color = 'black') +
  geom_point(data = df_seq_depth %>% filter(sample == 'mouse'), 
             aes(x = sum_count_avg, y = n_svs / n_gene_tested), color = '#8B0000') +
  geom_text(data = data.frame(x = 7000, y = 0.2), 
            aes(x = x, y = y), label = 'Mouse CBS', color = '#8B0000') +
  labs(x = 'Avg total UMI per spot', title = 'Prop of SVS genes', y = '') +
  # lims(y = c(0, 1.0)) +
  theme_classic() +
  theme(
    text=element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial')
  )

# sup, depth vs proportion of SVE
s1.2 <- ggplot(df_seq_depth %>% filter(sample != 'mouse'), 
      aes(x = sum_count_avg, y = n_sve / n_gene_tested)) +
  geom_smooth(method = 'lm', color = 'darkgray', fullrange=TRUE) +
  geom_point() +
  geom_text(data = data.frame(x = 3500, y = 0.7), 
           aes(x = x, y = y), label = 'Human DLPFC', color = 'black') +
  geom_point(data = df_seq_depth %>% filter(sample == 'mouse'), 
            aes(x = sum_count_avg, y = n_sve / n_gene_tested), color = '#8B0000') +
  geom_text(data = data.frame(x = 7000, y = 0.9), 
           aes(x = x, y = y), label = 'Mouse CBS', color = '#8B0000') +
  labs(x = 'Avg total UMI per spot', title = 'Prop of SVE genes', y = '') +
  lims(y = c(0, 1.0)) +
  theme_classic() +
  theme(
    text=element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial')
  )

s1 <- cowplot::plot_grid(s1.1, s1.2, ncol = 1, align = 'hv')
s1
ggsave(sprintf("%s/figures/sup_propsig_dataset.png", res_dir), s1, width = 2.5, height = 4.5)

## Fig 5d: Overlaps between human and mouse SV genes
df <- df_svs_pval %>% 
  mutate(occurrence = ifelse(occurrence > 4, '>4', as.character(occurrence)),
         occurrence = factor(occurrence, levels = c('1', '2', '3', '4', '>4'))) %>%
  count(is_mouse_svs, occurrence)
m3.1 <- ggplot(df, aes(x = occurrence, y = n, group = is_mouse_svs, fill = is_mouse_svs)) +
  geom_bar(stat="identity", position = "fill") + 
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'darkgray')) +
  theme_classic() + 
  labs(x = 'Occurrence', y = 'Proportion', title = 'Overlap with mouse SVS', fill = '') + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial'),
    legend.position = 'bottom'
  )
df <- df_sve_pval %>% 
  # mutate(occurrence = factor(occurrence, levels = c(1:12))) %>%
  mutate(occurrence = ifelse(occurrence < 7, '<7', as.character(occurrence)),
         occurrence = factor(occurrence, levels = c('<7', '7', '8', '9', '10', '11', '12'))) %>%
  count(is_mouse_sve, occurrence)
m3.2 <- ggplot(df, aes(x = occurrence, y = n, group = is_mouse_sve, fill = is_mouse_sve)) +
  geom_bar(stat="identity", position = "fill") + 
  geom_text(aes(label = n), position = position_fill(vjust = 0.5)) + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'darkgray')) +
  theme_classic() + 
  labs(x = 'Occurrence', y = 'Proportion', title = 'Overlap with mouse SVE', fill = '') + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial'),
    legend.position = 'bottom'
  ) + 
  guides(fill = guide_legend(nrow = 1, byrow = TRUE, ))
m3 <- cowplot::plot_grid(m3.1, m3.2, ncol = 2, align = 'hv', rel_widths = c(3, 4))
m3
ggsave(sprintf("%s/figures/main_sv_occurence.png", res_dir), m3, width = 7, height = 4)

## Fig S5C: SVS and SVE occurrence
# # merge dv_svs_pval and df_sve_pval occurrences size
# df_occur <- rbind(
#   df_svs_pval %>% group_by(occurrence) %>% summarize(size = n()) %>% mutate(type = 'SVS'),
#   df_sve_pval %>% group_by(occurrence) %>% summarize(size = n()) %>% mutate(type = 'SVE')
# ) %>%
#   # aggregate counts for occurrence > 8 as '8+'
#   # mutate(occurrence = ifelse(occurrence > 8, '8+', as.character(occurrence))) %>%
#   mutate(occurrence = factor(occurrence, levels = c(1:12))) %>%
#   group_by(occurrence, type) %>%
#   summarize(size = sum(size)) %>%
#   group_by(type) %>%
#   mutate(prop = size / sum(size))

# # occurrence vs n_genes
# s1 <- ggplot(df_occur, aes(x = occurrence, y = prop, fill = type)) +
#   geom_bar(stat = 'identity', position = 'dodge') +
#   geom_text(aes(label = size), position = position_dodge(width = 0.9), vjust = -0.5) +
#   scale_fill_manual(values = c('SVS'= '#FF9F1C', 'SVE'= '#2EC4B6')) +
#   labs(x = 'Occurrence', y = 'Proportion', fill = 'Gene set', title = 'Human DLPFC SV genes') +
#   lims(y = c(0, 0.65)) +
#   theme_classic() +
#   theme(
#     text = element_text(size = 12, family = 'Arial'),
#     plot.title = element_text(size = 12, family = 'Arial'),
#     legend.position = 'inside'
#   )
# s1
# ggsave(sprintf("%s/figures/sup_sv_ngenes.png", res_dir), s1, width = 6, height = 4)

# occurrence vs avg gene expr, SVS
df <- df_svs_pval %>% 
  mutate(occurrence = ifelse(occurrence > 4, '>4', as.character(occurrence)),
         occurrence = factor(occurrence, levels = c('1', '2', '3', '4', '>4')))
s2.1 <- ggplot(df, aes(x = occurrence, y = count_avg, group = interaction(is_mouse_svs, occurrence), 
                        fill = is_mouse_svs)) +
  geom_boxplot(position = 'dodge') + 
  stat_compare_means(label = 'p.signif') + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'darkgray')) +
  scale_y_log10() + 
  theme_classic() + 
  labs(x = 'Occurrence', y = 'Avg gene count per spot', fill = 'SVS in mouse', title = 'Human DLPFC SVS') + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial'),
    legend.position = 'bottom'
  )

# occurrence vs avg gene expr, SVE
df <- df_sve_pval %>%
  mutate(occurrence = ifelse(occurrence < 7, '<7', as.character(occurrence)),
         occurrence = factor(occurrence, levels = c('<7', '7', '8', '9', '10', '11', '12')))
s2.2 <- ggplot(df, aes(x = occurrence, y = count_avg, group = interaction(is_mouse_sve, occurrence), 
               fill = is_mouse_sve)) +
  geom_boxplot(position = 'dodge') + 
  stat_compare_means(label = 'p.signif') + 
  scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'darkgray')) +
  scale_y_log10() + 
  theme_classic() + 
  labs(x = 'Occurrence', y = 'Avg gene count per spot', fill = 'SVE in mouse', title = 'Human DLPFC SVE') + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(size = 14, family = 'Arial'),
    legend.position = 'bottom'
  )

s2 <- cowplot::plot_grid(s2.1, s2.2, ncol = 2, align = 'hv', rel_widths = c(3, 4))
s2
ggsave(sprintf("%s/figures/sup_sv_expr.png", res_dir), s2, width = 8, height = 4)

## Fig 6C: Test agreements between replicates
df_stat <- read.csv(paste0(res_dir, "figures/source_data/sv_agreements_stat_", date, ".csv")) %>%
  filter(stat != 'Perplexity') %>%
  mutate(
    comp = factor(comp, levels = c('replicates (n=6)', 'positions (n=12)', 'individuals (n=48)')),
    group = factor(group, levels = c('Low', 'Mid', 'High'), labels = c('(0, 0.5]', '(0.5, 1]', '(1, Inf)')),
    stat = factor(stat, levels = c('SPARK-X', 'HSIC-GC', 'HSIC-IC', 'HSIC-IR'))
  )

m4 <- ggplot(df_stat, aes(x = stat, y = spearmannr, group = interaction(group, stat), fill = group)) + 
  facet_wrap(~comp) + 
  # geom_point(position = position_dodge(width = 0.5), alpha = 0.5) + 
  stat_compare_means(label = 'p.signif', label.y = 0.95) +
  geom_boxplot() + 
  labs(x = '', y = 'Spearmann R', fill = 'Avg gene-level count per spot', 
       title = 'SV test agreements between sample pairs of different') + 
  scale_fill_manual(values = c('(0, 0.5]' = '#f4b45c', '(0.5, 1]' = '#e54b1d', '(1, Inf)' = '#8c2a2f')) +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    legend.position = 'bottom'
  )
m4
ggsave(sprintf("%s/figures/main_sv_agreements.png", res_dir), m4, width = 6, height = 4)
ggsave(sprintf("%s/figures/main_sv_agreements.pdf", res_dir), m4, width = 5, height = 4)

## Fig 6E: Evolutionary conservation of SVS and SVE peaks
# load conservation data
df_conservation = c()
for (group in c('svs_all', 'svens_all', 'neg_all', 'svs_shared', 'svens_shared')){
  for (stat in c('phylop', 'phastcons')){
    df <- read_table(
      paste0(res_dir, 'events/conservation/', group, '.', stat, '.tab'), 
      col_names = c('id', 'size', 'covered', 'sum', 'mean0', 'mean')
    ) %>%
      mutate(group = group, stat = stat)
    df_conservation <- rbind(df_conservation, df)
  }
}

df_conservation <- df_conservation %>%
  mutate(group = factor(
    group, 
    levels = c('neg_all', 'svens_all', 'svs_all', 'svens_shared', 'svs_shared'),
    labels = c('Negative controls', 'SVENS', 'SVS', 'SVENS shared with mouse', 'SVS shared with mouse')
  ),
    stat = factor(stat, levels = c('phylop', 'phastcons'))
  )

# boxplot, PhyloP score
m5.1 <- df_conservation %>% filter(stat == 'phylop') %>%
  ggplot(aes(x = group, y = mean, group = group, fill = group)) +
  stat_compare_means(label = 'p.signif', comparisons = list(
    c('Negative controls', 'SVENS'), c('Negative controls', 'SVS'), c('SVENS', 'SVS'),
    c('SVENS', 'SVENS shared with mouse'),
    c('SVS', 'SVS shared with mouse')
  )) +
  geom_boxplot() + 
  theme_classic() + 
  labs(x = '', y = 'Average PhyloP score', fill = 'Gene set', title = '') + 
  scale_fill_manual(
    # labels = c('Negative controls', 'SVENS', 'SVS', 'SVENS shared with mouse', 'SVS shared with mouse'),
    values = c('Negative controls' = 'darkgray', 'SVENS' = '#2EC4B6', 'SVS' = '#FF9F1C',
               'SVENS shared with mouse' = '#007B99', 'SVS shared with mouse' = '#cc6600')) +
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'none'
  )

# boxplot, PhastCons score
m5.2 <- df_conservation %>% filter(stat == 'phastcons') %>%
  ggplot(aes(x = group, y = mean, group = group, fill = group)) +
  stat_compare_means(label = 'p.signif', comparisons = list(
    c('Negative controls', 'SVENS'), c('Negative controls', 'SVS'), c('SVENS', 'SVS'),
    c('SVENS', 'SVENS shared with mouse'),
    c('SVS', 'SVS shared with mouse')
  )) +
  geom_boxplot() + 
  theme_classic() + 
  labs(x = '', y = 'Average PhastCons score', fill = 'Gene set', title = 'TREND region conservation') + 
  scale_fill_manual(
    # labels = c('Negative controls', 'SVENS', 'SVS', 'SVENS shared with mouse', 'SVS shared with mouse'),
    values = c('Negative controls' = 'darkgray', 'SVENS' = '#2EC4B6', 'SVS' = '#FF9F1C',
               'SVENS shared with mouse' = '#007B99', 'SVS shared with mouse' = '#cc6600')) +
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    axis.text.x = element_blank(),
    legend.position = 'right'
  )

# m5 <- cowplot::plot_grid(m5.1, m5.2, ncol = 2, align = 'hv', rel_widths = c(1, 1))
m5 <- cowplot::plot_grid(m5.1, m5.2, ncol = 2, align = 'h', rel_widths = c(1, 2.5))
m5
ggsave(sprintf("%s/figures/main_conservation.png", res_dir), m5, width = 6, height = 4)
ggsave(sprintf("%s/figures/main_conservation.pdf", res_dir), m5, width = 6, height = 4)

ggplot(df_conservation, aes(x = group, y = mean, group = group, fill = group)) +
  facet_wrap(~stat, scale = 'free') +
  stat_compare_means(label = 'p.signif', comparisons = list(
    c('Negative controls', 'SVS'), c('Negative controls', 'SVENS'), c('SVS', 'SVENS')
  )) + 
  geom_boxplot() + 
  theme_classic() + 
  labs(x = '', y = 'Average PhastCons score', fill = 'Gene set')

## Fig 6F and Fig S5D: SV functional enrichment
# load enrichment results of shared SV genes
df_kegg <- read_csv(sprintf("%s/figures/source_data/kegg_shared_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('SVENS shared with mouse', 'SVS shared with mouse'))
  ) %>%
  arrange(precision_diff, name) %>%
  # keep the top 8 terms with the largest precision difference
  dplyr::slice((n()-15):n())

df_reac <- read_csv(sprintf("%s/figures/source_data/reac_shared_top20.csv", res_dir)) %>%
  mutate(
    name = recode(
      name,
      'Diseases of signal transduction by growth factor receptors and second messengers' = 'Diseases of signal transduction by growth factor receptors',      
    ),
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('SVENS shared with mouse', 'SVS shared with mouse'))
    ) %>%
  arrange(precision_diff, name) %>%
  # keep the top 8 terms with the largest precision difference
  dplyr::slice((n()-15):n())

m6.1 <- ggplot(df_kegg, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
  geom_bar(stat='identity', position='dodge') + 
  labs(y = '', fill = 'Gene set', x = 'KEGG') + 
  coord_flip() +
  scale_fill_manual(values=c('SVENS shared with mouse' = '#007B99', 'SVS shared with mouse' = '#cc6600')) + 
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.position='none',
    legend.position.inside = c(0.3,0.85)
  )
m6.2 <- ggplot(df_reac, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = 'Precision (proportion of term genes)', fill = 'Gene set',
       x = 'Reactome') +
  coord_flip() +
  scale_fill_manual(values=c('SVENS shared with mouse' = '#007B99', 'SVS shared with mouse' = '#cc6600')) + 
  theme_classic() +
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.position='inside',
    legend.position.inside = c(-0.8, 1.25),
    legend.background = element_blank()
  )
m6 <- cowplot::plot_grid(m6.1, m6.2, nrow = 2, align = 'hv', rel_heights = c(1,1))
m6

ggsave(sprintf("%s/figures/main_enrich_new.png", res_dir), m6, width = 8.5, height = 4)
ggsave(sprintf("%s/figures/main_enrich_new.pdf", res_dir), m6, width = 8.5, height = 4)

## Fig S5D: all human sv genes
df_kegg <- read_csv(sprintf("%s/figures/source_data/kegg_all_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced', 
                                'Spatially variably spliced'))
  ) %>%
  arrange(precision_diff, name) %>%
  # keep the top 10 terms with the largest precision difference
  dplyr::slice((n()-19):n())

df_reac <- read_csv(sprintf("%s/figures/source_data/reac_all_top20.csv", res_dir)) %>%
  mutate(
    name_cat = factor(name, levels = unique(name[order(precision_diff, decreasing = FALSE)])),
    geneset = factor(geneset, levels = c('SVE', 'SVS'), 
                     labels = c('Spatially variably expressed but not spliced', 
                                'Spatially variably spliced'))
  ) %>%
  arrange(precision_diff, name) %>%
  # keep the top 10 terms with the largest precision difference
  dplyr::slice((n()-19):n())

s3.1 <- ggplot(df_kegg, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) + 
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
    legend.position='none',
    legend.position.inside = c(0.3,0.85)
  )
s3.2 <- ggplot(df_reac, aes(x = name_cat, y = precision, group = geneset, fill = geneset)) +
  geom_bar(stat='identity', position='dodge') +
  labs(y = 'Precision (proportion of term genes)', fill = 'Gene set',
       x = 'Reactome') +
  coord_flip() +
  # scale_fill_manual(values=c('SVS'= '#FF9F1C', 'SVENS'= '#2EC4B6')) + 
  scale_fill_manual(values=c('Spatially variably spliced'= '#FF9F1C', 
                             'Spatially variably expressed but not spliced'= '#2EC4B6')) + 
  theme_classic() +
  theme(
    text=element_text(size = 12, family = 'Arial'),
    axis.text=element_text(size = 12, family = 'Arial', color = 'black'),
    # axis.text.x=element_text(angle = 90, hjust = 1, vjust = 0.5),
    plot.title=element_text(hjust = 0.5, size = 16, family = 'Arial'),
    legend.position='inside',
    legend.position.inside = c(-1.3, 1.2),
    legend.background = element_blank()
  )
s3 <- cowplot::plot_grid(s3.1, s3.2, nrow = 2, align = 'hv', rel_heights = c(1,1))
s3

ggsave(sprintf("%s/figures/sup_enrich_sv_all_new.png", res_dir), s3, width = 8, height = 4.5)

## Fig 6G: RBP DU regulation
# load ranked DU results in human
df_du_rbp <- read.csv(paste0(res_dir, "trend_du_rbp_dlpfc_", date, ".csv")) %>%
  mutate(rank = row_number(),
         label = paste(gene, covariate, sep = ':')) 

# main, SEPTIN 8
m7 <- ggplot(df_du_rbp, aes(x = rank, y = pvalue_glmm)) + 
  geom_point(alpha = 0.1) + 
  geom_point(
    df_du_rbp %>% filter(gene == 'SEPTIN8') %>% top_n(7, -rank) %>%
      mutate(highlight = ifelse(covariate %in% c('QKI', 'CELF2', 'CELF3', 'G3BP2', 'TNRC6C'), 'TRUE', 'FALSE')),
    mapping = aes(x = rank, y = pvalue_glmm, color = highlight),
  ) +
  geom_text_repel(
    df_du_rbp %>% filter(gene == 'SEPTIN8') %>% top_n(7, -rank) %>%
      mutate(highlight = ifelse(covariate %in% c('QKI', 'CELF2', 'CELF3', 'G3BP2', 'TNRC6C'), 'TRUE', 'FALSE')),
    mapping = aes(x = rank, y = pvalue_glmm, label = label, color = highlight),
    fontface = 'italic'
  ) +
  labs(x = 'Rank (by occurrence)', y = 'Minimum GLMM p-value', 
       title = 'Potential SEPTIN8 regulators in human',
       color = 'Significant in mouse') +
  scale_y_continuous(transform = c('log10', 'reverse')) + 
  scale_x_log10() + 
  scale_color_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'black')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    legend.position = 'inside',
    legend.position.inside = c(0.7, 0.75),
    legend.background = element_blank()
  )
m7

ggsave(sprintf("%s/figures/main_du_septin8.png", res_dir), m7, width = 4, height = 3.5)
ggsave(sprintf("%s/figures/main_du_septin8.pdf", res_dir), m7, width = 4, height = 3.5)

# Fig S5F: MAP4
s4 <- ggplot(df_du_rbp, aes(x = rank, y = pvalue_glmm)) + 
  geom_point(alpha = 0.1) + 
  geom_point(
    df_du_rbp %>% filter(gene == 'MAP4') %>% top_n(7, -rank) %>%
      mutate(highlight = ifelse(covariate %in% c('CELF4', 'ARPP21', 'QKI', 'KHDRBS3', 'RBFOX1'), 'TRUE', 'FALSE')),
    mapping = aes(x = rank, y = pvalue_glmm, color = highlight),
  ) +
  geom_text_repel(
    df_du_rbp %>% filter(gene == 'MAP4') %>% top_n(7, -rank) %>%
      mutate(highlight = ifelse(covariate %in% c('CELF4', 'ARPP21', 'QKI', 'KHDRBS3', 'RBFOX1'), 'TRUE', 'FALSE')),
    mapping = aes(x = rank, y = pvalue_glmm, label = label, color = highlight),
    fontface = 'italic'
  ) +
  labs(x = 'Rank (by occurrence)', y = 'Minimum GLMM p-value', 
       title = 'Potential MAP4 regulators in human',
       color = 'Significant in mouse') +
  scale_y_continuous(transform = c('log10', 'reverse')) + 
  scale_x_log10() + 
  scale_color_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'black')) +
  theme_classic() + 
  theme(
    text=element_text(size = 12, family = 'Arial'),
    plot.title = element_text(size = 12, family = 'Arial'),
    legend.position = 'inside',
    legend.position.inside = c(0.25, 0.85),
    legend.background = element_blank()
  )
s4

ggsave(sprintf("%s/figures/sup_du_map4.png", res_dir), s4, width = 4, height = 3.5)

