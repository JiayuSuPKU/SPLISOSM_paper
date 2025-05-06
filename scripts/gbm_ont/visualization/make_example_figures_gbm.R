library(tidyverse)
library(GenomicFeatures)
library(Rsamtools)
library(GenomicAlignments)
library(rtracklayer)
library(plotgardener)
library(org.Hs.eg.db)
library(scales)

extrafont::loadfonts()

res_dir_ont <- "~/Projects/SPLISOSM_paper/results/gbm_ont/"
data_dir_trend <- "~/Projects/SPLISOSM_paper/data/gbm_visium_cell_24/"
res_dir_trend <- "/Users/jysumac/Projects/SPLISOSM_paper/results/gbm_visium/"
res_dir_dlpfc <- "~/Projects/SPLISOSM_paper/results/human_dlpfc/"
date_ont <- "0104"
date_vis <- "0310"
date_dlpfc <- "0123"

### Example (1): FTH, Fig 7I-M
## load ONT gtf and SR gtf
svs_ont <- import.gff(sprintf("%s/transcripts/DMG_3/DMG_3.svs_iso.oa.gtf", res_dir_ont))
# add transcript name to exons
txid2name <- svs_ont %>% 
  filter(type == 'transcript') %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>% 
  deframe()
svs_ont$transcript_name <- txid2name[svs_ont$transcript_id]
svs_ont$transcript_id_old <- svs_ont$transcript_id
svs_ont$transcript_id <- svs_ont$transcript_name

# add gene name and id to exons
svs_ont$gene_name <- svs_ont$transcript_name %>% strsplit("[_-]") %>% sapply(`[[`, 1)
# remove suffixes from gene_id
svs_ont$gene_id <- svs_ont$gene_id %>% str_remove("\\.\\d+$")
genename2id <- svs_ont %>% 
  as.data.frame() %>%
  dplyr::select(gene_name, gene_id) %>% 
  distinct() %>% 
  deframe()
svs_ont$gene_id <- genename2id[svs_ont$gene_name]
svs_ont$gene_symbol <- svs_ont$gene_name

## load short-read gtf
svs_sr <- import.bed(sprintf("%s/events/ZH8811Abulk/svs.exon.bed", res_dir_trend)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1)
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order)) # Create the new name

## load SR bam coverage
cov <- BamFile(sprintf("%s/bam/selected.ZH8811Abulk.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')

## load spatial expression data
expr_ont <- read.csv(
  sprintf('%s/figures/source_data/FTH1.ONT.DMG_3.csv', res_dir_ont)
)
expr_sr <- read.csv(
  sprintf('%s/figures/source_data/FTH1.SR.ZH8811Abulk.csv', res_dir_ont)
) %>%
  mutate(isoform_name = ifelse(isoform == 'FTH1:chr11:61964476-61964818:-1', 'FTH1-Event1', 'FTH1-Event2'))
# expression of hypoxia genes
hypoxia_ont <- read.csv(
  sprintf('%s/figures/source_data/FTH1.ONT.DMG_3.sup.csv', res_dir_ont)
)
# expression of hemoglobin genes
hb_sr <- read.csv(
  sprintf('%s/figures/source_data/HB.SR.ZH8811Abulk.csv', res_dir_ont)
)

## FTH1 transcripts
gene_to_plot <- 'FTH1'

# create ONT assembly object
subset_ont_gr <- svs_ont %>% 
  filter(gene_name == gene_to_plot)
subset_ont_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_ont_gr),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# create SR event bed
subset_sr_bed <- svs_sr %>% 
  filter(gene == gene_to_plot)

# genomic coordinates
chrom <- seqnames(subset_ont_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
chromstart <- min(start(subset_ont_gr))
chromend <- max(end(subset_ont_gr))
# chromstart <- min(subset_sr_bed$start)
# chromend <- max(subset_sr_bed$end) + 1000
strand <- strand(subset_ont_gr)[1] %>% as.character()
n_isoforms <- subset_ont_gr$transcript_name %>% unique() %>% length()

png(
  sprintf('%s/figures/FTH1.png', res_dir_ont), 
  width = 7.25, height = 3, units = 'in', res = 300
)
pdf(
  sprintf('%s/figures/FTH1.pdf', res_dir_ont), 
  width = 7.25, height = 3, useDingbats = FALSE
)

# create new page
pageCreate(width = 7.25, height = 3, default.units = "inches", showGuides = FALSE)

# plot and place ONT transcripts
height_ont = 0.3 * n_isoforms
plotTranscripts(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as,
  labels = 'transcript', fill = '#0066CC', 
  stroke = 0.5,
  x = 1, y = 0, width = 6, height = height_ont,
  just = c("left", "top"), default.units = "inches",
)

# plot ONT dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = 0.1, y1 = height_ont
)
plotText(
  label = 'DMG3 (ONT)', rot = 90,
  fontsize = 8,
  x = 0.75, y = (height_ont + 0.05) / 2,
  just = "center", default.units = "inches"
)

# plot SR bam coverage
plotSignal(
  cov,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, x = 1, y = height_ont + 0.2, 
  width = 6, height = 0.3,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#8E44AD', fill = '#8E44AD',
  # linecolor = '#E23D28', fill = '#E23D28',
  scale = T, fontsize = 8,
)

# plot and place SR events
plotRanges(
  subset_sr_bed,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, collapse = FALSE, fill = '#E23D28',
  x = 1, y = height_ont + 0.6, width = 6, height = 0.1,
)
for (peak in subset_sr_bed$peak_name){
  start = subset_sr_bed %>% filter(peak_name == peak) %>% pull(start)
  end = subset_sr_bed %>% filter(peak_name == peak) %>% pull(end)
  plotText(
    label = peak, fontsize = 8,
    x = 1 + (start - chromstart) / (chromend - chromstart) * 6,
    y = height_ont + 0.8,
    just = "left", default.units = "inches"
  )
}

# plot SR dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = height_ont + 0.1, y1 = height_ont + 0.9
)
plotText(
  label = 'ZH8811Abulk (SR)', rot = 90,
  fontsize = 8,
  x = 0.75, y = height_ont + 0.5,
  just = "center", default.units = "inches"
)

# Plot genome scale
plotGenomeLabel(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, fontsize = 8,
  x = 1, y = height_ont + 0.9,
  length = 6, default.units = "inches",
)

# Plot strand info
label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
label = paste0(gene_to_plot, ', ', label)
plotText(
  label = label, fontsize = 8,
  x = 1, y = 0.2,
  just = "left", default.units = "inches"
)

# pageGuideHide()
dev.off()

## Plot spatial region annotation
# hypoxia gene expression
covar_to_plot <- c('CA9', 'NDRG1')
p_list <- c()
for (i in seq_along(covar_to_plot)){
  p <- hypoxia_ont %>% filter(
    gene == covar_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(x = -array_row, y = -array_col, color = value)) + 
    geom_point(size = 0.2) +
    labs(color = 'LogExpr', title = sprintf("%s", covar_to_plot[i])) +
    scale_color_distiller(palette = 'Reds', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
m1.1 <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
m1.1

df <- expr_ont %>%
  filter(layer == 'ratios_obs', isoform == 'FTH1-211') %>%
  mutate(group = ifelse(value > 0.1, '>10%', '<=10%'))

df_prop <- df %>%
  group_by(group, region) %>%
  summarise(count = n()) %>%
  mutate(
    prop = count / sum(count)
  )

m1.2 <- df_prop %>%
  ggplot(aes(x = group, y = prop, fill = region)) +
  geom_bar(stat = 'identity') +
  geom_text(
    aes(label = region), position = position_stack(vjust = .5)
  ) +
  geom_text(
    data = data.frame(),
    mapping = aes(
      x = 1.5, y = 1.05, fill = NULL,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$region))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  # scale_fill_brewer(palette = 'Accent') +
  labs(x = 'FTH1-211', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text = element_text(family = "Arial", size = 12),
    # axis.text.x = element_text(angle = 90, hjust = 1),
    legend.position = "none"
  )

m1 <- cowplot::plot_grid(m1.1, m1.2, nrow = 1, rel_widths = c(1, 0.7))
m1

ggsave(
  sprintf('%s/figures/FTH1.ONT.anno.png', res_dir_ont), m1,
  width = 4.5, height = 3.5, units = 'in', dpi = 300
)
ggsave(
  sprintf('%s/figures/FTH1.ONT.anno.pdf', res_dir_ont), m1,
  width = 4, height = 3.5, units = 'in', dpi = 300
)

# SR region annotation
m2.1 <- ggplot(expr_sr, aes(x = -array_row, y = array_col, color = mp)) +
  geom_point(size = 0.8) + 
  theme_minimal() +
  # coord_fixed(ratio = 0.6) +
  labs(title = 'ZH8811Abulk (SR)', color = '') + 
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(family = "Arial", size = 12, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none"
  ) + 
  guides(color = guide_legend(
    override.aes = list(size = 3),
    nrow = 2
  ))

df <- expr_sr %>%
  filter(layer == 'ratios_obs', isoform_name == 'FTH1-Event2') %>%
  mutate(group = ifelse(value > 0.1, '>10%', '<=10%'))

df_prop <- df %>%
  group_by(group, mp) %>%
  summarise(count = n()) %>%
  mutate(
    prop = count / sum(count)
  )

# mp_order <- df_prop %>% pivot_wider(names_from = group, values_from = prop, id_cols = mp) %>%
#   mutate(diff = `>10%` - `<=10%`) %>%
#   arrange(desc(diff)) %>% pull(mp)
# 
# df_prop$mp <- factor(df_prop$mp, levels = mp_order)

m2.2 <- df_prop %>%
  ggplot(aes(x = group, y = prop, fill = mp)) +
  geom_bar(stat = 'identity') +
  geom_text(
    data = df_prop %>% filter(prop > 0.01),
    mapping = aes(label = mp), position = position_stack(vjust = .5)
  ) +
  geom_text(
    data = data.frame(),
    mapping = aes(
      x = 1.5, y = 1.05, fill = NULL,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$mp))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  labs(x = 'FTH1-Event2', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text = element_text(family = "Arial", size = 12),
    legend.position = "none"
  )

## Hemoglobin expression
m2.3 <- hb_sr %>% filter(layer == 'log1p') %>%
  # rename FTH1-Event2 to '>10%'
  mutate(group = ifelse(FTH1.group == 'FTH1-Event2', '>10%', '<=10%')) %>%
  ggplot(aes(x = group, y = value)) + 
  facet_wrap(~gene, ncol = 1, scale = 'free_y') + 
  geom_boxplot(aes(fill = group)) +
  # geom_violin(aes(fill = group)) + 
  stat_compare_means(
    # comparisons = list(c('>10%', '<=10%')),
    label = 'p.signif',
    method = 'wilcox.test', label.x.npc = 0.4, label.y.npc = 0.8) + 
  theme_classic() + 
  labs(x = 'FTH1-Event2', y = 'Log1p expression') + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text = element_text(family = "Arial", size = 12),
    strip.text = element_text(family = "Arial", size = 12, face = 'italic'),
    strip.background = element_blank(),
    legend.position = "none"
  )

m2 <- cowplot::plot_grid(m2.1, m2.2, m2.3, nrow = 1, rel_widths = c(1, 0.8, 0.8))
m2

ggsave(
  sprintf('%s/figures/FTH1.SR.anno.png', res_dir_ont), m2,
  width = 7, height = 3.5, units = 'in', dpi = 300
)

ggsave(
  sprintf('%s/figures/FTH1.SR.anno.pdf', res_dir_ont), m2,
  width = 7, height = 3.5, units = 'in', dpi = 300
)



## ONT spatial expression
# iso_to_plot <- c('FTH1-201', 'FTH1-211', 'FTH1-208', 'FTH1_Iso_1')
iso_to_plot <- c('FTH1-201', 'FTH1-211')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.5) +
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
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.5) +
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
  p_list2 <- c(p_list2, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 1, align = 'hv')
p

ggsave(
  sprintf('%s/figures/FTH1.ONT.spatial.png', res_dir_ont), p,
  width = 10, height = 2.5, units = 'in', dpi = 300
)
ggsave(
  sprintf('%s/figures/FTH1.ONT.spatial.pdf', res_dir_ont), p,
  width = 10, height = 2.5, units = 'in', dpi = 300
)

## SR spatial expression
iso_to_plot <- c('FTH1-Event1', 'FTH1-Event2')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 1, align = 'hv')
p

ggsave(
  sprintf('%s/figures/FTH1.SR.spatial.png', res_dir_ont), p,
  width = 10, height = 2.5, units = 'in', dpi = 300
)
ggsave(
  sprintf('%s/figures/FTH1.SR.spatial.pdf', res_dir_ont), p,
  width = 10, height = 2.5, units = 'in', dpi = 300
)

### Example (2): GFAP, Fig S6J
## load ONT gtf and SR gtf
svs_ont <- import.gff(sprintf("%s/transcripts/DMG_2/DMG_2.svs_iso.oa.gtf", res_dir_ont))
# add transcript name to exons
txid2name <- svs_ont %>% 
  filter(type == 'transcript') %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>% 
  deframe()
svs_ont$transcript_name <- txid2name[svs_ont$transcript_id]
svs_ont$transcript_id_old <- svs_ont$transcript_id
svs_ont$transcript_id <- svs_ont$transcript_name

# add gene name and id to exons
svs_ont$gene_name <- svs_ont$transcript_name %>% strsplit("[_-]") %>% sapply(`[[`, 1)
# remove suffixes from gene_id
svs_ont$gene_id <- svs_ont$gene_id %>% str_remove("\\.\\d+$")
genename2id <- svs_ont %>% 
  as.data.frame() %>%
  dplyr::select(gene_name, gene_id) %>% 
  distinct() %>% 
  deframe()
svs_ont$gene_id <- genename2id[svs_ont$gene_name]
svs_ont$gene_symbol <- svs_ont$gene_name

## load short-read gtf
svs_sr <- import.bed(sprintf("%s/events/ZH916bulk/svs.exon.bed", res_dir_trend)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1)
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order)) # Create the new name

## load SR bam coverage
cov <- BamFile(sprintf("%s/bam/selected.ZH916bulk.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')

## load spatial expression data
expr_ont <- read.csv(
  sprintf('%s/figures/source_data/GFAP.ONT.DMG_2.csv', res_dir_ont)
)
expr_sr <- read.csv(
  sprintf('%s/figures/source_data/GFAP.SR.ZH916bulk.csv', res_dir_ont)
) %>%
  mutate(isoform_name = ifelse(isoform == "GFAP:chr17:44905530-44905896:-1", 
                               'GFAP-Event1', 'GFAP-Event2'))

## GFAP transcripts
gene_to_plot <- 'GFAP'

# create ONT assembly object
subset_ont_gr <- svs_ont %>% 
  filter(gene_name == gene_to_plot) %>%
  filter(transcript_name %in% c('GFAP-230', 'GFAP-203', 'GFAP-214', 'GFAP-205', 'GFAP_Iso_5'))
         
subset_ont_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_ont_gr),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# create SR event bed
subset_sr_bed <- svs_sr %>% 
  filter(gene == gene_to_plot)

# genomic coordinates
chrom <- seqnames(subset_ont_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
chromstart <- min(start(subset_ont_gr))
chromend <- max(end(subset_ont_gr))
# chromstart <- min(subset_sr_bed$start)
# chromend <- max(subset_sr_bed$end) + 1000
strand <- strand(subset_ont_gr)[1] %>% as.character()
n_isoforms <- subset_ont_gr$transcript_name %>% unique() %>% length()

png(
  sprintf('%s/figures/GFAP.png', res_dir_ont), 
  width = 7.25, height = 3, units = 'in', res = 300
)

# create new page
pageCreate(width = 7.25, height = 3, default.units = "inches")

# plot and place ONT transcripts
height_ont = 0.3 * n_isoforms
plotTranscripts(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as,
  labels = 'transcript', fill = '#0066CC', 
  stroke = 0.5,
  x = 1, y = 0, width = 6, height = height_ont,
  just = c("left", "top"), default.units = "inches",
)

# plot ONT dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = 0.45, y1 = height_ont
)
plotText(
  label = 'DMG2 (ONT)', rot = 90,
  fontsize = 8,
  x = 0.75, y = (height_ont + 0.05) / 2,
  just = "center", default.units = "inches"
)

# plot SR bam coverage
plotSignal(
  cov,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, x = 1, y = height_ont + 0.2, 
  width = 6, height = 0.3,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#8E44AD', fill = '#8E44AD',
  # linecolor = '#E23D28', fill = '#E23D28',
  scale = T, fontsize = 8,
)

# plot and place SR events
plotRanges(
  subset_sr_bed,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, collapse = FALSE, fill = '#E23D28',
  x = 1, y = height_ont + 0.6, width = 6, height = 0.1,
)
for (peak in subset_sr_bed$peak_name){
  start = subset_sr_bed %>% filter(peak_name == peak) %>% pull(start)
  end = subset_sr_bed %>% filter(peak_name == peak) %>% pull(end)
  plotText(
    label = peak, fontsize = 8,
    x = 1 + (start - chromstart) / (chromend - chromstart) * 6,
    y = height_ont + 0.8,
    just = "left", default.units = "inches"
  )
}

# plot SR dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = height_ont + 0.1, y1 = height_ont + 0.9
)
plotText(
  label = 'ZH916bulk (SR)', rot = 90,
  fontsize = 8,
  x = 0.75, y = height_ont + 0.5,
  just = "center", default.units = "inches"
)

# Plot genome scale
plotGenomeLabel(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, fontsize = 8,
  x = 1, y = height_ont + 0.9,
  length = 6, default.units = "inches",
)

# Plot strand info
label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
label = paste0(gene_to_plot, ', ', label)
plotText(
  label = label, fontsize = 8,
  x = 1, y = 0.5,
  just = "left", default.units = "inches"
)

pageGuideHide()
dev.off()

## Plot spatial region annotation
# ONT region annotation
m3.1 <- ggplot(expr_ont, aes(x = -array_row, y = -array_col, color = anno_hzy)) +
  geom_point(size = 1.4) + 
  scale_color_brewer(palette = 'Accent') +
  theme_minimal() +
  coord_fixed(ratio = 0.7) +
  labs(title = 'DMG2 (ONT)', color = '') + 
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(family = "Arial", size = 12, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  ) + 
  guides(color = guide_legend(override.aes = list(size = 3)))

df <- expr_ont %>%
  filter(layer == 'ratios_obs', isoform == 'GFAP-230') %>%
  mutate(group = ifelse(value > 0.85, '>85%', '<=85%'))

m3.2 <- df %>%
  ggplot(aes(x = group)) +
  geom_bar(aes(fill = anno_hzy), position = 'fill', color = 'black') +
  geom_text(
    data = data.frame(),
    aes(
      x = 1.5, y = 1.05,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$anno_hzy))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  scale_fill_brewer(palette = 'Accent') +
  labs(x = 'GFAP-230', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text = element_text(family = "Arial", size = 12),
    legend.position = "none"
  )
m3 <- cowplot::plot_grid(m3.1, m3.2, nrow = 1, rel_widths = c(1, 0.6))
m3

ggsave(
  sprintf('%s/figures/GFAP.ONT.anno.png', res_dir_ont), m3,
  width = 5, height = 4, units = 'in', dpi = 300
)

# SR region annotation
m4.1 <- ggplot(expr_sr, aes(x = -array_row, y = array_col, color = mp)) +
  geom_point(size = 0.8) + 
  theme_minimal() +
  coord_fixed(ratio = 0.6) +
  labs(title = 'ZH916bulk (SR)', color = '') + 
  theme(
    text = element_text(family = "Arial", size = 12),
    plot.title = element_text(family = "Arial", size = 12, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "bottom"
  ) + 
  guides(color = guide_legend(
    override.aes = list(size = 3),
    nrow = 2
  ))

df <- expr_sr %>%
  filter(layer == 'ratios_obs', isoform_name == 'GFAP-Event1') %>%
  mutate(group = ifelse(value > 0.85, '>85%', '<=85%'))

m4.2 <- df %>%
  ggplot(aes(x = group)) +
  geom_bar(aes(fill = mp), position = 'fill', color = 'black') +
  geom_text(
    data = data.frame(),
    aes(
      x = 1.5, y = 1.05,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$mp))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  labs(x = 'GFAP-Event1', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text = element_text(family = "Arial", size = 12),
    legend.position = "none"
  )
m4 <- cowplot::plot_grid(m4.1, m4.2, nrow = 1, rel_widths = c(1, 0.5))
m4

ggsave(
  sprintf('%s/figures/GFAP.SR.anno.png', res_dir_ont), m4,
  width = 6, height = 4, units = 'in', dpi = 300
)

## ONT spatial expression
iso_to_plot <- c('GFAP-230', 'GFAP-214', 'GFAP_Iso_5', 'GFAP-203')

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.7) +
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
  p_list2 <- c(p_list2, list(p))
}

p <- cowplot::plot_grid(plotlist = p_list2, nrow = 1, align = 'hv')
p

ggsave(
  sprintf('%s/figures/GFAP.ONT.spatial.png', res_dir_ont), p,
  width = 10, height = 2.5, units = 'in', dpi = 300
)

## SR spatial expression
iso_to_plot <- c('GFAP-Event1', 'GFAP-Event2')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 1, align = 'hv')
p

ggsave(
  sprintf('%s/figures/GFAP.SR.spatial.png', res_dir_ont), p,
  width = 10, height = 2.5, units = 'in', dpi = 300
)


### Example (3): RPL5, Fig S7B
## load ONT gtf and SR gtf
svs_ont <- import.gff(sprintf("%s/transcripts/GBM_3/GBM_3.svs_iso.oa.gtf", res_dir_ont))
# add transcript name to exons
txid2name <- svs_ont %>% 
  filter(type == 'transcript') %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>% 
  deframe()
svs_ont$transcript_name <- txid2name[svs_ont$transcript_id]
svs_ont$transcript_id_old <- svs_ont$transcript_id
svs_ont$transcript_id <- svs_ont$transcript_name

# add gene name and id to exons
svs_ont$gene_name <- svs_ont$transcript_name %>% strsplit("[_-]") %>% sapply(`[[`, 1)
# remove suffixes from gene_id
svs_ont$gene_id <- svs_ont$gene_id %>% str_remove("\\.\\d+$")
genename2id <- svs_ont %>% 
  as.data.frame() %>%
  dplyr::select(gene_name, gene_id) %>% 
  distinct() %>% 
  deframe()
svs_ont$gene_id <- genename2id[svs_ont$gene_name]
svs_ont$gene_symbol <- svs_ont$gene_name

## load short-read gtf
svs_sr <- import.bed(sprintf("%s/events/ZH8811Abulk/svs.exon.bed", res_dir_trend)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1)
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order)) # Create the new name

## load SR bam coverage
cov <- BamFile(sprintf("%s/bam/selected.ZH8811Abulk.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')

## load spatial expression data
expr_ont <- read.csv(
  sprintf('%s/figures/source_data/RPL5.ONT.GBM_3.csv', res_dir_ont)
)
expr_sr <- read.csv(
  sprintf('%s/figures/source_data/RPL5.SR.ZH8811Abulk.csv', res_dir_ont)
) %>%
  mutate(isoform_name = ifelse(isoform == "RPL5:chr1:92836272-92837576:1", 
                               'RPL5-Event1', 'RPL5-Event2'))

## RPL5 transcripts
gene_to_plot <- 'RPL5'

# create ONT assembly object
subset_ont_gr <- svs_ont %>% 
  filter(gene_name == gene_to_plot) %>%
  filter(transcript_name %in% c('RPL5-202', 'RPL5_Iso_1', 'RPL5_Iso_2'))

subset_ont_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_ont_gr),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# create SR event bed
subset_sr_bed <- svs_sr %>% 
  filter(gene == gene_to_plot)

# genomic coordinates
chrom <- seqnames(subset_ont_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
chromstart <- min(start(subset_ont_gr))
chromend <- max(end(subset_ont_gr))
# chromstart <- min(subset_sr_bed$start)
# chromend <- max(subset_sr_bed$end) + 1000
strand <- strand(subset_ont_gr)[1] %>% as.character()
n_isoforms <- subset_ont_gr$transcript_name %>% unique() %>% length()

png(
  sprintf('%s/figures/RPL5.png', res_dir_ont), 
  width = 7.25, height = 3, units = 'in', res = 300
)

# create new page
pageCreate(width = 7.25, height = 3, default.units = "inches")

# plot and place ONT transcripts
height_ont = 0.3 * n_isoforms
plotTranscripts(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as,
  labels = 'transcript', fill = '#0066CC', 
  stroke = 0.5,
  x = 1, y = 0, width = 6, height = height_ont,
  just = c("left", "top"), default.units = "inches",
)

# plot ONT dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = 0.1, y1 = height_ont
)
plotText(
  label = 'GBM3 (ONT)', rot = 90,
  fontsize = 8,
  x = 0.75, y = (height_ont + 0.05) / 2,
  just = "center", default.units = "inches"
)

# plot SR bam coverage
plotSignal(
  cov,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, x = 1, y = height_ont + 0.2, 
  width = 6, height = 0.3,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#8E44AD', fill = '#8E44AD',
  # linecolor = '#E23D28', fill = '#E23D28',
  scale = T, fontsize = 8,
)

# plot and place SR events
plotRanges(
  subset_sr_bed,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, collapse = FALSE, fill = '#E23D28',
  x = 1, y = height_ont + 0.6, width = 6, height = 0.2,
)
for (peak in subset_sr_bed$peak_name){
  start = subset_sr_bed %>% filter(peak_name == peak) %>% pull(start)
  end = subset_sr_bed %>% filter(peak_name == peak) %>% pull(end)
  plotText(
    label = peak, fontsize = 8,
    x = 1 + (start - chromstart) / (chromend - chromstart) * 6,
    y = height_ont + 0.9,
    just = "left", default.units = "inches"
  )
}

# plot SR dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = height_ont + 0.1, y1 = height_ont + 1
)
plotText(
  label = 'ZH8811Abulk (SR)', rot = 90,
  fontsize = 8,
  x = 0.75, y = height_ont + 0.55,
  just = "center", default.units = "inches"
)

# Plot genome scale
plotGenomeLabel(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, fontsize = 8,
  x = 1, y = height_ont + 1.0,
  length = 6, default.units = "inches",
)

# Plot strand info
label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
label = paste0(gene_to_plot, ', ', label)
plotText(
  label = label, fontsize = 8,
  x = 1, y = 0.15,
  just = "left", default.units = "inches"
)

pageGuideHide()
dev.off()

## Plot spatial region annotation
# ONT region annotation
m5.1 <- ggplot(expr_ont, aes(x = -array_row, y = -array_col, color = anno_hzy)) +
  geom_point(size = 1) + 
  scale_color_brewer(palette = 'Set1') +
  theme_minimal() +
  coord_fixed(ratio = 0.5) +
  labs(title = 'GBM3 (ONT)', color = '') + 
  theme(
    text = element_text(family = "Arial", size = 16),
    plot.title = element_text(family = "Arial", size = 16, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) + 
  guides(color = guide_legend(override.aes = list(size = 3)))

df <- expr_ont %>%
  filter(layer == 'ratios_obs', isoform == 'RPL5_Iso_2') %>%
  mutate(group = ifelse(value > 0.5, '>50%', '<=50%'))

m5.2 <- df %>%
  ggplot(aes(x = group)) +
  geom_bar(aes(fill = anno_hzy), position = 'fill', color = 'black') +
  geom_text(
    data = data.frame(),
    aes(
      x = 1.5, y = 1.05,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$anno_hzy))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  scale_fill_brewer(palette = 'Set1') +
  labs(x = 'RPL5_Iso_2', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.text = element_text(family = "Arial", size = 14),
    legend.position = "none"
  )
m5 <- cowplot::plot_grid(m5.1, m5.2, nrow = 1, rel_widths = c(1, 0.4), align = 'h')
m5

ggsave(
  sprintf('%s/figures/RPL5.ONT.anno.png', res_dir_ont), m5,
  width = 8.5, height = 4, units = 'in', dpi = 300
)

# SR region annotation
m6.1 <- ggplot(expr_sr, aes(x = -array_row, y = array_col, color = mp)) +
  geom_point(size = 0.8) + 
  theme_minimal() +
  coord_fixed(ratio = 0.6) +
  labs(title = 'ZH8811Abulk (SR)', color = '') + 
  theme(
    text = element_text(family = "Arial", size = 16),
    plot.title = element_text(family = "Arial", size = 16, hjust = 0.5),
    axis.text = element_blank(),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) + 
  guides(color = guide_legend(
    override.aes = list(size = 3),
    ncol = 2
  ))

df <- expr_sr %>%
  filter(layer == 'ratios_obs', isoform_name == 'RPL5-Event1') %>%
  mutate(group = ifelse(value > 0.4, '>40%', '<=40%'))

m6.2 <- df %>%
  ggplot(aes(x = group)) +
  geom_bar(aes(fill = mp), position = 'fill', color = 'black') +
  geom_text(
    data = data.frame(),
    aes(
      x = 1.5, y = 1.05,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$mp))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  labs(x = 'RPL5-Event1', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 14),
    axis.text = element_text(family = "Arial", size = 14),
    legend.position = "none"
  )
m6 <- cowplot::plot_grid(m6.1, m6.2, nrow = 1, rel_widths = c(1, 0.4))
m6

ggsave(
  sprintf('%s/figures/RPL5.SR.anno.png', res_dir_ont), m6,
  width = 8.5, height = 4, units = 'in', dpi = 300
)

## ONT spatial expression
iso_to_plot <- c('RPL5-202', 'RPL5_Iso_1', 'RPL5_Iso_2')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i],
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) +
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
    theme_void() +
    theme(
      text = element_text(size = 14, family = 'Arial'),
      strip.text = element_text(size = 14, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 14, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
    theme_void() + 
    theme(
      text = element_text(size = 14, family = 'Arial'),
      strip.text = element_text(size = 14, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 14, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list2 <- c(p_list2, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 2, align = 'hv')
p

ggsave(
  sprintf('%s/figures/RPL5.ONT.spatial.png', res_dir_ont), p,
  width = 7.5, height = 5, units = 'in', dpi = 300
)

## SR spatial expression
iso_to_plot <- c('RPL5-Event1', 'RPL5-Event2')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
    theme_void() + 
    theme(
      text = element_text(size = 14, family = 'Arial'),
      strip.text = element_text(size = 14, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 14, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
    theme_void() + 
    theme(
      text = element_text(size = 14, family = 'Arial'),
      strip.text = element_text(size = 14, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 14, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list1 <- c(p_list1, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 2, align = 'hv')
p

ggsave(
  sprintf('%s/figures/RPL5.SR.spatial.png', res_dir_ont), p,
  width = 5, height = 5, units = 'in', dpi = 300
)


### Example (4): RPS8, Fig S7A
## load ONT gtf and SR gtf
svs_ont <- import.gff(sprintf("%s/transcripts/GBM_3/GBM_3.svs_iso.oa.gtf", res_dir_ont))
# add transcript name to exons
txid2name <- svs_ont %>% 
  filter(type == 'transcript') %>%
  as.data.frame() %>%
  dplyr::select(transcript_id, transcript_name) %>% 
  distinct() %>% 
  deframe()
svs_ont$transcript_name <- txid2name[svs_ont$transcript_id]
svs_ont$transcript_id_old <- svs_ont$transcript_id
svs_ont$transcript_id <- svs_ont$transcript_name

# add gene name and id to exons
svs_ont$gene_name <- svs_ont$transcript_name %>% strsplit("[_-]") %>% sapply(`[[`, 1)
# remove suffixes from gene_id
svs_ont$gene_id <- svs_ont$gene_id %>% str_remove("\\.\\d+$")
genename2id <- svs_ont %>% 
  as.data.frame() %>%
  dplyr::select(gene_name, gene_id) %>% 
  distinct() %>% 
  deframe()
svs_ont$gene_id <- genename2id[svs_ont$gene_name]
svs_ont$gene_symbol <- svs_ont$gene_name

## load short-read gtf
svs_sr <- import.bed(sprintf("%s/events/ZH8811Abulk/svs.exon.bed", res_dir_trend)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1)
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order)) # Create the new name

## load SR bam coverage
cov_gbm <- BamFile(sprintf("%s/bam/RPS8.ZH8811Abulk.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')

## load healthy control
cov_dlpfc <- BamFile(sprintf("%s/bam/RPS8.151675.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')
expr_dlpfc <- read.csv(
  sprintf('%s/figures/source_data/RPS8.DLPFC.151675.csv', res_dir_ont)
)

## load spatial expression data
expr_ont <- read.csv(
  sprintf('%s/figures/source_data/RPS8.ONT.GBM_3.csv', res_dir_ont)
)
expr_sr <- read.csv(
  sprintf('%s/figures/source_data/RPS8.SR.ZH8811Abulk.csv', res_dir_ont)
) %>%
  mutate(isoform_name = ifelse(isoform == "RPS8:chr1:44776530-44778779:1", 
                               'RPS8-Event1', 'RPS8-Event2'))

## RPS8 transcripts
gene_to_plot <- 'RPS8'

# create ONT assembly object
subset_ont_gr <- svs_ont %>% 
  filter(gene_name == gene_to_plot) %>%
  filter(transcript_name %in% c('RPS8-202', 'RPS8-207', 'RPS8_Iso_1', 'RPS8_Iso_2'))

subset_ont_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_ont_gr),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# create SR event bed
subset_sr_bed <- svs_sr %>% 
  filter(gene == gene_to_plot)

# genomic coordinates
chrom <- seqnames(subset_ont_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
chromstart <- min(start(subset_ont_gr))
chromend <- max(end(subset_ont_gr))
# chromstart <- min(subset_sr_bed$start)
# chromend <- max(subset_sr_bed$end) + 1000
strand <- strand(subset_ont_gr)[1] %>% as.character()
n_isoforms <- subset_ont_gr$transcript_name %>% unique() %>% length()

png(
  sprintf('%s/figures/RPS8.png', res_dir_ont), 
  width = 7.25, height = 3, units = 'in', res = 300
)

# create new page
pageCreate(width = 7.25, height = 3, default.units = "inches")

# plot and place ONT transcripts
height_ont = 0.3 * n_isoforms
plotTranscripts(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as,
  labels = 'transcript', fill = '#0066CC', 
  stroke = 0.5,
  x = 1, y = 0, width = 6, height = height_ont,
  just = c("left", "top"), default.units = "inches",
)

# plot ONT dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = 0.1, y1 = height_ont
)
plotText(
  label = 'GBM3 (ONT)', rot = 90,
  fontsize = 8,
  x = 0.75, y = (height_ont + 0.05) / 2,
  just = "center", default.units = "inches"
)

# plot SR bam coverage
plotSignal(
  cov_gbm,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, x = 1, y = height_ont + 0.2, 
  width = 6, height = 0.3,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#8E44AD', fill = '#8E44AD',
  # linecolor = '#E23D28', fill = '#E23D28',
  scale = T, fontsize = 8,
)

# plot and place SR events
plotRanges(
  subset_sr_bed,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, collapse = FALSE, fill = '#E23D28',
  x = 1, y = height_ont + 0.6, width = 6, height = 0.2,
)
for (peak in subset_sr_bed$peak_name){
  plotText(
    label = peak, fontsize = 8,
    x = 2,
    y = height_ont + ifelse(peak == 'RPS8-Event1', 0.75, 0.65),
    just = "left", default.units = "inches"
  )
}

# plot SR dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = height_ont + 0.1, y1 = height_ont + 0.9
)
plotText(
  label = 'ZH8811Abulk (SR)', rot = 90,
  fontsize = 8,
  x = 0.75, y = height_ont + 0.5,
  just = "center", default.units = "inches"
)

# plot healthy DLPFC bam coverage
plotSignal(
  cov_dlpfc,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, x = 1, y = height_ont + 1, 
  width = 6, height = 0.3,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#1e6510', fill = '#1e6510',
  scale = T, fontsize = 8,
)


# plot SR dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = height_ont + 1, y1 = height_ont + 1.4
)
plotText(
  label = '151675 (SR, healthy DLPFC)',
  fontsize = 8,
  x = 1, y = height_ont + 1.2,
  just = "left", default.units = "inches"
)


# Plot genome scale
plotGenomeLabel(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ont_as, fontsize = 8,
  x = 1, y = height_ont + 1.4,
  length = 6, default.units = "inches",
)

# Plot strand info
label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
label = paste0(gene_to_plot, ', ', label)
plotText(
  label = label, fontsize = 8,
  x = 1, y = 0.15,
  just = "left", default.units = "inches"
)

pageGuideHide()
dev.off()

## ONT spatial expression
iso_to_plot <- c('RPS8-202', 'RPS8-207', 'RPS8_Iso_1', 'RPS8_Iso_2')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i],
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) +
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_ont %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = -array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (ONT)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
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
  p_list2 <- c(p_list2, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 2, align = 'hv')
p

ggsave(
  sprintf('%s/figures/RPS8.ONT.spatial.png', res_dir_ont), p,
  width = 10, height = 5, units = 'in', dpi = 300
)

## SR spatial expression
iso_to_plot <- c('RPS8-Event1', 'RPS8-Event2')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    isoform_name == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 2, align = 'hv')
p

ggsave(
  sprintf('%s/figures/RPS8.SR.spatial.png', res_dir_ont), p,
  width = 5, height = 5, units = 'in', dpi = 300
)

## DLPFC spatial expression
iso_to_plot <- c("RPS8:chr1:44777629-44778779:1", "RPS8:chr1:44776092-44777647:1")

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_dlpfc %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(y = - array_row, x = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s", str_replace(iso_to_plot[i], 'RPS8:', ''))) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 1.7) +
    theme_void() + 
    theme(
      text = element_text(size = 14, family = 'Arial'),
      strip.text = element_text(size = 14, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 14, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratio
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_dlpfc %>% filter(
    isoform == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(y = - array_row, x = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s", str_replace(iso_to_plot[i], 'RPS8:', ''))) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 1.7) +
    theme_void() + 
    theme(
      text = element_text(size = 14, family = 'Arial'),
      strip.text = element_text(size = 14, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
      legend.position = 'right',
      legend.text = element_text(size = 14, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list1 <- c(p_list1, list(p))
}

p <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 2, align = 'hv')
p

ggsave(
  sprintf('%s/figures/RPS8.DLPFC.spatial.png', res_dir_ont), p,
  width = 6, height = 6, units = 'in', dpi = 300
)


### Example (5): TPM3
## load reference gtf and SR gtf
hg38_gtf <- import.gff("~/reference/hg38/gencode.v47.annotation.gtf.gz")
hg38_gtf$gene_id <- sapply(strsplit(hg38_gtf$gene_id, "\\."), `[`, 1)

## load short-read gtf
svs_sr <- import.bed(sprintf("%s/events/ZH916bulk/svs.exon.bed", res_dir_trend)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1)
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order)) # Create the new name

## load SR bam coverage
cov_gbm_916bulk <- BamFile(sprintf("%s/bam/TPM3.ZH916bulk.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')

## load healthy control
cov_gbm_8811Abulk <- BamFile(sprintf("%s/bam/TPM3.ZH8811Abulk.bam", data_dir_trend)) %>%
  readGAlignments() %>% coverage() %>% 
  as(., 'GRanges')

## load spatial expression data
expr_sr <- read.csv(
  sprintf('%s/figures/source_data/TPM3.SR.ZH916bulk.csv', res_dir_ont)
) %>%
  left_join(
    svs_sr %>% dplyr::select(name, peak_name),
    by = c('isoform' = 'name')
  )

## TPM3 transcripts
gene_to_plot <- 'TPM3'

subset_ref_gr <- hg38_gtf %>% 
  filter(
    gene_name == gene_to_plot,
    transcript_type == 'protein_coding',
    transcript_name %in% c('TPM3-202', 'TPM3-204', 'TPM3-206', 
                           'TPM3-211', 'TPM3-212', 'TPM3-223')
  ) %>%
  mutate(transcript_id = transcript_name)

subset_ref_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_ref_gr),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# # create SR event bed
# subset_sr_bed <- svs_sr %>% 
#   filter(gene == gene_to_plot)

# create SR event gtf
# SR gtf
subset_sr_gtf <- import.gff(sprintf("%s/events/ZH916bulk/svs.exon.gtf", res_dir_trend)) %>%
  as.data.frame() %>%
  # convert gene_symbol to gene_id
  mutate(gene_name = gene_id, .keep = 'unused') %>%
  filter(gene_name == gene_to_plot) %>%
  left_join(
    subset_ref_gr %>% as.data.frame() %>% dplyr::select(gene_name, gene_id) %>% unique(),
    by = 'gene_name'
  ) %>%
  # convert transcript_id to name
  left_join(
    svs_sr %>% dplyr::select(name, peak_name),
    by = c('transcript_id' = 'name')
  ) %>%
  mutate(
    transcript_name = transcript_id,
    transcript_id = peak_name
  ) %>%
  GRanges()

subset_sr_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_sr_gtf),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# genomic coordinates
chrom <- seqnames(subset_ref_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
chromstart <- min(start(subset_ref_gr))
chromend <- max(end(subset_ref_gr))
strand <- strand(subset_ref_gr)[1] %>% as.character()
n_isoforms <- subset_ref_gr$transcript_name %>% unique() %>% length()

png(
  sprintf('%s/figures/TPM3.png', res_dir_ont), 
  width = 7.25, height = 3.5, units = 'in', res = 300
)
pdf(
  sprintf('%s/figures/TPM3.pdf', res_dir_ont), 
  width = 7.25, height = 3.5, useDingbats = FALSE
)

# create new page
pageCreate(width = 7.25, height = 3.5, default.units = "inches", showGuides = FALSE)

# plot and place reference transcripts
height_ont = 0.3 * n_isoforms
plotTranscripts(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ref_as,
  labels = 'transcript', fill = '#0066CC', 
  stroke = 0.5,
  x = 1, y = 0, width = 6, height = height_ont,
  just = c("left", "top"), default.units = "inches",
)

# plot reference label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = 0.1, y1 = height_ont
)
plotText(
  label = 'Reference isoforms (coding)', rot = 90,
  fontsize = 8,
  x = 0.75, y = (height_ont + 0.05) / 2,
  just = "center", default.units = "inches"
)


# plot SR bam coverages
plotSignal(
  cov_gbm_916bulk,
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ref_as, x = 1, y = height_ont + 0.2, 
  width = 6, height = 0.3,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#8E44AD', fill = '#8E44AD',
  # linecolor = '#E23D28', fill = '#E23D28',
  scale = T, fontsize = 8,
)
# plotText(
#   label = 'ZH916bulk (SR)',
#   fontsize = 8,
#   x = 2.5, y = height_ont + 0.3,
#   just = "left", default.units = "inches"
# )

# plotSignal(
#   cov_gbm_8811Abulk,
#   chrom = chrom, chromstart = chromstart, chromend = chromend,
#   assembly = subset_ref_as, x = 1, y = height_ont + 0.6, 
#   width = 6, height = 0.3,
#   just = c("left", "top"), default.units = "inches", 
#   linecolor = '#8E44AD', fill = '#8E44AD',
#   # linecolor = '#E23D28', fill = '#E23D28',
#   scale = T, fontsize = 8,
# )
# plotText(
#   label = 'ZH8811Abulk (SR)',
#   fontsize = 8,
#   x = 2.5, y = height_ont + 0.7,
#   just = "left", default.units = "inches"
# )

# Plot TREND events
plotTranscripts(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_sr_as,
  labels = 'transcript', fill = '#E23D28', 
  stroke = 0.5,
  x = 1, y = height_ont + 0.4, width = 6, height = 0.6,
  just = c("left", "top"), default.units = "inches",
)

# plot SR dataset label
plotSegments(
  x0 = 0.85, x1 = 0.85, y0 = height_ont + 0.1, y1 = height_ont + 1.2
)
plotText(
  label = 'ZH916bulk (SR)', rot = 90,
  fontsize = 8,
  x = 0.75, y = (height_ont + 0.1 + height_ont + 1.2) / 2,
  just = "center", default.units = "inches"
)

# Plot genome scale
plotGenomeLabel(
  chrom = chrom, chromstart = chromstart, chromend = chromend,
  assembly = subset_ref_as, fontsize = 8,
  x = 1, y = height_ont + 1.2,
  length = 6, default.units = "inches",
)

# Plot strand info
label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
label = paste0(gene_to_plot, ', ', label)
plotText(
  label = label, fontsize = 8,
  x = 1, y = 0.15,
  just = "left", default.units = "inches"
)

pageGuideHide()
dev.off()

## SR spatial expression
iso_to_plot <- c('TPM3-Event1', 'TPM3-Event4')

# log1p
p_list1 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    peak_name == iso_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Purples', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
  p_list1 <- c(p_list1, list(p))
}

# smoothed ratios
p_list2 <- c()
for (i in seq_along(iso_to_plot)){
  p <- expr_sr %>% filter(
    peak_name == iso_to_plot[i], 
    layer == 'ratios_smoothed'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'Ratio', title = sprintf("%s (SR)", iso_to_plot[i])) +
    scale_color_distiller(palette = 'Spectral') +
    coord_fixed(ratio = 0.6) +
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
  p_list2 <- c(p_list2, list(p))
}


m7.1 <- cowplot::plot_grid(plotlist = c(p_list1, p_list2), nrow = 2, align = 'hv')
m7.1


## SR region annotation
df <- expr_sr %>%
  filter(layer == 'ratios_obs', peak_name == 'TPM3-Event1') %>%
  mutate(group = ifelse(value > 0.5, '>50%', '<=50%'))
table(df$group)

df_prop <- df %>%
  group_by(group, mp) %>%
  summarise(count = n()) %>%
  mutate(
    prop = count / sum(count)
  )

m7.2 <- df_prop %>%
  ggplot(aes(x = group, y = prop, fill = mp)) +
  geom_bar(stat = 'identity') +
  geom_text(
    data = df_prop %>% filter(prop > 0.01),
    mapping = aes(label = mp), position = position_stack(vjust = .5)
  ) +
  geom_text(
    data = data.frame(),
    mapping = aes(
      x = 1.5, y = 1.05, fill = NULL,
      label = sprintf('p = %.2e', chisq.test(table(df$group, df$mp))$p.value)
    ), vjust = 0.5, hjust = 0.5
  ) +
  labs(x = 'TPM3-Event1', y = 'Proportion', fill = '') +
  theme_classic() + 
  theme(
    text = element_text(family = "Arial", size = 12),
    axis.text = element_text(family = "Arial", size = 12),
    legend.position = "none"
  )
m7.2

## Marker expression
sup_expr_sr <- read.csv(
  sprintf('%s/figures/source_data/TPM3.SR.ZH916bulk.sup.csv', res_dir_ont)
)

covar_to_plot <- c('SNAP25', 'ACTB', 'IGKC', 'NDRG1')
p_list3 <- c()
for (i in seq_along(covar_to_plot)){
  p <- sup_expr_sr %>% filter(
    gene == covar_to_plot[i], 
    layer == 'log1p'
  ) %>%
    ggplot(aes(-array_row, y = array_col, color = value)) + 
    geom_point(size = 0.5) +
    labs(color = 'LogExpr', title = sprintf("%s", covar_to_plot[i])) +
    scale_color_distiller(palette = 'Reds', direction = 1) +
    coord_fixed(ratio = 0.6) +
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
  p_list3 <- c(p_list3, list(p))
}
m7.3 <- cowplot::plot_grid(plotlist = p_list3, nrow = 2, align = 'hv')
m7.3

m7 <- cowplot::plot_grid(
  plotlist = list(m7.1, m7.2, m7.3), nrow = 1, rel_widths = c(1, 0.5, 1)
)
m7
ggsave(
  sprintf('%s/figures/TPM3.SR.spatial.png', res_dir_ont), m7,
  width = 12, height = 5, units = 'in', dpi = 300
)
ggsave(
  sprintf('%s/figures/TPM3.SR.spatial.pdf', res_dir_ont), m7,
  width = 12, height = 5, units = 'in', dpi = 300
)
