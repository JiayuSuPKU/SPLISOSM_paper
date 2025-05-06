library(tidyverse)
library(scales)
library(ggpubr)

extrafont::loadfonts()

data_dir <- "~/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir <- "~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
date <- "1119"

source("~/Projects/SPLISOSM_paper/scripts/visium_mouse_cbs/visualization/plotgardener_peak.R")

## make figure directory if doesn't exist
dir.create(sprintf("%s/figures/main_example/", res_dir), recursive = TRUE, showWarnings = FALSE)
dir.create(sprintf("%s/figures/svs_rbp/", res_dir), recursive = TRUE, showWarnings = FALSE)

### Visualize example TREND events
## Load data
# Load cbs SVS peak bed
svs_bed <- import.bed(sprintf('%s/events/cbs_svs.peak.bed', res_dir)) %>%
  as.data.frame() %>%
  mutate(
    gene = sapply(strsplit(name, ":"), `[`, 1),
    gene = str_replace(gene, 'Sept', 'Septin') # rename Sept to Septin
  ) %>%
  arrange(seqnames, start, end, gene) %>%
  group_by(gene) %>%  # Group by gene
  mutate(order = row_number()) %>%  # Assign relative order within each group
  ungroup() %>%  # Ungroup to avoid grouped output
  mutate(peak_name = paste0(gene, "-Event", order))  # Create the new name

# Load isoform ratio
df_svs_ratio <- read_csv(sprintf("%s/figures/source_data/svs_ratio.csv", res_dir)) %>%
  left_join(svs_bed %>% dplyr::select(name, peak_name), by = c('isoform' = 'name')) %>%
  mutate(gene = str_replace(gene, 'Sept', 'Septin'))  # rename Sept to Septin

# load RBP expression
df_rbp_expr <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_dir)) %>%
  filter(array_col > 15)

# load reference transcripts annotation
ensembl_gtf <- import.gff('~/reference/mm10/Mus_musculus.GRCm38.102.gtf')
ensembl_gtf <- renameSeqlevels(ensembl_gtf, value = paste0('chr', seqlevels(ensembl_gtf)))
# refseq does not have transcript_name column and uses gene symbol as gene_id
refseq_gtf <- import.gff('~/reference/mm10/mm10.ncbiRefSeq.gtf')
refseq_gtf$transcript_name <- refseq_gtf$transcript_id

# # For gencode, need to remove the version suffix in gene_id
# gencode_gtf <- import.gff('~/reference/mm10/gencode.vM10.annotation.gtf')
# gencode_gtf$gene_id <- sapply(strsplit(gencode_gtf$gene_id, "\\."), `[`, 1)

# load FIMO RBP binding sites
fimo_dir <- sprintf('%s/events/fimo.out/', res_dir)
fimo_rbp_dirs <- list.files(fimo_dir, pattern = 'cbs_svs_', full.names = FALSE)
rbs_fimo_gtf <- lapply(fimo_rbp_dirs, function(rbp_dir){
  rbp_gtf <- import.gff(sprintf('%s/%s/fimo.gff', fimo_dir, rbp_dir))
  rbp_gtf$RBP_Name <- toupper(gsub('cbs_svs_', '', rbp_dir))
  rbp_gtf$RBP_Name <- paste(rbp_gtf$RBP_Name, '(motif)')
  return(rbp_gtf)
}) %>% do.call(c, .)

# load POSTAR RBP binding sites
rbs_clip_gtf <- import.gff('~/reference/POSTAR3/mouse.gtf')
rbs_clip_gtf$RBP_Name <- paste(rbs_clip_gtf$gene_id, '(CLIP)')

# combine all RBP binding sites
rbs_all_gtf <- c(rbs_fimo_gtf, rbs_clip_gtf)

## Example (1): Arpp21, Fig 4H
gene_to_plot = 'Arpp21'
transcript_name_list = c(
  'Arpp21-223', 'Arpp21-231', 'Arpp21-230', 'Arpp21-203', 'Arpp21-229',
  'Arpp21-225', 'Arpp21-212', 'Arpp21-214'
)
# peak_to_plot <- c('Arpp21-Event1', 'Arpp21-Event2', 'Arpp21-Event4')
peak_to_plot <- c('Arpp21-Event1', 'Arpp21-Event2')
# transcript structure
png(
  sprintf("%s/figures/main_example/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('RBM24 (motif)', 'CELF4 (CLIP)'),
  transcript_name_list = transcript_name_list,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F, label_peak = F,
  peak_padding = c(500, 1000)
)
# add Mir128-2 transcript
plotRanges(
  ensembl_gtf %>% filter(gene_name == 'Mir128-2', type == 'exon'),
  chrom = 'chr9', chromstart = 112065181 - 500, chromend = 112185786 + 1000,
  assembly = 'mm10', collapse = TRUE, fill = '#8B0000',
  x = 0.5, y = 0.5, width = 6, height = 0.8, strike = 10
)
plotText(
  label = 'Mir128-2', fontsize = 8, fontcolor = '#8B0000',
  x = 2.7, y = 0.5,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.2) +
    # coord_fixed(ratio = 1.5) + 
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') +
    scale_color_distiller(palette = 'Blues', direction = 1) +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/main_example/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)
ggsave(sprintf("%s/figures/log1p_%s.pdf", res_dir, gene_to_plot), p, width = 4, height = 2)

## Example (2): Gnao1, Figure 4G
gene_to_plot <- 'Gnao1'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# peak_to_plot <- c('Gnao1-Event1', 'Gnao1-Event4', 'Gnao1-Event5')
peak_to_plot <- c('Gnao1-Event1', 'Gnao1-Event4')
# transcript structure
png(
  sprintf("%s/figures/main_example/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 4, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  fill = '#8B0000',
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(9000, 100),
)
# Highlight exon
plotRect(
  x = 2.1, y = 0.6, width = 2.2, height = 0.9,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKD – WT) in motor neuron = 0.18', fontsize = 8,
  x = 2 , y = 0.5,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.2) +
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') +
    scale_color_distiller(palette = 'Reds', direction = 1) +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(
  sprintf("%s/figures/main_example/ratio_%s.png", res_dir, gene_to_plot), 
  p, width = 7.5, height = 2.5
)
ggsave(
  sprintf("%s/figures/log1p_%s.pdf", res_dir, gene_to_plot), 
  p, width = 4, height = 2
)

# spatial ratio vs rbp expression
df <- df_svs_ratio %>% filter(
  # peak_name %in% peak_to_plot,
  peak_name == 'Gnao1-Event1',
  layer == 'ratios_obs'
) %>%
  left_join(
    df_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  # facet_wrap(~ peak_name, scales = 'free', nrow = 1) +
  stat_density_2d(
    data = df %>% filter(ratio > 0), 
    mapping = aes(x = expression, y = ratio, fill = ..level..),
    bins = 5, geom = "polygon"
  ) +
  geom_point(alpha = 0.1, color = 'gray', size = 0.5) +
  scale_fill_distiller(palette = 'Reds', direction = 1) +
  stat_smooth(method = 'lm') +
  labs(x = 'Rbfox3 log-normalized expression', y = 'Observed ratio',
       title = 'Gnao1 - Event1') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

ggsave(
  sprintf("%s/figures/main_example/%s_vs_Rbfox3.png", res_dir, gene_to_plot), 
  p2, width = 3, height = 3
)

# Plot Gnao1 isoform sequence alignment
library(msa)
setwd("~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/figures")

# load fasta
gnao1_fa <- readAAStringSet(sprintf('%s/figures/source_data/Gnao1.fa', res_dir)) %>%
  `names<-`(c('Gnao1-201', 'Gnao1-207'))

# run alignment  
gnao1_aln <- msa(gnao1_fa, method = "Muscle", type = "protein") 

# save as pdf
msaPrettyPrint(
  gnao1_aln, y=c(249, 354), output="pdf",
  paperWidth = 6, paperHeight = 3,
  # showConsensus = 'none',
  # shadingColors = 'grays',
  shadingMode = 'functional',
  # psFonts = TRUE,
  # shadingModeArg = 'structure',
  shadingModeArg = 'accessible area',
  furtherCode = '\\setfamily{residues}{sf}',
  showLogo="none", showLegend = TRUE,
  askForOverwrite = FALSE
)


## Example (3): Pcbp2, Fig S3E
gene_to_plot <- 'Pcbp2'
transcript_list_to_plot <- refseq_gtf %>% 
  filter(
    gene_id %in% gene_to_plot,
    type == 'transcript'
  ) %>% as.data.frame() %>%
  # keep two transcript with the same end position
  group_by(end) %>%
  slice_sample(n = 2) %>% pull(transcript_name)
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)

# transcript structure of all events
png(
  sprintf("%s/figures/main_example/peak_Pcbp2_mouse.png", res_dir), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = refseq_gtf,
  ref_gtf_gene_id = 'SYMBOL',
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('PCBP2 (motif)'),
  transcript_name_list = transcript_list_to_plot,
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(3800, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage of all events
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i],
    layer == 'log1p',
    # layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    # scale_color_distiller(palette = 'Spectral') + 
    scale_color_distiller(palette = 'Blues', direction = 1) +
    labs(color = '', title = peak_to_plot[i]) + 
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
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/main_example/ratio_%s_mouse.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)

# spatial ratio vs rbp expression
df <- df_svs_ratio %>% filter(
  # peak_name %in% peak_to_plot,
  peak_name == 'Pcbp2-Event1',
  # layer == 'ratios_obs'
  layer == 'counts'
) %>%
  left_join(
    df_rbp_expr %>% filter(layer == 'log1p', gene == 'Pcbp2'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  # facet_wrap(~ peak_name, scales = 'free', nrow = 1) +
  geom_point(alpha = 0.1, color = 'gray') + 
  stat_smooth(method = 'lm') + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'Observed ratio',
       title = 'Pcbp2-Event1') + 
  theme_classic() + 
  theme(
    text = element_text(size = 16, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

df <- df_svs_ratio %>%
  filter(
    peak_name %in% c('Pcbp2-Event1', 'Pcbp2-Event2'),
    layer == 'log1p'
  ) %>%
  pivot_wider(names_from = peak_name, values_from = ratio, 
              id_cols = c(array_row, array_col, barcode))

p2 <- ggplot(df, aes(x = `Pcbp2-Event2`, y = `Pcbp2-Event1`)) +
  geom_point(alpha = 0.2, color = 'gray') + 
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  stat_smooth(method = 'lm') + 
  stat_cor() +
  labs(x = 'Pcbp2-Event2', y = 'Pcbp2-Event1', title = 'Log-normalized expression') + 
  theme_classic() + 
  theme(
    text = element_text(size = 16, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

ggsave(
  sprintf("%s/figures/main_example/%s_event_logexpr.png", res_dir, gene_to_plot), 
  p2, width = 3, height = 3
)

# ratio vs gene expression
df <- df_svs_ratio %>%
  filter(
    peak_name %in% c('Pcbp2-Event1'),
    layer == 'ratios_obs'
  ) %>%
  left_join(
    df_rbp_expr %>% filter(layer == 'log1p', gene == 'Pcbp2'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )

p3 <- ggplot(df, aes(x = expression, y = ratio)) +
  geom_point(alpha = 0.2, color = 'gray') + 
  stat_smooth(
    method = 'lm'
  ) + 
  stat_cor() +
  labs(x = 'Pcbp2 log-norm expression', y = 'Event1 observed ratio', 
       title = 'Potential self-regulation') + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    strip.text = element_text(size = 14, family = 'Arial'),
    strip.background = element_blank()
  )
p3

ggsave(
  sprintf("%s/figures/main_example/%s_event_ratio.png", res_dir, gene_to_plot), 
  p3, width = 3, height = 3
)

## Human PCBP2 eCLIP
# load gencode hg38 annotation
library(org.Hs.eg.db)
hg38_gtf <- import.gff("~/reference/hg38/gencode.v47.annotation.gtf.gz")
hg38_gtf$gene_id <- sapply(strsplit(hg38_gtf$gene_id, "\\."), `[`, 1)

# focus on PCBP2
subset_gr <- hg38_gtf %>% 
  filter(gene_name == 'PCBP2') %>%
  filter(transcript_name %in% c(
    'PCBP2-201', 'PCBP2-202', 'PCBP2-203', 'PCBP2-205', 'PCBP2-207', 'PCBP2-225'
  ))
subset_as <- assembly(
  Genome = 'custom',
  TxDb = makeTxDbFromGRanges(subset_gr),
  OrgDb = org.Hs.eg.db,
  gene.id.column = 'ENSEMBL',
  display.column = 'SYMBOL'
)

# load ENCODE eCLIP data
eclip_rep1 <- import.bw(sprintf('%s/external_data/PCBP2_eCLIP.hg38.ENCFF679YWY.bigWig', data_dir))
eclip_rep2 <- import.bw(sprintf('%s/external_data/PCBP2_eCLIP.hg38.ENCFF842LPQ.bigWig', data_dir))

# visualize the human PCBP2 region structure
png(
  sprintf("%s/figures/main_example/peak_Pcbp2_human.png", res_dir), 
  width = 7.25, height = 3, units = "in", res = 300
)
pageCreate(width = 7.25, height = 3, default.units = "inches")
# plot transcripts
plotTranscripts(
  chrom = 'chr12', chromstart = 53471000, chromend = 53481500,
  assembly = subset_as, 
  labels = NULL, fill = '#8E44AD', 
  stroke = 0.5,
  x = 1, y = 0, width = 6, height = 0.9,
  just = c("left", "top"), default.units = "inches",
)
# add eCLIP tracks
plotSignal(
  eclip_rep1, chrom = 'chr12', chromstart = 53471000, chromend = 53481500,
  assembly = subset_as, x = 1, y = 1, width = 6, height = 0.15,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#E23D28', fill = '#E23D28',
  scale = F, fontsize = 8,
)
plotText(
  label = 'PCBP2 eCLIP (rep1)', fontsize = 8,
  x = 0.3, y = 1.075,
  just = "left", default.units = "inches"
)
plotSignal(
  eclip_rep2, chrom = 'chr12', chromstart = 53471000, chromend = 53481500,
  assembly = subset_as, x = 1, y = 1.2, width = 6, height = 0.15,
  just = c("left", "top"), default.units = "inches", 
  linecolor = '#E23D28', fill = '#E23D28',
  scale = F, fontsize = 8,
)
plotText(
  label = 'PCBP2 eCLIP (rep2)', fontsize = 8,
  x = 0.3, y = 1.275,
  just = "left", default.units = "inches"
)
# plot genome scale
plotGenomeLabel(
  chrom = 'chr12', chromstart = 53471000, chromend = 53481500,
  assembly = subset_as, fontsize = 8,
  x = 1, y = 1.4,
  length = 6, default.units = "inches",
)
# plot strand info
plotText(
  label = "PCBP2 (human), + strand (5'>3')", fontsize = 8,
  x = 0.3, y = 0.2,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

## GTEX PCBP2 expression
# Load necessary libraries
library(jsonlite)
library(httr)
library(dplyr)

# download and parse gene expression data from GTEx
url <- "https://gtexportal.org/api/v2/expression/clusteredMedianGeneExpression?gencodeId=ENSG00000197111.16&datasetId=gtex_v10"
pcbp2_gene_expr <- content(GET(url), as = "text") %>% fromJSON() %>% .$medianGeneExpression

# download and parse exon expression data from GTEx
url <- "https://gtexportal.org/api/v2/expression/clusteredMedianExonExpression?gencodeId=ENSG00000197111.16&datasetId=gtex_v10"
pcbp2_exon_expr <- content(GET(url), as = "text") %>% fromJSON() %>% .$medianExonExpression %>%
  mutate(exonName = paste('Exon', gsub('ENSG00000197111.16_', '', exonId), sep = '_'))

# download and parse junction expression data from GTEx
url <- "https://gtexportal.org/api/v2/expression/clusteredMedianJunctionExpression?gencodeId=ENSG00000197111.16&datasetId=gtex_v10"
pcbp2_junction_expr <- content(GET(url), as = "text") %>% fromJSON() %>% .$medianJunctionExpression

# exon visualization
df_exon <- merge(
  pcbp2_gene_expr, pcbp2_exon_expr, 
  by = c('tissueSiteDetailId','ontologyId', 'datasetId', 'gencodeId', 'geneSymbol'),
  suffixes = c('.gene', '.exon')
) %>% 
  mutate(IsBrainTissue = grepl('Brain', tissueSiteDetailId)) %>%
  dplyr::select(-exonId) %>%
  pivot_wider(names_from = exonName, values_from = median.exon)

p1.1 <- df_exon %>%
  ggplot(aes(x = Exon_16, y = Exon_15)) + 
  geom_point(aes(color = IsBrainTissue)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = 'Exon 16', y = 'Exon 15', title = 'GTEx median exon expression') +
  scale_color_manual(values = c('FALSE' = 'gray', 'TRUE' = 'red')) +
  theme_classic() +
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    legend.position = 'none'
  )

p1.2 <- df_exon %>%
  ggplot(aes(x = Exon_15 + Exon_16, y = Exon_15 / (Exon_15 + Exon_16))) + 
  geom_point(aes(color = IsBrainTissue)) +
  geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(x = 'Exon 15 + 16', y = 'Exon 15 / (Exon 15 + 16)', title = 'GTEx median exon expression') +
  scale_color_manual(values = c('FALSE' = 'gray', 'TRUE' = 'red')) +
  theme_classic() +
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    legend.position = 'none'
  )

# junction
df_junction <- merge(
  pcbp2_gene_expr, pcbp2_junction_expr, 
  by = c('tissueSiteDetailId','ontologyId', 'datasetId', 'gencodeId', 'geneSymbol'),
  suffixes = c('.gene', '.junc')
) %>% 
  mutate(IsBrainTissue = grepl('Brain', tissueSiteDetailId)) %>%
  filter(junctionId %in% c("chr12:53468833-53471637:+", "chr12:53471808-53474868:+", "chr12:53471808-53479405:+")) %>%
  pivot_wider(names_from = junctionId, values_from = median.junc) %>%
  dplyr::rename(
    junction_16 = `chr12:53468833-53471637:+`,
    junction_17 = `chr12:53471808-53474868:+`, 
    junction_18 = `chr12:53471808-53479405:+`)

p1.3 <- df_junction %>%
  ggplot(aes(x = junction_18, y = junction_17)) + 
  geom_point(aes(color = IsBrainTissue)) + 
  geom_smooth(method = 'lm') +
  labs(x = 'Junction 18 (Exon 16)', y = 'Junction 17 (Exon 15)', 
       title = 'GTEx median junction expression') +
  scale_color_manual(values = c('FALSE' = 'gray', 'TRUE' = 'red')) +
  stat_cor() +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    legend.position = 'none'
  )

p1.4 <- df_junction %>%
  ggplot(aes(x = junction_17 + junction_18, y = junction_17 / (junction_17 + junction_18))) + 
  geom_point(aes(color = IsBrainTissue)) + 
  geom_smooth(method = 'lm') +
  labs(x = 'Junction 17 + 18', y = 'Junction 17 / (Junction 17 + 18)', 
       title = 'GTEx median junction expression') +
  scale_color_manual(values = c('FALSE' = 'gray', 'TRUE' = 'red')) +
  stat_cor() +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    legend.position = 'inside',
    legend.position.inside = c(0.7, 0.7),
    legend.background = element_blank()
  )

p1 <- cowplot::plot_grid(p1.1, p1.2, p1.3, p1.4, nrow = 1, align = 'hv')
p1
ggsave(sprintf("%s/figures/main_example/GTEx_PCBP2_expr.png", res_dir), p1, width = 14, height = 3.5)

## SupFig: Celf2
gene_to_plot <- 'Celf2'
peak_to_plot <- c('Celf2-Event1', 'Celf2-Event2', 'Celf2-Event4')
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('RBMS3 (motif)', 'CELF4 (CLIP)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(300, 500),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Fig 6H: Septin8
gene_to_plot <- 'Septin8'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Septin8-Event1', 'Septin8-Event2', 'Septin8-Event4')

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('CELF2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(1500, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.2) +
    labs(color = '', title = peak_to_plot[i]) + 
    # scale_color_distiller(palette = 'Spectral') +
    scale_color_distiller(palette = 'Purples', direction = 1) +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)
ggsave(sprintf("%s/figures/log1p_%s.pdf", res_dir, gene_to_plot), p, width = 6, height = 2)


## SupFig: Septin11
gene_to_plot <- 'Septin11'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('PABPC1 (CLIP)', 'CELF2 (CLIP)', 'CELF4 (CLIP)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(5000, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## SupFig: Hsp90aa1
gene_to_plot <- 'Hsp90aa1'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Hsp90aa1-Event1', 'Hsp90aa1-Event4', 'Hsp90aa1-Event5')

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c("KHDRBS2 (motif)"),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 1000),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Fig S3C: Hsp90ab1
gene_to_plot <- 'Hsp90ab1'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Hsp90ab1-Event1', 'Hsp90ab1-Event3', 'Hsp90ab1-Event4')

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 5, height = 3, units = "in", res = 300
)
plot_peaks_short(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('CPSF6 (CLIP)', 'PCBP2 (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 1000),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
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
      legend.key.height = unit(0.3, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)

## Fig S3C: Olfm1
gene_to_plot <- 'Olfm1'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Olfm1-Event1', 'Olfm1-Event2', 'Olfm1-Event4')

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 5, height = 3, units = "in", res = 300
)
plot_peaks_short(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('FUS (CLIP)', 'PCBP2 (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(500, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
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
      legend.key.height = unit(0.3, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)

## Fig S3C: Rexo2
gene_to_plot <- 'Rexo2'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 5.25, height = 3, units = "in", res = 300
)
plot_peaks_short(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('PUM1 (motif)', 'QK (motif)', 'RBFOX3 (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 1000),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
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
    scale_color_distiller(palette = 'Reds', direction = 1) +
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
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Extended fig: Kalrn
gene_to_plot <- 'Kalrn'
peak_to_plot <- c('Kalrn-Event2', 'Kalrn-Event3')

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('RBM47 (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(1000, 3000),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)

## Extended fig: Ap1s2
gene_to_plot <- 'Ap1s2'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('RBM47 (motif)', 'CELF4 (CLIP)', 'RBFOX3 (CLIP)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(3000, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## Extended fig: Ppp2r5c
gene_to_plot <- 'Ppp2r5c'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Ppp2r5c-Event1', 'Ppp2r5c-Event2', 'Ppp2r5c-Event4')

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('QK (motif)', 'PABPC1 (CLIP)', 'CELF4 (CLIP)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(1000, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Extended fig: Lgi3
gene_to_plot <- 'Lgi3'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)

# transcript structure
png(
  sprintf("%s/figures/svs_cds/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  rbs_gtf = rbs_all_gtf,
  rbp_name_list = c('SNRNP70 (motif)', 'QK (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(1000, 100),
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_cds/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


### All other RBPs with varying length UTRs
## Celf1
gene_to_plot <- 'Celf1'
peak_to_plot <- c('Celf1-Event1', 'Celf1-Event2')
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(3000, 100),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## Celf4
gene_to_plot <- 'Celf4'
peak_to_plot <- c('Celf4-Event1', 'Celf4-Event2', 'Celf4-Event3', 'Celf4-Event4')
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 1000), 
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 10, height = 2.5)


## Elavl2
gene_to_plot <- 'Elavl2'
peak_to_plot <- c('Elavl2-Event1', 'Elavl2-Event2', 'Elavl2-Event3')
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 3000), 
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Elavl3
gene_to_plot <- 'Elavl3'
peak_to_plot <- c('Elavl3-Event1', 'Elavl3-Event4', 'Elavl3-Event5')
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 1000), 
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Qk
gene_to_plot <- 'Qk'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(4500, 4000),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## Rbfox1
gene_to_plot <- 'Rbfox1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(10000, 100),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 7.5, height = 2.5)


## Pabpn1
gene_to_plot <- 'Pabpn1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(3400, 100),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## Pabpc1
gene_to_plot <- 'Pabpc1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 7000),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## Hnrnpa2b1
gene_to_plot <- 'Hnrnpa2b1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 1000),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


# Hnrnpk
gene_to_plot <- 'Hnrnpk'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 2000),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


## Zcrb1
gene_to_plot <- 'Zcrb1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/svs_rbp/peak_%s.png", res_dir, gene_to_plot), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 3000),
)
pageGuideHide()
dev.off()
# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'ratios_smoothed',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p
ggsave(sprintf("%s/figures/svs_rbp/ratio_%s.png", res_dir, gene_to_plot), p, width = 5, height = 2.5)


### Plot expression of selected RBPs
df_rbp_expr <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_dir)) %>%
  filter(array_col > 15)

rbp_list <- c(
  # 'Arpp21', 'Celf2', 'Pcbp2',
  'Celf1', 'Celf4', 'Elavl2', 'Elavl3', 'Qk', 'Rbfox1', 'Hnrnpa2b1', 'Hnrnpk'
)
p_list <- c()
for (rbp in rbp_list){
  p <- df_rbp_expr %>% filter(
    gene == rbp, layer == 'log1p', array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = expression)) + 
    geom_point(size = 0.5) + 
    labs(color = '', title = rbp) + 
    scale_color_distiller(palette = 'Spectral') +
    # scale_color_gradient2(low = 'white', high = '#0066CC') +
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
p <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p

ggsave(sprintf("%s/figures/svs_rbp/sub_rbp_expr.png", res_dir), p, width = 8*2.5, height = 2.5)


### Top regulators of SVS CDS and RBP genes
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

svs_to_plot <- c(
  'Sept8', 'Sept11', 'Hsp90aa1', 'Hsp90ab1', 'Olfm1', 
  'Kalrn', 'Rexo2', 'Ppp2r5c', 'Ap1s2', 'Lgi3'
) %>% sort()
rbp_to_include <- union(
  union(rbp_with_clip$rbp, rbp_with_motifs$rbp),
  c('Celf5', 'Arpp21', 'Adar', 'Adarb1', 'Adarb2')
) %>% sort()
rbp_to_highlight <- list(
  'Sept8' = c('Celf2', 'Rbfox3', 'Celf4'),
  'Sept11' = c('Pabpc1', 'Celf3'),
  'Hsp90aa1' = c('Khdrbs2'),
  'Hsp90ab1' = c('Cpsf6', 'Pcbp4'),
  'Olfm1' = c('Fus', 'Pcbp3'),
  'Kalrn' = c('Rbm47'),
  'Rexo2' = c('Pum1', 'Qk', 'Rbfox3'),
  'Ppp2r5c' = c('Qk', 'Pabpc1', 'Celf3'),
  'Ap1s2' = c('Rbm47', 'Celf4', 'Rbfox3'),
  'Lgi3' = c('Snrnp70', 'Qk')
)

p_list <- c()
for (svs in svs_to_plot){
  p <- df_du_rbp %>% 
    filter(gene == svs, covariate %in% rbp_to_include) %>%
    mutate(highlight = covariate %in% rbp_to_highlight[[svs]]) %>%
    slice_min(`pvalue_hsic-gp`, n = 5) %>%
    mutate(covariate = fct_reorder(covariate, `pvalue_hsic-gp`, .desc = F)) %>%
    ggplot(aes(x = covariate, y = `pvalue_hsic-gp`)) +
    geom_bar(stat = 'identity', aes(fill = highlight)) +
    geom_hline(yintercept = 0.01, linetype = 'dashed', color = 'black') +
    labs(title = svs, x = '', y = 'DU test pvalue (HSIC-GP)', fill = 'Has binding') + 
    scale_fill_manual(values = c('TRUE' = '#8B0000', 'FALSE' = 'gray')) +
    scale_y_continuous(trans = c('log10', 'reverse')) +
    theme_classic() + 
    theme(
      text = element_text(size = 12, family = 'Arial'),
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = 'none',
    )
  p_list <- c(p_list, list(p))
}

p <- cowplot::plot_grid(plotlist = p_list, nrow = 1)
p
ggsave(sprintf("%s/figures/svs_cds/du_top5.png", res_dir), p, width = 20, height = 3)
