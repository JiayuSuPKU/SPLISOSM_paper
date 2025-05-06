library(tidyverse)
library(scales)

extrafont::loadfonts()

data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/human_dlpfc/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/human_dlpfc/"
date = '0123'

source("~/Projects/SPLISOSM_paper/scripts/dlpfc_visium/visualization/plotgardener_peak.R")

### Visualize example TREND events
## Load data
# Load SVS peak bed
svs_bed <- import.bed(sprintf('%s/events/svs_all.exon.bed', res_dir)) %>%
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

# load reference transcripts annotation
ensembl_gtf <- import.gff('~/reference/hg38/gencode.v47.annotation.sorted.gtf.gz')

# load POSTAR RBP binding sites
rbs_clip_gtf <- import.gff('~/reference/POSTAR3/human.sorted.gtf.gz')
rbs_clip_gtf$RBP_Name <- paste(rbs_clip_gtf$gene_id, '(CLIP)')

## Main figure example: SEPTIN8
gene_to_plot <- 'SEPTIN8'

# rename SEPTIN8 human events to match with mouse
# (human) SEPTIN8-Event1 -> (mouse) SEPTIN8-Event4
# (human) SEPTIN8-Event3 -> (mouse) SEPTIN8-Event2
# (human) SEPTIN8-Event4 -> (mouse) SEPTIN8-Event1
peak_bed <- svs_bed %>% filter(gene == gene_to_plot) %>%
  mutate(
    peak_name = case_when(
      peak_name == 'SEPTIN8-Event1' ~ 'SEPTIN8-Event4',
      peak_name == 'SEPTIN8-Event2' ~ 'SEPTIN8-Event3',
      peak_name == 'SEPTIN8-Event3' ~ 'SEPTIN8-Event2',
      peak_name == 'SEPTIN8-Event4' ~ 'SEPTIN8-Event1',
      TRUE ~ peak_name
    )
  )
df_sept8_ratio <- df_svs_ratio %>% filter(gene == gene_to_plot) %>%
  mutate(
    peak_name = case_when(
      peak_name == 'SEPTIN8-Event1' ~ 'SEPTIN8-Event4',
      peak_name == 'SEPTIN8-Event2' ~ 'SEPTIN8-Event3',
      peak_name == 'SEPTIN8-Event3' ~ 'SEPTIN8-Event2',
      peak_name == 'SEPTIN8-Event4' ~ 'SEPTIN8-Event1',
      TRUE ~ peak_name
    )
  )

peak_to_plot <- c('SEPTIN8-Event1', 'SEPTIN8-Event2', 'SEPTIN8-Event4')

# Fig 6H: SEPTIN8
png(
  sprintf("%s/figures/septin8_struct.png", res_dir), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf,
  ref_gtf_gene_id = 'SYMBOL',
  transcript_name_list = c('SEPTIN8-205', 'SEPTIN8-206', 'SEPTIN8-201', 'SEPTIN8-207', 'SEPTIN8-212', 'SEPTIN8-208'), 
  rbs_gtf = rbs_clip_gtf,
  rbp_name_list = c('QKI (CLIP)', 'CELF2 (CLIP)'), # no peaks for 'SRRM4 (CLIP)', 'YWHAG (CLIP)', 'TNRC6C (CLIP)'
  peak_bed = peak_bed,
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 2000),
  # fill = '#1a4e66'
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_sept8_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    # layer == 'ratios_smoothed'
    layer == 'log1p'
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.15) +
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
ggsave(
  sprintf("%s/figures/septin8_spatial.png", res_dir), 
  p, width = 7.5, height = 2.5, units = "in", dpi = 300
)
ggsave(
  sprintf("%s/figures/septin8_spatial.pdf", res_dir), 
  p, width = 6, height = 2, units = "in", dpi = 300
)


## Fig S5G: MAP4
gene_to_plot <- 'MAP4'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)

# transcript structure of all events
png(
  sprintf("%s/figures/map4_struct.png", res_dir), 
  width = 7.25, height = 3, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot,
  ref_gtf = ensembl_gtf,
  ref_gtf_gene_id = 'SYMBOL',
  rbs_gtf = rbs_clip_gtf,
  transcript_name_list = c('MAP4-209', 'MAP4-210', 'MAP4-216', 'MAP4-201', 'MAP4-202', 'MAP4-204', 'MAP4-208'),
  rbp_name_list = c('QKI (CLIP)', 'CELF2 (CLIP)', 'RBFOX1 (CLIP)'),
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 4000),
  fill = '#1a4e66'
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i],
    # layer == 'ratios_smoothed'
    layer == 'log1p'
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
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
  sprintf("%s/figures/map4_spatial.png", res_dir), 
  p, width = 7.5, height = 2.5, units = "in", dpi = 300
)

## Fig S5E, RBP expression
# load RBP expression
df_rbp_expr <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_dir))

# visualization
rbp_to_plot <- c('NOVA2', 'PTBP2', 'QKI', 'CELF2')
p_list <- c()
for (i in seq_along(rbp_to_plot)){
  p <- df_rbp_expr %>% filter(
    gene == rbp_to_plot[i],
    layer == 'log1p'
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = expression)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = rbp_to_plot[i]) + 
    scale_color_distiller(palette = 'Spectral') +
    theme_void() + 
    theme(
      text = element_text(size = 16, family = 'Arial'),
      strip.text = element_text(size = 16, family = 'Arial'),
      plot.title = element_text(hjust = 0.5, size = 16, family = 'Arial', face = 'italic'),
      legend.position = 'right',
      legend.text = element_text(size = 16, family = 'Arial'),
      legend.key.width = unit(0.1, "in"),
      legend.key.height = unit(0.3, "in")
    )
  p_list <- c(p_list, list(p))
}
p <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
p
ggsave(
  sprintf("%s/figures/sup_rbp_spatial.png", res_dir), 
  p, width = 5, height = 5, units = "in", dpi = 300
)
