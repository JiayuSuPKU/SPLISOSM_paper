library(tidyverse)
library(scales)
library(ggrepel)
library(ggpubr)
library(rtracklayer)

extrafont::loadfonts()

res_sc_dir <- "~/Projects/SPLISOSM_paper/results/sc_tasic_nature_18"
res_ont_dir <- "~/Projects/SPLISOSM_paper/results/sit_nar_23"
res_trend_dir <- "~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
ont_date <- "1107"
trend_date <- "1119"

### Cooperative splicing regulation figures
## make figure directory if doesn't exist
dir.create(sprintf("%s/figures/", res_sc_dir), recursive = TRUE, showWarnings = FALSE)

## Load POSTAR RBP binding sites
rbs_clip_gtf <- import.gff('~/reference/POSTAR3/mouse.gtf')
rbs_clip_gtf$RBP_Name <- paste(rbs_clip_gtf$gene_id, '(CLIP)')

## Load exon usage per cell type
ct_exon_psi <- read_csv(sprintf("%s/svs_ratio.csv", res_sc_dir))

### ONT examples
source("~/Projects/SPLISOSM_paper/scripts/sit_data_analysis/visualization/plotgardener_isoform_structure.R")

## Load RBP expression and isoform usage
ont_rbp_expr <- read_csv(sprintf("%s/figures/cbs/source_data/rbp_expr.csv", res_ont_dir))
ont_svs_ratio <- read_csv(sprintf("%s/figures/cbs/source_data/svs_ratio.csv", res_ont_dir))

## Load necessary data for visualization
fimo_dir <- sprintf('%s/transcripts/cbs/fimo.out/', res_ont_dir)

# Load ONT-CBS SVS transcript GTF
svs_gtf <- import.gff(sprintf('%s/transcripts/cbs/cbs_svs.txid.gtf', res_ont_dir))
svs_gtf$gene_id <- gsub("\\..*", "", svs_gtf$gene_id)

# Find and load FIMO motif scanning results
fimo_rbp_dirs <- list.files(fimo_dir, pattern = 'cbs_svs_', full.names = FALSE)
rbs_fimo_gtf <- lapply(fimo_rbp_dirs, function(rbp_dir){
  rbp_gtf <- import.gff(sprintf('%s/%s/fimo.gff', fimo_dir, rbp_dir))
  rbp_gtf$RBP_Name <- toupper(gsub('cbs_svs_', '', rbp_dir))
  rbp_gtf$RBP_Name <- paste(rbp_gtf$RBP_Name, '(motif)')
  return(rbp_gtf)
}) %>% do.call(c, .)

# Combine all RBP binding sites
rbs_ont_gtf <- c(rbs_fimo_gtf, rbs_clip_gtf)

## Fig 5A: Clta
gene_to_plot <- 'Clta'
df <- ont_svs_ratio %>% filter(gene == gene_to_plot) %>% 
  dplyr::select(transcript_id, transcript_name) %>%
  distinct() %>%
  filter(transcript_name %in% c('Clta-203', 'Clta-204', 'Clta-205', 'Clta-206'))
transcript_to_plot <- setNames(df$transcript_id, df$transcript_name)

png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot),
  width = 7.25, height = 4, units = "in", res = 300
)
# transcript structure
plot_transcripts(
  transcript_ids = transcript_to_plot,
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  rbs_gtf = rbs_ont_gtf,
)
# highlight exon
plotRect(
  x = 6.15, y = 0.2, width = 0.15, height = 1.1,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'Exon 5 dPSI (Rbfox tKO – WT) = 0.21 in motor neurons', fontsize = 8,
  x = 4 , y = 0.2,
  just = "left", default.units = "inches"
)
plotRect(
  x = 6.5, y = 0.2, width = 0.15, height = 1.1,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'Exon 5+6 dPSI (Rbfox tKO – WT) = 0.12 in motor neurons', fontsize = 8,
  x = 4.2 , y = 1.25,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial expression
p_list <- c()
for (i in seq_along(transcript_to_plot)){
  p <- ont_svs_ratio %>% filter(
    transcript_id == transcript_to_plot[i],
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.2) +
    labs(color = '', title = names(transcript_to_plot)[i]) + 
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
p1 <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
p1
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p1, width = 4.5, height = 4
)
ggsave(
  sprintf("%s/figures/%s_sp_log1p.pdf", res_sc_dir, gene_to_plot), 
  p1, width = 4, height = 4
)

# ratio of Exon 5 + 6 vs Rbfox3 expression in the spatial data
df_clta_sp <- ont_svs_ratio %>% filter(
  transcript_name %in% c('Clta-206', 'Clta-205'),
  layer == 'ratios_smoothed'
) %>% group_by(barcode) %>%
  mutate(ratio_adj = ratio / sum(ratio)) %>%
  filter(transcript_name == 'Clta-205') %>%
  left_join(
    ont_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )

p2.1 <- ggplot(df_clta_sp, aes(x = expression, y = ratio_adj)) + 
  geom_point(alpha = 0.2, color = 'gray') + 
  stat_smooth(method = 'lm', color = '#1976d2') + 
  labs(x = '', y = 'Smoothed ratio', 
       title = 'Clta-205 / (205 + 206) (sp-ONT)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
    axis.title.x = element_blank()
  )

# ratio of Clta exons vs Rbfox3 expression in the single-cell data
# 12757//Clta,chr4,44025448,44031434,CA-12757-15918-20572-20626-21757[INC][5/1]
# 12757//Clta,chr4,44025448,44032844,CA-12757-15918-20572-20626-22758[INC][7/201][DNT]
df_clta_sc <- ct_exon_psi %>% filter(
  gene == 'Clta',
  event %in% c(
    # 'CA-12757-15918-20572-20626-21757[INC][5/1]', 
    'CA-12757-15918-20572-20626-22758[INC][7/201][DNT]'
  )
) %>%
  mutate(
    is_neuron = (!NP %in% c('CR_Meis2', 'Non Neuronal')),
    event_name = factor(
      event, levels = c(
        # 'CA-12757-15918-20572-20626-21757[INC][5/1]', 
        'CA-12757-15918-20572-20626-22758[INC][7/201][DNT]'
      ),
      labels = c(
        # 'Exon 5', 
        'Exon 5 and 6'
      )
    )
  )
p2.2 <- ggplot(df_clta_sc, aes(x = Rbfox3, y = psi)) + 
  # facet_wrap(~ event_name, scales = 'free', nrow = 2) +
  geom_point(aes(color = NP), size = 1) + 
  stat_smooth(method = 'lm', color = '#1976d2') +
  stat_smooth(
    data = df_clta_sc %>% filter(is_neuron),
    # aes(color = NP), 
    color = '#ff4d4d',
    method = 'lm',
    fullrange = TRUE
  ) + 
  
  labs(x = '', y = 'Exon inclusion ratio', 
       color = 'Cell type', title = 'Exon 5 + 6 PSI (single cell)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'none',
  )
p2.3 <- ggplot(df_clta_sc, aes(x = Rbfox3, y = psi, color = NP)) + 
  geom_point() + 
  labs(color = 'Cell type') +
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'bottom',
  ) + 
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)))
p2.3 <- ggpubr::get_legend(p2.3)

p <- cowplot::plot_grid(p2.1, p2.2, nrow = 1, align = 'hv')
p <- cowplot::add_sub(p, "Rbfox3 normalized expression", hjust = 0.5)
p <- cowplot::plot_grid(p, p2.3, nrow = 2, rel_heights = c(1, 0.2))
p

ggsave(
  sprintf("%s/figures/Clta_vs_Rbfox3.png", res_sc_dir), p, width = 6, height = 4
)
ggsave(
  sprintf("%s/figures/Clta_vs_Rbfox3.pdf", res_sc_dir), p, width = 5.5, height = 4
)


## Fig S4C: Myl6
gene_to_plot <- 'Myl6'
df <- ont_svs_ratio %>% filter(gene == gene_to_plot) %>% 
  dplyr::select(transcript_id, transcript_name) %>%
  distinct() %>%
  filter(transcript_name %in% c('Myl6-201', 'Myl6-206'))
transcript_to_plot <- setNames(df$transcript_id, df$transcript_name)

png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot),
  width = 5.25, height = 4, units = "in", res = 300
)
# transcript structure
plot_transcripts_short(
  transcript_ids = transcript_to_plot,
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  rbs_gtf = rbs_ont_gtf,
)
# Highlight exon
plotRect(
  x = 2.1, y = 0.2, width = 0.2, height = 0.4,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox2 KD – WT) = -0.26 in myoblasts', fontsize = 8,
  x = 2 , y = 0.1,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial expression
p_list <- c()
for (i in seq_along(transcript_to_plot)){
  p <- ont_svs_ratio %>% filter(
    transcript_id == transcript_to_plot[i],
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = names(transcript_to_plot)[i]) + 
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
p1 <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p1
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p1, width = 5, height = 2.5
)

# Myl6 ratio vs Rbfox3 expression in the spatial data
df_myl6_sp <- ont_svs_ratio %>% filter(
  transcript_name %in% c('Myl6-201', 'Myl6-206'),
  layer == 'ratios_smoothed'
) %>% group_by(barcode) %>%
  mutate(ratio_adj = ratio / sum(ratio)) %>%
  filter(transcript_name == 'Myl6-206') %>%
  left_join(
    ont_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )

p2.1 <- ggplot(df_myl6_sp, aes(x = expression, y = ratio_adj)) + 
  geom_point(alpha = 0.2, color = 'gray') + 
  stat_smooth(method = 'lm', color = '#1976d2') + 
  labs(x = '', y = 'Smoothed ratio', 
       title = 'Myl6-206 / (201 + 206) (sp)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    axis.title.x = element_blank()
  )

# ratio of Myl6 exons vs Rbfox3 expression in the single-cell data
# 17904//Myl6,chr10,128490860,128492135,CA-17904-4767-5061-5106-5792[INC][169/447][DNT]
# 17904//Myl6,chr10,128490860,128492135,CA-17904-4767-5478-5536-5792[INC][5/447][DNT]
# 17904//Myl6,chr10,128490860,128492135,CA-17904-4767-5478-5540-5792[INC][1/447][DNT]
df_myl6_sc <- ct_exon_psi %>% filter(
  gene == 'Myl6',
  event %in% c('CA-17904-4767-5478-5540-5792[INC][1/447][DNT]')
) %>%
  mutate(
    is_neuron = (!NP %in% c('CR_Meis2', 'Non Neuronal')),
    event_name = factor(
      event, levels = c('CA-17904-4767-5478-5540-5792[INC][1/447][DNT]'),
      labels = c('Exon 6')
    )
  )
p2.2 <- ggplot(df_myl6_sc, aes(x = Rbfox3, y = psi)) + 
  geom_point(aes(color = NP), size = 1) + 
  stat_smooth(method = 'lm', color = '#1976d2') +
  stat_smooth(
    data = df_myl6_sc %>% filter(is_neuron),
    color = '#ff4d4d',
    method = 'lm',
    fullrange = TRUE
  ) + 
  labs(x = '', y = 'Exon inclusion ratio', 
       color = 'Cell type', title = 'Exon 6 PSI (sc)') + 
  ylim(0, 1) + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'none',
  )
p2.3 <- ggplot(df_myl6_sc, aes(x = Rbfox3, y = psi, color = NP)) + 
  geom_point() + 
  labs(color = 'Cell type') +
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    strip.text = element_text(size = 14, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'bottom',
  ) + 
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)))
p2.3 <- ggpubr::get_legend(p2.3)

p <- cowplot::plot_grid(p2.1, p2.2, nrow = 1, align = 'hv')
p <- cowplot::add_sub(p, "Rbfox3 normalized expression", hjust = 0.5)
p <- cowplot::plot_grid(p, p2.3, nrow = 2, rel_heights = c(1, 0.2))
p

ggsave(
  sprintf("%s/figures/Myl6_vs_Rbfox3.png", res_sc_dir), p, width = 6, height = 4
)

## Fig S4C: Cltb
gene_to_plot <- 'Cltb'
df <- ont_svs_ratio %>% filter(gene == gene_to_plot) %>% 
  dplyr::select(transcript_id, transcript_name) %>%
  distinct()
transcript_to_plot <- setNames(df$transcript_id, df$transcript_name)

png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot),
  width = 5.25, height = 4, units = "in", res = 300
)
# transcript structure
plot_transcripts_short(
  transcript_ids = transcript_to_plot,
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  rbs_gtf = rbs_ont_gtf,
)
# Highlight exon
plotRect(
  x = 1.9, y = 0.2, width = 0.2, height = 0.4,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKD – WT) = 0.15 in motor neurons', fontsize = 8,
  x = 1.8 , y = 0.1,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial expression
p_list <- c()
for (i in seq_along(transcript_to_plot)){
  p <- ont_svs_ratio %>% filter(
    transcript_id == transcript_to_plot[i],
    # layer == 'ratios_smoothed',
    layer == 'log1p',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = names(transcript_to_plot)[i]) + 
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
p1 <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
p1
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p1, width = 5, height = 2.5
)

# Cltb ratio vs Rbfox3 expression in the spatial data
df_cltb_sp <- ont_svs_ratio %>% filter(
  transcript_name %in% c('Cltb-201', 'Cltb-202'),
  layer == 'ratios_smoothed'
) %>% group_by(barcode) %>%
  mutate(ratio_adj = ratio / sum(ratio)) %>%
  filter(transcript_name == 'Cltb-201') %>%
  left_join(
    ont_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )

p2.1 <- ggplot(df_cltb_sp, aes(x = expression, y = ratio_adj)) + 
  geom_point(alpha = 0.2, color = 'gray') + 
  stat_smooth(method = 'lm', color = '#1976d2') + 
  labs(x = '', y = 'Smoothed ratio', 
       title = 'Cltb-201 / (201 + 202) (sp)') + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    axis.title.x = element_blank()
  )

# ratio of Cltb exons vs Rbfox3 expression in the single-cell data
# 74325//Cltb,chr13,54593584,54598832,CA-74325-15551-17193-17247-20230[INC][11/56][DNT]
df_cltb_sc <- ct_exon_psi %>% filter(
  gene == 'Cltb',
  event %in% c('CA-74325-15551-17193-17247-20230[INC][11/56][DNT]')
) %>%
  mutate(
    is_neuron = (!NP %in% c('CR_Meis2', 'Non Neuronal')),
    event_name = factor(
      event, levels = c('CA-74325-15551-17193-17247-20230[INC][11/56][DNT]'),
      labels = c('Exon 5')
    )
  )
p2.2 <- ggplot(df_cltb_sc, aes(x = Rbfox3, y = psi)) + 
  geom_point(aes(color = NP), size = 1) + 
  stat_smooth(method = 'lm', color = '#1976d2') +
  stat_smooth(
    data = df_cltb_sc %>% filter(NP == "GABA"),
    color = '#007d66',
    method = 'lm',
    fullrange = TRUE
  ) + 
  labs(x = '', y = 'Exon inclusion ratio', 
       color = 'Cell type', title = 'Exon 5 PSI (sc)') + 
  ylim(0, 1) + 
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    axis.title.x = element_blank(),
    strip.text = element_text(size = 14, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'none',
  )
p2.3 <- ggplot(df_cltb_sc, aes(x = Rbfox3, y = psi, color = NP)) + 
  geom_point() + 
  labs(color = 'Cell type') +
  theme_classic() + 
  theme(
    text = element_text(size = 14, family = 'Arial'),
    plot.title = element_text(hjust = 0.5, size = 14, family = 'Arial'),
    strip.text = element_text(size = 14, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'bottom',
  ) + 
  guides(color = guide_legend(nrow = 2, override.aes = list(size = 3)))
p2.3 <- ggpubr::get_legend(p2.3)

p <- cowplot::plot_grid(p2.1, p2.2, nrow = 1, align = 'hv')
p <- cowplot::add_sub(p, "Rbfox3 normalized expression", hjust = 0.5)
p <- cowplot::plot_grid(p, p2.3, nrow = 2, rel_heights = c(1, 0.2))
p

ggsave(
  sprintf("%s/figures/Cltb_vs_Rbfox3.png", res_sc_dir), p, width = 6, height = 4
)

## Nnat
gene_to_plot <- 'Nnat'
df <- ont_svs_ratio %>% filter(gene == gene_to_plot) %>% 
  dplyr::select(transcript_id, transcript_name) %>%
  distinct() %>%
  filter(transcript_name %in% c('Nnat-202', 'Nnat-205'))
transcript_to_plot <- setNames(df$transcript_id, df$transcript_name)

png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot),
  width = 7.25, height = 4, units = "in", res = 300
)
# transcript structure
plot_transcripts(
  transcript_ids = transcript_to_plot,
  ref_gtf = svs_gtf,
  rbp_names = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  rbs_gtf = rbs_ont_gtf,
)
# Highlight exon
plotRect(
  x = 3.75, y = 0.2, width = 0.3, height = 0.4,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKD – WT) in motor neuron = 0.11', fontsize = 8,
  x = 3 , y = 0.1,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial expression
p_list <- c()
for (i in seq_along(transcript_to_plot)){
  p <- ont_svs_ratio %>% filter(
    transcript_id == transcript_to_plot[i],
    layer == 'ratios_smoothed',
    array_row > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = names(transcript_to_plot)[i]) + 
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
p1 <- cowplot::plot_grid(plotlist = p_list, nrow = 2, align = 'hv')
p1
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p1, width = 2.5, height = 4
)

# spatial ratio vs rbp expression
df <- ont_svs_ratio %>% filter(
  transcript_id %in% transcript_to_plot,
  layer == 'ratios_smoothed'
) %>%
  # add transcript names according to iso_to_plot (names -> values)
  mutate(isoform = factor(
    transcript_id, levels = transcript_to_plot, 
    labels = names(transcript_to_plot)
  )) %>%
  left_join(
    ont_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  facet_wrap(~ isoform, scales = 'free', nrow = 2) +
  geom_point(alpha = 0.2, color = 'gray') + 
  stat_smooth(method = 'lm') + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'Smoothed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

ggsave(
  sprintf("%s/figures/%s_vs_Rbfox3.png", res_sc_dir, gene_to_plot), 
  p2, width = 2.5, height = 4.5
)

# ratio vs Rbfox3 expression in single-cell data
# 18111//Nnat,chr2,157560109,157562106,CA-18111-3151-4103-4179-4432[INC][7/61][UPT][DNT]
# 18111//Nnat,chr2,157560109,157562106,CA-18111-3151-4103-4130-4432[INC][3/61][UPT][DNT]
# 18111//Nnat,chr2,157560109,157562106,CA-18111-3151-4103-4184-4432[INC][126/61][UPT][DNT]
# 18111//Nnat,chr2,157560109,157562106,CA-18111-3151-3338-3426-4432[INC][2/61][UPT][DNT]
df <- ct_exon_psi %>% filter(
  gene == gene_to_plot,
  event %in% c('CA-18111-3151-4103-4184-4432[INC][126/61][UPT][DNT]')
) %>%
  mutate(
    is_neuron = (!NP %in% c('CR_Meis2', 'Non Neuronal')),
    event_name = factor(
      event, levels = c('CA-18111-3151-4103-4184-4432[INC][126/61][UPT][DNT]'),
      labels = c('Exon 2')
    )
  )
p3 <- ggplot(df, aes(x = Rbfox3, y = psi)) + 
  facet_wrap(~ event_name, scales = 'free', nrow = 2) +
  geom_point(aes(color = NP), size = 0.5) + 
  stat_smooth(method = 'lm', color = 'blue') +
  stat_smooth(
    data = df %>% filter(is_neuron), method = 'lm', color = '#8B0000'
  ) + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'PSI (inclusion ratio)', color = 'Cell type') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'right',
  )
p3
ggsave(
  sprintf("%s/figures/%s_tasic_Rbfox3.png", res_sc_dir, gene_to_plot), 
  p3, width = 4, height = 3
)


### TREND examples
source("~/Projects/SPLISOSM_paper/scripts/visium_mouse_cbs/visualization/plotgardener_peak.R")

## Load RBP expression and event ratio
trend_rbp_expr <- read_csv(sprintf("%s/figures/source_data/rbp_expr.csv", res_trend_dir))

# Load cbs SVS peak bed
svs_bed <- import.bed(sprintf('%s/events/cbs_svs.peak.bed', res_trend_dir)) %>%
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
trend_svs_ratio <- read_csv(sprintf("%s/figures/source_data/svs_ratio.csv", res_trend_dir)) %>%
  left_join(svs_bed %>% dplyr::select(name, peak_name), by = c('isoform' = 'name')) %>%
  mutate(gene = str_replace(gene, 'Sept', 'Septin'))  # rename Sept to Septin

## Load transcript and RBS annotations
# load reference transcripts annotation
ensembl_gtf <- import.gff('~/reference/mm10/Mus_musculus.GRCm38.102.gtf')
ensembl_gtf <- renameSeqlevels(ensembl_gtf, value = paste0('chr', seqlevels(ensembl_gtf)))

# load FIMO RBP binding sites
fimo_dir <- sprintf('%s/events/fimo.out/', res_trend_dir)
fimo_rbp_dirs <- list.files(fimo_dir, pattern = 'cbs_svs_', full.names = FALSE)
rbs_fimo_gtf <- lapply(fimo_rbp_dirs, function(rbp_dir){
  rbp_gtf <- import.gff(sprintf('%s/%s/fimo.gff', fimo_dir, rbp_dir))
  rbp_gtf$RBP_Name <- toupper(gsub('cbs_svs_', '', rbp_dir))
  rbp_gtf$RBP_Name <- paste(rbp_gtf$RBP_Name, '(motif)')
  return(rbp_gtf)
}) %>% do.call(c, .)

# Combine all RBP binding sites
rbs_trend_gtf <- c(rbs_fimo_gtf, rbs_clip_gtf)


## Gnao1
gene_to_plot <- 'Gnao1'
# peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
peak_to_plot <- c('Gnao1-Event1', 'Gnao1-Event4', 'Gnao1-Event5')
# transcript structure
png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot), 
  width = 7.25, height = 4, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  rbs_gtf = rbs_trend_gtf,
  rbp_name_list = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                    'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(10000, 100),
)
# Highlight exon
plotRect(
  x = 2.1, y = 0.8, width = 2.2, height = 0.9,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKD – WT) in motor neuron = 0.18', fontsize = 8,
  x = 2 , y = 0.7,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- trend_svs_ratio %>% filter(
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
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p, width = 7.5, height = 2.5
)

# spatial ratio vs rbp expression
df <- trend_svs_ratio %>% filter(
  peak_name %in% peak_to_plot,
  layer == 'ratios_obs'
) %>%
  left_join(
    trend_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  facet_wrap(~ peak_name, scales = 'free', nrow = 1) +
  geom_point(alpha = 0.1, color = 'gray') + 
  stat_smooth(method = 'lm') + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'Observed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

ggsave(
  sprintf("%s/figures/%s_vs_Rbfox3.png", res_sc_dir, gene_to_plot), 
  p2, width = 7.5, height = 2.5
)


## Tln1
gene_to_plot <- 'Tln1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot), 
  width = 7.25, height = 4, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  rbs_gtf = rbs_trend_gtf,
  rbp_name_list = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                    'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(100, 3000),
)
# Highlight exon
plotRect(
  x = 6.3, y = 0.8, width = 0.2, height = 0.5,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKD – WT) in motor neuron = 0.20', fontsize = 8,
  x = 4, y = 0.7,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- trend_svs_ratio %>% filter(
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
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p, width = 5, height = 2.5
)

# spatial ratio vs rbp expression
df <- trend_svs_ratio %>% filter(
  peak_name %in% peak_to_plot,
  layer == 'ratios_obs'
) %>%
  left_join(
    trend_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  facet_wrap(~ peak_name, scales = 'free', nrow = 1) +
  geom_point(alpha = 0.1, color = 'gray') + 
  stat_smooth(method = 'lm') + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'Observed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

ggsave(
  sprintf("%s/figures/%s_vs_Rbfox3.png", res_sc_dir, gene_to_plot), 
  p2, width = 5, height = 2.5
)


## Klc1
gene_to_plot <- 'Klc1'
peak_to_plot <- svs_bed %>% filter(gene == gene_to_plot) %>% pull(peak_name)
# transcript structure
png(
  sprintf("%s/figures/%s_event_struct.png", res_sc_dir, gene_to_plot), 
  width = 7.25, height = 4, units = "in", res = 300
)
plot_peaks(
  gene_name = gene_to_plot, 
  ref_gtf = ensembl_gtf, 
  peak_bed = svs_bed %>% filter(gene == gene_to_plot),
  rbs_gtf = rbs_trend_gtf,
  rbp_name_list = c('RBFOX2 (CLIP)', 'RBFOX3 (CLIP)', 'CELF4 (CLIP)', 
                    'RBFOX3 (motif)', 'CELF4 (motif)', 'QK (motif)'),
  peak_name_list = peak_to_plot,
  zoom_in_peak = T, collapse_peak = F,
  peak_padding = c(6000, 100),
)
# Highlight exon
plotRect(
  x = 2.6, y = 0.7, width = 0.6, height = 0.9,
  just = c("left", "top"), default.units = "inches",
  lwd = 2, fill = '#FFD700', alpha = 0.4, linecolor = NA
)
plotText(
  label = 'dPSI (Rbfox tKD – WT) in motor neuron = 0.12', fontsize = 8,
  x = 2 , y = 0.65,
  just = "left", default.units = "inches"
)
pageGuideHide()
dev.off()

# spatial peak usage
p_list <- c()
for (i in seq_along(peak_to_plot)){
  p <- trend_svs_ratio %>% filter(
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
ggsave(
  sprintf("%s/figures/%s_sp_ratio.png", res_sc_dir, gene_to_plot), 
  p, width = 5, height = 2.5
)

# spatial ratio vs rbp expression
df <- trend_svs_ratio %>% filter(
  peak_name %in% peak_to_plot,
  layer == 'ratios_obs'
) %>%
  left_join(
    trend_rbp_expr %>% filter(layer == 'log1p', gene == 'Rbfox3'), 
    by = c('array_row' = 'array_row', 'array_col' = 'array_col', 'barcode' = 'barcode')
  )
p2 <- ggplot(df, aes(x = expression, y = ratio)) + 
  facet_wrap(~ peak_name, scales = 'free', nrow = 1) +
  geom_point(alpha = 0.1, color = 'gray') + 
  stat_smooth(method = 'lm') + 
  labs(x = 'Rbfox3 log-normalized expression', y = 'Observed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank()
  )
p2

ggsave(
  sprintf("%s/figures/%s_vs_Rbfox3.png", res_sc_dir, gene_to_plot), 
  p2, width = 5, height = 2.5
)
