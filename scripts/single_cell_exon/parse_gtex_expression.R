# Load necessary libraries
library(jsonlite)
library(httr)
library(dplyr)
library(ggpubr)

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
  mutate(isBrainTissue = grepl('Brain', tissueSiteDetailId)) %>%
  dplyr::select(-exonId) %>%
  pivot_wider(names_from = exonName, values_from = median.exon)

df_exon %>%
  ggplot(aes(x = median.gene, y = Exon_16)) + 
  geom_point(aes(color = isBrainTissue)) + 
  geom_smooth(method = 'lm') +
  stat_cor() +
  theme_classic()

# junction
df_junction <- merge(
  pcbp2_gene_expr, pcbp2_junction_expr, 
  by = c('tissueSiteDetailId','ontologyId', 'datasetId', 'gencodeId', 'geneSymbol'),
  suffixes = c('.gene', '.junc')
) %>% 
  mutate(isBrainTissue = grepl('Brain', tissueSiteDetailId)) %>%
  filter(junctionId %in% c("chr12:53468833-53471637:+", "chr12:53471808-53474868:+", "chr12:53471808-53479405:+")) %>%
  pivot_wider(names_from = junctionId, values_from = median.junc) %>%
  dplyr::rename(
    junction_16 = `chr12:53468833-53471637:+`,
    junction_17 = `chr12:53471808-53474868:+`, 
    junction_18 = `chr12:53471808-53479405:+`)

df_junction %>%
  ggplot(aes(x = median.gene, y = junction_17)) + 
  geom_point(aes(color = isBrainTissue)) + 
  geom_smooth(method = 'lm') +
  stat_cor(label.y = 0.06) +
  theme_classic()


df_junction %>%
  ggplot(aes(x = median.gene, y = junction_2 / (junction_1 + junction_2))) + 
  geom_point() + 
  geom_smooth(method = 'lm') +
  stat_cor(label.y = 1) +
  theme_classic()


data_dir <- "~/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir <- "~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
date <- "1119"

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

gene_to_plot <- 'Pcbp2'
peak_to_plot <- c('Pcbp2-Event1', 'Pcbp2-Event2')
p_list <- c()

# spatial usage of selected peaks
for (i in seq_along(peak_to_plot)){
  p <- df_svs_ratio %>% filter(
    peak_name == peak_to_plot[i], 
    layer == 'counts',
    array_col > 15
  ) %>%
    ggplot(aes(x = array_col, y = - array_row, color = ratio)) + 
    geom_point(size = 0.5) +
    labs(color = '', title = peak_to_plot[i]) + 
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
  if (i == 2){
    p <- p + scale_color_gradient2(
      high = '#8B0000', mid = 'white', low = '#0066CC',
      # limits = c(0.6, 1), oob = scales::squish,
      # high = '#E23D28', mid = 'white', low = '#3F7CA4', 
      midpoint = 0.75)
  } else{
    p <- p + scale_color_gradient(low = 'white', high = '#0066CC') #'#3F7CA4'
  }
  p_list <- c(p_list, list(p))
}
m9.1 <- cowplot::plot_grid(plotlist = p_list, nrow = 1, align = 'hv')
m9.1


# Pcbp2 expression vs peak ratio
df <- df_svs_ratio %>% filter(
  gene == gene_to_plot,
  # layer == 'ratios_obs',
  layer == 'log1p',
  peak_name %in% peak_to_plot
) %>% dplyr::select(-isoform) %>%
  pivot_wider(names_from = peak_name, values_from = ratio)

ggplot(df, aes(x = `Pcbp2-Event2`, 
                       y = `Pcbp2-Event1` )) + 
  geom_point(alpha = 0.1) +
  geom_smooth(method = 'lm') + 
  labs(x = 'Pcbp2 log-normalized expression', y = 'Observed ratio') + 
  theme_classic() + 
  theme(
    text = element_text(size = 12, family = 'Arial'),
    strip.text = element_text(size = 12, family = 'Arial'),
    strip.background = element_blank(),
    legend.position = 'none'
  )



