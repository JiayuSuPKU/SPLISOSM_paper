library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(tidyverse)
library(ggpubr)
library(ggrepel)
library(UpSetR)

project_dir <- "~/Projects/SPLISOSM_paper/"

# Function to read BED file and parse gene information and return a GRanges object
read_xenium_bed <- function(file_path) {
  # Read the BED file, skipping the track header
  bed_data <- read_tsv(file_path, 
                       col_names = c("chrom", "chromStart", "chromEnd", "name", 
                                     "score", "strand", "thickStart", "thickEnd", 
                                     "itemRgb", "blockCount", "blockSizes", "blockStarts", 
                                     "description"),
                       comment = "track",
                       show_col_types = FALSE)
  
  # Extract gene name and codeword from the name column
  bed_data <- bed_data %>%
    mutate(
      gene_name = sub("\\s+\\(.*", "", name),
      codeword = as.numeric(sub(".*\\|\\s+(\\d+)\\s+\\|.*", "\\1", name)),
      probeset = as.numeric(sub(".*\\|\\s+(\\d+)\\s*$", "\\1", name))
    )
  
  # Find genes with more than one codewords
  gene_codeword_counts <- bed_data %>%
    group_by(gene_name) %>%
    summarise(
      unique_codewords = n_distinct(codeword),
      total_regions = n(),
      .groups = 'drop'
    ) %>%
    filter(unique_codewords >= 2)
  
  cat(sprintf("Found %d genes with more than one codeword\n", 
              length(gene_codeword_counts$gene_name %>% unique)))
  
  # Filter Xenium data to genes with >=2 codewords
  filtered_xenium <- bed_data %>%
    filter(gene_name %in% gene_codeword_counts$gene_name)
  
  # Convert to GRanges object
  xenium_gr <- GRanges(
    seqnames = filtered_xenium$chrom,
    ranges = IRanges(start = filtered_xenium$chromStart + 1, # BED is 0-based, GRanges is 1-based
                     end = filtered_xenium$chromEnd),
    strand = filtered_xenium$strand,
    gene_name = filtered_xenium$gene_name,
    codeword = filtered_xenium$codeword,
    probeset = filtered_xenium$probeset
  )
  
  return(xenium_gr)
}

# Function to convert the multi-exon BED12 format to exon GRanges
read_exons_from_bed12 <- function(bed12_path) {
  # Import the BED12 file
  bed12_gr <- import(bed12_path, format = "BED")

  # Initialize list to store all exons
  all_exons <- list()
  
  for (i in seq_along(bed12_gr)) {
    # Get the current feature
    feature <- bed12_gr[i]
    
    # Get the blocks (exons) - these are already parsed by rtracklayer
    if (!is.null(feature$blocks) && length(feature$blocks[[1]]) > 0) {
      exon_ranges <- feature$blocks[[1]]
      
      # The blocks are relative to the feature start, so we need to shift them
      feature_start <- start(feature) - 1  # Convert to 0-based to match BED format
      
      # Create absolute exon coordinates
      exon_starts <- start(exon_ranges) + feature_start
      exon_ends <- end(exon_ranges) + feature_start
      
      # Get feature information
      chr <- as.character(seqnames(feature))
      strand_val <- as.character(strand(feature))
      feature_name <- if (!is.null(feature$name) && !is.na(feature$name)) {
        as.character(feature$name)
      } else {
        paste0("feature_", i)
      }
      
      # Create GRanges object for exons
      exons <- GRanges(
        seqnames = rep(chr, length(exon_ranges)),
        ranges = IRanges(start = exon_starts, end = exon_ends),
        strand = rep(strand_val, length(exon_ranges)),
        name = rep(feature_name, length(exon_ranges)),
        parent_feature = rep(i, length(exon_ranges)),
        exon_number = seq_along(exon_ranges)
      )
      
      all_exons[[i]] <- exons
    }
  }

  cat(sprintf("Extracted a total of %d exons from %d features\n", 
              sum(sapply(all_exons, length)), length(bed12_gr)))
  
  # Combine all exons into a single GRanges object
  if (length(all_exons) > 0) {
    return(do.call(c, all_exons))
  } else {
    return(GRanges())
  }
}


## (1) Check event overlaps between Xenium, ONT, and SR datasets
### load bed files
gr_ont <- import(
  sprintf("%s/results/sit_nar_23/transcripts/cbs/events/cbs_svs.suppa.all_exon.bed", project_dir),
  format = "BED"
)
gr_sr <- read_exons_from_bed12(
  sprintf(
    "%s/results/visium_mouse_cbs/events/cbs_svs.exon.bed", 
    project_dir
  )
)
gr_xenium <- read_xenium_bed(
  sprintf(
    "%s/data/xenium_5k_mouse_brain/xenium_mouse_5K_gene_expression_panel_probe_locations.bed",
    project_dir
  )
)

### Compare event overlaps and visualize as barplots
# ONT events compared to SR and Xenium events
p1.1 <- data.frame(
  group = c("Total ONT-variable", "Overlapped with SR", "In Xenium 5K panel"),
  count = c(length(gr_ont),
            length(unique(queryHits(findOverlaps(gr_ont, gr_sr)))),
            length(unique(queryHits(findOverlaps(gr_ont, gr_xenium, ignore.strand = TRUE))))
  )
) %>%
  ggplot(aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = 0.5) +
  scale_x_discrete(limits = c("Total ONT-variable", "Overlapped with SR", "In Xenium 5K panel")) +
  labs(
    title = "ONT exon detectability",
    x = "Group",
    y = "Number of exons",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    # legend.position = "none",
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.8),
    legend.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank()
  )

# SR events compared to ONT and Xenium events
p1.2 <- data.frame(
  group = c("Total SR-variable", "Overlapped with ONT", "In Xenium 5K panel"),
  count = c(length(gr_sr),
            length(unique(queryHits(findOverlaps(gr_sr, gr_ont)))),
            length(unique(queryHits(findOverlaps(gr_sr, gr_xenium, ignore.strand = TRUE))))
  )
) %>%
  ggplot(aes(x = group, y = count, fill = group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = count), vjust = 0.5) +
  scale_x_discrete(limits = c("Total SR-variable", "Overlapped with ONT", "In Xenium 5K panel")) +
  labs(
    title = "SR exon detectability",
    x = "Group",
    y = "Number of exons",
    fill = ""
  ) +
  theme_classic() +
  theme(
    text = element_text(size = 12),
    # legend.position = "none",
    # axis.text.x = element_text(angle = 45, hjust = 1, size = 12)
    legend.position = "inside",
    legend.position.inside = c(0.7, 0.8),
    legend.background = element_blank(),
    axis.text = element_text(size = 12),
    axis.text.x = element_blank()
  )

p1 <- cowplot::plot_grid(p1.1, p1.2, nrow = 1, align = 'hv')
p1

## (2) Compare HSIC-IR p-values across technologies
### load SV results from ONT, SR, and Xenium datasets
# load short-reads-based TREND results
df_sv_sr <- read_csv("~/Projects/SPLISOSM_paper/results/visium_mouse_cbs/sv_results/cbs_peak_sv_combined_1119.csv") %>%
  mutate(
    is_sr_svs = `padj_hsic-ir` < 0.01,
    is_sr_gene = TRUE,
  )
# load and combine ONT results
df_sv_ont1 <- read_csv("~/Projects/SPLISOSM_paper/results/sit_nar_23/sv_results/cbs1_sv_combined_1107.csv") %>%
  mutate(
    is_ont_svs = `padj_hsic-ir` < 0.05,
    is_ont_gene = TRUE,
  )
df_sv_ont2 <- read_csv("~/Projects/SPLISOSM_paper/results/sit_nar_23/sv_results/cbs2_sv_combined_1107.csv") %>%
  mutate(
    is_ont_svs = `padj_hsic-ir` < 0.05,
    is_ont_gene = TRUE
  )
df_sv_ont <- inner_join(df_sv_ont1, df_sv_ont2, by = 'gene', suffix = c('.1', '.2')) %>%
  mutate(
    is_ont_svs = is_ont_svs.1 & is_ont_svs.2,
    is_ont_gene = TRUE,
  )
# load Xenium results
df_sv_xenium <- read_csv("~/Projects/SPLISOSM_paper/results/xenium_5k_mouse_brain/sv_results/all_multi_codeword_genes.csv") %>%
  mutate(
    # is_xenium_svs = `padj_hsic-ir` < 0.01,
    is_xenium_gene = TRUE,
  )
# set the pval threshold to be the min pval-adj of negative controls and intergenic genes
# which starts with "Neg" and "Intergenic_Region"
padj_threshold <- min(
  df_sv_xenium$`padj_hsic-ir`[grepl("^NegControl", df_sv_xenium$gene)],
  df_sv_xenium$`padj_hsic-ir`[grepl("^Intergenic_Region", df_sv_xenium$gene)]
)
df_sv_xenium <- df_sv_xenium %>%
  mutate(
    is_xenium_svs = `padj_hsic-ir` < padj_threshold,
    is_xenium_gene = TRUE
  )


### Boxplot of p-values per gene groups
# Plot the SR HSIC-IR pvalue distribution of different ONT groups
p2.1 <- df_sv_sr %>%
  filter(gene %in% df_sv_ont$gene) %>%
  mutate(
    is_ont_svs = gene %in% df_sv_ont$gene[df_sv_ont$is_ont_svs],
  ) %>%
  mutate(
    # add obs number to is_ont_svs label
    x = ifelse(is_ont_svs, 
               paste0("T (n=", sum(is_ont_svs), ")"), 
               paste0("F (n=", sum(!is_ont_svs), ")"))
  ) %>%
  ggplot(aes(x = x, y = `pvalue_hsic-ir`, fill = is_ont_svs)) +
  geom_boxplot() +
  stat_compare_means() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("grey", "red")) +
  labs(
    title = "ONT vs SR",
    x = "Is ONT SVP",
    y = "HSIC-IR p-value (SR)",
    fill = "Is ONT SVP"
  ) +
  theme_classic() +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    # axis.title.x = element_blank()
  )

# plot the Xenium HSIC-IR pvalue distribution of different ONT groups
p2.2 <- df_sv_xenium %>%
  filter(gene %in% df_sv_ont$gene) %>%
  mutate(
    is_ont_svs = gene %in% df_sv_ont$gene[df_sv_ont$is_ont_svs],
  ) %>%
  mutate(
    # add obs number to is_ont_svs label
    x = ifelse(is_ont_svs, 
               paste0("T (n=", sum(is_ont_svs), ")"), 
               paste0("F (n=", sum(!is_ont_svs), ")"))
  ) %>%
  ggplot(aes(x = x, y = `pvalue_hsic-ir`, fill = is_ont_svs)) +
  geom_boxplot() +
  stat_compare_means() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("grey", "red")) +
  labs(
    title = "ONT vs Xenium",
    x = "Is ONT SVP",
    y = "HSIC-IR p-value (Xenium Prime 5K)",
    fill = "Is ONT SVP"
  ) +
  theme_classic() +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    # axis.title.x = element_blank()
  )

# plot the Xenium HSIC-IR pvalue distribution of different SR groups
p2.3 <- df_sv_xenium %>%
  filter(gene %in% df_sv_sr$gene) %>%
  mutate(
    is_sr_svs = gene %in% df_sv_sr$gene[df_sv_sr$is_sr_svs],
  ) %>%
  mutate(
    # add obs number to is_sr_svs label
    x = ifelse(is_sr_svs, 
               paste0("T (n=", sum(is_sr_svs), ")"), 
               paste0("F (n=", sum(!is_sr_svs), ")"))
  ) %>%
  ggplot(aes(x = x, y = `pvalue_hsic-ir`, fill = is_sr_svs)) +
  geom_boxplot() +
  stat_compare_means() + 
  scale_y_log10() + 
  scale_fill_manual(values = c("grey", "red")) +
  labs(
    title = "SR vs Xenium",
    x = "Is SR SVP",
    y = "HSIC-IR p-value (Xenium Prime 5K)",
    fill = "Is SR SVP"
  ) +
  theme_classic() +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    # axis.title.x = element_blank()
  )

p2 <- cowplot::plot_grid(p2.1, p2.2, p2.3, nrow = 1)
p2

p <- cowplot::plot_grid(p1, NULL, p2, nrow = 1, rel_widths = c(1, 0.1, 1.5))
p

# save to PNG
ggsave(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/tech_comparison.png", project_dir), 
  p, width = 16, height = 3.5
)

ggsave(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/tech_comparison.pdf", project_dir), 
  p, width = 16, height = 3.5
)

###
# Boxplots of Xenium NegControl and Intergenic_Region genes
p3 <- df_sv_xenium %>%
  # group genes into NegControl, Intergenic_Region, and others
  mutate(
    gene_type = case_when(
      grepl("^NegControl", gene) ~ "NegControl",
      grepl("^Intergenic_Region", gene) ~ "Intergenic",
      TRUE ~ "Others"
    ),
    gene_type = factor(gene_type, levels = c("Others", "Intergenic", "NegControl"))
  ) %>%
  ggplot(aes(x = gene_type, y = `padj_hsic-ir`, fill = gene_type)) +
  geom_boxplot() +
  stat_compare_means() + 
  geom_hline(yintercept = 0.01, linetype = "dashed", color = "blue") +
  geom_hline(yintercept = padj_threshold, linetype = "dashed", color = "red") +
  # add labels to the hlines
  geom_text(
    data = data.frame(
      x = c(3, 2.5), y = c(0.01/10, padj_threshold/10), 
      label = c("p-adj = 0.01", sprintf("p-adj = %.2g", padj_threshold)),
      gene_type= c(NA, NA)  # to match the data frame structure
    ),
    mapping = aes(x = x, y = y, label = label), size = 4, color = "black",
  ) +
  # scale_y_log10(limits = c(1e-40, 1)) + 
  scale_y_log10() + 
  coord_cartesian(ylim = c(1e-40, 1)) +
  labs(
    title = "Xenium artifacts",
    x = "Gene type",
    y = "HSIC-IR adjusted p-value (Xenium Prime)",
    fill = "Gene type"
  ) +
  theme_classic() +
  theme(
    # legend.position = "bottom",
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    # axis.title.x = element_blank()
  )
p3
# save to PNG
ggsave(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/xenium_neg_intergenic.png", project_dir), 
  p3, width = 3.5, height = 3.5
)
ggsave(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/xenium_neg_intergenic.pdf", project_dir), 
  p3, width = 3.5, height = 3.5
)

### Scatter plots of p-values per pair of technologies
# ONT vs SR
df_sv_sr %>%
  filter(gene %in% df_sv_ont$gene) %>%
  mutate(
    is_ont_svs = gene %in% df_sv_ont$gene[df_sv_ont$is_ont_svs]
  ) %>%
  ggplot(aes(x = `padj_hsic-ir`, y = df_sv_ont$`padj_hsic-ir.2`[match(gene, df_sv_ont$gene)], color = is_ont_svs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "gray") +
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = 0.01, linetype = "dotted", color = "red") +
  scale_x_continuous(trans = c('log10', 'reverse')) +
  scale_y_continuous(trans = c('log10', 'reverse')) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    title = sprintf(
      "ONT-SR Spearman's Rho = %.2f",
      cor(df_sv_sr$`padj_hsic-ir`, df_sv_ont$`padj_hsic-ir.2`[match(df_sv_sr$gene, df_sv_ont$gene)], 
          use = "complete.obs", method = "spearman")
    ),
    x = "HSIC-IR adjusted p-value (SR)",
    y = "HSIC-IR adjusted p-value (ONT, CBS2)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12)
  )

# ONT vs Xenium
p4.1 <- df_sv_xenium %>%
  filter(gene %in% df_sv_ont$gene) %>%
  mutate(
    is_ont_svs = gene %in% df_sv_ont$gene[df_sv_ont$is_ont_svs]
  ) %>%
  ggplot(aes(x = `padj_hsic-ir` + 1e-200, y = df_sv_ont$`padj_hsic-ir.2`[match(gene, df_sv_ont$gene)], color = is_ont_svs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "gray") +
  geom_point() +
  # show the name of top 10 genes with the lowest p-values
  geom_text_repel(
    data = df_sv_ont %>%
      inner_join(df_sv_xenium, by = "gene") %>%
      filter(is_ont_svs) %>%
      arrange(`padj_hsic-ir.2`) %>%
      head(10),
    aes(label = gene),
    vjust = -1, hjust = 0.5, size = 4, color = "black", fontface = "italic"
  ) +
  # geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = padj_threshold, linetype = "dashed", color = "red") +
  scale_x_continuous(trans = c('log10', 'reverse')) +
  scale_y_continuous(trans = c('log10', 'reverse')) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    title = sprintf(
      "Spearman's Rho = %.2f",
      cor(df_sv_xenium$`padj_hsic-ir`, df_sv_ont$`padj_hsic-ir.2`[match(df_sv_xenium$gene, df_sv_ont$gene)], 
          use = "complete.obs", method = "spearman")
    ),
    x = "HSIC-IR adjusted p-value (Xenium Prime)",
    y = "HSIC-IR adjusted p-value (ONT, CBS2)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 12)
  )  

# SR vs Xenium
p4.2 <- df_sv_xenium %>%
  filter(gene %in% df_sv_sr$gene) %>%
  mutate(
    is_sr_svs = gene %in% df_sv_sr$gene[df_sv_sr$is_sr_svs],
    xenium_padj = `padj_hsic-ir` + 1e-200  # Avoid log(0) issues
  ) %>%
  ggplot(aes(x = xenium_padj, y = df_sv_sr$`padj_hsic-ir`[match(gene, df_sv_sr$gene)], color = is_sr_svs)) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "gray") +
  geom_point() +
  geom_text_repel(
    data = df_sv_sr %>%
      inner_join(df_sv_xenium, by = "gene", suffix = c('.sr', '.xenium')) %>%
      mutate(
        `padj_hsic-ir` = `padj_hsic-ir.sr` + 1e-200,  # Avoid log(0) issues
        xenium_padj = `padj_hsic-ir.xenium` + 1e-200
      ) %>%
      filter(is_sr_svs) %>%
      arrange(`padj_hsic-ir`) %>%
      head(10),
    aes(label = gene),
    vjust = -1, hjust = 0.5, size = 4, color = "black", fontface = "italic"
  ) +
  # geom_hline(yintercept = 0.05, linetype = "dotted", color = "red") +
  geom_vline(xintercept = padj_threshold, linetype = "dashed", color = "red") +
  scale_x_continuous(trans = c('log10', 'reverse')) +
  scale_y_continuous(trans = c('log10', 'reverse')) +
  scale_color_manual(values = c("grey", "red")) +
  labs(
    title = sprintf(
      "Spearman's Rho = %.2f",
      cor(df_sv_xenium$`padj_hsic-ir`, df_sv_sr$`padj_hsic-ir`[match(df_sv_xenium$gene, df_sv_sr$gene)], 
          use = "complete.obs", method = "spearman")
    ),
    x = "HSIC-IR adjusted p-value (Xenium Prime)",
    y = "HSIC-IR adjusted p-value (SR)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    text = element_text(size = 12),
    axis.text = element_text(size = 12),
    axis.text.x = element_text(angle = 270, hjust = 0, vjust = 0.5, size = 12)
  ) 

p4 <- cowplot::plot_grid(p4.1, p4.2, nrow = 1)
p4

# save to PNG
ggsave(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/pval_scatter.png", project_dir), 
  p4, width = 8, height = 4
)
ggsave(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/pval_scatter.pdf", project_dir), 
  p4, width = 8, height = 4
)

### UpSet plot for SVS genes across technologies
# Combine all results
df_combined <- full_join(df_sv_sr %>% select(gene, is_sr_svs, is_sr_gene), 
                         df_sv_ont %>% select(gene, is_ont_svs, is_ont_gene), by = 'gene') %>%
  full_join(df_sv_xenium %>% select(gene, is_xenium_svs, is_xenium_gene), by = 'gene') %>%
  mutate(
    # For genes profiled by each technology
    is_sr_gene = !is.na(is_sr_gene),
    is_ont_gene = !is.na(is_ont_gene),
    is_xenium_gene = !is.na(is_xenium_gene),
    # For SVS genes
    is_sr_svs = !is.na(is_sr_svs) & is_sr_svs,
    is_ont_svs = !is.na(is_ont_svs) & is_ont_svs,
    is_xenium_svs = !is.na(is_xenium_svs) & is_xenium_svs
  )

# Create binary matrix for UpSetR
upset_data <- list(
  `SR All` = df_combined %>% filter(is_sr_gene) %>% pull(gene),
  `ONT All` = df_combined %>% filter(is_ont_gene) %>% pull(gene),
  `Xenium All` = df_combined %>% filter(is_xenium_gene) %>% pull(gene),
  `SR SVP` = df_combined %>% filter(is_sr_svs) %>% pull(gene),
  `ONT SVP` = df_combined %>% filter(is_ont_svs) %>% pull(gene),
  `Xenium SVP` = df_combined %>% filter(is_xenium_svs) %>% pull(gene)
)

# Create the UpSet plot
p5 <- upset(
  fromList(upset_data),
  # sets = c("ONT All", "SR All", "Xenium All", "ONT SVS", "SR SVS", "Xenium SVS"),
  sets = c("ONT SVP", "SR SVP", "Xenium SVP"),
  order.by = c("freq", "degree"),
  keep.order = TRUE,
  text.scale = 2,
  mb.ratio = c(0.6, 0.4),
  # group.by = "sets",
  # query.legend = "top",
  queries = list(
    list(
      query = intersects,
      params = list("ONT SVP", "SR SVP", "Xenium SVP"), color = "red", active = T,
      query.name = "Detected by more than one technology"
    ),
    list(
      query = intersects,
      params = list("ONT SVP", "SR SVP"), color = "red", active = T,
      query.name = "Detected by more than one technology"
    ),
    list(
      query = intersects,
      params = list("SR SVP", "Xenium SVP"), color = "red", active = T,
      query.name = "Detected by more than one technology"
    ),
    list(
      query = intersects,
      params = list("ONT SVP","Xenium SVP"), color = "red", active = T,
      query.name = "Detected by more than one technology"
    )
  )
)
p5
# save to PNG
png(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/upset_svs.png", project_dir), 
  width = 8, height = 4, units = "in", res = 300
)
show(p5)
dev.off()

pdf(
  sprintf("%s/results/xenium_5k_mouse_brain/figures/upset_svs.pdf", project_dir), 
  width = 8, height = 4
)
show(p5)
dev.off()

