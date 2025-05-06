#!/usr/bin/env Rscript

# Description: Extracts the entire region of genes from a reference GTF file based on a list of gene symbols.
# Usage: Rscript extract_region_by_symbol.R <gene_symbol_file> <gtf_file> <output_bed_file>
# Arguments:
#   <gene_symbol_file>: Path to a file containing a list of gene symbols. Single column, no header.
#   <gtf_file>: Path to the reference GTF file.
#   <output_bed_file>: Path to the output BED file.
# Output:
#   <output_bed_file>: A BED file containing the entire region of the genes.

suppressMessages(library(tidyverse))
suppressMessages(library(rtracklayer))

# Parse the command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
  stop("Usage: extract_region_by_symbol.R <gene_symbol_file> <gtf_file> <output_bed_file>")
}

gene_symbol_file <- args[1]
gtf_file <- args[2]
output_bed_file <- args[3]

# Load the reference GTF file
ref_gtf <- import.gff(gtf_file)
ref_gtf <- renameSeqlevels(ref_gtf, value = paste0('chr', seqlevels(ref_gtf)))

# Load the SVS gene list
svs_genes <- read.table(
  gene_symbol_file,
  col.names = "gene_symbol",
  header = FALSE,
  stringsAsFactors = FALSE
)

# Save the entire region to BED file
svs_region_bed <- ref_gtf %>% 
  as.data.frame() %>%
  dplyr::filter(gene_name %in% (svs_genes %>% pull(gene_symbol)), type == "transcript") %>%
  group_by(gene_id) %>%
  summarise(seqnames = seqnames[1], start = min(start), end = max(end), strand = strand[1]) %>%
  mutate(score = ".") %>%
  dplyr::select(seqnames, start, end, gene_id, score, strand)

write.table(
  svs_region_bed,
  output_bed_file,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)