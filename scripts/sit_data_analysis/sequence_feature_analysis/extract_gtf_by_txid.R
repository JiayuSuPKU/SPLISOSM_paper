#!/usr/bin/env Rscript

# Description: Extracts GTF entries from a reference GTF file based on a list of transcript IDs.
# Usage: Rscript extract_gtf_by_txid.R <txid_csv> [reference_gtf]
# Arguments:
#   <txid_csv>: Path to a CSV file with a column 'transcript_id' containing transcript IDs.
#   <reference_gtf>: Path to the reference GTF file.
# Output: (in the same directory as the input CSV file)
#   A GTF file containing the entries from the reference GTF file that match the transcript IDs in the CSV file.

suppressMessages(library(rtracklayer))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: extract_gtf_by_txid.R <txid_csv> <reference_gtf>")
}

txid_dir <- args[1] # path to txid_csv with column 'transcript_id'
ref_gtf_dir <- args[2] # path to reference GTF file
out_gtf_dir <- sprintf("%s.gtf", tools::file_path_sans_ext(txid_dir))

message("=== Extracting GTF entries by transcript ID ===")
message("TXID CSV: ", txid_dir)
message("Reference GTF: ", ref_gtf_dir)
message("Output GTF: ", out_gtf_dir)

# Load data
ref_gtf <- import(ref_gtf_dir,format = 'gtf')
txid <- read.csv(txid_dir, header = TRUE, stringsAsFactors = FALSE)

# Switch chromosome naming convention (1 -> chr1)
ref_gtf <- renameSeqlevels(ref_gtf, value = paste0('chr', seqlevels(ref_gtf)))

# Remove transcript ID suffix (version)
ref_gtf$transcript_id_clean <- sub("\\..*$", "", ref_gtf$transcript_id)
txid$transcript_id_clean <- sub("\\..*$", "", txid$transcript_id)

# Find matching transcript IDs and write output GTF
out_gtf <- ref_gtf[ref_gtf$transcript_id_clean %in% txid$transcript_id_clean, ]
out_gtf$transcript_id_clean <- NULL
export(out_gtf, out_gtf_dir, format = "gtf")
