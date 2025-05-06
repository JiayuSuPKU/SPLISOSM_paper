#!/usr/bin/env R

# ****** NMD detector based on the position of the premature termination codon (PTC)
# https://github.com/TheJacksonLaboratory/BRCA-LRseq-pipeline/blob/main/2_ORF_annotation/2.6-BRCA_isoforms.R

suppressPackageStartupMessages({
	library(tidyverse)
	library(rtracklayer)
	library(IRanges)
	library(parallel)
	library(optparse)
})

option_list <- list(
	make_option(
		c("-i", "--input"), 
		type="character",
		help="Path to the GTF file with CDS predicted by ORFanage to annotate"),
	make_option(
		c("-o", "--output"), 
		type="character",
		help="Path to save the GTF with NMD annotation")
)

# Parse arguments
opt <- parse_args(OptionParser(
	option_list=option_list,
	description="Detect nonsense-mediated decay (NMD) events using the 55nt rule."
))

# Transcript gff files with CDS annotation predicted by ORFanage
input_gtf <- opt$input

# Output annotated gtf file
output_gtf <- opt$output

# Load CDS coordinates
all_gtf <- import.gff(input_gtf)
cds_gtf <- all_gtf[all_gtf$type == "CDS"]

# Keep transcripts with complete ORF
t_ids <- cds_gtf$transcript_id %>% unique()
exon_gtf <- all_gtf[(all_gtf$transcript_id %in% t_ids) & (all_gtf$type == "exon")]

# Group GRange records by transcript
exon_by_transcript <- split(exon_gtf, exon_gtf$transcript_id)
cds_by_transcript <- split(cds_gtf, cds_gtf$transcript_id)

# Predict whether each ORF is productive or not
# NMD: yes if stop codon locates > 55nt upstream the 3' most junction
res <- mclapply(t_ids, function(id){
	# exon locations for one ORF
	exon <- exon_by_transcript[[id]]

	# CDS locations for one ORF, split by exons
	cds <- cds_by_transcript[[id]]
	
	# find last junction and stop codon positions
	if (as.character(strand(exon))[1] == '-'){
		# negative strand
		last_junc_pos <- start(exon)[if(length(exon) > 1) 2 else 1]
		stop_codon_pos <- start(cds)[1]
		dist_to_junc <- stop_codon_pos - last_junc_pos
	} else{
		# positive strand
		last_junc_pos <- end(exon)[if(length(exon) > 1) length(exon) - 1 else 1]
		stop_codon_pos <- end(cds)[length(cds)]
		dist_to_junc <- last_junc_pos - stop_codon_pos
	}
	
	# NMD: yes if stop codon locates > 55nt upstream the 3' most junction
	if (dist_to_junc > 55){
		NMD <- "True"
	} else{
		NMD <- "False"
	}

	data.frame(
		iso_id = id,
		NMD_pred = NMD,
		last_junc_pos = last_junc_pos,
		stop_codon_pos = stop_codon_pos
	) %>% return()

}, mc.cores = detectCores()) %>% do.call(what = 'rbind')

# Check warnings
warnings()

# Add NMD annotation to the GTF
# Match NMD predictions back to the GTF by transcript_id
nmd_map <- res %>% select(iso_id, NMD_pred) %>% deframe()
all_gtf$NMD_prediction <- nmd_map[as.character(all_gtf$transcript_id)]

# Save the updated GTF with NMD annotation
export.gff(all_gtf, con = output_gtf)

# Print summary statistics
cat(sprintf("Number of transcripts processed: %d\n", nrow(res)))
cat("NMD predictions:\n")
print(table(res$NMD_pred))