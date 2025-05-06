library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)
library(plotgardener)
library(org.Mm.eg.db)

### Visualization
# Function to plot transcripts and RBP binding sites
# given a list of transcript IDs, visualize them using plotgardener
plot_transcripts <- function(
  transcript_ids, ref_gtf, 
  rbp_names = NULL, rbs_gtf = NULL, 
  fill = '#008080'
) {
  # Remove suffix from transcript IDs if present
  transcript_ids <- gsub("\\..*$", "", transcript_ids)
  
  # Select a subset of transcripts
  subset_gr <- ref_gtf[ref_gtf$transcript_id %in% transcript_ids]
  subset_gr$transcript_id <- subset_gr$transcript_name
  n_isoforms <- length(unique(subset_gr$transcript_id))
  
  # convert GRanges to genome assembly
  subset_as <- assembly(
    Genome = 'custom',
    TxDb = makeTxDbFromGRanges(subset_gr),
    OrgDb = org.Mm.eg.db,
    gene.id.column = "ENSEMBL",
    display.column = 'SYMBOL'
  )
  
  # Extract chrom location
  chrom <- seqnames(subset_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
  chromstart <- min(start(subset_gr))
  chromend <- max(end(subset_gr))
  strand <- strand(subset_gr)[1] %>% as.character()
  
  # Create page
  pageCreate(width = 7.25, height = 2, default.units = "inches")
  
  # Plot and place transcripts
  plotTranscripts(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as, 
    labels = 'transcript', fill = fill, 
    stroke = 0.5,
    x = 1, y = 0, width = 6, height = 0.3 * n_isoforms,
    just = c("left", "top"), default.units = "inches",
  )
  
  # Plot RBP binding sites
  if (!is.null(rbs_gtf)){
    if (is.null(rbp_names)){
      rbp_names <- unique(rbs_gtf$RBP_Name) # Plot all RBPs if not specified
    }
    
    # Loop over each RBP and plot the binding sites
    for (i in seq_along(rbp_names)){
      # Plot RBS
      plotRanges(
        data = rbs_gtf[rbs_gtf$RBP_Name == rbp_names[i]],
        chrom = chrom, chromstart = chromstart, chromend = chromend,
        assembly = subset_as,
        collapse = TRUE,
        fill = '#8E44AD', linecolor = 'fill',
        x = 1, y =  0.3 * n_isoforms + 0.15 * i + 0.05, 
        width = 6, height = 0.10, 
        just = c("left", "top"), default.units = "inches",
      )
      # Plot RBP name
      plotText(
        label = rbp_names[i], fontsize = 8,
        x = 0.2, y = 0.3 * n_isoforms + 0.15 * i + 0.075,
        just = "left", default.units = "inches"
      )
      # Plot horizontal line
      plotSegments(
        x0 = 1, x1 = 7, 
        y0 = 0.3 * n_isoforms + 0.15 * (i + 1), y1 = 0.3 * n_isoforms + 0.15 * (i + 1),
        linecolor = "gray"
      )
    }  
  }
  n_rbps <- length(rbp_names) # 0 if rbs_gtf is NULL
  
  # Plot genome scale
  plotGenomeLabel(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as, fontsize = 8,
    x = 1, y = 0.3 * n_isoforms + 0.15 * (n_rbps + 1) + 0.05,
    length = 6, default.units = "inches",
  )
  
  # Plot strand info
  plotText(
    label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')"), fontsize = 8,
    x = 0.2, y = 0.1,
    just = "left", default.units = "inches"
  )
}

### Plot with reduced width
plot_transcripts_short <- function(
    transcript_ids, ref_gtf, 
    rbp_names = NULL, rbs_gtf = NULL, 
    fill = '#008080'
) {
  # Remove suffix from transcript IDs if present
  transcript_ids <- gsub("\\..*$", "", transcript_ids)
  
  # Select a subset of transcripts
  subset_gr <- ref_gtf[ref_gtf$transcript_id %in% transcript_ids]
  subset_gr$transcript_id <- subset_gr$transcript_name
  n_isoforms <- length(unique(subset_gr$transcript_id))
  
  # convert GRanges to genome assembly
  subset_as <- assembly(
    Genome = 'custom',
    TxDb = makeTxDbFromGRanges(subset_gr),
    OrgDb = org.Mm.eg.db,
    gene.id.column = "ENSEMBL",
    display.column = 'SYMBOL'
  )
  
  # Extract chrom location
  chrom <- seqnames(subset_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
  chromstart <- min(start(subset_gr))
  chromend <- max(end(subset_gr))
  strand <- strand(subset_gr)[1] %>% as.character()
  
  # Create page
  pageCreate(width = 5.25, height = 2, default.units = "inches")
  
  # Plot and place transcripts
  plotTranscripts(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as, 
    labels = 'transcript', fill = fill, 
    stroke = 0.5,
    x = 1, y = 0, width = 4, height = 0.3 * n_isoforms,
    just = c("left", "top"), default.units = "inches",
  )
  
  # Plot RBP binding sites
  if (!is.null(rbs_gtf)){
    if (is.null(rbp_names)){
      rbp_names <- unique(rbs_gtf$RBP_Name) # Plot all RBPs if not specified
    }
    
    # Loop over each RBP and plot the binding sites
    for (i in seq_along(rbp_names)){
      # Plot RBS
      plotRanges(
        data = rbs_gtf[rbs_gtf$RBP_Name == rbp_names[i]],
        chrom = chrom, chromstart = chromstart, chromend = chromend,
        assembly = subset_as,
        collapse = TRUE,
        fill = '#8E44AD', linecolor = 'fill',
        x = 1, y =  0.3 * n_isoforms + 0.15 * i + 0.05, 
        width = 4, height = 0.10, 
        just = c("left", "top"), default.units = "inches",
      )
      # Plot RBP name
      plotText(
        label = rbp_names[i], fontsize = 8,
        x = 0.2, y = 0.3 * n_isoforms + 0.15 * i + 0.075,
        just = "left", default.units = "inches"
      )
      # Plot horizontal line
      plotSegments(
        x0 = 1, x1 = 5, 
        y0 = 0.3 * n_isoforms + 0.15 * (i + 1), y1 = 0.3 * n_isoforms + 0.15 * (i + 1),
        linecolor = "gray"
      )
    }  
  }
  n_rbps <- length(rbp_names) # 0 if rbs_gtf is NULL
  
  # Plot genome scale
  plotGenomeLabel(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as, fontsize = 8,
    x = 1, y = 0.3 * n_isoforms + 0.15 * (n_rbps + 1) + 0.05,
    length = 4, default.units = "inches",
  )
  
  # Plot strand info
  plotText(
    label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')"), fontsize = 8,
    x = 0.2, y = 0.1,
    just = "left", default.units = "inches"
  )
}
