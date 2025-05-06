library(tidyverse)
library(GenomicFeatures)
library(rtracklayer)
library(plotgardener)
library(org.Mm.eg.db)

### Visualization
# Function to plot transcripts and RBP binding sites
# given a list of transcript IDs, visualize them using plotgardener
plot_peaks <- function(
    gene_name,
    ref_gtf,
    peak_bed,
    rbs_gtf = NULL,
    rbp_name_list = NULL,
    ref_gtf_gene_id = 'ENSEMBL',
    transcript_name_list = NULL,
    peak_name_list = NULL,
    zoom_in_peak = TRUE,
    collapse_peak = FALSE,
    label_peak = TRUE,
    peak_padding = c(1000, 1000),
    fill = '#0066CC'
){
  # extract transcripts for the gene
  subset_gr <- ref_gtf[ref_gtf$gene_name == gene_name]
  if (!is.null(transcript_name_list)){
    subset_gr <- subset_gr[subset_gr$transcript_name %in% transcript_name_list]
  }
  subset_gr$transcript_id <- subset_gr$transcript_name
  subset_gr$phase[subset_gr$type == 'stop_codon'] <- NA

  # focus on selected peaks
  if (!is.null(peak_name_list)){
    peak_bed <- peak_bed %>% filter(peak_bed$peak_name %in% peak_name_list)
  }
  peak_name_list <- peak_bed$peak_name
  n_peaks <- length(unique(peak_bed$peak_name))
  peakstart <- min(peak_bed$start) - peak_padding[1]
  peakend <- max(peak_bed$end) + peak_padding[2]
  
  # focus on the peak region
  if (zoom_in_peak){
    # remove transcripts that do not overlap with the peak
    t_keep <- subset_gr[
      (subset_gr$type == 'transcript') & (end(subset_gr) >= peakstart) & (start(subset_gr) <= peakend)
    ]$transcript_id
    subset_gr <- subset_gr[subset_gr$transcript_id %in% t_keep]
  }
  n_isoforms <- length(unique(subset_gr$transcript_id))

  # make sure there are transcripts and peaks to visualize
  stopifnot(n_isoforms > 0, n_peaks > 0)
  
  # create assembly object
  subset_as <- assembly(
    Genome = 'custom',
    TxDb = makeTxDbFromGRanges(subset_gr),
    OrgDb = org.Mm.eg.db,
    gene.id.column = ref_gtf_gene_id,
    display.column = 'SYMBOL'
  )
  
  # genomic coordinates
  chrom <- seqnames(subset_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
  chromstart <- min(start(subset_gr))
  chromend <- max(end(subset_gr))
  strand <- strand(subset_gr)[1] %>% as.character()

  if (zoom_in_peak){
    chromstart <- peakstart
    chromend <- peakend
  }

  # create new page
  pageCreate(width = 7.25, height = 3, default.units = "inches")
  
  # plot and place transcripts
  plotTranscripts(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as,
    labels = NULL, fill = fill, 
    stroke = 0.5,
    x = 1, y = 0, width = 6, height = 0.15 * n_isoforms,
    just = c("left", "top"), default.units = "inches",
  )
  y_interval = 0.15 * n_isoforms
  
  # plot RNA peaks
  if (collapse_peak){ # all peaks in one track
    n_peaks <- 1
    plotRanges(
      peak_bed,
      chrom = chrom, chromstart = chromstart, chromend = chromend,
      assembly = subset_as, collapse = TRUE, fill = '#E23D28',
      x = 1, y = y_interval + 0.1, width = 6, height = 0.1,
    )
    plotText(
      label = 'Events', fontsize = 8,
      x = 1, y = y_interval + 0.15,
      just = "right", default.units = "inches"
    )
    # add label under the peaks
    if (label_peak){
      for (peak in peak_bed$peak_name){
        start = peak_bed %>% filter(peak_name == peak) %>% pull(start)
        end = peak_bed %>% filter(peak_name == peak) %>% pull(end)
        plotText(
          label = peak, fontsize = 8,
          x = 1 + (start - chromstart) / (chromend - chromstart) * 6,
          y = y_interval + 0.25,
          just = "left", default.units = "inches"
        )
      }
      y_interval <- y_interval + 0.05
    }
  } else{ # peaks in separate tracks
    for (i in seq_along(peak_name_list)){
      # plot peak
      plotRanges(
        peak_bed %>% filter(peak_name == peak_name_list[i]),
        chrom = chrom, chromstart = chromstart, chromend = chromend,
        assembly = subset_as, collapse = TRUE, fill = '#E23D28',
        x = 1, y = y_interval + 0.15*i + 0.05, width = 6, height = 0.1,
      )
      # plot peak name
      plotText(
        label = str_split_1(peak_name_list[i], '-')[2], fontsize = 8,
        x = 1, y = y_interval + 0.15*i + 0.075,
        just = "right", default.units = "inches"
      )
      # plot horizontal line
      plotSegments(
        x0 = 1, x1 = 7, 
        y0 = y_interval + 0.15*(i + 1), y1 = y_interval + 0.15*(i + 1),
        linecolor = "gray"
      )
    }
  }
  y_interval <- y_interval + 0.15 * n_peaks
  
  # Plot RBP binding sites
  if (!is.null(rbs_gtf)){
    if (is.null(rbp_name_list)){
      rbp_name_list <- unique(rbs_gtf$RBP_Name) # Plot all RBPs if not specified
    }

    # Loop over each RBP and plot the binding sites
    for (i in seq_along(rbp_name_list)){
      # Plot RBS
      plotRanges(
        data = rbs_gtf[rbs_gtf$RBP_Name == rbp_name_list[i]],
        chrom = chrom, chromstart = chromstart, chromend = chromend,
        assembly = subset_as, collapse = TRUE,
        fill = '#8E44AD', linecolor = 'fill',
        x = 1, y = y_interval + 0.15*i + 0.05, 
        width = 6, height = 0.10, 
        just = c("left", "top"), default.units = "inches",
      )
      # Plot RBP name
      plotText(
        label = rbp_name_list[i], fontsize = 8,
        x = 1, y = y_interval + 0.15*i + 0.075,
        just = "right", default.units = "inches"
      )
      # Plot horizontal line
      plotSegments(
        x0 = 1, x1 = 7, 
        y0 = y_interval + 0.15*(i + 1), 
        y1 = y_interval + 0.15*(i + 1),
        linecolor = "gray"
      )
    }
    y_interval <- y_interval + 0.15*length(rbp_name_list)
  }

  y_interval <- y_interval + 0.2

  # Plot genome scale
  plotGenomeLabel(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as, fontsize = 8,
    x = 1, y = y_interval,
    length = 6, default.units = "inches",
  )
  # Plot strand info
  label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
  label = paste0(gene_name, ', ', label)
  plotText(
    label = label, fontsize = 8,
    x = 0.3, y = 0.3,
    just = "left", default.units = "inches"
  )
}

## Reduce plot width
plot_peaks_short <- function(
    gene_name,
    ref_gtf,
    peak_bed,
    rbs_gtf = NULL,
    rbp_name_list = NULL,
    ref_gtf_gene_id = 'ENSEMBL',
    transcript_name_list = NULL,
    peak_name_list = NULL,
    zoom_in_peak = TRUE,
    collapse_peak = FALSE,
    label_peak = TRUE,
    peak_padding = c(1000, 1000),
    fill = '#0066CC'
){
  # extract transcripts for the gene
  subset_gr <- ref_gtf[ref_gtf$gene_name == gene_name]
  if (!is.null(transcript_name_list)){
    subset_gr <- subset_gr[subset_gr$transcript_name %in% transcript_name_list]
  }
  subset_gr$transcript_id <- subset_gr$transcript_name
  subset_gr$phase[subset_gr$type == 'stop_codon'] <- NA
  
  # focus on selected peaks
  if (!is.null(peak_name_list)){
    peak_bed <- peak_bed %>% filter(peak_bed$peak_name %in% peak_name_list)
  }
  peak_name_list <- peak_bed$peak_name
  n_peaks <- length(unique(peak_bed$peak_name))
  peakstart <- min(peak_bed$start) - peak_padding[1]
  peakend <- max(peak_bed$end) + peak_padding[2]
  
  # focus on the peak region
  if (zoom_in_peak){
    # remove transcripts that do not overlap with the peak
    t_keep <- subset_gr[
      (subset_gr$type == 'transcript') & (end(subset_gr) >= peakstart) & (start(subset_gr) <= peakend)
    ]$transcript_id
    subset_gr <- subset_gr[subset_gr$transcript_id %in% t_keep]
  }
  n_isoforms <- length(unique(subset_gr$transcript_id))
  
  # make sure there are transcripts and peaks to visualize
  stopifnot(n_isoforms > 0, n_peaks > 0)
  
  # create assembly object
  subset_as <- assembly(
    Genome = 'custom',
    TxDb = makeTxDbFromGRanges(subset_gr),
    OrgDb = org.Mm.eg.db,
    gene.id.column = ref_gtf_gene_id,
    display.column = 'SYMBOL'
  )
  
  # genomic coordinates
  chrom <- seqnames(subset_gr)[1] %>% as.character() # assuming all transcripts are on the same chromosome
  chromstart <- min(start(subset_gr))
  chromend <- max(end(subset_gr))
  strand <- strand(subset_gr)[1] %>% as.character()
  
  if (zoom_in_peak){
    chromstart <- peakstart
    chromend <- peakend
  }
  
  # create new page
  pageCreate(width = 5.25, height = 3, default.units = "inches")
  
  # plot and place transcripts
  plotTranscripts(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as,
    labels = NULL, fill = fill, 
    stroke = 0.5,
    x = 1, y = 0, width = 4, height = 0.15 * n_isoforms,
    just = c("left", "top"), default.units = "inches",
  )
  y_interval = 0.15 * n_isoforms
  
  # plot RNA peaks
  if (collapse_peak){ # all peaks in one track
    n_peaks <- 1
    plotRanges(
      peak_bed,
      chrom = chrom, chromstart = chromstart, chromend = chromend,
      assembly = subset_as, collapse = TRUE, fill = '#E23D28',
      x = 1, y = y_interval + 0.1, width = 4, height = 0.1,
    )
    plotText(
      label = 'Events', fontsize = 8,
      x = 1, y = y_interval + 0.15,
      just = "right", default.units = "inches"
    )
    # add label under the peaks
    if (label_peak){
      for (peak in peak_bed$peak_name){
        start = peak_bed %>% filter(peak_name == peak) %>% pull(start)
        end = peak_bed %>% filter(peak_name == peak) %>% pull(end)
        plotText(
          label = peak, fontsize = 8,
          x = 1 + (start - chromstart) / (chromend - chromstart) * 6,
          y = y_interval + 0.25,
          just = "left", default.units = "inches"
        )
      }
      y_interval <- y_interval + 0.05
    }
  } else{ # peaks in separate tracks
    for (i in seq_along(peak_name_list)){
      # plot peak
      plotRanges(
        peak_bed %>% filter(peak_name == peak_name_list[i]),
        chrom = chrom, chromstart = chromstart, chromend = chromend,
        assembly = subset_as, collapse = TRUE, fill = '#E23D28',
        x = 1, y = y_interval + 0.15*i + 0.05, width = 4, height = 0.1,
      )
      # plot peak name
      plotText(
        label = str_split_1(peak_name_list[i], '-')[2], fontsize = 8,
        x = 1, y = y_interval + 0.15*i + 0.075,
        just = "right", default.units = "inches"
      )
      # plot horizontal line
      plotSegments(
        x0 = 1, x1 = 5, 
        y0 = y_interval + 0.15*(i + 1), y1 = y_interval + 0.15*(i + 1),
        linecolor = "gray"
      )
    }
  }
  y_interval <- y_interval + 0.15 * n_peaks
  
  # Plot RBP binding sites
  if (!is.null(rbs_gtf)){
    if (is.null(rbp_name_list)){
      rbp_name_list <- unique(rbs_gtf$RBP_Name) # Plot all RBPs if not specified
    }
    
    # Loop over each RBP and plot the binding sites
    for (i in seq_along(rbp_name_list)){
      # Plot RBS
      plotRanges(
        data = rbs_gtf[rbs_gtf$RBP_Name == rbp_name_list[i]],
        chrom = chrom, chromstart = chromstart, chromend = chromend,
        assembly = subset_as, collapse = TRUE,
        fill = '#8E44AD', linecolor = 'fill',
        x = 1, y = y_interval + 0.15*i + 0.05, 
        width = 4, height = 0.10, 
        just = c("left", "top"), default.units = "inches",
      )
      # Plot RBP name
      plotText(
        label = rbp_name_list[i], fontsize = 8,
        x = 1, y = y_interval + 0.15*i + 0.075,
        just = "right", default.units = "inches"
      )
      # Plot horizontal line
      plotSegments(
        x0 = 1, x1 = 5, 
        y0 = y_interval + 0.15*(i + 1), 
        y1 = y_interval + 0.15*(i + 1),
        linecolor = "gray"
      )
    }
    y_interval <- y_interval + 0.15*length(rbp_name_list)
  }
  
  y_interval <- y_interval + 0.2
  
  # Plot genome scale
  plotGenomeLabel(
    chrom = chrom, chromstart = chromstart, chromend = chromend,
    assembly = subset_as, fontsize = 8,
    x = 1, y = y_interval,
    length = 4, default.units = "inches",
  )
  # Plot strand info
  label = ifelse(strand == '+', "+ strand (5'>3')", "- strand (3'<5')")
  label = paste0(gene_name, ', ', label)
  plotText(
    label = label, fontsize = 8,
    x = 0.3, y = 0.3,
    just = "left", default.units = "inches"
  )
}
