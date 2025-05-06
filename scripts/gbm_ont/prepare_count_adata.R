#! /usr/bin/env Rscript

# Convert Seurat objects to anndata format for downstream analysis
# Dataset: 6 GBM and 5 DMG ONT samples from https://www.nature.com/articles/s41467-023-36707-6

setwd("~/Projects/SPLISOSM_paper/data/gbm_ont_nc_23")
suppressMessages({
  library(tidyverse)
  library(Seurat)
  library(SeuratDisk)
  library(rtracklayer)
  library(org.Hs.eg.db)
})

# load the reference transcriptome (hg38 gencode v32)
hg38_gtf <- import.gff("~/reference/hg38/gencode.v32.annotation.gtf.gz")

# helper function to annotate ONT features
merge_gtf_iso_meta <- function(sample_id, ont_transcript_id_list){
  # match the gtf naming convention
  sample_id = paste(strsplit(sample_id, '_')[[1]], collapse = '')

  # load the gtf file for novel transcripts
  novel_iso_gtf <- readGFF(sprintf("./new_isoform_gtf/%s_NewIsoform.gtf", sample_id))
  novel_iso <- novel_iso_gtf %>% 
    as.data.frame() %>%
    filter(type == 'transcript') %>%
    distinct(transcript_id, .keep_all=TRUE) %>%
    mutate(
      chr = seqid,
      id_to_match = transcript_id,
      gene_name = NA
    ) %>%
    dplyr::select(c(source, chr, start, end, strand, gene_id, gene_name, id_to_match))
  
  # novel transcripts are named as 'chr:start-end:strand',
  # additionally, transcripts of new genes have the prefix 'transcript:' in transcript_id
  # table((novel_iso%>% filter(str_detect(id_to_match, '^transcript')))$source)
  # table((novel_iso%>% filter(!str_detect(id_to_match, '^transcript')))$source)

  # process reference isoforms
  ref_iso <- hg38_gtf %>% 
    as.data.frame() %>%
    filter(type == 'transcript') %>%
    distinct(transcript_name, .keep_all=TRUE) %>%
    mutate(
      chr = seqnames,
      id_to_match = transcript_name
    ) %>%
    dplyr::select(c(source, chr, start, end, strand, gene_id, gene_name, id_to_match))
  
  # merge reference and novel isoforms
  all_iso_meta <- rbind(ref_iso, novel_iso)
  
  # seperate ONT features into reference, novel transcript, and novel genes
  df_iso_meta <- data.frame(
    transcript_id = ont_transcript_id_list
  ) %>% mutate(
    is_novel = str_detect(transcript_id, ':'),
    # for novel transcripts, find the corresponding transcript_id in the novel_iso
    has_matched_transcript = ifelse(is_novel, transcript_id %in% novel_iso$id_to_match, NA),
    # for reference transcripts, find the corresponding transcript_name in hg38_gtf
    has_matched_transcript = ifelse(
      is_novel, 
      has_matched_transcript,
      transcript_id %in% hg38_gtf$transcript_name
    ),
    # for novel genes, find the corresponding gene_id in the novel_iso
    has_matched_new_gene = ifelse(has_matched_transcript, NA, transcript_id %in% novel_iso$gene_id)
  ) %>% mutate(
    group = ifelse(is_novel, 'not matched', 'reference transcript'),
    group = ifelse(
      is_novel & (has_matched_transcript),
      'novel transcript',
      group
    ),
    group = ifelse(
      is_novel & (! is.na(has_matched_new_gene)) & has_matched_new_gene, 
      'novel gene',
      group
    )
  ) %>% mutate(
    # for novel genes, append 'transcript:' to transcript_id
    id_to_match = ifelse(
      group == 'novel gene',
      paste0('transcript:', transcript_id),
      transcript_id
    )
  )

  # merge with GTF info
  df_iso_meta <- df_iso_meta %>% left_join(all_iso_meta, by = 'id_to_match')

  # prepare ensembl-id-to-symbol mapping
  x <- unlist(
    df_iso_meta %>% 
      filter(group == 'novel transcript') %>% 
      drop_na(gene_id) %>%
      dplyr::select(gene_id) %>%
      unique()
  )
  id2sym <- AnnotationDbi::select(
    org.Hs.eg.db, 
    keys= x, 
    columns="SYMBOL", keytype="ENSEMBL"
  ) %>%
    distinct(ENSEMBL, .keep_all = TRUE) %>%
    column_to_rownames(var = 'ENSEMBL')
  
  # add gene symbol for novel transcripts
  df_iso_meta <- df_iso_meta %>% mutate(
    gene_name = ifelse(
      group != 'novel transcript',
      gene_name,
      id2sym[gene_id, 'SYMBOL']
    ),
  )
  
  # make sure feature order is the same
  stopifnot(nrow(df_iso_meta) == length(ont_transcript_id_list))
  df_iso_meta <- df_iso_meta %>% arrange(match(transcript_id, ont_transcript_id_list))

  return(df_iso_meta)
}

# load ONT counts of all spots
ont <- readRDS("./Isoforms_filter_Allsample_exp.RDS")
# split spot names to sample_id and barcode_id
split_names <- strsplit(colnames(ont), "_")
ont_sample_id_list <- sapply(split_names, function(x) paste(x[1], x[2], sep = "_"))
ont_barcode_id_list <- sapply(split_names, `[`, 3)

# load visium data with images and spatial coordinates
vis_list <- readRDS("./NGS_Allsample_seu.RDS")

# create folders to store processed anndata and spatial meta
dir.create("./all_spots/anndata", recursive = TRUE)
dir.create("./all_spots/visium_spatial_meta", recursive = TRUE)

# iterate through samples
# sample_id <- "DMG_2"
for (sample_id in sort(names(vis_list))){
  cat(sprintf("Preparing anndata for %s ...\n", sample_id))

  # prepare Visium slide seurat obj
  vis_s <- vis_list[[sample_id]]
  vis_s <- RenameCells(vis_s, new.names = str_split_i(colnames(vis_s), '_', 1))

  # prepare ONT slide seurat obj
  ont_s <- ont[, ont_sample_id_list == sample_id]
  colnames(ont_s) <- ont_barcode_id_list[ont_sample_id_list == sample_id]
  ont_s <- CreateSeuratObject(
    counts = ont_s,
    assay = 'ONT',
    min.cells = 1,
    min.features = 0
  )  

  # keep only shared spots
  df_spot_meta <- vis_s@meta.data %>% 
    rownames_to_column(var = 'barcode') %>% 
    inner_join((ont_s@meta.data %>% rownames_to_column(var = 'barcode')), by = c('barcode')) %>%
    column_to_rownames(var = 'barcode')
  
  cat(sprintf("Number of spots: Visium/ONT/shared %d/%d/%d \n", 
      ncol(vis_s), ncol(ont_s), nrow(df_spot_meta)))
  
  # converting factor columns to characters
  j <- sapply(df_spot_meta, is.factor)
  df_spot_meta[j] <- lapply(df_spot_meta[j], as.character)

  # subset Visium data and only keep the counts slot
  vis_s <- vis_s[, rownames(df_spot_meta)]
  vis_s_obj <- CreateSeuratObject(
    counts = GetAssayData(vis_s, slot = "counts"), 
    meta.data = df_spot_meta,
    assay = 'Visium',
  )

  # subset ONT data
  ont_s_obj = ont_s[, rownames(df_spot_meta)]
  ont_s_obj@meta.data <- df_spot_meta

  # prepare isoform meta
  ont_iso_meta <- merge_gtf_iso_meta(sample_id, rownames(ont_s_obj))
  cat(sprintf("Number of transcripts: all/annotated %d/%d \n", 
              nrow(ont_s_obj), sum(ont_iso_meta$group != 'not matched')))
  print(table(ont_iso_meta$group))

  # save seurat objects and convert to anndata
  dir.create(sprintf("./all_spots/anndata/%s", sample_id), recursive = TRUE)
  
  # save isoform meta info
  write_csv(
    ont_iso_meta,
    sprintf("./all_spots/anndata/%s/ont_features.csv", sample_id)
  )

  # ONT data
  SaveH5Seurat(
    ont_s_obj,
    sprintf("./all_spots/anndata/%s/ont.h5Seurat", sample_id),
    overwrite = TRUE
  )
  Convert(
    sprintf("./all_spots/anndata/%s/ont.h5Seurat", sample_id),
    dest = "h5ad", overwrite = TRUE
  )
  file.remove(sprintf("./all_spots/anndata/%s/ont.h5Seurat", sample_id))
  
  # Visium data
  SaveH5Seurat(
    vis_s_obj,
    sprintf("./all_spots/anndata/%s/visium.h5Seurat", sample_id),
    overwrite = TRUE
  )
  Convert(
    sprintf("./all_spots/anndata/%s/visium.h5Seurat", sample_id),
    dest = "h5ad", overwrite = TRUE
  )
  file.remove(sprintf("./all_spots/anndata/%s/visium.h5Seurat", sample_id))
  
  
  # save spatial information
  dir.create(sprintf("./all_spots/visium_spatial_meta/%s", sample_id), recursive = TRUE)
  # sp_meta <- ont_s@images[[sample_id]] # missing spots in ONT images
  sp_meta <- vis_s@images[[1]]
  
  # write png image
  png::writePNG(
    sp_meta@image,
    sprintf("./all_spots/visium_spatial_meta/%s/tissue_lowres_image.png", sample_id)
  )
  # write scale factor
  sf <- sp_meta@scale.factors %>%
    `class<-`("list") %>%
    `names<-`(c(
      "spot_diameter_fullres", "fiducial_diameter_fullres",
      "tissue_hires_scalef", "tissue_lowres_scalef"
    ))
  # 'spot_diameter_fullres' is incorrect, overwrite with 'fiducial_diameter_fullres'
  sf$spot_diameter_fullres <- sf$fiducial_diameter_fullres * 0.619
  jsonlite::write_json(
    sf,
    sprintf("./all_spots/visium_spatial_meta/%s/scalefactors_json.json", sample_id),
    auto_unbox = T
  )
  # write tissue positions
  # pos <- ont@images[[sample_id]]@coordinates %>% # missing spots in ONT images
  pos <- sp_meta@coordinates %>%
    `colnames<-`(
      c("in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
    ) %>%
    rownames_to_column(var = 'barcode')
  
  write_csv(
    pos,
    sprintf("./all_spots/visium_spatial_meta/%s/tissue_positions.csv", sample_id)
  )
}
