# Convert SiT data (Seurat object) to h5ad format for later use
setwd("~/Projects/SPLISOSM_paper/")
library(tidyverse)
library(Seurat)
library(SeuratDisk)

# loop through the three samples
for (sample in c('mob', 'cbs1', 'cbs2')){
  rds_data <- readRDS(sprintf('data/sit_nar_23/GSE153859_%s.rds', toupper(sample)))

  # create data directories
  dir.create(sprintf('data/sit_nar_23/%s/anndata', sample), recursive = TRUE)
  dir.create(sprintf('data/sit_nar_23/%s/visium_spatial_meta', sample), recursive = TRUE)

  # separate and save data by assay
  for (assay in c('Spatial', 'ISO', 'ISOG')){
    seu_obj = CreateSeuratObject(
      counts = GetAssayData(rds_data, slot = 'counts', assay = assay),
      meta.data = data.frame(
        'orig.ident' = rds_data@meta.data$orig.ident,
        'region' = as.character(rds_data@meta.data$ClusterName),
        'AtoI_total' = rds_data@meta.data$AtoI_total,
        'AtoI_ratio' = rds_data@meta.data$AtoI_ratio,
        'AtoI_detected' = rds_data@meta.data$AtoI_detected,
        row.names = rownames(rds_data@meta.data))
    )
    SaveH5Seurat(
      seu_obj, 
      sprintf('data/sit_nar_23/%s/anndata/%s_%s.h5Seurat', sample, sample, tolower(assay)), 
      overwrite = TRUE
    )
    Convert(
      sprintf('data/sit_nar_23/%s/anndata/%s_%s.h5Seurat', sample, sample, tolower(assay)),
      dest = "h5ad", overwrite = T
    )
    file.remove(
      sprintf('data/sit_nar_23/%s/anndata/%s_%s.h5Seurat', sample, sample, tolower(assay))
    )
  }

  # separate and save AtoI editing data
  Acounts = GetAssayData(rds_data, slot = 'counts', assay = 'AtoI')
  rownames(Acounts) <- paste0(rownames(Acounts), '..A')
  Gcounts = GetAssayData(rds_data, slot = 'data', assay = 'AtoI')
  rownames(Gcounts) <- paste0(rownames(Gcounts), '..G')
  AtoIcounts = rbind(Acounts, Gcounts)

  seu_obj = CreateSeuratObject(
    counts = AtoIcounts,
    meta.data = data.frame(
      'orig.ident' = rds_data@meta.data$orig.ident,
      'region' = as.character(rds_data@meta.data$ClusterName),
      'AtoI_total' = rds_data@meta.data$AtoI_total,
      'AtoI_ratio' = rds_data@meta.data$AtoI_ratio,
      'AtoI_detected' = rds_data@meta.data$AtoI_detected,
      row.names = rownames(rds_data@meta.data))
  )
  SaveH5Seurat(
    seu_obj, 
    sprintf('data/sit_nar_23/%s/anndata/%s_atoi.h5Seurat', sample, sample), 
    overwrite = TRUE
  )
  Convert(
    sprintf('data/sit_nar_23/%s/anndata/%s_atoi.h5Seurat', sample, sample),
    dest = "h5ad", overwrite = TRUE
  )
  file.remove(
    sprintf('data/sit_nar_23/%s/anndata/%s_atoi.h5Seurat', sample, sample)
  )

  # save spatial information
  sp_meta = rds_data@images[[1]]
  # write png image
  png::writePNG(
    sp_meta@image, 
    sprintf('data/sit_nar_23/%s/visium_spatial_meta/tissue_lowres_image.png', sample)
  )
  # write scale factor
  sf = sp_meta@scale.factors %>% `class<-`('list') %>% 
    `names<-`(c('spot_diameter_fullres', 'fiducial_diameter_fullres', 
                'tissue_hires_scalef', 'tissue_lowres_scalef'))
  # 'spot_diameter_fullres' is incorrect, overwrite with 'fiducial_diameter_fullres'
  sf$spot_diameter_fullres = sf$fiducial_diameter_fullres * 0.619
  jsonlite::write_json(
    sf, 
    sprintf('data/sit_nar_23/%s/visium_spatial_meta/scalefactors_json.json', sample), 
    auto_unbox = T
  )
  # write tissue positions
  pos <- sp_meta@coordinates %>% `colnames<-`(
    c('in_tissue', 'array_row', 'array_col', 'pxl_row_in_fullres', 'pxl_col_in_fullres')
  ) %>% rownames_to_column('barcode')
  write_csv(pos, sprintf('data/sit_nar_23/%s/visium_spatial_meta/tissue_positions.csv', sample))

}

