#!/bin/bash

# Specify raw data paths
data_dir="/gpfs/commons/home/jsu/data/visium_ted/mouse_cbs"
bam_file=$data_dir"/Visium_Fresh_Frozen_Adult_Mouse_Brain_possorted_genome_bam.bam"
junc_file=$data_dir"/possorted_genome_bam.junc.bed"

# Extract junction reads from bam
module load regtools
regtools junctions extract $bam_file -o $junc_file -s 0

# Prepare Sierra configs and inputs
gtf_file="/gpfs/commons/home/jsu/reference/annotations/gencode.vM10.annotation.gtf"
res_dir="/gpfs/commons/home/jsu/results/visium_ted/mouse_cbs"
peak_file=$res_dir"/peak_no_cutoff.txt"
peak_annot_file=$res_dir"/peak_no_cutoff.annot.txt"
peak_annot_overlap_file=$res_dir"/peak_no_cutoff.annot.overlap.txt"
whitelist_file=$res_dir"/barcodes.txt"
output_dir=$res_dir"/counts_no_cutoff"

# convert 10x Visium spatial/tissue_positions.csv to the required whitelist file
awk 'NR > 1' $data_dir/spatial/tissue_positions.csv | sed 's/,/\t/g' > $whitelist_file

mamba activate r-4.4

# Run FindPeaks in R
Rscript -e "
library(Sierra)  # Load Sierra package
FindPeaks(
  output.file = '${peak_file}', 
  gtf.file = '${gtf_file}', 
  bamfile = '${bam_file}', 
  junctions.file = '${junc_file}', 
  ncores = 8, 
  min.jcutoff.prop = 0.0, 
  min.cov.prop = 0.0, 
  min.peak.prop = 0.0
)
"

# Run CountPeaks in R
Rscript -e "
library(Sierra)  # Load Sierra package
CountPeaks(
  peak.sites.file = '${peak_file}', 
  gtf.file = '${gtf_file}', 
  bamfile = '${bam_file}', 
  whitelist.file = '${whitelist_file}', 
  output.dir = '${output_dir}', 
  ncores = 8
)
"

# Run AnnotatePeaksFromGTF in R
Rscript -e "
library(Sierra)  # Load Sierra package
AnnotatePeaksFromGTF(
  peak.sites.file = '${peak_file}', 
  gtf.file = '${gtf_file}', 
  output.file = '${peak_annot_file}',
  transcriptDetails = TRUE,
)
"

# Run AnnotatePeaksFromGTF in R, without 3'UTR overwriting other annotations
Rscript -e "
library(Sierra)  # Load Sierra package
AnnotatePeaksFromGTF(
  peak.sites.file = '${peak_file}', 
  gtf.file = '${gtf_file}', 
  output.file = '${peak_annot_overlap_file}',
  transcriptDetails = TRUE,
  annotation_correction = FALSE
)
"

echo "Analysis complete. Outputs saved in ${output_dir}."
