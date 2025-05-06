#!/bin/bash

# Specify raw data paths
data_dir="/gpfs/commons/home/jsu/data/human_cortex_visium"
gtf_file="/gpfs/commons/home/jsu/reference/cellranger/GRCh38-2024-A.gtf"
# gtf_file="/gpfs/commons/home/jsu/reference/annotations/gencode.v32.annotation.gtf"
res_dir="/gpfs/commons/home/jsu/results/visium_ted/human_dlpfc"

# Create the results directory
mkdir -p $res_dir

# Extract all subdirectories in the data directory
samples=$(find $data_dir -mindepth 1 -maxdepth 1 -type d -print0 | sort -z | xargs -0 -n 1 basename)

module load regtools
mamba activate r-4.4

# Loop through each sample, extract junction reads and run Sierra FindPeaks
for sample in $samples; do
    echo "Processing sample $sample ..."
    # Specify the input and output files
    bam_file=${data_dir}/${sample}/${sample}_mRNA.bam
    new_bam_file=${data_dir}/${sample}/${sample}.renamed.bam
    junc_file=${data_dir}/${sample}/${sample}.junc.bed
    peak_file=${res_dir}/${sample}/peak_default.txt

    # Rename bam chromosome names to match gtf (from 1, 2, 3, ... to chr1, chr2, chr3, ...)
    echo "Renaming chromosome names for sample $sample ..."
    samtools view -H "$bam_file" | sed -e 's/SN:\([0-9XY]\)/SN:chr\1/' -e 's/SN:MT/SN:chrM/' | samtools reheader - "$bam_file" > "$new_bam_file"
    samtools index -@ 32 "$new_bam_file"
    
    # Extract junction reads from bam
    echo "Extracting junction reads for sample $sample ..."
    regtools junctions extract "$new_bam_file" -o "$junc_file" -s 0

    mkdir -p "$res_dir/$sample"

    # Run Sierra FindPeaks
    echo "Running Sierra FindPeaks for sample $sample ..."
    Rscript -e "
    suppressMessages(library(Sierra))  # Load Sierra package
    FindPeaks(
      output.file = '${peak_file}', 
      gtf.file = '${gtf_file}', 
      bamfile = '${new_bam_file}', 
      junctions.file = '${junc_file}', 
      ncores = 8
      # min.jcutoff.prop = 0.0, 
      # min.cov.prop = 0.0, 
      # min.peak.prop = 0.0
    )
    "
done

# Run Sierra mergePeaks to combine peaks from all samples
merged_peak_file=$res_dir"/peak_default.merged.txt"
Rscript -e "
# split the samples into a vector
sample_list <- '$samples'
sample_list <- strsplit(sample_list, '\n')[[1]]

# create a data frame with the peak files and sample identifiers
peak.output.file <- c(
  paste0('${res_dir}/', 
  sample_list,
  '/peak_default.txt')
)
peak.dataset.table = data.frame(
  Peak_file = peak.output.file,
  Identifier = sample_list, 
  stringsAsFactors = FALSE
)
head(peak.dataset.table)

suppressMessages(library(Sierra))  # Load Sierra package

# merge the peaks
peak.merge.output.file = '${merged_peak_file}'
MergePeakCoordinates(peak.dataset.table, output.file = peak.merge.output.file, ncores = 8)
"

# Run AnnotatePeaksFromGTF to annotate the called peaks
annot_peak_file=$res_dir"/peak_default.merged.annot.txt"
Rscript -e "
suppressMessages(library(Sierra))  # Load Sierra package
AnnotatePeaksFromGTF(
  peak.sites.file = '${merged_peak_file}', 
  gtf.file = '${gtf_file}', 
  output.file = '${annot_peak_file}',
  transcriptDetails = TRUE,
)
"

# Loop through each sample, run Sierra CountPeaks
for sample in $samples; do
    echo "Counting peaks for sample $sample ..."
    # Specify the input and output files
    new_bam_file=${data_dir}/${sample}/${sample}.renamed.bam
    output_dir=${res_dir}/${sample}/counts_default

    # convert 10x Visium spatial/tissue_positions.csv to the required whitelist file
    whitelist_file=${res_dir}/${sample}/barcodes.tsv
    sed 's/,/\t/g' ${data_dir}/${sample}/tissue_positions_list.txt > "$whitelist_file"

    # Run CountPeaks in R
    Rscript -e "
    suppressMessages(library(Sierra))  # Load Sierra package
    CountPeaks(
      peak.sites.file = '${merged_peak_file}', 
      gtf.file = '${gtf_file}', 
      bamfile = '${new_bam_file}', 
      whitelist.file = '${whitelist_file}', 
      output.dir = '${output_dir}', 
      ncores = 8
    )
    "
done
