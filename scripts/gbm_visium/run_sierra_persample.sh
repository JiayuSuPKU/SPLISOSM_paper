#!/bin/bash
#SBATCH --job-name=Sierra-GBM                  # Job name
#SBATCH --partition=cpu                      # Partition Name
#SBATCH --mail-type=END,FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jsu@nygenome.org      # Where to send mail  
#SBATCH --mem=64G                            # Job memory request. Different units can be specified using the suffix [K|M|G|T]
#SBATCH --cpus-per-task=8                     # Number of CPU cores per task
#SBATCH --time=24:00:00                       # Time limit 24 hours
#SBATCH --output=stdout_%j.log               # Standard output and error log

# Specify raw data paths
data_dir="/gpfs/commons/home/jsu/data/visium_ted/gbm_visium_cell_24/count"
gtf_file="/gpfs/commons/home/jsu/reference/cellranger/GRCh38-2024-A.gtf"
# gtf_file="/gpfs/commons/home/jsu/reference/annotations/gencode.v32.annotation.gtf"
res_dir="/gpfs/commons/home/jsu/projects/SPLISOSM_paper/data/gbm_visium_cell_24/sierra_peaks_default_individual"

# Create the results directory
mkdir -p $res_dir

# Extract all subdirectories in the data directory
samples=$(find $data_dir -mindepth 1 -maxdepth 1 -type d -print0 | sort -z | xargs -0 -n 1 basename)

module load regtools

conda_dir="/gpfs/commons/home/jsu/softwares/miniforge3/"
source "$conda_dir/etc/profile.d/conda.sh"
conda activate r-4.4

# Loop through each sample, extract junction reads and run Sierra FindPeaks
for sample in $samples; do
    echo "Processing sample $sample ..."

    # Specify the input and output files
    bam_file=${data_dir}/${sample}/outs/possorted_genome_bam.bam
    junc_file=${data_dir}/${sample}/outs/possorted_genome_bam.junc.bed
    peak_file=${res_dir}/${sample}/peak_default.txt
    # peak_file=${res_dir}/${sample}/peak_no_cutoff.txt
    peak_annot_file=${res_dir}/${sample}/peak_default.annot.no_correction.txt
    output_dir=${res_dir}/${sample}/counts_default

    # Create the results directory if it does not exist
    mkdir -p "$res_dir/$sample"

    # Extract junction reads from bam
    echo "Extracting junction reads for sample $sample ..."
    regtools junctions extract "$bam_file" -o "$junc_file" -s 0

    # convert 10x Visium spatial/tissue_positions.csv to the required whitelist file
    whitelist_file=${res_dir}/${sample}/barcodes.tsv
    sed 's/,/\t/g' ${data_dir}/${sample}/outs/spatial/tissue_positions.csv > "$whitelist_file"


    # Run Sierra FindPeaks
    echo "Running Sierra FindPeaks for sample $sample ..."
    Rscript -e "
    suppressMessages(library(Sierra))  # Load Sierra package
    FindPeaks(
      output.file = '${peak_file}', 
      gtf.file = '${gtf_file}', 
      bamfile = '${bam_file}', 
      junctions.file = '${junc_file}', 
      ncores = 8
      # min.jcutoff.prop = 0.0, 
      # min.cov.prop = 0.0, 
      # min.peak.prop = 0.0
    )
    "

    # Run AnnotatePeaksFromGTF in R, without 3'UTR overwriting other annotations
    echo "Running Sierra AnnotatePeaksFromGTF for sample $sample ..."
    Rscript -e "
    library(Sierra)  # Load Sierra package
    AnnotatePeaksFromGTF(
      peak.sites.file = '${peak_file}', 
      gtf.file = '${gtf_file}', 
      output.file = '${peak_annot_file}',
      transcriptDetails = TRUE,
      annotation_correction = FALSE
    )
    "

    # Run CountPeaks in R
    echo "Running Sierra CountPeaks for sample $sample ..."
    Rscript -e "
    suppressMessages(library(Sierra))  # Load Sierra package
    CountPeaks(
      peak.sites.file = '${peak_file}', 
      gtf.file = '${gtf_file}', 
      bamfile = '${bam_file}', 
      whitelist.file = '${whitelist_file}', 
      output.dir = '${output_dir}', 
      ncores = 8
    )
    "
done
