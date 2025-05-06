#! /bin/bash

# This script runs MEME on the sequences extracted from the reference genome
# to find the motifs in the sequences

MEME_DIR="/Users/jysumac/Projects/Packages/meme"
CISBP_ALL_MOTIFS="/Users/jysumac/reference/cisbp-rna/human_pm/meme_all_motifs.cleaned.meme" # converted by pwm2meme.sh
REF_GENOME_FASTA="/Users/jysumac/reference/hg38/hg38.genome.fa"

res_dir="/Users/jysumac/Projects/SPLISOSM_paper/results/gbm_visium/events/recurrent"

# conda activate bioinfo

cd $res_dir || exit


### De novo motif discovery using XSTREME
mkdir -p $res_dir/xstreme.out/

for dataset in rec_all rec_gbm
do
  echo "=== Running XSTREME for $dataset ..."

  mkdir -p $res_dir/xstreme.out/$dataset

  # Extract fasta sequence of exon skipping events
  for group in svs_"$dataset" svens_"$dataset"
  do
    # # Filter out super long regions (> 10,000bp) for MEME analysis
    # # split $11 by ',' and sum to get the total length of the exons
    # awk 'BEGIN {FS="\t"; OFS="\t"} {
    #     split($11, exon_lengths, ",");                # Split exon lengths into an array
    #     total_length = 0;                             # Initialize total exon length
    #     for (i in exon_lengths) total_length += exon_lengths[i]; # Sum up exon lengths
    #     if (total_length <= 10000) print $0;          # Print only rows with total_length <= 10000
    # }' "$res_dir/$group.exon.bed" > "$res_dir/$group.exon.trunc.bed"
    # echo "Extracting sequences for $(wc -l < "$res_dir/$group.exon.trunc.bed") out of $(wc -l < "$res_dir/$group.exon.bed") exons ..."

    echo "Extracting sequences for $(wc -l < "$res_dir/$group.exon.bed") exons ..."

    # Extract the RNA sequences from the reference genome
    bedtools getfasta -name -s -split -fi $REF_GENOME_FASTA \
      -bed "$res_dir/$group.exon.bed" \
      -fo "$res_dir/$group.exon.fa"
    sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/$group.exon.fa" # Mac OS sed
  done

  # Run XSTREME to find the enriched motifs
  $MEME_DIR/bin/xstreme -oc $res_dir/xstreme.out/$dataset \
    -rna --evt 0.05 --fimo-skip \
    -meme-p 4 \
    -m $CISBP_ALL_MOTIFS \
    -n "$res_dir"/svens_"$dataset".exon.fa \
    -p "$res_dir"/svs_"$dataset".exon.fa

  echo "=== Done with XSTREME for $dataset ..."

done