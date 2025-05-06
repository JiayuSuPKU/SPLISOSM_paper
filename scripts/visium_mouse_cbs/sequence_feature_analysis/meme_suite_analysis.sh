#! /bin/bash

# This script runs MEME on the sequences extracted from the reference GTF file
# to find the motifs in the sequences

MEME_DIR="/Users/jysumac/Projects/Packages/meme"
CISBP_MEME_DIR="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_by_rbp.cleaned" # converted by pwm2meme.sh
CISBP_ALL_MOTIFS="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme" # converted by pwm2meme.sh
REF_GENOME_FASTA="/Users/jysumac/reference/mm10/mm10.genome.fa"
REF_ANNO_GTF="/Users/jysumac/reference/mm10/Mus_musculus.GRCm38.102.gtf"

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RSCRIPT_PATH="$SCRIPT_DIR/extract_region_by_symbol.R"

res_dir="/Users/jysumac/Projects/SPLISOSM_paper/results/visium_mouse_cbs/events"
dataset="cbs"

# conda activate bioinfo

cd $res_dir || exit

### FIMO for scanning RBP with known motifs in the SVS sequences
mkdir -p $res_dir/fimo.out

### Prepare the sequences for MEME analysis
group="$dataset"_svs # Use the SVS sequences for motif analysis

echo "=== Extracting gene regions for $group ..."

# For every SVS gene, extract the start and end coordinates from the GTF file
/usr/local/bin/Rscript "$RSCRIPT_PATH" \
  $res_dir/$group.gene_symbol.txt $REF_ANNO_GTF $res_dir/fimo.out/$group.gene_region.bed

# Then extract the sequences from the reference genome
bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
  -bed "$res_dir/fimo.out/$group.gene_region.bed" \
  -fo "$res_dir/fimo.out/$group.gene_region.fa"
sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/fimo.out/$group.gene_region.fa" # Mac OS sed

echo "=== Running FIMO for $group ..."

# Prepare the background model for FIMO
$MEME_DIR/libexec/meme-5.5.7/fasta-get-markov -rna -norc \
  "$res_dir/fimo.out/$group.gene_region.fa" > "$res_dir/fimo.out/$group.gene_region.meme.bg"

# Run FIMO to find the motifs
for rbp in Rbfox3 Celf4 Qk Pcbp2 Rbm24 Rbms3 Rbm47 Khdrbs2 Pum1 Enox1 Snrnp70 
do
  $MEME_DIR/bin/fimo --oc $res_dir/fimo.out/"$group"_"$rbp" \
    --thresh 1e-3  --norc \
    --bfile "$res_dir/fimo.out/$group.gene_region.meme.bg" \
    "$CISBP_MEME_DIR/$rbp.meme" "$res_dir/fimo.out/$group.gene_region.fa"
done

echo "=== Done with FIMO for $group ..."

### De novo motif discovery using XSTREME
mkdir -p $res_dir/xstreme.out/

echo "=== Running XSTREME for $dataset ..."

# Extract fasta sequence of exon skipping events
for group in "$dataset"_svs "$dataset"_svens
do
  # # Filter out super long regions (> 100000bp) for MEME analysis
  # awk 'BEGIN {FS="\t"; OFS="\t"} {if ($3 - $2 < 100000) print $0}' "$res_dir/$group.peak.dedup.bed" \
  #   > "$res_dir/$group.peak.dedup.trunc.bed"

  # # Extract the RNA sequences from the reference genome
  # bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
  #   -bed "$res_dir/$group.peak.dedup.trunc.bed" \
  #   -fo "$res_dir/$group.peak.dedup.fa"
  # sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/$group.peak.dedup.fa" # Mac OS sed


  # Extract the RNA sequences from the reference genome
  bedtools getfasta -name -s -split -fi $REF_GENOME_FASTA \
    -bed "$res_dir/$group.exon.bed" \
    -fo "$res_dir/$group.exon.fa"
  sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/$group.exon.fa" # Mac OS sed
done

# # Run XSTREME to find the enriched motifs
# $MEME_DIR/bin/xstreme -oc $res_dir/xstreme.out \
#   -rna --evt 0.05 --fimo-skip \
#   -meme-p 4 \
#   -m $CISBP_ALL_MOTIFS \
#   -n "$res_dir"/"$dataset"_svens.peak.dedup.fa \
#   -p "$res_dir"/"$dataset"_svs.peak.dedup.fa

# Run XSTREME to find the enriched motifs
$MEME_DIR/bin/xstreme -oc $res_dir/xstreme.out \
  -rna --evt 0.05 --fimo-skip \
  -meme-p 4 \
  -m $CISBP_ALL_MOTIFS \
  -n "$res_dir"/"$dataset"_svens.exon.fa \
  -p "$res_dir"/"$dataset"_svs.exon.fa

echo "=== Done with XSTREME for $dataset ..."
