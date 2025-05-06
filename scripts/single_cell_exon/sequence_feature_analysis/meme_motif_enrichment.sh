#!/bin/bash


SUPPA_DIR="/Users/jysumac/Projects/Packages/SUPPA/"
MEME_DIR="/Users/jysumac/Projects/Packages/meme"
CISBP_ALL_MOTIFS="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme" # converted by pwm2meme.sh
CISBP_MEME_DIR="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_by_rbp.cleaned" # converted by pwm2meme.sh
REF_GENOME_FASTA="/Users/jysumac/reference/mm10/mm10.genome.fa"

res_dir="/Users/jysumac/Projects/SPLISOSM_paper/results/sc_tasic_nature_18/events"

# conda activate bioinfo

### (1) SVS RBFOX targets from the ONT CBS2 dataset
if [ ! -d "$res_dir/ont" ]; then
    mkdir -p "$res_dir/ont"
fi
sample_id=ont_rbfox_targets

## First run SUPPA on the ONT RBFOX SVS targets to get exon skipping events
if [ ! -d "$res_dir/ont/suppa.out/$sample_id" ]; then
    mkdir -p "$res_dir/ont/suppa.out/$sample_id"
fi
python $SUPPA_DIR/suppa.py generateEvents -i $res_dir/ont/$sample_id.gtf \
  -o "$res_dir/ont/suppa.out/$sample_id/" \
  -f ioe -e SE SS MX RI FL

## Then extract sequences for the exon events
# Aggregate exons in all SE and MX events
true > "$res_dir/ont/$sample_id.exon.gtf"
for event_type in SE MX; do
  tail -n +2 "$res_dir/ont/suppa.out/$sample_id/_${event_type}_strict.gtf" >> "$res_dir/ont/$sample_id.exon.gtf"
done

# Remove duplicated exons and add exon_id to the gtf file
sort -u -k1,1 -k4,4n -k5,5n "$res_dir/ont/$sample_id.exon.gtf" | \
  awk '{print $0" exon_id \"exon_"NR"\";"}' > "$res_dir/ont/$sample_id.exon.gtf.tmp"
mv "$res_dir/ont/$sample_id.exon.gtf.tmp" "$res_dir/ont/$sample_id.exon.gtf"

# Convert the transcript GTF file to BED format by extract the corresponding columns
awk '{print $1"\t"$4"\t"$5"\t""exon_"NR"\t0\t"$7}' "$res_dir/ont/$sample_id.exon.gtf" \
  > "$res_dir/ont/$sample_id.exon.bed"

# Extract the RNA sequences of exons from the reference genome
bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
  -bed "$res_dir/ont/$sample_id.exon.bed" \
  -fo "$res_dir/ont/$sample_id.exon.fa"
sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/ont/$sample_id.exon.fa" # Mac OS sed

# Extract the RNA sequences of exons with flanking 150bp on each side
bedtools slop -i "$res_dir/ont/$sample_id.exon.bed" \
  -g $REF_GENOME_FASTA.fai -b 150 > "$res_dir/ont/$sample_id.exon.slop.bed"
bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
  -bed "$res_dir/ont/$sample_id.exon.slop.bed" \
  -fo "$res_dir/ont/$sample_id.exon.slop.fa"
sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/ont/$sample_id.exon.slop.fa" # Mac OS sed


## Run XSTREME for motif discovery
mkdir -p $res_dir/ont/xstreme.out

# XSTREME for exon only
$MEME_DIR/bin/xstreme -oc $res_dir/ont/xstreme.out/exon \
  -rna --evt 0.05 --fimo-skip \
  --meme-p 4 \
  -m $CISBP_ALL_MOTIFS \
  -p "$res_dir/ont/ont_rbfox_targets.exon.fa"

# XSTREME for exon + flanking 150bp on each side
$MEME_DIR/bin/xstreme -oc $res_dir/ont/xstreme.out/slop \
  -rna --evt 0.05 --fimo-skip \
  --meme-p 4 \
  -m $CISBP_ALL_MOTIFS \
  -p "$res_dir/ont/ont_rbfox_targets.exon.slop.fa"

## Run SEA motif enrichment analysis
mkdir -p $res_dir/ont/sea.out

# SEA for exon only
$MEME_DIR/bin/sea --oc $res_dir/ont/sea.out/exon \
  -p "$res_dir/ont/ont_rbfox_targets.exon.fa" \
  -m $CISBP_MEME_DIR/Rbfox3.meme -m $CISBP_MEME_DIR/Celf4.meme \
  -m $CISBP_MEME_DIR/Qk.meme -m $CISBP_MEME_DIR/Khdrbs3.meme \
  -m $CISBP_MEME_DIR/Matr3.meme -m $CISBP_MEME_DIR/Ralyl.meme \
  -m $CISBP_MEME_DIR/Hnrnpk.meme -m $CISBP_MEME_DIR/Elavl2.meme

# SEA for exon + flanking 150bp on each side
$MEME_DIR/bin/sea --oc $res_dir/ont/sea.out/slop \
  -p "$res_dir/ont/ont_rbfox_targets.exon.slop.fa" \
  -m $CISBP_MEME_DIR/Rbfox3.meme -m $CISBP_MEME_DIR/Celf4.meme \
  -m $CISBP_MEME_DIR/Qk.meme -m $CISBP_MEME_DIR/Khdrbs3.meme \
  -m $CISBP_MEME_DIR/Matr3.meme -m $CISBP_MEME_DIR/Ralyl.meme \
  -m $CISBP_MEME_DIR/Hnrnpk.meme -m $CISBP_MEME_DIR/Elavl2.meme


### (2) RBFOX-associated SVS TREND events
if [ ! -d "$res_dir/trend" ]; then
    mkdir -p "$res_dir/trend"
fi
sample_id=trend_rbfox_targets

## Extract sequences for the TREND events
# Remove potential duplicated peaks
sort -u -k1,1 -k6,6 -k2,2n -k3,3n "$res_dir/trend/$sample_id.bed" > "$res_dir/trend/$sample_id.dedup.bed"
# Filter out super long regions (> 100000bp) for MEME analysis
awk 'BEGIN {FS="\t"; OFS="\t"} {if ($3 - $2 < 100000) print $0}' "$res_dir/trend/$sample_id.dedup.bed" \
  > "$res_dir/trend/$sample_id.dedup.bed.tmp"
mv "$res_dir/trend/$sample_id.dedup.bed.tmp" "$res_dir/trend/$sample_id.dedup.bed"

# Extract exon sequences only
bedtools getfasta -name -s -split -fi $REF_GENOME_FASTA \
  -bed "$res_dir/trend/$sample_id.dedup.bed" \
  -fo "$res_dir/trend/$sample_id.exon.fa"
sed -i '' 's/T/U/g; s/t/u/g' "$res_dir/trend/$sample_id.exon.fa" # Mac OS sed

# Extract exon + 300bp upstream sequences
bedtools slop -i "$res_dir/trend/$sample_id.dedup.bed" \
  -g $REF_GENOME_FASTA.fai -s -l 300 -r 0 > "$res_dir/trend/$sample_id.slop.bed"
bedtools getfasta -name -s -split -fi $REF_GENOME_FASTA \
  -bed "$res_dir/trend/$sample_id.slop.bed" \
  -fo "$res_dir/trend/$sample_id.slop.fa"

## Run XSTREME for motif discovery
mkdir -p $res_dir/trend/xstreme.out

# XSTREME for exon only
$MEME_DIR/bin/xstreme -oc $res_dir/trend/xstreme.out/exon \
  -rna --evt 0.05 --fimo-skip \
  --meme-p 4 \
  -m $CISBP_ALL_MOTIFS \
  -p "$res_dir/trend/trend_rbfox_targets.exon.fa"

# XSTREME for exon + flanking 300bp upstream
$MEME_DIR/bin/xstreme -oc $res_dir/trend/xstreme.out/slop \
  -rna --evt 0.05 --fimo-skip \
  --meme-p 4 \
  -m $CISBP_ALL_MOTIFS \
  -p "$res_dir/trend/trend_rbfox_targets.slop.fa"

## Run SEA motif enrichment analysis
mkdir -p $res_dir/trend/sea.out

# SEA for exon only
$MEME_DIR/bin/sea --oc $res_dir/trend/sea.out/exon \
  -p "$res_dir/trend/trend_rbfox_targets.exon.fa" \
  -m $CISBP_MEME_DIR/Rbfox3.meme -m $CISBP_MEME_DIR/Celf4.meme \
  -m $CISBP_MEME_DIR/Qk.meme -m $CISBP_MEME_DIR/Khdrbs3.meme \
  -m $CISBP_MEME_DIR/Matr3.meme -m $CISBP_MEME_DIR/G3bp2.meme \
  -m $CISBP_MEME_DIR/Hnrnph1.meme -m $CISBP_MEME_DIR/Elavl2.meme

# SEA for exon + flanking 300bp upstream
$MEME_DIR/bin/sea --oc $res_dir/trend/sea.out/slop \
  -p "$res_dir/trend/trend_rbfox_targets.slop.fa" \
  -m $CISBP_MEME_DIR/Rbfox3.meme -m $CISBP_MEME_DIR/Celf4.meme \
  -m $CISBP_MEME_DIR/Qk.meme -m $CISBP_MEME_DIR/Khdrbs3.meme \
  -m $CISBP_MEME_DIR/Matr3.meme -m $CISBP_MEME_DIR/G3bp2.meme \
  -m $CISBP_MEME_DIR/Hnrnph1.meme -m $CISBP_MEME_DIR/Elavl2.meme