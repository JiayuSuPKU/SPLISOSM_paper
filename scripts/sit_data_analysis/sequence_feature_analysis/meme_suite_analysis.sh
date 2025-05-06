#! /bin/bash

# This script runs MEME on the sequences extracted from the reference GTF file
# to find the motifs in the sequences

MEME_DIR="/Users/jysumac/Projects/Packages/meme"
CISBP_MEME_DIR="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_by_rbp.cleaned" # converted by pwm2meme.sh
CISBP_ALL_MOTIFS="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme" # converted by pwm2meme.sh
REF_GENOME_FASTA="/Users/jysumac/reference/mm10/mm10.genome.fa"

res_dir="/Users/jysumac/Projects/SPLISOSM_paper/results/sit_nar_23/transcripts"

# conda activate bioinfo

# Iterate over the samples
for dataset in mob cbs
do
  echo "=== Processing $dataset ..."

  # Create the result directory if it does not exist
  res_sample_dir="$res_dir/$dataset"
  if [ ! -d "$res_sample_dir" ]; then
      mkdir -p "$res_sample_dir"
  fi

  cd $res_sample_dir || exit

  ### (1) FIMO for scanning RBP with known motifs in the SVS sequences
  mkdir -p $res_sample_dir/fimo.out

  ### Prepare the sequences for MEME analysis
  group="$dataset"_svs # Use the SVS sequences for motif analysis
  # Convert the transcript GTF file to BED format
  gxf2bed -i "$res_sample_dir/$group.txid.gtf" -o "$res_sample_dir/fimo.out/$group.txid.bed"

  # Extract the RNA sequences (exon + intron) from the reference genome
  bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
    -bed "$res_sample_dir/fimo.out/$group.txid.bed" \
    -fo "$res_sample_dir/fimo.out/$group.txid.fa"
  sed -i '' 's/T/U/g; s/t/u/g' "$res_sample_dir/fimo.out/$group.txid.fa" # Mac OS sed

  # Prepare the background model for FIMO
  $MEME_DIR/libexec/meme-5.5.7/fasta-get-markov -rna -norc \
    "$res_sample_dir/fimo.out/$group.txid.fa" > "$res_sample_dir/fimo.out/$group.txid.meme.bg"

  # Run FIMO to find the motifs
  for rbp in Rbfox3 Celf4 Qk Cirbp Fmr1 Zfp36l1 Pcbp2 Khdrbs3 Elavl2
  do
    $MEME_DIR/bin/fimo --oc $res_sample_dir/fimo.out/"$group"_"$rbp" \
      --thresh 1e-3  --norc \
      --bfile "$res_sample_dir/fimo.out/$group.txid.meme.bg" \
      "$CISBP_MEME_DIR/$rbp.meme" "$res_sample_dir/fimo.out/$group.txid.fa"
  done

  ### (2) De novo motif discovery using XSTREME
  mkdir -p $res_sample_dir/xstreme.out/
  mkdir -p $res_sample_dir/xstreme.out/exon/ # to store results for exon sequence only
  mkdir -p $res_sample_dir/xstreme.out/slop/ # to store results for exon + flanking 300bp on each side
  mkdir -p $res_sample_dir/events/ # to store genomic coordinates of the events

  # Extract fasta sequence of exon skipping events
  for group in "$dataset"_svs "$dataset"_svens
  do
    # Aggregate exons in all SE and MX events
    true > "$res_sample_dir/events/$group.suppa.exon.gtf"
    for event_type in SE MX; do
      tail -n +2 "$res_sample_dir/suppa.out/$group/_${event_type}_strict.gtf" >> "$res_sample_dir/events/$group.suppa.exon.gtf"
    done

    # Remove duplicated exons and add exon_id to the gtf file
    sort -u -k1,1 -k4,4n -k5,5n "$res_sample_dir/events/$group.suppa.exon.gtf" | \
      awk '{print $0" exon_id \"exon_"NR"\";"}' >"$res_sample_dir/events/$group.suppa.exon.gtf.tmp"
    mv "$res_sample_dir/events/$group.suppa.exon.gtf.tmp" "$res_sample_dir/events/$group.suppa.exon.gtf"

    # Convert the transcript GTF file to BED format
    awk '{print $1"\t"$4"\t"$5"\t""exon_"NR"\t0\t"$7}' "$res_sample_dir/events/$group.suppa.exon.gtf" \
      > "$res_sample_dir/events/$group.suppa.exon.bed"

    # Extract the RNA sequences of exons from the reference genome
    bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
      -bed "$res_sample_dir/events/$group.suppa.exon.bed" \
      -fo "$res_sample_dir/events/$group.suppa.exon.fa"
    sed -i '' 's/T/U/g; s/t/u/g' "$res_sample_dir/events/$group.suppa.exon.fa" # Mac OS sed

    # Increase AS exons with flanking 300 bp on each side and extract the RNA sequences
    bedtools slop -i "$res_sample_dir/events/$group.suppa.exon.bed" \
      -g $REF_GENOME_FASTA.fai -b 300 > "$res_sample_dir/events/$group.suppa.exon.slop.bed"

    bedtools getfasta -name -s -fi $REF_GENOME_FASTA \
      -bed "$res_sample_dir/events/$group.suppa.exon.slop.bed" \
      -fo "$res_sample_dir/events/$group.suppa.exon.slop.fa"
    sed -i '' 's/T/U/g; s/t/u/g' "$res_sample_dir/events/$group.suppa.exon.slop.fa" # Mac OS sed
  done

  # Run XSTREME to find the enriched motifs in exon sequences
  $MEME_DIR/bin/xstreme -oc $res_sample_dir/xstreme.out/exon \
    -rna --evt 0.05 --fimo-skip \
    --meme-p 4 \
    -m $CISBP_ALL_MOTIFS \
    -n "$res_sample_dir"/events/"$dataset"_svens.suppa.exon.fa \
    -p "$res_sample_dir"/events/"$dataset"_svs.suppa.exon.fa

  # Run XSTREME to find the enriched motifs in exon + flanking intron sequences
  $MEME_DIR/bin/xstreme -oc $res_sample_dir/xstreme.out/slop \
    -rna --evt 0.05 --fimo-skip \
    --meme-p 4 \
    -m $CISBP_ALL_MOTIFS \
    -n "$res_sample_dir"/events/"$dataset"_svens.suppa.exon.slop.fa \
    -p "$res_sample_dir"/events/"$dataset"_svs.suppa.exon.slop.fa

done