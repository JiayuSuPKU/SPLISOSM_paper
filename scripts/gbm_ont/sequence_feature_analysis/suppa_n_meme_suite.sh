#! /bin/bash

# This script runs MEME on the sequences extracted from the reference GTF file
# to find the motifs in the sequences

# MEME_DIR="/Users/jysumac/Projects/Packages/meme"
# CISBP_MEME_DIR="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_by_rbp.cleaned" # converted by pwm2meme.sh
# CISBP_ALL_MOTIFS="/Users/jysumac/reference/cisbp-rna/mouse_pm/meme_all_motifs.cleaned.meme" # converted by pwm2meme.sh
REF_GENOME_FASTA="/Users/jysumac/reference/hg38/hg38.genome.fa"
# REF_GENCODE_GTF="/Users/jysumac/reference/hg38/MANE.GRCh38.v0.6.select_ensembl_genomic.gtf"
REF_GENCODE_GTF="/Users/jysumac/reference/hg38/gencode.v32.annotation.gtf"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
# GTF2BED="$SCRIPT_DIR/gtf2bed.py"
NMD_CALLER="$SCRIPT_DIR/nmd_calling_orfanage.R"

SUPPA_DIR="/Users/jysumac/Projects/Packages/SUPPA/"


res_dir="/Users/jysumac/Projects/SPLISOSM_paper/results/gbm_ont/transcripts"

# conda activate bioinfo

# Iterate over the samples
for sample_id in DMG_1 DMG_2 DMG_3 DMG_4 DMG_5 GBM_1 GBM_2 GBM_3 GBM_4 GBM_5 GBM_6
do
  echo "=== Processing $sample_id ..."

  # Create the result directory if it does not exist
  res_sample_dir="$res_dir/$sample_id"
  if [ ! -d "$res_sample_dir" ]; then
      mkdir -p "$res_sample_dir"
  fi

  # cd $res_sample_dir || exit

  for group in svs_iso svens_iso
  do
    ## (1) Run SUPPA to find the splicing events
    if [ ! -d "$res_sample_dir/suppa.out/$group" ]; then
        mkdir -p "$res_sample_dir/suppa.out/$group"
    fi

    python $SUPPA_DIR/suppa.py generateEvents -i "$res_sample_dir/$sample_id.$group.gtf" \
      -o "$res_sample_dir/suppa.out/$group/" -f ioe -e SE SS MX RI FL

    # Empty the overall stats file
    true > "$res_sample_dir/suppa.out/$group/n_events.txt"

    # Count number of events in each category and store in a file
    for event in A3 A5 AF AL MX RI SE
    do
        file=$res_sample_dir/suppa.out/$group/_"$event"_strict.ioe
        n_events=$(( $(wc -l < "$file") - 1 ))
        echo $event $n_events >> "$res_sample_dir/suppa.out/$group/n_events.txt"
    done

    ## (2) Extract RNA sequences from the reference genome
    # Convert the transcript GTF file to BED format, one transcript per line
    # python $GTF2BED --gtf "$res_sample_dir/$sample_id.$group.gtf" \
    #   --bed "$res_sample_dir/$sample_id.$group.bed" \
    #   -t "transcript_name"
    gxf2bed -i "$res_sample_dir/$sample_id.$group.gtf" -o "$res_sample_dir/$sample_id.$group.bed" \
      --parent "transcript" --child "exon" --feature "transcript_name"

    # Extract the RNA sequences from the reference genome
    bedtools getfasta -name -s -split -fi $REF_GENOME_FASTA \
      -bed "$res_sample_dir/$sample_id.$group.bed" \
      -fo "$res_sample_dir/$sample_id.$group.rna.fa"
    sed -i '' 's/T/U/g; s/t/u/g' "$res_sample_dir/$sample_id.$group.rna.fa" # Mac OS sed

  done

  ## (3) Call ORF for novel transcripts using orfanage
  orfanage --reference $REF_GENOME_FASTA \
    --keep_all_cds \
    --output "$res_sample_dir/$sample_id.svs_iso.oa.gtf" \
    --query "$res_sample_dir/$sample_id.svs_iso.gtf" \
    --minlen 300 \
    --stats "$res_sample_dir/svs_iso.oa.stats" \
    --cleant \
    $REF_GENCODE_GTF
  
  ## (4) Call nonsense-mediated decay (NMD)
  /usr/local/bin/Rscript "$NMD_CALLER" -i $res_sample_dir/$sample_id.svs_iso.oa.gtf -o $res_sample_dir/$sample_id.svs_iso.oa.tmp.gtf
  mv $res_sample_dir/$sample_id.svs_iso.oa.tmp.gtf $res_sample_dir/$sample_id.svs_iso.oa.gtf

done