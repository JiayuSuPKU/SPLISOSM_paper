#!/bin/bash

# This script calculates the conservation score for the sequences extracted from the reference tracks

bigWigAverageOverBed="/Users/jysumac/Projects/Packages/bigWigAverageOverBed"
phyloP100way="/Users/jysumac/reference/evo_conservation/hg38.phyloP100way.bw"
phastCons100way="/Users/jysumac/reference/evo_conservation/hg38.phastCons20way.bw"
res_dir="/Users/jysumac/Projects/SPLISOSM_paper/results/human_dlpfc/events"


mkdir -p $res_dir/conservation

# loop over five groups [svs_all, svs_cbs, svs_cbs_up, svs_cbs_down, svs_cbs_exon]
for group in svs_all svens_all svs_shared svens_shared neg_all
do
  echo "=== Calculating conservation score for $group ..."
  # PhyloP score
  $bigWigAverageOverBed $phyloP100way "$res_dir/$group.exon.bed" \
    "$res_dir/conservation/$group.phylop.tab"

  # PhastCons score
  $bigWigAverageOverBed $phastCons100way "$res_dir/$group.exon.bed" \
    "$res_dir/conservation/$group.phastcons.tab"
done
