#!/bin/bash

VAST_DIR="$HOME/softwares/vast-tools-2.5.1"
MATT_DIR="$HOME/softwares/matt"

data_dir='/gpfs/commons/home/jsu/data/rbp_ko/rbfox_tko_mouse_neuron'
ref_gtf='/gpfs/commons/home/jsu/reference/annotations/Mus_musculus.GRCm38.102.gtf'
ref_fa='/gpfs/commons/home/jsu/reference/mm10/mm10.fa'

conda activate bioinfo
module load bowtie

cd $data_dir || exit

## Step 1: Download RNA-seq data
echo "SRR6441294 rbfox_tko_d10_a
SRR6441279 rbfox_tko_d10_b
SRR6441280 rbfox_tko_d10_c
SRR6441289 rbfox_wt_d10_a
SRR6441290 rbfox_wt_d10_b
SRR6441293 rbfox_wt_d10_c
" > accession_numbers.txt

"$MATT_DIR"/matt retr_rnaseq accession_numbers.txt -fasta -keepsra -o rnaseq_data -p 16

## Step 2: Run VAST-TOOLS to get exon skipping events
# Step 2.1: aligning the RNA-seq data
for file in rbfox_tko_d10_a rbfox_tko_d10_b rbfox_tko_d10_c rbfox_wt_d10_a rbfox_wt_d10_b rbfox_wt_d10_c
do
  $VAST_DIR/vast-tools align rnaseq_data/${file}_1.fasta.gz rnaseq_data/${file}_2.fasta.gz \
    -sp Mm2 -c 32
done

# Step 2.2: Merging the three replicates into one sample to increase the read coverage.
echo -e "rbfox_tko_d10_a_1\trbfox_tko_d10
rbfox_tko_d10_b_1\trbfox_tko_d10
rbfox_tko_d10_c_1\trbfox_tko_d10
rbfox_wt_d10_a_1\trbfox_wt_d10
rbfox_wt_d10_b_1\trbfox_wt_d10
rbfox_wt_d10_c_1\trbfox_wt_d10
" > merge_samples.txt

"$VAST_DIR"/vast-tools merge -g merge_samples.txt -sp Mm2

# Step 2.3: Estimate inclusion levels and generate final output table, which will be vast_out/INCLUSION_LEVELS_FULL-Mmu8-mm9.tab.
# This table contains the estimated inclusion levels for each single sample and for the merged samples.
"$VAST_DIR"/vast-tools combine -sp Mm2 --cores 4

## Step 3: Run MATT to extract exon features
mkdir -p $data_dir/matt_rbp_map
cd $data_dir/matt_rbp_map || exit

# The final output table of vast-tools with inclusion levels genome-wide for AS events of type SE, IR, ALT3, ALT5 is:
vts_file=$data_dir/vast_out/INCLUSION_LEVELS_FULL-mm10-8.tab

# Step 3.0: Inspecting the column names of this table reveals it's contents and column structure. It contains, besides other
# pieces of information, for each AS event its estimated inclusion level for all original samples 
# (3 replicates per condition) and for the two merged samples.
"$MATT_DIR"/matt get_colnms $vts_file
# inclusion levels are in the range of [0;100])
"$MATT_DIR"/matt prnt_tab $vts_file -W 10 | less -S -

# Step 3.1: Extract from vast-tools output table skipped exons
# get_vast is a utility command for extracting splicing events from vast-tools output tables
# - complex : restricts retrieved events; vast-tools encodes skipped-exon events with IDs S, C1, C2, C3, MIC
# -a nova2kd : is the merged Nova2 KD sample
# -b nova2wt : is the merged Nova2 WT sample
# -minqab VLOW : vast-tools reports for each estimated PSI value a quality (N,VLOW,LOW,OK,SOK) and we restrict here the output to events to those for which the nova2kd and the nova2wt PSI could be estimated with at least VLOW quality.
# -minqglob N : similar to -minqab but applies to all other samples in the table beside the samples defined with arguments -a and -b.
# Setting it to N prevents the exclusion of events, which have VLOW in -a and -b samples, but worse quality in the other samples.
# -gtf and -f : gene IDs will be extracted automatically
"$MATT_DIR"/matt get_vast $vts_file -complex S,C1,C2,C3,MIC -a rbfox_tko_d10 -b rbfox_wt_d10 \
  -minqab VLOW -minqglob N -gtf $ref_gtf -f gene_id > exons.tab

# Check the number of events (all are skipped exons) obtained.
"$MATT_DIR"/matt col_uniq exons.tab COMPLEX

# Check column names of table with retrieved skipped exons
# Applying get_vast as split the information mixed in some columns of vast-tool's output file into distinct columns and has added columns like
# STRAND, flanking borders of flanking exons, GENEID, DPSI_GRPA_MINUS_GRPB with deltaPSIs nova2kd vs. nova2wt, and the number of inclusion- and exclusion reads.
"$MATT_DIR"/matt get_colnms exons.tab

# Step 3.2: Definition of groups of skipped exons to be compared in subsequent analyses.
# The command def_cats outputs a table with one column named GROUP containing the group IDs.
# We specify three groups:
# 'silenced=DPSI_GRPA_MINUS_GRPB[25,100]' : exons silenced by Nova2 have a dPSI in [25,100], i.e., upon Nova2 KD they are more included meaning, under normal conditions, Nova2 would silence them.
# 'enhanced=DPSI_GRPA_MINUS_GRPB[-100,-25]' : exons enhanced by Nova2 have a dPSI in [-100,-25]
# 'unregulated=DPSI_GRPA_MINUS_GRPB[-1,1] PSI_nova2wt[10,90]' : unregulated exons have a dPSI in [-1,1] and a avg. PSI in the WT sample in [10,90], i.e., they are truly alternatively skipped-exons not constitutive.
# All events not falling in any of these groups will get the group ID NA = Not Available.
# Note: Here, piping is applied and the command col_uniq on the output table is used to directly count the number of events per group.
# By this, the user might change the conditions defined and immediately check the number of events per group.

# Exons with strong RBFOX KO effect
"$MATT_DIR"/matt def_cats exons.tab STRONG 'silenced=DPSI_GRPA_MINUS_GRPB[25,100]' 'enhanced=DPSI_GRPA_MINUS_GRPB[-100,-25]' 'unregulated=DPSI_GRPA_MINUS_GRPB[-1,1] PSI_rbfox_wt_d10[5,95]' | \
  "$MATT_DIR"/matt add_cols exons.tab -

# Exons with weak RBFOX KO effect (dPSI in 10% - 25%)
"$MATT_DIR"/matt def_cats exons.tab WEAK 'silenced=DPSI_GRPA_MINUS_GRPB[10,25]' 'enhanced=DPSI_GRPA_MINUS_GRPB[-25,-10]' 'unregulated=DPSI_GRPA_MINUS_GRPB[-1,1] PSI_rbfox_wt_d10[5,95]' | \
  "$MATT_DIR"/matt add_cols exons.tab -

# Final check of number events per group in exons.tab
"$MATT_DIR"/matt col_uniq exons.tab STRONG
"$MATT_DIR"/matt col_uniq exons.tab WEAK

# Step 3.3: Generate motif RNA-maps for all CISBP-RNA consensus motifs comparing the three exon groups
# A FASTA file with chromosome sequences of mouse (Mmu09) must be given as sequences of the exons and flanking introns 
# need to be retrieved for generating a motif RNA-map. Sequence IDs therein must match sequence IDs in exons.tab.
"$MATT_DIR"/matt rna_maps_cisbp exons.tab UPSTRM_EX_BORDER START END DOSTRM_EX_BORDER SCAFFOLD STRAND STRONG[silenced,enhanced,unregulated] 49 50 150 $ref_fa cisbprna_regexps -d cisbprna_maps_strong
"$MATT_DIR"/matt rna_maps_cisbp exons.tab UPSTRM_EX_BORDER START END DOSTRM_EX_BORDER SCAFFOLD STRAND WEAK[silenced,enhanced,unregulated] 49 50 150 $ref_fa cisbprna_regexps -d cisbprna_maps_weak
# Infos: the output folder contains a
# 1.) summary with all motif RNA-maps (0_all_motif_rna_maps.pdf); each motif RNA-map is printed twice: with and without a data coverage plot.
# The maps are ordered according to differences of the enrichment scores between groups
# (silenced vs. unregulated and enhanced vs unregulated) from most positive differences to most negative differences. 
# Hence, maps with most differences come first and last.
# 2.) all motif RNA-maps as PDFs
# More infos: This analysis can be done with all input tables describing exons, wherever they come from, if they contain the necessary pieces of information:
# 1.) exon start coordinate
# 2.) exon end
# 3.) scaffold/chromosome of exon
# 4.) strand
# 5.) flanking borders of flanking exons
# 6.) group IDs

# Step 3.4: Extract 150 nt long sequences upstream of exons
# For Extracting sequences Matt offers the utility command get_seqs, which can be invoked in three different modi.
# Sequences from strand - will be reverse complemented automatically.
fa_renamed=mm10_renamed.fasta
seqkit replace -p "chrM" -r "chrMT" $ref_fa -o $fa_renamed
# # Piping allows to add the generated two columns with sequences (SEQ_UP and SEQ_DOWN) directly to the exon table.
# "$MATT_DIR"/matt get_seqs exons.tab START END 50 50 SCAFFOLD STRAND $fa_renamed | \
#   "$MATT_DIR"/matt rn_cols - SEQ_UP:SEQ_UP50 SEQ_DOWN:SEQ_DOWN50 | "$MATT_DIR"/matt add_cols exons.tab -
# "$MATT_DIR"/matt get_seqs exons.tab START END 100 100 SCAFFOLD STRAND $fa_renamed | \
#   "$MATT_DIR"/matt rn_cols - SEQ_UP:SEQ_UP100 SEQ_DOWN:SEQ_DOWN100 | "$MATT_DIR"/matt add_cols exons.tab -

# extract the exon sequences
"$MATT_DIR"/matt get_seqs exons.tab START END SCAFFOLD STRAND $fa_renamed | \
  "$MATT_DIR"/matt rn_cols - SEQ:SEQ_EXON | "$MATT_DIR"/matt add_cols exons.tab -

# 150 upstream and downstream of the exons
"$MATT_DIR"/matt get_seqs exons.tab START END 150 150 SCAFFOLD STRAND $fa_renamed | \
  "$MATT_DIR"/matt rn_cols - SEQ_UP:SEQ_UP150 SEQ_DOWN:SEQ_DOWN150 | "$MATT_DIR"/matt add_cols exons.tab -

# Checking extended table with sequences
"$MATT_DIR"/matt get_colnms exons.tab

# Step 3.5: Enrichment analysis in extracted 50 nt long upstream sequences with all CISBP-RNA motifs
# This command outputs results in form of a table to STDOUT which will be redirected into a file.
# output_file=cisbprna_enrichment_test.tab
# Comparisons done are specified by GROUP[silenced,unregulated,enhanced,unregulated], meaning silenced vs. unregulated and enhanced vs. unregulated

# # Strongly regulated exons
# # 50nt enriched motifs
# "$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_UP50 GROUP25[silenced,unregulated,enhanced,unregulated] 1000 \
#   -p fox,celf,qki,khdrbs > motif_enrichment.dpsi25.up50.tab
# "$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_DOWN50 GROUP25[silenced,unregulated,enhanced,unregulated] 1000 \
#   -p fox,celf,qki,khdrbs > motif_enrichment.dpsi25.down50.tab

# # Strongly regulated exons
# # 100 nt enriched motifs
# "$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_UP100 GROUP25[silenced,unregulated,enhanced,unregulated] 1000 \
#   -p fox,celf,qki,khdrbs > motif_enrichment.dpsi25.up100.tab
# "$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_DOWN100 GROUP25[silenced,unregulated,enhanced,unregulated] 1000 \
#   -p fox,celf,qki,khdrbs > motif_enrichment.dpsi25.down100.tab

# Strongly regulated exons
# exon and up/downstream 150 nt enriched motifs
"$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_EXON STRONG[silenced,unregulated,enhanced,unregulated] 2000 \
  -p fox,celf,qki,khdrbs,elavl > motif_enrichment.strong.exon.tab
"$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_UP150 STRONG[silenced,unregulated,enhanced,unregulated] 2000 \
  -p fox,celf,qki,khdrbs,elavl > motif_enrichment.strong.up150.tab
"$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_DOWN150 STRONG[silenced,unregulated,enhanced,unregulated] 2000 \
  -p fox,celf,qki,khdrbs,elavl > motif_enrichment.strong.down150.tab

# Weakly regulated exons
# exon and up/downstream 150 nt enriched motifs
"$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_EXON WEAK[silenced,unregulated,enhanced,unregulated] 2000 \
  -p fox,celf,qki,khdrb,elavl > motif_enrichment.weak.exon.tab
"$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_UP150 WEAK[silenced,unregulated,enhanced,unregulated] 2000 \
  -p fox,celf,qki,khdrbs,elavl > motif_enrichment.weak.up150.tab
"$MATT_DIR"/matt test_cisbp_enrich exons.tab SEQ_DOWN150 WEAK[silenced,unregulated,enhanced,unregulated] 2000 \
  -p fox,celf,qki,khdrbs,elavl > motif_enrichment.weak.down150.tab


