#! /bin/bash

# This script runs Matt on AS exons to extract and compare the sequence features
# Matt: https://gitlab.com/aghr/matt
# Matt exon features: https://gitlab.com/aghr/matt/-/wikis/sec_examples/Comparing_sets_of_exons_wrt_their_exon_features

MATT_DIR="$HOME/softwares/matt/"
REF_FA="$HOME/reference/mm10/mm10.fa"
REF_GTF="$HOME/reference/annotations/Mus_musculus.GRCm38.102.gtf"

# Function to process a single GTF and event table pair
process_suppa_gtf() {
    local gtf_file=$1
    local ioe_file=$2
    local group_label=$3
    
    awk -F'\t' -v group="$group_label" '
        # Read the event-to-gene table into an associative array
        NR == FNR {
            if (FNR == 1) next; # Skip the header of the first file
            split($3, event_parts, ";");
            event_id = event_parts[2];   # Extract the event_id (second part of column 3)
            gene_id = $2;               # Extract the gene_id (column 2)
            event_to_gene[event_id] = gene_id; # Map event_id to gene_id
            next;
        }

        # Process the GTF file
        FNR == 1 {
            next; # Skip the header of the second file
        }
        NR > 1 {
            # Extract start, end, scaffold, strand, and gene_id
            start = $4;
            end = $5;
            scaffold = $1;
            strand = $7;

            # Extract the event_id from the gene_id field in the GTF file
            split($9, gene_parts, "\"");
            full_event_id = gene_parts[2]; # Extract the full gene_id (e.g., SE:chr1:...)
            ensembl_gene_id = event_to_gene[full_event_id]; # Lookup the gene_id from the mapping

            # Print the output with Group and ENSEMBL_GENEID as the last columns
            print start "\t" end "\t" scaffold "\t" strand "\t" ensembl_gene_id "\t" group;
        }
    ' "$ioe_file" "$gtf_file"
}

# Specify the output directory and files
res_dir="$HOME/projects/SPLISOSM_paper/results/sit_nar_23/transcripts/"

for dataset in mob cbs
do
    echo "=== Processing $dataset ..."
    res_sample_dir="$res_dir/$dataset"
    output_dir="$res_sample_dir/matt.out"
    exon_tab=$output_dir/exon.matt.tab
    exon_w_efeatures_tab=$output_dir/exon.with_efeatures.matt.tab

    mkdir -p "$output_dir"

    # Convert SUPPA SE events to MATT tab format
    # Tab header: START END SCAFFOLD STRAND ENSEMBL_GENEID
    printf "START\tEND\tSCAFFOLD\tSTRAND\tENSEMBL_GENEID\tGROUP\n" > "$exon_tab"

    # Process both GTF/Event table pairs and combine the outputs
    {
        # Process the SE events of SVS genes
        process_suppa_gtf "$res_sample_dir"/suppa.out/"$dataset"_svs/_SE_strict.gtf \
            "$res_sample_dir"/suppa.out/"$dataset"_svs/_SE_strict.ioe "$dataset"_svs

        # Process the MX events of SVS genes
        process_suppa_gtf "$res_sample_dir"/suppa.out/"$dataset"_svs/_MX_strict.gtf \
            "$res_sample_dir"/suppa.out/"$dataset"_svs/_MX_strict.ioe "$dataset"_svs

        # Process the SE events of SVENS genes
        process_suppa_gtf "$res_sample_dir"/suppa.out/"$dataset"_svens/_SE_strict.gtf \
            "$res_sample_dir"/suppa.out/"$dataset"_svens/_SE_strict.ioe "$dataset"_svens

        # Process the MX events of SVENS genes
        process_suppa_gtf "$res_sample_dir"/suppa.out/"$dataset"_svens/_MX_strict.gtf \
            "$res_sample_dir"/suppa.out/"$dataset"_svens/_MX_strict.ioe "$dataset"_svens

    } >> "$exon_tab"

    # Keep unique exons per group
    head -n 1 "$exon_tab" > "$exon_tab.tmp"
    tail -n +2 "$exon_tab" | sort -u -k6,6 -k3,3 -k4,4 -k1,1n -k2,2n >> "$exon_tab.tmp"
    mv "$exon_tab.tmp" "$exon_tab"

    # Run MATT to extract exon features
    cp "$exon_tab" "$exon_w_efeatures_tab"
    $MATT_DIR/matt get_efeatures "$exon_w_efeatures_tab" START END SCAFFOLD STRAND ENSEMBL_GENEID \
    "$REF_GTF" "$REF_FA" Mmus 150 -notrbts | $MATT_DIR/matt add_cols "$exon_w_efeatures_tab" -

    # Check number of exons in each group
    $MATT_DIR/matt col_uniq "$exon_tab" GROUP

    # Compare exon features between the two groups
    $MATT_DIR/matt cmpr_exons "$exon_tab" START END SCAFFOLD STRAND ENSEMBL_GENEID \
    "$REF_GTF" "$REF_FA" Mmus 150 GROUP["$dataset"_svs,"$dataset"_svens] "$output_dir/exon.compare" -notrbts

done