#! /bin/bash

# This script runs SUPPA on the transcripts extracted from the reference GTF file
# to find the splicing events in the transcripts

SUPPA_DIR="/Users/jysumac/Projects/Packages/SUPPA/"
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
RSCRIPT_PATH="$SCRIPT_DIR/extract_gtf_by_txid.R"

ref_gtf_file="/Users/jysumac/reference/mm10/Mus_musculus.GRCm38.102.gtf"
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

    # Extract AS events for the two groups of transcripts
    for tx_id_file in "$dataset"_svs.txid.csv "$dataset"_svens.txid.csv
    do
        sample_id=$(basename $tx_id_file .txid.csv)
    
        # Extract the transcripts from the reference GTF file
        /usr/local/bin/Rscript "$RSCRIPT_PATH" $res_sample_dir/$tx_id_file $ref_gtf_file
        sub_gtf=$(basename $tx_id_file .csv).gtf

        # Run SUPPA to find the splicing events
        if [ ! -d "$res_sample_dir/suppa.out/$sample_id" ]; then
            mkdir -p "$res_sample_dir/suppa.out/$sample_id"
        fi

        python $SUPPA_DIR/suppa.py generateEvents -i $sub_gtf -o "$res_sample_dir/suppa.out/$sample_id/" -f ioe -e SE SS MX RI FL

        # Empty the file
        true > "$res_sample_dir/suppa.out/$sample_id/n_events.txt"

        # Count number of events in each category and store in a file
        for event in A3 A5 AF AL MX RI SE
        do
            file=$res_sample_dir/suppa.out/$sample_id/_"$event"_strict.ioe
            n_events=$(( $(wc -l < "$file") - 1 ))
            echo $event $n_events >> "$res_sample_dir/suppa.out/$sample_id/n_events.txt"
        done
    done

done