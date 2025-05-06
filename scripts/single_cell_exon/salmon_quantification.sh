#!/bin/bash

### Quantify isoform expression using salmon (transcriptomic alignment)

# Set paths to input and output directories
data_dir='/gpfs/commons/home/jsu/data/rbp_ko/rbfox_tko_mouse_neuron'
salmon_ind='/gpfs/commons/home/jsu/reference/salmon_ind/mm10/default'
salmon_dir='/gpfs/commons/home/jsu/softwares/salmon-latest_linux_x86_64/bin/'

mkdir -p $data_dir/salmon_quant

# Quantify using salmon
for sample in 'rbfox_tko_d10_a' 'rbfox_tko_d10_b' 'rbfox_tko_d10_c' 'rbfox_wt_d10_a' 'rbfox_wt_d10_b' 'rbfox_wt_d10_c'; do
  echo $sample
  $salmon_dir/salmon quant \
    -i $salmon_ind -l A -p 32 \
    -1 "$data_dir/rnaseq_data/$sample"_1.fasta.gz \
    -2 "$data_dir/rnaseq_data/$sample"_2.fasta.gz \
    --validateMappings --softclip \
    -o $data_dir/salmon_quant/$sample
done

# merge salmon results (estimated number of reads)
$salmon_dir/salmon quantmerge \
  --quants $(ls -d $data_dir/salmon_quant/*/) \
  --names $(ls -d $data_dir/salmon_quant/*/ | xargs -I{} basename {}) \
  --column numreads \
  -o $data_dir/salmon_quant/all_samples.numreads.tsv

# merge salmon results (TPM)
$salmon_dir/salmon quantmerge \
  --quants $(ls -d $data_dir/salmon_quant/*/) \
  --names $(ls -d $data_dir/salmon_quant/*/ | xargs -I{} basename {}) \
  --column TPM \
  -o $data_dir/salmon_quant/all_samples.tpm.tsv