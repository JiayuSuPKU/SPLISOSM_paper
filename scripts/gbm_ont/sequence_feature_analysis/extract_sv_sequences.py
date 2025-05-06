### Extract isoform-level annotations for SV genes per sample
import os
import pyranges as pr
import pandas as pd
import itertools
import scanpy as sc


# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/gbm_ont_nc_23/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/gbm_ont/"
date = '0104'

if __name__ == '__main__':

    list_sample_ids = [
        'DMG_1', 'DMG_2', 'DMG_3', 'DMG_4', 'DMG_5',
        'GBM_1', 'GBM_2', 'GBM_3', 'GBM_4', 'GBM_5', 'GBM_6'
    ]

    # load reference GTF annotation
    reference_gtf = pr.read_gtf("/Users/jysumac/reference/hg38/gencode.v32.annotation.gtf")
    # keep only transcript and exon-related features
    reference_gtf = reference_gtf[~reference_gtf.transcript_id.isna()]
    reference_gtf = reference_gtf[reference_gtf.Feature.isin([
        'transcript', 'exon', 'CDS', 'UTR', 'start_codon', 'stop_codon'
    ])]

    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/transcripts"]:
        os.makedirs(_dir, exist_ok=True)

    for sample_id in list_sample_ids:
        ## create the results subdirectory
        os.makedirs(f"{res_dir}/transcripts/{sample_id}", exist_ok = True)

        ## (1) load the SV test results
        # load SV test results
        df_sv_pval = pd.read_csv(f"{res_dir}/sv_results/{sample_id}.sv_combined_{date}.csv", index_col=0)
        df_sv_pval['is_ont_svs'] = df_sv_pval['padj_hsic-ir'] < 0.05
        df_sv_pval['is_ont_sve'] = df_sv_pval['padj_hsic-gc'] < 0.05

        # load the ONT count data
        adata_ont = sc.read_h5ad(f"{res_dir}/anndata/ont/{sample_id}.ont_filtered_{date}.h5ad")

        ## (2) load GTF annotations for reference and novel transcripts
        # load GTF of novel transcripts
        novel_gtf = pr.read_gtf(f"{data_dir}/new_isoform_gtf/{sample_id.replace("_", "")}_NewIsoform.gtf")
        novel_gtf = novel_gtf[~novel_gtf.transcript_id.isna()]
        # collapse transcripts with the same transcript_id
        novel_gtf = novel_gtf.merge(by = ['Feature', 'Source', 'transcript_id'], strand = True)

        # merge reference and novel GTF
        merged_df = pr.concat([reference_gtf, novel_gtf]).df
        merged_df['id_to_match'] = merged_df['transcript_name'].fillna(merged_df['transcript_id'])

        ## (3) extract isoform-level annotations for SVS genes
        # focus on all SVS genes
        _iso_ind = adata_ont.var['gene_symbol'].isin(df_sv_pval.index[df_sv_pval['is_ont_svs']])
        iso_svs = adata_ont[:, _iso_ind].var.reset_index()[
            ['transcript_name', 'gene_symbol', 'transcript_id', 'n_cells', 'is_novel', 'group', 'id_to_match', 
                'source', 'chr', 'start', 'end', 'strand', 'gene_id', 'gene_name']
        ]

        # save isoform info for SVS genes
        iso_svs.to_csv(f"{res_dir}/transcripts/{sample_id}/{sample_id}.svs_iso.csv", header = True, index = False)
        
        print(f"({sample_id}) Saving {iso_svs['transcript_id'].nunique()} isoforms of {iso_svs['gene_symbol'].nunique()} SVS genes.")

        # extract GTF for isoforms of SVS genes
        svs_df = merged_df[merged_df['id_to_match'].isin(iso_svs['id_to_match'])].reset_index(drop = True)
        svs_df = iso_svs.set_index('id_to_match').join(svs_df.set_index('id_to_match'), how = 'inner', rsuffix = '_r')
        svs_df['Score'] = svs_df['Score'].fillna('.')
        svs_df['Frame'] = svs_df['Frame'].fillna('.')

        print(f"({sample_id}) Found {svs_df['transcript_id'].nunique()}/{iso_svs['transcript_id'].nunique()} isoforms "
                f"of {svs_df['gene_symbol'].nunique()}/{iso_svs['gene_symbol'].nunique()} SVS genes in the GTF.")

        # ensure all columns are not Categorical
        svs_df = svs_df.astype({col: str for col in svs_df.columns if svs_df[col].dtype.name == 'category'})
        # save GTF for isoforms of SVS genes
        svs_gtf = pr.PyRanges(svs_df[
            ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
                'Frame', 'gene_symbol', 'gene_id', 'gene_name', 'transcript_id', 'transcript_name'
            ]
        ])
        svs_gtf.to_gtf(f"{res_dir}/transcripts/{sample_id}/{sample_id}.svs_iso.gtf")

        ## (4) extract isoform-level annotations for SVE but not SVS genes
        # focus on all SVE but not SVS genes
        _iso_ind = adata_ont.var['gene_symbol'].isin(df_sv_pval.index[
            df_sv_pval['is_ont_sve'] & (~df_sv_pval['is_ont_svs'])
        ])
        iso_svens = adata_ont[:, _iso_ind].var.reset_index()[
            ['transcript_name', 'gene_symbol', 'transcript_id', 'n_cells', 'is_novel', 'group', 'id_to_match', 
                'source', 'chr', 'start', 'end', 'strand', 'gene_id', 'gene_name']
        ]
        
        # save isoform info for SVS genes
        iso_svens.to_csv(f"{res_dir}/transcripts/{sample_id}/{sample_id}.svens_iso.csv", header = True, index = False)
        
        print(f"({sample_id}) Saving {iso_svens['transcript_id'].nunique()} isoforms of {iso_svens['gene_symbol'].nunique()} SVENS genes.")

        # extract GTF for isoforms of SVS genes
        svens_df = merged_df[merged_df['id_to_match'].isin(iso_svens['id_to_match'])].reset_index(drop = True)
        svens_df = iso_svens.set_index('id_to_match').join(svens_df.set_index('id_to_match'), how = 'inner', rsuffix = '_r')
        svens_df['Score'] = svens_df['Score'].fillna('.')
        svens_df['Frame'] = svens_df['Frame'].fillna('.')

        print(f"({sample_id}) Found {svens_df['transcript_id'].nunique()}/{iso_svens['transcript_id'].nunique()} isoforms "
                f"of {svens_df['gene_symbol'].nunique()}/{iso_svens['gene_symbol'].nunique()} SVENS genes in the GTF.")

        # ensure all columns are not Categorical
        svens_df = svens_df.astype({col: str for col in svens_df.columns if svens_df[col].dtype.name == 'category'})
        # save GTF for isoforms of SVS genes
        svens_gtf = pr.PyRanges(svens_df[
            ['Chromosome', 'Source', 'Feature', 'Start', 'End', 'Score', 'Strand',
                'Frame', 'gene_symbol', 'gene_id', 'gene_name', 'transcript_id', 'transcript_name'
            ]
        ])
        svens_gtf.to_gtf(f"{res_dir}/transcripts/{sample_id}/{sample_id}.svens_iso.gtf")

                
        
