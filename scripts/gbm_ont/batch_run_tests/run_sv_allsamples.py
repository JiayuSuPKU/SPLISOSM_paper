### Run spatial variability tests on all GBM/DMG ONT samples
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios, extract_gene_level_statistics
from isosde.hyptest_np import SplisosmNP

# data preprocessing functions
def load_ont_adata(data_dir, sample_id):
    """Load and filter ONT anndata to remove lowly expressed genes and single isoform genes"""
    # load the ONT count data
    adata = sc.read_h5ad(f"{data_dir}/anndata/{sample_id}/ont.h5ad") # ONT isoform counts

    # load the isoform-level annotation and add to adata.var
    df_iso_meta = pd.read_csv(
        f"{data_dir}/anndata/{sample_id}/ont_features.csv", 
        index_col='transcript_id',
        na_values=['', 'NA'],
        dtype = str
    )
    adata.var = pd.merge(adata.var, df_iso_meta, left_index=True, right_on='transcript_id')

    # save raw data
    adata.raw = adata.copy() # raw counts data
    adata.layers['counts'] = adata.X.copy()

    # keep only reference isoforms and novel isoforms from known genes
    adata = adata[:, adata.var['group'].isin(['reference transcript', 'novel transcript'])]
    # filter out multi-mapping isoforms
    adata = adata[:, ~adata.var['gene_name'].isna()]
    adata.var['gene_symbol'] = adata.var['gene_name'].astype(str)

    # log normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

    # load the Visium spatial metadata
    adata = load_visium_sp_meta(adata, f"{data_dir}/visium_spatial_meta/{sample_id}", library_id='adata_ont')

    print(f"Number of spots: {adata.shape[0]}")
    print(f"Number of genes before QC: {adata.var['gene_symbol'].nunique()}")
    print(f"Number of isoforms before QC: {adata.shape[1]}")
    print(f"Average number of isoforms per gene before QC: {adata.shape[1] / adata.var['gene_symbol'].nunique()}")

    # first filter out lowly expressed isoforms
    sc.pp.filter_genes(adata, min_cells=0.05 * adata.shape[0])

    # prepare gene-level metadata
    df_iso_meta = adata.var.copy()
    df_gene_meta = df_iso_meta.groupby('gene_symbol').size().reset_index(name='n_iso')
    df_gene_meta = df_gene_meta.set_index('gene_symbol')

    # calculate the total counts per gene
    mapping_matrix = pd.get_dummies(df_iso_meta['gene_symbol'])
    mapping_matrix = mapping_matrix.loc[df_iso_meta.index, df_gene_meta.index]
    isog_counts = adata[:, mapping_matrix.index].layers['counts'] @ mapping_matrix

    # calculate mean and sd of total gene counts
    df_gene_meta['pct_spot_on'] = (isog_counts > 0).mean(axis = 0)
    df_gene_meta['count_avg'] = isog_counts.mean(axis = 0)
    df_gene_meta['count_std'] = isog_counts.std(axis = 0)

    # filter out lowly expressed genes by total gene-level counts
    _gene_keep = df_gene_meta['pct_spot_on'] > 0.05
    # _gene_keep = (df_gene_meta['count_avg'] > 0.5) & _gene_keep

    # filter out genes with single isoform
    _gene_keep = (df_gene_meta['n_iso'] > 1) & _gene_keep

    # keep only isoforms from genes that pass QC
    _iso_keep = df_iso_meta['gene_symbol'].isin(df_gene_meta.index[_gene_keep])

    # update feature meta
    df_gene_meta = df_gene_meta.loc[_gene_keep, :]
    adata = adata[:, _iso_keep]
    # fill NaN values with empty string to avoid errors when saving the anndata
    adata.var = df_iso_meta.loc[_iso_keep, :].copy().fillna('')

    print(f"Number of genes after QC: {sum(_gene_keep)}")
    print(f"Number of isoforms after QC: {sum(_iso_keep)}")
    print(f"Average number of isoforms per gene after QC: {sum(_iso_keep) / sum(_gene_keep)}")

    return adata, df_gene_meta


# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/gbm_ont_nc_23/all_spots/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/gbm_ont/"
date = '0104'

if __name__ == '__main__':
    
    list_sample_ids = [
        'DMG_1', 'DMG_2', 'DMG_3', 'DMG_4', 'DMG_5',
        'GBM_1', 'GBM_2', 'GBM_3', 'GBM_4', 'GBM_5', 'GBM_6'
    ]

    # create the results directory if it does not exist
    for _dir in [f"{res_dir}/anndata/visium", f"{res_dir}/anndata/ont", f"{res_dir}/sv_results"]:
        os.makedirs(_dir, exist_ok=True)

    # loop over all samples
    for _sample_id in list_sample_ids:
        print(f"=== Running spatial variability tests for {_sample_id} ...")
        
        # load filtered ONT anndata
        adata_ont, _ = load_ont_adata(data_dir, _sample_id)

        # rename novel isoform to gene_symbol_Iso_N
        df_iso_meta = adata_ont.var.reset_index()
        df_iso_meta['transcript_name'] = df_iso_meta.groupby('gene_symbol', observed = True).cumcount() + 1
        df_iso_meta['transcript_name'] = df_iso_meta['gene_symbol'].astype(str) + "_Iso_" + \
            df_iso_meta['transcript_name'].astype(str)
        df_iso_meta.loc[df_iso_meta['is_novel'] == 'FALSE', 'transcript_name'] = \
            df_iso_meta.loc[df_iso_meta['is_novel'] == 'FALSE', 'transcript_id']
        
        # reorganize columns
        df_iso_meta = df_iso_meta[
            ['transcript_name', 'gene_symbol', 'transcript_id', 'n_cells', 'is_novel', 'group', 'id_to_match', 'source',
             'chr', 'start', 'end', 'strand', 'gene_id', 'gene_name']
        ]
        
        adata_ont.var = df_iso_meta.set_index('transcript_name')

        # extract lists of isoform counts and ratios
        counts_list, _, gene_name_list, ratio_obs = extract_counts_n_ratios(
            adata_ont, layer = 'counts', group_iso_by = 'gene_symbol'
        )
        adata_ont.layers['ratios_obs'] = ratio_obs.copy()

        # save anndata
        adata_ont.write(f"{res_dir}/anndata/ont/{_sample_id}.ont_filtered_{date}.h5ad")

        # add major ratio avg to df_gene_meta
        df_gene_meta = extract_gene_level_statistics(adata_ont, layer = 'counts', group_iso_by = 'gene_symbol')

        # spatial coordinates
        coords = adata_ont.obs.loc[:, ['array_row', 'array_col']]

        # non-parametric testings
        model_np = SplisosmNP()
        model_np.setup_data(data = counts_list, coordinates = coords, gene_names = gene_name_list)

        # run all SV tests
        df_sv_res = {}
        for _test_method in ['hsic-ir', 'hsic-ic', 'hsic-gc', 'spark-x']:
            model_np.test_spatial_variability(
                method = _test_method,
                ratio_transformation = 'none', # only applicable to 'hsic-ir'
                nan_filling = 'mean', # how to fill NaN values in the data, can be 'mean' (global mean), 'none' (ignoring NaN spots)
                return_results = False, 
                print_progress = True
            )
            df_sv_res[_test_method] = model_np.get_formatted_test_results(test_type = 'sv') # per gene test statistics

        # merge SV test results
        df_sv_pval = df_gene_meta.join(pd.DataFrame({
            'pvalue_hsic-ir': df_sv_res['hsic-ir']['pvalue'].values,
            'padj_hsic-ir': df_sv_res['hsic-ir']['pvalue_adj'].values,
            'pvalue_hsic-ic': df_sv_res['hsic-ic']['pvalue'].values,
            'padj_hsic-ic': df_sv_res['hsic-ic']['pvalue_adj'].values,
            'pvalue_hsic-gc': df_sv_res['hsic-gc']['pvalue'].values,
            'padj_hsic-gc': df_sv_res['hsic-gc']['pvalue_adj'].values,
            'pvalue_spark-x': df_sv_res['spark-x']['pvalue'].values,
            'padj_spark-x': df_sv_res['spark-x']['pvalue_adj'].values,
            }, index=gene_name_list))
        df_sv_pval = df_sv_pval.sort_values('pvalue_hsic-ir', ascending=True)

        # save results
        df_sv_pval.to_csv(f"{res_dir}/sv_results/{_sample_id}.sv_combined_{date}.csv", index=True)

