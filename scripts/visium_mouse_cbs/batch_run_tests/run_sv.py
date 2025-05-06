### Run spatial variability tests on transcriptome 3'end diversity (TREND) events in the mouse CBS sample
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios, extract_gene_level_statistics
from isosde.hyptest_np import SplisosmNP

# data preprocessing functions
def _prepare_anndata(sierra_out_dir, sp_meta_dir):
    # load the Sierra outputs
    adata = sc.read(
        f"{sierra_out_dir}/matrix.mtx.gz",
        cache_compression='cache_compression',
    ).T

    # load peak metadata
    peaks = pd.read_csv(
        f"{sierra_out_dir}/sitenames.tsv.gz",
        header=None,
        sep="\t",
    )
    df_var = peaks[0].str.split(':', expand = True)
    df_var.columns = ['gene_symbol', 'chr', 'position', 'strand']
    df_var.index = peaks[0].values

    # load barcode metadata
    barcodes = pd.read_csv(f"{sierra_out_dir}/barcodes.tsv.gz", header=None)
    adata.var_names = peaks[0].values
    adata.obs_names = barcodes[0].values
    adata.var = df_var
    adata.var['gene_id'] = adata.var['gene_symbol']

    # load Visium spatial metadata
    adata = load_visium_sp_meta(adata, f"{sp_meta_dir}/", library_id='adata_peak')
    adata = adata[adata.obs['in_tissue'].astype(bool), :].copy()

    # store raw counts
    adata.raw = adata.copy() # raw counts data
    adata.layers['counts'] = adata.X.copy()

    # log normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

    # filter out lowly expressed peaks
    sc.pp.filter_genes(adata, min_cells=0.01 * adata.shape[0])

    # extract gene symbols and peak ids
    df_iso_meta = adata.var.copy() # gene_symbol, chr, position, strand, gene_id
    df_iso_meta['peak_id'] = adata.var_names

    # prepare gene-level metadata
    df_gene_meta = df_iso_meta.groupby('gene_symbol').size().reset_index(name='n_peak')
    df_gene_meta = df_gene_meta.set_index('gene_symbol')

    print(f"Number of spots: {adata.shape[0]}")
    print(f"Number of genes before QC: {df_gene_meta.shape[0]}")
    print(f"Number of peaks before QC: {adata.shape[1]}")
    print(f"Average number of peaks per gene before QC: {adata.shape[1] / df_gene_meta.shape[0]}")

    # calculate the total counts per gene
    mapping_matrix = pd.get_dummies(df_iso_meta['gene_symbol'])
    mapping_matrix = mapping_matrix.loc[df_iso_meta.index, df_gene_meta.index]
    isog_counts = adata[:, mapping_matrix.index].layers['counts'] @ mapping_matrix

    # calculate mean and sd of total gene counts
    df_gene_meta['pct_spot_on'] = (isog_counts > 0).mean(axis = 0)
    df_gene_meta['count_avg'] = isog_counts.mean(axis = 0)
    df_gene_meta['count_std'] = isog_counts.std(axis = 0)

    # filter out lowly expressed genes
    _gene_keep = df_gene_meta['pct_spot_on'] > 0.01
    # _gene_keep = (df_gene_meta['count_avg'] > 0.5) & _gene_keep

    # filter out genes with single isoform
    _gene_keep = (df_gene_meta['n_peak'] > 1) & _gene_keep

    # filter for isoforms
    _iso_keep = df_iso_meta['gene_symbol'].isin(df_gene_meta.index[_gene_keep])

    # update feature meta
    df_gene_meta = df_gene_meta.loc[_gene_keep, :]
    adata = adata[:, _iso_keep]
    adata.var = df_iso_meta.loc[_iso_keep, :].copy()

    print(f"Number of genes after QC: {sum(_gene_keep)}")
    print(f"Number of peaks after QC: {sum(_iso_keep)}")
    print(f"Average number of peaks per gene after QC: {sum(_iso_keep) / sum(_gene_keep)}")

    return adata, df_gene_meta


# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
date = '1119'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/anndata", f"{res_dir}/sv_results"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    # load and prepare peak anndata
    adata, _ = _prepare_anndata(f'{data_dir}/counts_no_cutoff', f'{data_dir}/visium_spatial_meta')

    # extract lists of isoform counts and ratios
    counts_list, ratios_list, gene_name_list, ratio_obs = extract_counts_n_ratios(
        adata, layer = 'counts', group_iso_by = 'gene_symbol'
    )
    adata.layers['ratios_obs'] = ratio_obs.copy()

    # extract gene-level statistics
    df_gene_meta = extract_gene_level_statistics(adata, layer = 'counts', group_iso_by = 'gene_symbol')

    # spatial coordinates
    coords = adata.obs.loc[:, ['array_row', 'array_col']]

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
    df_sv_pval.to_csv(f"{res_dir}/sv_results/cbs_peak_sv_combined_{date}.csv", index=True)
    
    # save anndata
    adata.write(f"{res_dir}/anndata/cbs_peak_filtered_{date}.h5ad")
