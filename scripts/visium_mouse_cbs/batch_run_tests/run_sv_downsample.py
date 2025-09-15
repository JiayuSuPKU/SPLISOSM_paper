### Run spatial variability tests with downsampled data
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from splisosm.utils import load_visium_sp_meta, extract_counts_n_ratios, extract_gene_level_statistics
from splisosm.hyptest_np import SplisosmNP

# data preprocessing functions
def _load_raw_anndata(sierra_out_dir, sp_meta_dir):
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
    df_var = peaks[0].str.split(':', expand=True)
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

    return adata

def _downsample_anndata(adata, downsample_factor=0.1):
    """Downsample the anndata object by a given factor."""
    if not (0 < downsample_factor <= 1):
        raise ValueError("Downsample factor must be between 0 and 1.")

    # sample a dropout rate per spot around the downsample factor
    dropout_rates = np.random.uniform(
        low=downsample_factor * 0.8, high=downsample_factor * 1.2, size=adata.shape[0]
    )
    dropout_rates = np.clip(dropout_rates, 0, 1)  # ensure rates are between 0 and 1
    counts_per_spot_after = np.asarray(adata.X.sum(1)) * dropout_rates[:, np.newaxis]
    counts_per_spot_after = np.round(counts_per_spot_after).astype(int).squeeze()

    # downsample the counts per spot
    adata_downsampled = sc.pp.downsample_counts(
        adata, counts_per_cell=counts_per_spot_after, replace=False, copy=True
    )
    adata_downsampled.layers['counts'] = adata_downsampled.X.copy()

    # filter out lowly expressed peaks
    sc.pp.filter_genes(adata_downsampled, min_cells=0.01 * adata_downsampled.shape[0])

    # extract gene symbols and peak ids
    df_iso_meta = adata_downsampled.var.copy() # gene_symbol, chr, position, strand, gene_id
    df_iso_meta['peak_id'] = adata_downsampled.var_names

    # prepare gene-level metadata
    df_gene_meta = df_iso_meta.groupby('gene_symbol').size().reset_index(name='n_peak')
    df_gene_meta = df_gene_meta.set_index('gene_symbol')

    print(f"Number of spots: {adata_downsampled.shape[0]}")
    print(f"Number of genes before QC: {df_gene_meta.shape[0]}")
    print(f"Number of peaks before QC: {adata_downsampled.shape[1]}")
    print(f"Average number of peaks per gene before QC: {adata_downsampled.shape[1] / df_gene_meta.shape[0]}")

    # calculate the total counts per gene
    mapping_matrix = pd.get_dummies(df_iso_meta['gene_symbol'])
    mapping_matrix = mapping_matrix.loc[df_iso_meta.index, df_gene_meta.index]
    isog_counts = adata_downsampled[:, mapping_matrix.index].layers['counts'] @ mapping_matrix

    # calculate mean and sd of total gene counts
    df_gene_meta['pct_spot_on'] = (isog_counts > 0).mean(axis = 0)
    df_gene_meta['count_avg'] = isog_counts.mean(axis = 0)
    df_gene_meta['count_std'] = isog_counts.std(axis = 0)

    # filter out lowly expressed genes
    _gene_keep = df_gene_meta['pct_spot_on'] > 0.01

    # filter out genes with single isoform
    _gene_keep = (df_gene_meta['n_peak'] > 1) & _gene_keep

    # filter for isoforms
    _iso_keep = df_iso_meta['gene_symbol'].isin(df_gene_meta.index[_gene_keep])

    # update feature meta
    df_gene_meta = df_gene_meta.loc[_gene_keep, :]
    adata_downsampled = adata_downsampled[:, _iso_keep]
    adata_downsampled.var = df_iso_meta.loc[_iso_keep, :].copy()

    print(f"Number of genes after QC: {sum(_gene_keep)}")
    print(f"Number of peaks after QC: {sum(_iso_keep)}")
    print(f"Average number of peaks per gene after QC: {sum(_iso_keep) / sum(_gene_keep)}")
    
    return adata_downsampled, df_gene_meta

def _run_sv_tests(adata):
    # extract lists of isoform counts and ratios
    counts_list, ratios_list, gene_name_list, ratio_obs = extract_counts_n_ratios(
        adata, layer = 'counts', group_iso_by = 'gene_symbol'
    )

    # extract gene-level statistics
    df_gene_meta = extract_gene_level_statistics(adata, layer = 'counts', group_iso_by = 'gene_symbol')

    # spatial coordinates
    coords = adata.obs.loc[:, ['array_row', 'array_col']]

    # non-parametric testings
    model_np = SplisosmNP()
    model_np.setup_data(data = counts_list, coordinates = coords, gene_names = gene_name_list)

    # run all SV tests
    df_sv_res = {}
    for _test_method in ['hsic-ir', 'hsic-gc']:
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
        'pvalue_hsic-gc': df_sv_res['hsic-gc']['pvalue'].values,
        'padj_hsic-gc': df_sv_res['hsic-gc']['pvalue_adj'].values,
        }, index=gene_name_list))
    df_sv_pval = df_sv_pval.sort_values('pvalue_hsic-ir', ascending=True)

    return df_sv_pval

# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"

if __name__ == '__main__':

    # load raw anndata
    adata = _load_raw_anndata(
        f'{data_dir}/counts_no_cutoff', f'{data_dir}/visium_spatial_meta'
    )

    # set random seed for reproducibility
    np.random.seed(42)

    # run downsampling at multiple levels
    df_nsig_summary = []
    for downsample_factor in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8]:
        # for rep in range(1, 3):
        for rep in range(1, 2): # more factors, less replicates
            # downsample the anndata
            print(f"\n=== Downsample factor: {downsample_factor}, replicate: {rep} ===")
            adata_downsampled, _ = _downsample_anndata(
                adata, downsample_factor=downsample_factor
            )
            # run SV tests
            df_sv_pval = _run_sv_tests(adata_downsampled)

            # count number of significant genes and save results
            for method in ['hsic-ir', 'hsic-gc']:
                n_sig = sum(df_sv_pval[f'padj_{method}'] < 0.01)
                df_nsig_summary.append({
                    'downsample_factor': downsample_factor,
                    'replicate': rep,
                    'method': method,
                    'n_significant_genes': n_sig,
                    'n_tested_genes': df_sv_pval.shape[0],
                    'sum_count_avg': df_sv_pval['count_avg'].sum(),
                    'median_spot_on': df_sv_pval['pct_spot_on'].median(),
                })
                print(f"Method: {method}, number of significant genes: {n_sig}.")
    
    df_nsig_summary = pd.DataFrame(df_nsig_summary)
    df_nsig_summary.to_csv(f"{res_dir}/sv_results/downsampled_nsig_summary.csv", index=False)
    print("\nSummary of number of significant genes saved.")

    


