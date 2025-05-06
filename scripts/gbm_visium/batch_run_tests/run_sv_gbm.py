### Run spatial variability tests on transcriptome 3'end diversity (TREND) events for all 13 GBM samples
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios, extract_gene_level_statistics
from isosde.hyptest_np import SplisosmNP

# data preprocessing functions
def _prepare_anndata(sample_id, data_dir):
    ## Prepare the Visium adata
    print(f"=== Processing sample {sample_id} (gene-level Visium) ===")
    # load space ranger outputs
    adata_visium = sc.read_visium(f"{data_dir}/tirosh_inputs/general/GBM_data/{sample_id}/outs/")

    # load the spot annotations
    meta = pd.read_csv(f'{data_dir}/tirosh_inputs/general/visium_metadata.csv')
    meta = meta[meta['sample'] == sample_id]
    print(f"Number of spots: {adata_visium.shape[0]}")
    print(f"Number of spots with annotations: {meta.shape[0]}")

    # remove spots without annotations
    adata_visium = adata_visium[adata_visium.obs.index.isin(meta['spot_id'])]
    adata_visium.obs = adata_visium.obs.join(meta.set_index('spot_id'), how = 'left')

    # Preprocess the Visium adata
    adata_visium.raw = adata_visium.copy() # raw counts data
    adata_visium.var['gene_symbol'] = adata_visium.var_names
    adata_visium.var_names_make_unique()

    # filter out barely expressed genes
    sc.pp.filter_genes(adata_visium, min_cells=0.01 * adata_visium.shape[0])
    adata_visium.layers['counts'] = adata_visium.X.copy()

    # log normalize
    sc.pp.normalize_total(adata_visium, target_sum=1e4)
    sc.pp.log1p(adata_visium)
    adata_visium.layers['log1p'] = adata_visium.X.copy()
    
    ## (2) Prepare the peak-level adata
    print(f"=== Processing sample {sample_id} (peak-level TREND) ===")
    # load the conversion table from GEO IDs to sample IDs
    geo2sample = pd.read_table(
        f'{data_dir}/GEO_sample_list.tsv', sep = '\t', header = None, 
        names = ['sample', 'GEO', 'group', 'description']
    )
    geo_id = geo2sample[geo2sample['sample'] == sample_id]['GEO'].values[0]
    # sierra_out_dir = f'{data_dir}/sierra_peaks_default/{geo_id}/counts_default' # date = '0104'
    # sierra_out_dir = f'{data_dir}/sierra_peaks_individual_no_cutoff/{geo_id}/counts_no_cutoff' # date = '0110'
    sierra_out_dir = f'{data_dir}/sierra_peaks_default_individual/{geo_id}/counts_default' # date = '0310'

    # load Sierra peak counts
    adata_peak = sc.read(
        f"{sierra_out_dir}/matrix.mtx.gz", 
        cache_compression='cache_compression',
    ).T

    # load peak metadata and add to adata.var
    peaks = pd.read_csv(
        f"{sierra_out_dir}/sitenames.tsv.gz",
        header=None,
        sep="\t",
    )
    df_var = peaks[0].str.split(':', expand = True)
    df_var.columns = ['gene_symbol', 'chr', 'position', 'strand']
    df_var.index = peaks[0].values
    adata_peak.var_names = peaks[0].values
    adata_peak.var = df_var
    adata_peak.var['gene_id'] = adata_peak.var['gene_symbol']

    # # keep only peak detected in more than one sample aka "Merged" # date = '0104'
    # merged_peaks = pd.read_table(f"{data_dir}/sierra_peaks_default/merged.txt")
    # adata_peak = adata_peak[:, adata_peak.var_names.isin(
    #     merged_peaks.loc[merged_peaks['PeakClass'] == 'Merged', 'polyA_ID']
    # )]
    # adata_peak.var = adata_peak.var.join(merged_peaks.set_index('polyA_ID'), how = 'left')

    # keep all peaks per individual
    # peaks = pd.read_table(f"{data_dir}/sierra_peaks_individual_no_cutoff/{geo_id}/peak_no_cutoff.txt") # date = '0110'
    peaks = pd.read_table(f"{data_dir}/sierra_peaks_default_individual/{geo_id}/peak_default.txt") # date = '0310'
    adata_peak.var = adata_peak.var.join(peaks.set_index('polyA_ID'), how = 'left')

    # filter out peaks for novel genes aka starting with "ENSG"
    adata_peak = adata_peak[:, ~adata_peak.var['gene_symbol'].str.startswith('ENSG')]

    # load barcode metadata
    barcodes = pd.read_csv(f"{sierra_out_dir}/barcodes.tsv.gz", header=None, sep="\t")
    adata_peak.obs_names = barcodes[0].values
    # remove spots without annotations
    adata_peak = adata_peak[adata_peak.obs_names.isin(meta['spot_id']), :]
    # reorder spots to match Visium anndata
    adata_peak = adata_peak[adata_visium.obs.index, :]

    # add Visium spatial metadata to peak-level adata
    adata_peak.obs = adata_peak.obs.join(adata_visium.obs, how = 'left')
    adata_peak.uns = adata_visium.uns
    adata_peak.obsm = adata_visium.obsm

    # Preprocess the peak-level adata
    adata_peak.raw = adata_peak.copy() # raw counts data
    adata_peak.layers['counts'] = adata_peak.X.copy()

    # log normalize
    sc.pp.normalize_total(adata_peak, target_sum=1e4)
    sc.pp.log1p(adata_peak)
    adata_peak.layers['log1p'] = adata_peak.X.copy()

    # filter out lowly expressed peaks
    sc.pp.filter_genes(adata_peak, min_cells=0.01 * adata_peak.shape[0])

    # extract gene symbols and peak ids
    df_iso_meta = adata_peak.var.copy() # gene_symbol, chr, position, strand, gene_id
    df_iso_meta['peak_id'] = adata_peak.var_names

    # prepare gene-level metadata
    df_gene_meta = df_iso_meta.groupby('gene_symbol').size().reset_index(name='n_peak')
    df_gene_meta = df_gene_meta.set_index('gene_symbol')

    print(f"Number of spots: {adata_peak.shape[0]}")
    print(f"Number of genes before QC: {df_gene_meta.shape[0]}")
    print(f"Number of peaks before QC: {adata_peak.shape[1]}")
    print(f"Average number of peaks per gene before QC: {adata_peak.shape[1] / df_gene_meta.shape[0]}")

    # calculate the total counts per gene
    mapping_matrix = pd.get_dummies(df_iso_meta['gene_symbol'])
    mapping_matrix = mapping_matrix.loc[df_iso_meta.index, df_gene_meta.index]
    isog_counts = adata_peak[:, mapping_matrix.index].layers['counts'] @ mapping_matrix

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
    adata_peak = adata_peak[:, _iso_keep]
    adata_peak.var = df_iso_meta.loc[_iso_keep, :].copy()

    print(f"Number of genes after QC: {sum(_gene_keep)}")
    print(f"Number of peaks after QC: {sum(_iso_keep)}")
    print(f"Average number of peaks per gene after QC: {sum(_iso_keep) / sum(_gene_keep)}")

    return adata_visium, adata_peak, df_gene_meta


# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/gbm_visium_cell_24/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/gbm_visium/"
# date = '0104'
# date = '0110'
date = '0310'

if __name__ == '__main__':
    # create the results directory
    for _dir in [f"{res_dir}/anndata/visium", f"{res_dir}/anndata/sierra", f"{res_dir}/sv_results"]:
        os.makedirs(_dir, exist_ok=True)

    # process all 13 GBM samples
    df_sample = pd.read_table(
        f'{data_dir}/GEO_sample_list.tsv', sep = '\t', header = None,
        names = ['sample', 'GEO', 'group', 'description']
    )
    df_sample = df_sample[df_sample['group'] == 'GBM']

    for sample_id in df_sample['sample']:
        # prepare gene-level and peak-level adata
        adata_visium, adata_peak, _ = _prepare_anndata(sample_id, data_dir)

        # extract lists of isoform counts and ratios
        counts_list, ratios_list, gene_name_list, ratio_obs = extract_counts_n_ratios(
            adata_peak, layer = 'counts', group_iso_by = 'gene_symbol'
        )
        adata_peak.layers['ratios_obs'] = ratio_obs.copy()

        # save anndata
        adata_visium.write(f"{res_dir}/anndata/visium/{sample_id}.visium_{date}.h5ad")
        adata_peak.write(f"{res_dir}/anndata/sierra/{sample_id}.peak_filtered_{date}.h5ad")

        # extract gene-level statistics
        df_gene_meta = extract_gene_level_statistics(adata_peak, layer = 'counts', group_iso_by = 'gene_symbol')

        # spatial coordinates
        coords = adata_peak.obs.loc[:, ['array_row', 'array_col']]

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
        df_sv_pval.to_csv(f"{res_dir}/sv_results/{sample_id}.sv_combined_{date}.csv", index=True)
