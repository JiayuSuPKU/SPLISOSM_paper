### Run spatial variability tests on one mouse olfactory bulb (MOB) and two coronal brain section (CBS) samples
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios, extract_gene_level_statistics
from isosde.hyptest_np import SplisosmNP

def parse_gtf(gtf_file):
    """Parse the GTF file to build a dictionary mapping transcript IDs to transcript names."""
    transcript_dict = {}

    # open and read the GTF file
    with open(gtf_file, "r") as file:
        for line in file:
            # skip header lines
            if line.startswith("#"):
                continue

            # split the GTF line into columns
            fields = line.strip().split("\t")
            if fields[2] == "transcript":  # Only process transcript entries
                attributes = fields[8]  # The attributes column
                transcript_id = None
                transcript_name = None

                # extract transcript_id and transcript_name from attributes
                for attribute in attributes.split(";"):
                    attribute = attribute.strip()
                    if attribute.startswith("transcript_id"):
                        transcript_id = attribute.split('"')[1].split(".")[0]  # Strip version
                    elif attribute.startswith("transcript_name"):
                        transcript_name = attribute.split('"')[1]

                # add to dictionary if both values are found
                if transcript_id and transcript_name:
                    transcript_dict[transcript_id] = transcript_name

    return transcript_dict

# data preprocessing functions
def _prepare_anndata(data_dir, sample_id):
    # load the ONT count data
    adata = sc.read_h5ad(f"{data_dir}/anndata/{sample_id}_iso.h5ad") # long-read isoform counts
    adata.raw = adata.copy() # raw counts data
    adata.layers['counts'] = adata.X.copy()

    # log normalize
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    adata.layers['log1p'] = adata.X.copy()

    # load the Visium spatial metadata
    adata = load_visium_sp_meta(adata, f"{data_dir}/visium_spatial_meta", library_id='adata_ont')

    # filter out lowly expressed isoforms
    # sc.pp.filter_genes(adata_ont, min_counts=10)
    sc.pp.filter_genes(adata, min_cells=0.01 * adata.shape[0])

    # extract gene symbols and isoform ids
    df_iso_meta = adata.var.copy()
    df_iso_meta['gene_symbol'] = df_iso_meta['features'].str.split('\\.\\.').str[0]
    df_iso_meta['transcript_id'] = df_iso_meta['features'].str.split('\\.\\.').str[1]

    # convert ensembl transcript ids to transcript names
    id2name = parse_gtf("/Users/jysumac/reference/mm10/Mus_musculus.GRCm38.102.gtf")
    id_list = df_iso_meta['transcript_id'].str.split('.').str[0]
    df_iso_meta['transcript_name'] = [id2name.get(_id, 'N/A') for _id in id_list]

    # prepare gene-level metadata
    df_gene_meta = df_iso_meta.groupby('gene_symbol').size().reset_index(name='n_iso')
    df_gene_meta = df_gene_meta.set_index('gene_symbol')

    print(f"Number of spots: {adata.shape[0]}")
    print(f"Number of genes before QC: {df_gene_meta.shape[0]}")
    print(f"Number of isoforms before QC: {adata.shape[1]}")
    print(f"Average number of isoforms per gene before QC: {adata.shape[1] / df_gene_meta.shape[0]}")

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
    _gene_keep = (df_gene_meta['n_iso'] > 1) & _gene_keep

    # filter for isoforms
    _iso_keep = df_iso_meta['gene_symbol'].isin(df_gene_meta.index[_gene_keep])

    # update feature meta
    df_gene_meta = df_gene_meta.loc[_gene_keep, :]
    adata = adata[:, _iso_keep]
    adata.var = df_iso_meta.loc[_iso_keep, :].copy()

    print(f"Number of genes after QC: {sum(_gene_keep)}")
    print(f"Number of isoforms after QC: {sum(_iso_keep)}")
    print(f"Average number of isoforms per gene after QC: {sum(_iso_keep) / sum(_gene_keep)}")

    return adata, df_gene_meta

# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/sit_nar_23/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/sit_nar_23/"
date = '1107'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/anndata", f"{res_dir}/sv_results"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    for _sample_id in ['mob', 'cbs1', 'cbs2']:
        # load and prepare ONT anndata
        adata_ont, _ = _prepare_anndata(f'{data_dir}/{_sample_id}', _sample_id)

        # extract lists of isoform counts and ratios
        counts_list, _, gene_name_list, ratio_obs = extract_counts_n_ratios(
            adata_ont, layer = 'counts', group_iso_by = 'gene_symbol'
        )
        adata_ont.layers['ratios_obs'] = ratio_obs.copy()

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
        df_sv_pval.to_csv(f"{res_dir}/sv_results/{_sample_id}_sv_combined_{date}.csv", index=True)
        
        # save anndata
        adata_ont.write(f"{res_dir}/anndata/{_sample_id}_ont_filtered_{date}.h5ad")
