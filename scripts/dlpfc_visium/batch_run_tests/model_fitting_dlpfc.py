### Fit the parametric GLM and GLMM models on filtered TREND anndata
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios
from isosde.hyptest_glmm import SplisosmGLMM

# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/human_dlpfc/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/human_dlpfc/"
date = '0123'

if __name__ == '__main__':

    # process all 12 DLPFC samples
    sample_list = [
        '151507', '151508', '151509', '151510', 
        '151669', '151670', '151671', '151672', 
        '151673', '151674', '151675', '151676'
    ]

    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/model"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    for sample_id in sample_list:
        ## load the filtered peak-level anndata and SV test results
        adata_peak = sc.read(f"{res_dir}/anndata/sierra/{sample_id}.peak_filtered_{date}.h5ad")
        coords = adata_peak.obs.loc[:, ['array_row', 'array_col']]
        df_sv_pval = pd.read_csv(f"{res_dir}/sv_results/{sample_id}.sv_combined_{date}.csv", index_col=0)

        # focus on SVS genes only
        _is_sig_svs = (df_sv_pval['padj_hsic-ir'] < 0.05) & (df_sv_pval['count_avg'] > 0.5)
        _iso_ind = adata_peak.var['gene_symbol'].isin(df_sv_pval.index[_is_sig_svs])
        adata_peak_sv = adata_peak[:, _iso_ind].copy()

        print(f"Fitting {adata_peak_sv.var['gene_symbol'].nunique()} genes for {sample_id}")
        print(f"Average number of isoforms per gene: {adata_peak_sv.shape[1] / adata_peak_sv.var['gene_symbol'].nunique()}")

        # extract lists of isoform counts and ratios
        counts_list, _, gene_name_list, _ = extract_counts_n_ratios(adata_peak_sv, layer = 'counts', group_iso_by = 'gene_symbol')

        # set random seed
        torch.manual_seed(0)

        ## fit the GLM and GLMM models
        # set up and fit the GLM model
        glm = SplisosmGLMM(model_type = 'glm')
        glm.setup_data(
            data = counts_list, # list of length n_genes, each element is a (n_spots, n_isos) tensor
            coordinates = coords, # (n_spots, 2), 2D array/tensor of spatial coordinates
            design_mtx = None, # (n_spots, n_covariates), 2D array/tensor of covariates
            gene_names = gene_name_list, # list of length n_genes, gene names
            covariate_names = None, # list of length n_covariates, covariate names
            group_gene_by_n_iso = False, # whether to group genes by the number of isoforms for batch processing
        )
        glm.fit(
            n_jobs = 2, # number of cores to use
            batch_size = 1, # number of genes with the same number of isoforms to process in parallel per core
            quiet=True, print_progress=True,
            with_design_mtx = False, 
            random_seed = 0
        )
        glm.save(f"{res_dir}/model/{sample_id}.glm_{date}.pkl")

        # set up and fit the GLMM model
        glmm = SplisosmGLMM(model_type = 'glmm-full')
        glmm.setup_data(
            data = counts_list, # list of length n_genes, each element is a (n_spots, n_isos) tensor
            coordinates = coords, # (n_spots, 2), 2D array/tensor of spatial coordinates
            design_mtx = None, # (n_spots, n_covariates), 2D array/tensor of covariates
            gene_names = gene_name_list, # list of length n_genes, gene names
            covariate_names = None, # list of length n_covariates, covariate names
            group_gene_by_n_iso = True, # whether to group genes by the number of isoforms for batch processing
        )
        glmm.fit(
            n_jobs = 2, # number of cores to use
            batch_size = 20, # number of genes with the same number of isoforms to process in parallel per core
            quiet=True, print_progress=True,
            with_design_mtx = False, 
            from_null = False, refit_null = True, # for LLR test
            random_seed = 0
        )
        glmm.save(f"{res_dir}/model/{sample_id}.glmm_{date}.pkl")
