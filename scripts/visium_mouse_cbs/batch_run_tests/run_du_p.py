### Run the parametric differential usage tests on the mouse CBS visium sample
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios
from isosde.hyptest_glmm import SplisosmGLMM


# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
date = '1119'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/model", f"{res_dir}/du_results"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    ## (1) load the filtered peak counts anndata and SV test results
    adata_peak = sc.read(f"{res_dir}/anndata/cbs_peak_filtered_{date}.h5ad")
    df_sv_pval = pd.read_csv(f"{res_dir}/sv_results/cbs_peak_sv_combined_{date}.csv", index_col=0)

    # focus on SVS genes only
    _is_sig_svs = (df_sv_pval['padj_hsic-ir'] < 0.01) & (df_sv_pval['count_avg'] > 0.5) & \
        (~ df_sv_pval.index.str.startswith('Gm')) & \
        (~ df_sv_pval.index.str.endswith('Rik')) # exclude unannotated genes
    _iso_ind = adata_peak.var['gene_symbol'].isin(df_sv_pval.index[_is_sig_svs])
    adata_peak_sv = adata_peak[:, _iso_ind].copy()

    # extract lists of isoform counts and ratios
    counts_list, _, gene_name_list, _ = extract_counts_n_ratios(adata_peak_sv, layer = 'counts', group_iso_by = 'gene_symbol')
    coords = adata_peak_sv.obs.loc[:, ['array_row', 'array_col']]

    ## (2) fit the GLM and GLMM models
    print(f"Fitting {adata_peak_sv.var['gene_symbol'].nunique()} SVS genes")
    print(f"Average number of isoforms per gene: {adata_peak_sv.shape[1] / adata_peak_sv.var['gene_symbol'].nunique()}")

    # set random seed
    torch.manual_seed(0)

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
        n_jobs = 4, # number of cores to use
        batch_size = 1, # number of genes with the same number of isoforms to process in parallel per core
        quiet=True, print_progress=True,
        with_design_mtx = False, 
        random_seed = 0
    )
    glm.save(f"{res_dir}/model/cbs_peak_glm_{date}.pkl")

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
        n_jobs = 4, # number of cores to use
        batch_size = 20, # number of genes with the same number of isoforms to process in parallel per core
        quiet=True, print_progress=True,
        with_design_mtx = False, 
        from_null = False, refit_null = False, # only fit the full model
        random_seed = 0
    )
    glmm.save(f"{res_dir}/model/cbs_peak_glmm_{date}.pkl")

    ## (3) run the differential usage tests
    # load processed RBP visium data
    adata_visium_rbp = sc.read(f"{res_dir}/anndata/cbs_visium_rbp_{date}.h5ad")
    # make sure spots are in the same order
    covariates = adata_visium_rbp[adata_peak_sv.obs.index, adata_visium_rbp.var['is_visium_sve']].layers['log1p'].toarray()
    covariate_names = adata_visium_rbp.var.loc[adata_visium_rbp.var['is_visium_sve'], 'gene_symbols']
    design_mtx = torch.from_numpy(covariates).float()

    # load the fitted GLM and GLMM models
    glm = torch.load(f"{res_dir}/model/cbs_peak_glm_{date}.pkl")
    glmm = torch.load(f"{res_dir}/model/cbs_peak_glmm_{date}.pkl")

    # set random seed
    torch.manual_seed(0)

    # GLM-score test
    glm.setup_data(
        data = counts_list, 
        coordinates = coords, 
        design_mtx = design_mtx,
        gene_names = gene_name_list, 
        covariate_names = covariate_names,
        group_gene_by_n_iso = False
    )
    glm.test_differential_usage(method = 'score', print_progress = True, return_results = False)
    df_du_glm = glm.get_formatted_test_results(test_type = 'du')

    # GLMM-score test
    glmm.setup_data(
        data = counts_list, 
        coordinates = coords, 
        design_mtx = design_mtx,
        gene_names = gene_name_list, 
        covariate_names = covariate_names,
        group_gene_by_n_iso = False
    )
    glmm.test_differential_usage(method = 'score', print_progress = True, return_results = False)
    df_du_glmm = glmm.get_formatted_test_results(test_type = 'du')

    # merge all test results
    df_rbp_pval = []
    for _test_method, _res in zip(['glm', 'glmm'], [df_du_glm, df_du_glmm]):
        _df = _res[['gene', 'covariate', 'pvalue']].copy()
        _df.rename(columns = {'pvalue': f'pvalue_{_test_method}'}, inplace=True)
        _df.set_index(['gene', 'covariate'], inplace=True)
        df_rbp_pval.append(_df)

    df_rbp_pval = df_rbp_pval[0].join(df_rbp_pval[1:]).reset_index()

    # save results
    df_rbp_pval.to_csv(f"{res_dir}/du_results/cbs_rbp_du_p_combined_{date}.csv", index=False)
