### Run the parametric differential usage tests on one mouse olfactory bulb (MOB) and two coronal brain section (CBS) samples
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import load_visium_sp_meta, extract_counts_n_ratios
from isosde.hyptest_glmm import SplisosmGLMM

# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/sit_nar_23/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/sit_nar_23/"
date = '1107'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/model", f"{res_dir}/du_p_results"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    for _sample_id in ['mob', 'cbs1', 'cbs2']:
        ## (0) load the filtered ONT anndata and SV test results
        adata_ont = sc.read(f"{res_dir}/anndata/{_sample_id}_ont_filtered_{date}.h5ad")
        coords = adata_ont.obs.loc[:, ['array_row', 'array_col']]
        df_sv_pval = pd.read_csv(f"{res_dir}/sv_results/{_sample_id}_sv_combined_{date}.csv", index_col=0)

        # focus on SVS genes only
        if _sample_id == 'mob':
            _is_sig_svs = (df_sv_pval['pvalue_hsic-ir'] < 0.05) & (df_sv_pval['count_avg'] > 0.5)
        else:
            _is_sig_svs = (df_sv_pval['padj_hsic-ir'] < 0.05) & (df_sv_pval['count_avg'] > 0.5)

        _iso_ind = adata_ont.var['gene_symbol'].isin(df_sv_pval.index[_is_sig_svs])
        adata_ont_sv = adata_ont[:, _iso_ind].copy()

        print(f"Fitting {adata_ont_sv.var['gene_symbol'].nunique()} genes for {_sample_id}")
        print(f"Average number of isoforms per gene: {adata_ont_sv.shape[1] / adata_ont_sv.var['gene_symbol'].nunique()}")

        # extract lists of isoform counts and ratios
        counts_list, _, gene_name_list, _ = extract_counts_n_ratios(adata_ont_sv, layer = 'counts', group_iso_by = 'gene_symbol')

        # set random seed
        torch.manual_seed(0)

        ## (1) fit the GLM and GLMM models
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
        glm.save(f"{res_dir}/model/{_sample_id}_glm_{date}.pkl")

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
        glmm.save(f"{res_dir}/model/{_sample_id}_glmm_{date}.pkl")

        ## (2) run the differential usage tests
        # load processed RBP visium data
        adata_visium_rbp = sc.read(f"{res_dir}/anndata/{_sample_id}_visium_rbp_{date}.h5ad")
        covariates = adata_visium_rbp[:, adata_visium_rbp.var['is_visium_sve']].layers['log1p'].toarray()
        covariate_names = adata_visium_rbp.var.loc[adata_visium_rbp.var['is_visium_sve'], 'features']
        design_mtx = torch.from_numpy(covariates).float()

        # # load the fitted GLM and GLMM models
        # glm = torch.load(f"{res_dir}/model/{_sample_id}_glm_{date}.pkl")
        # glmm = torch.load(f"{res_dir}/model/{_sample_id}_glmm_{date}.pkl")

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
        df_rbp_pval.to_csv(f"{res_dir}/du_p_results/{_sample_id}_rbp_du_p_combined_{date}.csv", index=False)
