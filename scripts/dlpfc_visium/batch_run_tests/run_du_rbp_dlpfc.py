### Run differential usage tests against RBP expression for all DLPFC samples
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import get_cov_sp, load_visium_sp_meta, false_discovery_control, extract_counts_n_ratios
from isosde.likelihood import liu_sf
from isosde.hyptest_glmm import SplisosmGLMM
from isosde.hyptest_np import SplisosmNP


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
    for _dir in [res_dir, f"{res_dir}/du_rbp_results"]:
        os.makedirs(_dir, exist_ok=True)

    for sample_id in sample_list:
        ## load the filtered peak-level anndata and SV test results
        adata_peak = sc.read(f"{res_dir}/anndata/sierra/{sample_id}.peak_filtered_{date}.h5ad")
        coords = adata_peak.obs.loc[:, ['array_row', 'array_col']]
        df_sv_pval = pd.read_csv(f"{res_dir}/sv_results/{sample_id}.sv_combined_{date}.csv", index_col=0)

        # focus on SVS genes only
        _is_sig_svs = (df_sv_pval['padj_hsic-ir'] < 0.05) & (df_sv_pval['count_avg'] > 0.5)
        _iso_ind = adata_peak.var['gene_symbol'].isin(df_sv_pval.index[_is_sig_svs])
        adata_peak_sv = adata_peak[:, _iso_ind].copy()
        
        # load adata of RBP genes
        adata_visium_rbp = sc.read_h5ad(f"{res_dir}/anndata/visium/{sample_id}.rbp_{date}.h5ad")

        # extract lists of isoform counts and ratios
        counts_list, _, gene_name_list, _ = extract_counts_n_ratios(
            adata_peak_sv, layer = 'counts', group_iso_by = 'gene_symbol'
        )

        # extract RBP expression data
        covariates = adata_visium_rbp[:, adata_visium_rbp.var['is_visium_sve']].layers['log1p'].toarray()
        covariate_names = adata_visium_rbp.var_names[adata_visium_rbp.var['is_visium_sve']]
        design_mtx = torch.from_numpy(covariates).float()

        print(f"=== Running differential usage tests for {sample_id} ...")
        print(f"Number of SVS genes: {len(gene_name_list)}")
        print(f"Number of SVE RBP genes: {design_mtx.shape[1]}")

        # run the DU tests for RBP expression
        df_du_res = {}

        ## parametric testings (GLM and GLMM)
        # load the fitted GLM and GLMM models
        glm = torch.load(f"{res_dir}/model/{sample_id}.glm_{date}.pkl")
        glmm = torch.load(f"{res_dir}/model/{sample_id}.glmm_{date}.pkl")

        for method, model_p in zip(['glm', 'glmm'], [glm, glmm]):
            torch.manual_seed(0) # for reproducibility
            model_p.setup_data(
                data = counts_list, 
                coordinates = coords, 
                design_mtx = design_mtx,
                gene_names = gene_name_list, 
                covariate_names = covariate_names,
                group_gene_by_n_iso = False
            )
            model_p.test_differential_usage(method = 'score', print_progress = True, return_results = False)
            df_du_res[method] = model_p.get_formatted_test_results(test_type = 'du')

        # non-parametric testings
        model_np = SplisosmNP()
        model_np.setup_data(
            counts_list, coords, design_mtx = design_mtx, 
            gene_names = gene_name_list, covariate_names = covariate_names
        )

        for _test_method in ['hsic', 'hsic-gp']:
            torch.manual_seed(0) # for reproducibility
            model_np.test_differential_usage(
                method = _test_method,
                ratio_transformation = 'none', nan_filling = 'mean', # same as above
                hsic_eps = None, # regularization parameter kernel regression, only applicable to 'hsic'. If set to None, will be the unconditional HSIC test.
                gp_configs = None, # dictionary of configs for the Gaussian process regression, only applicable to 'hsic-gp'
                print_progress = True, return_results = False
            )
            df_du_res[_test_method] = model_np.get_formatted_test_results(test_type = 'du') # per gene-factor pair test statistics

        # merge all test results
        df_rbp_pval = []
        for _test_method, _res in df_du_res.items():
            _df = _res[['gene', 'covariate', 'pvalue']].copy()
            _df.rename(columns = {'pvalue': f'pvalue_{_test_method}'}, inplace=True)
            _df.set_index(['gene', 'covariate'], inplace=True)
            df_rbp_pval.append(_df)

        df_rbp_pval = df_rbp_pval[0].join(df_rbp_pval[1:]).reset_index()

        # save results
        df_rbp_pval.to_csv(f"{res_dir}/du_rbp_results/{sample_id}.du_combined_{date}.csv", index=False)


