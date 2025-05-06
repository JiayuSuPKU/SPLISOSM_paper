### Run parametric differential usage tests on simulated data
import os
import torch
import pickle
import numpy as np
import pandas as pd
import itertools
from isosde.hyptest_glmm import SplisosmGLMM

# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/simulation_data/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/benchmark/"
date = '1205'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/model", f"{res_dir}/du_results"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    ### (1) run DU tests on the general six scenarios
    ## (1.1) fit the GLM and GLMM models
    for group_gene, group_iso in itertools.product(
        ['none', 'donut'],
        ['none', 'mvn', 'donut']
    ):
        print(f"=== Fitting GLM and GLMM models on sample G={group_gene} I={group_iso} ===")
        # load data and params
        with open(f"{data_dir}/general_six_scenarios/data_gene-{group_gene}_iso-{group_iso}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/general_six_scenarios/param_gene-{group_gene}_iso-{group_iso}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']

        # set random seed
        torch.manual_seed(0)

        ## fit the GLM and GLMM models without covariates (using score-test for DU)
        # set up and fit the GLM model
        glm = SplisosmGLMM(model_type = 'glm')
        glm.setup_data(
            data = counts_list,
            coordinates = coords,
            design_mtx = None,
            group_gene_by_n_iso = False, # whether to group genes by the number of isoforms for batch processing
        )
        glm.fit(
            n_jobs = 2, # number of cores to use
            batch_size = 1, # number of genes with the same number of isoforms to process in parallel per core
            quiet=True, print_progress=True,
            with_design_mtx = False, 
            random_seed = 0
        )
        glm.save(f"{res_dir}/model/glm_gene-{group_gene}_iso-{group_iso}_{date}.pkl")

        # set up and fit the GLMM model
        glmm = SplisosmGLMM(model_type = 'glmm-full')
        glmm.setup_data(
            data = counts_list,
            coordinates = coords,
            design_mtx = None,
            group_gene_by_n_iso = True, # whether to group genes by the number of isoforms for batch processing
        )
        glmm.fit(
            n_jobs = 2, # number of cores to use
            batch_size = 50, # number of genes with the same number of isoforms to process in parallel per core
            quiet=True, print_progress=True,
            with_design_mtx = False, 
            from_null = False,
            random_seed = 0
        )
        glmm.save(f"{res_dir}/model/glmm_gene-{group_gene}_iso-{group_iso}_{date}.pkl")

    ## (1.2) run the differential usage tests
    df_du_list = []
    for group_gene, group_iso in itertools.product(
        ['none', 'donut'],
        ['none', 'mvn', 'donut']
    ):
        print(f"=== Running DU tests on sample G={group_gene} I={group_iso} ===")
        # load data and params
        with open(f"{data_dir}/general_six_scenarios/data_gene-{group_gene}_iso-{group_iso}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/general_six_scenarios/param_gene-{group_gene}_iso-{group_iso}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']
        design_mtx = data['design_mtx']
        covariate_names = ['C1', 'C2', 'B1', 'B2']  

        # load the fitted GLM and GLMM models
        glm = torch.load(f"{res_dir}/model/glm_gene-{group_gene}_iso-{group_iso}_{date}.pkl")
        glmm = torch.load(f"{res_dir}/model/glmm_gene-{group_gene}_iso-{group_iso}_{date}.pkl")

        # set random seed
        torch.manual_seed(0)

        # GLM-score test
        glm.setup_data(
            data = counts_list, 
            coordinates = coords, 
            design_mtx = design_mtx,
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
            covariate_names = covariate_names,
            group_gene_by_n_iso = False
        )
        glmm.test_differential_usage(method = 'score', print_progress = True, return_results = False)
        df_du_glmm = glmm.get_formatted_test_results(test_type = 'du')

        # merge GLM and GLMM results
        df_du_pval_sample = []
        for _test_method, _res in zip(['glm', 'glmm'], [df_du_glm, df_du_glmm]):
            _df = _res[['gene', 'covariate', 'pvalue']].copy()
            # set covar_type by covariate name
            _df['covar_type'] = _df['covariate'].apply(lambda x: 'continuous' if x.startswith('C') else 'binary')
            _df.rename(columns = {'pvalue': f'pvalue_{_test_method}'}, inplace=True)
            _df.set_index(['gene', 'covariate', 'covar_type'], inplace=True)
            df_du_pval_sample.append(_df)

        df_du_pval_sample = df_du_pval_sample[0].join(df_du_pval_sample[1:]).reset_index()
        df_du_pval_sample['group_gene'] = group_gene
        df_du_pval_sample['group_iso'] = group_iso
        df_du_list.append(df_du_pval_sample)
        
    # save results
    df_du_pval = pd.concat(df_du_list)
    df_du_pval.to_csv(f"{res_dir}/du_results/du_p_general_six_{date}.csv", index=False)

    ### (2) Power analysis of DU tests on dataset with varying degrees of effect sizes
    ## (2.1) fit the GLM and GLMM models
    for beta_scale in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        print(f"=== Fitting GLM and GLMM models on sample G=donut I=donut (beta_scale={beta_scale}) ===")

        # load data and params
        with open(f"{data_dir}/du_donut_gene/data_beta-{beta_scale}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/du_donut_gene/param_beta-{beta_scale}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']
        design_mtx = data['design_mtx']
        covariate_names = ['C1', 'C2', 'B1', 'B2']

        # set random seed
        torch.manual_seed(0)

        ## fit the GLM and GLMM models with and without covariates (using score-test for DU)
        # set up and fit the GLM model
        for _covar in [False, True]:
            glm = SplisosmGLMM(model_type = 'glm')
            glm.setup_data(
                data = counts_list,
                coordinates = coords,
                design_mtx = design_mtx,
                covariate_names = covariate_names,
                group_gene_by_n_iso = False, # whether to group genes by the number of isoforms for batch processing
            )
            glm.fit(
                n_jobs = 2, # number of cores to use
                batch_size = 1, # number of genes with the same number of isoforms to process in parallel per core
                quiet=True, print_progress=True,
                with_design_mtx = _covar, 
                random_seed = 0
            )
            glm.save(f"{res_dir}/model/glm_covar-{_covar}_beta-{beta_scale}_{date}.pkl")
            
        # set up and fit the GLMM model
        for _covar in [False, True]:
            glmm = SplisosmGLMM(model_type = 'glmm-full')
            glmm.setup_data(
                data = counts_list,
                coordinates = coords,
                design_mtx = design_mtx,
                covariate_names = covariate_names,
                group_gene_by_n_iso = True, # whether to group genes by the number of isoforms for batch processing
            )
            glmm.fit(
                n_jobs = 2, # number of cores to use
                batch_size = 50, # number of genes with the same number of isoforms to process in parallel per core
                quiet=True, print_progress=True,
                with_design_mtx = _covar, 
                from_null = False,
                random_seed = 0
            )
            glmm.save(f"{res_dir}/model/glmm_covar-{_covar}_beta-{beta_scale}_{date}.pkl")

    ## (2.2) run the differential usage tests
    df_du_power_list = []
    for beta_scale in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        print(f"=== Fitting GLM and GLMM models on sample G=donut I=donut (beta_scale={beta_scale}) ===")

        # load data and params
        with open(f"{data_dir}/du_donut_gene/data_beta-{beta_scale}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/du_donut_gene/param_beta-{beta_scale}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']
        design_mtx = data['design_mtx']
        covariate_names = ['C1', 'C2', 'B1', 'B2']  

        # load the fitted GLM and GLMM models
        glm = torch.load(f"{res_dir}/model/glm_covar-False_beta-{beta_scale}_{date}.pkl")
        glmm = torch.load(f"{res_dir}/model/glmm_covar-False_beta-{beta_scale}_{date}.pkl")
        glm_w_covar = torch.load(f"{res_dir}/model/glm_covar-True_beta-{beta_scale}_{date}.pkl")
        glmm_w_covar = torch.load(f"{res_dir}/model/glmm_covar-True_beta-{beta_scale}_{date}.pkl")

        # set random seed
        torch.manual_seed(0)

        df_du_power_res = {}

        # GLM-score test
        glm.setup_data(
            data = counts_list, 
            coordinates = coords, 
            design_mtx = design_mtx,
            covariate_names = covariate_names,
            group_gene_by_n_iso = False
        )
        glm.test_differential_usage(method = 'score', print_progress = True, return_results = False)
        df_du_power_res['glm-score'] = glm.get_formatted_test_results(test_type = 'du')
    
        # GLMM-score test
        glmm.setup_data(
            data = counts_list, 
            coordinates = coords, 
            design_mtx = design_mtx,
            covariate_names = covariate_names,
            group_gene_by_n_iso = False
        )
        glmm.test_differential_usage(method = 'score', print_progress = True, return_results = False)
        df_du_power_res['glmm-score'] = glmm.get_formatted_test_results(test_type = 'du')

        # GLM-Wald test
        glm_w_covar.test_differential_usage(method = 'wald', print_progress = True, return_results = False)
        df_du_power_res['glm-wald'] = glm_w_covar.get_formatted_test_results(test_type = 'du')

        # GLMM-Wald test
        glmm_w_covar.test_differential_usage(method = 'wald', print_progress = True, return_results = False)
        df_du_power_res['glmm-wald'] = glmm_w_covar.get_formatted_test_results(test_type = 'du')

        # merge all results
        df_du_pval_sample = []
        for _test_method, _res in df_du_power_res.items():
            _df = _res[['gene', 'covariate', 'pvalue']].copy()
            _df['covar_type'] = _df['covariate'].apply(lambda x: 'continuous' if x.startswith('C') else 'binary')
            _df.rename(columns = {'pvalue': f'pvalue_{_test_method}'}, inplace=True)
            _df.set_index(['gene', 'covariate', 'covar_type'], inplace=True)
            df_du_pval_sample.append(_df)

        df_du_pval_sample = df_du_pval_sample[0].join(df_du_pval_sample[1:]).reset_index()
        df_du_pval_sample['beta_scale'] = beta_scale
        df_du_power_list.append(df_du_pval_sample)
        
    # save results
    df_du_pval = pd.concat(df_du_power_list)
    df_du_pval.to_csv(f"{res_dir}/du_results/du_p_power_{date}.csv", index=False)
