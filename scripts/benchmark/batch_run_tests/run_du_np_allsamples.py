### Run non-parametric differential usage tests on simulated data
import os
import torch
import pickle
import numpy as np
import pandas as pd
import itertools
from isosde.hyptest_np import SplisosmNP

# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/simulation_data/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/benchmark/"
date = '1205'

if __name__ == '__main__':
    # create the results directory if it does not exist
    if not os.path.exists(f"{res_dir}/du_results"):
        os.makedirs(f"{res_dir}/du_results")

    ### (1) run DU tests on the general six scenarios
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

        # run DU tests for continuous covariates
        model_np_cont = SplisosmNP()
        model_np_cont.setup_data(
            data = counts_list, coordinates = coords,
            design_mtx = design_mtx[:, :2], covariate_names = ['C1', 'C2']
        )
        df_du_cont_res = {}
        for _test_method in ['hsic', 'hsic-gp']:
            model_np_cont.test_differential_usage(
                method = _test_method,
                ratio_transformation = 'none',
                nan_filling = 'mean',
                hsic_eps = None,
                return_results = False,
                print_progress = True
            )
            df_du_cont_res[_test_method] = model_np_cont.get_formatted_test_results(test_type = 'du')

        # run DU tests for binary covariates
        model_np_bi = SplisosmNP()
        model_np_bi.setup_data(
            data = counts_list, coordinates = coords,
            design_mtx = design_mtx[:, 2:], covariate_names = ['B1', 'B2']
        )
        df_du_bi_res = {}
        for _test_method in ['hsic', 'hsic-gp', 't-fisher']:
            model_np_bi.test_differential_usage(
                method = _test_method,
                ratio_transformation = 'none',
                nan_filling = 'mean',
                hsic_eps = None,
                return_results = False,
                print_progress = True
            )
            df_du_bi_res[_test_method] = model_np_bi.get_formatted_test_results(test_type = 'du')

        # merge DU test results
        df_du_pval_sample = []
        for dict_res, covar_type in zip([df_du_cont_res, df_du_bi_res], ['continuous', 'binary']):
            # merge results within covariate type first
            _df_du_pval = []
            for _test_method, _res in dict_res.items():
                _df = _res[['gene', 'covariate', 'pvalue']].copy()
                _df['covar_type'] = covar_type
                _df.rename(columns = {'pvalue': f'pvalue_{_test_method}'}, inplace=True)
                _df.set_index(['gene', 'covariate', 'covar_type'], inplace=True)
                _df_du_pval.append(_df)
            _df_du_pval = _df_du_pval[0].join(_df_du_pval[1:]).reset_index()
            df_du_pval_sample.append(_df_du_pval)

        # merge continuous and binary covariate results
        df_du_pval_sample = pd.concat(df_du_pval_sample)
        df_du_pval_sample['group_gene'] = group_gene
        df_du_pval_sample['group_iso'] = group_iso
        df_du_list.append(df_du_pval_sample)

    # save results
    df_du_pval = pd.concat(df_du_list)
    df_du_pval.to_csv(f"{res_dir}/du_results/du_np_general_six_{date}.csv", index=False)

    ### (2) Power analysis of DU tests on dataset with varying degrees of effect sizes
    df_du_power_list = []
    for beta_scale in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        print(f"=== Running DU tests on sample G=donut I=donut (beta_scale={beta_scale}) ===")

        # load data and params
        with open(f"{data_dir}/du_donut_gene/data_beta-{beta_scale}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/du_donut_gene/param_beta-{beta_scale}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']
        design_mtx = data['design_mtx']

        # run DU tests for continuous covariates
        model_np_cont = SplisosmNP()
        model_np_cont.setup_data(
            data = counts_list, coordinates = coords,
            design_mtx = design_mtx[:, :2], covariate_names = ['C1', 'C2']
        )
        df_du_cont_res = {}
        for _test_method in ['hsic', 'hsic-knn', 'hsic-eps']:
            model_np_cont.test_differential_usage(
                method = _test_method if _test_method != 'hsic-eps' else 'hsic',
                ratio_transformation = 'none',
                nan_filling = 'mean',
                hsic_eps = None if _test_method != 'hsic-eps' else 1e-3,
                return_results = False,
                print_progress = True
            )
            df_du_cont_res[_test_method] = model_np_cont.get_formatted_test_results(test_type = 'du')

        # run DU tests for binary covariates
        model_np_bi = SplisosmNP()
        model_np_bi.setup_data(
            data = counts_list, coordinates = coords,
            design_mtx = design_mtx[:, 2:], covariate_names = ['B1', 'B2']
        )
        df_du_bi_res = {}
        for _test_method in ['hsic', 'hsic-knn', 'hsic-eps', 't-fisher', 't-tippett']:
            model_np_bi.test_differential_usage(
                method = _test_method if _test_method != 'hsic-eps' else 'hsic',
                ratio_transformation = 'none',
                nan_filling = 'mean',
                hsic_eps = None if _test_method != 'hsic-eps' else 1e-3,
                return_results = False,
                print_progress = True
            )
            df_du_bi_res[_test_method] = model_np_bi.get_formatted_test_results(test_type = 'du')

        # merge DU test results
        df_du_pval_sample = []
        for dict_res, covar_type in zip([df_du_cont_res, df_du_bi_res], ['continuous', 'binary']):
            # merge results within covariate type first
            _df_du_pval = []
            for _test_method, _res in dict_res.items():
                _df = _res[['gene', 'covariate', 'pvalue']].copy()
                _df['covar_type'] = covar_type
                _df.rename(columns = {'pvalue': f'pvalue_{_test_method}'}, inplace=True)
                _df.set_index(['gene', 'covariate', 'covar_type'], inplace=True)
                _df_du_pval.append(_df)
            _df_du_pval = _df_du_pval[0].join(_df_du_pval[1:]).reset_index()
            df_du_pval_sample.append(_df_du_pval)

        # merge continuous and binary covariate results
        df_du_pval_sample = pd.concat(df_du_pval_sample)
        df_du_pval_sample['beta_scale'] = beta_scale

        # add hsic-gp results
        model_np = SplisosmNP()
        model_np.setup_data(
            data = counts_list, coordinates = coords,
            design_mtx = design_mtx, covariate_names = ['C1', 'C2', 'B1', 'B2']
        )
        model_np.test_differential_usage(
            method = 'hsic-gp',
            ratio_transformation = 'none',
            nan_filling = 'mean',
            hsic_eps = None,
            return_results = False,
            print_progress = True
        )
        df_du_gp = model_np.get_formatted_test_results(test_type = 'du')
        df_du_gp = df_du_gp[['gene', 'covariate', 'pvalue']].rename(columns = {'pvalue': 'pvalue_hsic-gp'})
        df_du_gp['covar_type'] = df_du_gp['covariate'].apply(lambda x: 'continuous' if x.startswith('C') else 'binary')

        # add as a new column to other test results
        df_du_pval_sample = df_du_pval_sample.merge(df_du_gp, on = ['gene', 'covariate', 'covar_type'], how = 'left')

        df_du_power_list.append(df_du_pval_sample)

    # save results
    df_du_power = pd.concat(df_du_power_list)
    df_du_power.to_csv(f"{res_dir}/du_results/du_np_power_{date}.csv", index=False)

