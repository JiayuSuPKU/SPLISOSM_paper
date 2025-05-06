### Run spatial variability tests on simulated data
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
    if not os.path.exists(f"{res_dir}/sv_results"):
        os.makedirs(f"{res_dir}/sv_results")

    ### (1) run SV tests on the general six scenarios
    df_sv_list = []
    for group_gene, group_iso in itertools.product(
        ['none', 'donut'],
        ['none', 'mvn', 'donut']
    ):
        print(f"=== Running SV tests on sample G={group_gene} I={group_iso} ===")
        # load data and params
        with open(f"{data_dir}/general_six_scenarios/data_gene-{group_gene}_iso-{group_iso}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/general_six_scenarios/param_gene-{group_gene}_iso-{group_iso}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']

        # non-parametric testings
        model_np = SplisosmNP()
        model_np.setup_data(data = counts_list, coordinates = coords)

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
        df_sv_pval = pd.DataFrame({
            'pvalue_hsic-ir': df_sv_res['hsic-ir']['pvalue'].values,
            'pvalue_hsic-ic': df_sv_res['hsic-ic']['pvalue'].values,
            'pvalue_hsic-gc': df_sv_res['hsic-gc']['pvalue'].values,
            'pvalue_spark-x': df_sv_res['spark-x']['pvalue'].values,
        })
        df_sv_pval['group_gene'] = group_gene
        df_sv_pval['group_iso'] = group_iso
        df_sv_list.append(df_sv_pval)

    # save results
    df_sv_pval = pd.concat(df_sv_list)
    df_sv_pval.to_csv(f"{res_dir}/sv_results/sv_general_six_{date}.csv", index=False)

    ### (2) Power analysis of HSIC-IR on dataset with varying degrees of spatial variability
    df_sv_power_list = []
    # iterate over simulated data
    for prop_var_sp in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
        print(f"=== Running SV tests on sample G=donut I=mvn (prop_var_sp={prop_var_sp}) ===")

        # load data and params
        with open(f"{data_dir}/sv_donut_gene/data_prop_var_sp-{prop_var_sp}.pkl", "rb") as f:
            data = pickle.load(f)

        with open(f"{data_dir}/sv_donut_gene/param_prop_var_sp-{prop_var_sp}.pkl", "rb") as f:
            params = pickle.load(f)

        n_genes = 1000
        counts_list = data['counts'][:n_genes]
        coords = data['coords']

        # non-parametric testings
        model_np = SplisosmNP()
        model_np.setup_data(data = counts_list, coordinates = coords)

        df_sv_res = {}
        for _trans, _na_fill in itertools.product(
            ['none', 'alr', 'radial'], ['mean', 'none']
        ):
            model_np.test_spatial_variability(
                method = 'hsic-ir',
                ratio_transformation = _trans,
                nan_filling = _na_fill,
                return_results = False,
                print_progress = True
            )
            _key = f"hsic-ir_{_trans}_{_na_fill}"
            df_sv_res[_key] = model_np.get_formatted_test_results(test_type = 'sv') # per gene test statistics

        # merge SV test results
        df_sv_pval = pd.DataFrame({
            _key: df_sv_res[_key]['pvalue'].values
            for _key in df_sv_res.keys()
        })
        df_sv_pval['prop_var_sp'] = prop_var_sp
        df_sv_power_list.append(df_sv_pval)

    # save results
    df_sv_power = pd.concat(df_sv_power_list)
    df_sv_power.to_csv(f"{res_dir}/sv_results/sv_power_{date}.csv", index=False)
