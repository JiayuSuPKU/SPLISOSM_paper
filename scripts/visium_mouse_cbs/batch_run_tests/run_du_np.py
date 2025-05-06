### Run the non-parametric differential usage tests on the mouse CBS visium sample
import os
import torch
import numpy as np
import pandas as pd
import itertools
import scanpy as sc
from isosde.utils import get_cov_sp, load_visium_sp_meta, false_discovery_control, extract_counts_n_ratios
from isosde.likelihood import liu_sf
from isosde.hyptest_np import SplisosmNP

def _load_rbp_list():
    # load the RBP list (downloaded from http://eurbpdb.gzsys.org.cn/index.php)
    df_rbp = pd.read_table("/Users/jysumac/Projects/SPLISOSM_paper/data/Mus_musculus.RBP.txt", header=None, sep='\t')
    df_rbp.columns = ['ensembl_id', 'gene_symbol', 'entrez_id', 'alias', 'full_name', 'family']
    return df_rbp

def _prepare_visium_rbp_anndata(data_dir):
    # load the count data
    adata_visium = sc.read_10x_h5(f"{data_dir}/Visium_Fresh_Frozen_Adult_Mouse_Brain_filtered_feature_bc_matrix.h5") # short-read gene counts
    adata_visium.raw = adata_visium.copy() # raw counts data
    adata_visium.var['gene_symbols'] = adata_visium.var_names
    adata_visium.var_names_make_unique()

    # filter out barely expressed genes
    sc.pp.filter_genes(adata_visium, min_cells=0.01 * adata_visium.shape[0])
    adata_visium.layers['counts'] = adata_visium.X.copy()

    # log normalize
    sc.pp.normalize_total(adata_visium, target_sum=1e4)
    sc.pp.log1p(adata_visium)
    adata_visium.layers['log1p'] = adata_visium.X.copy()

    # load the Visium spatial metadata 
    adata_visium = load_visium_sp_meta(adata_visium, f"{data_dir}/visium_spatial_meta", library_id='adata_visium')

    # load the RBP list and keep only non-ribosomal canonical RBPs
    df_rbp = _load_rbp_list()
    _ind_rbp_remove = df_rbp['family'].str.contains('ribosomal', case=False) | \
        df_rbp['full_name'].str.contains('ribosomal', case=False) | \
        (df_rbp['family'] == 'Non-canonical_RBP')

    # subset to RBP genes
    adata_visium.var['is_rbp'] = adata_visium.var['gene_symbols'].isin(
        df_rbp.loc[~ _ind_rbp_remove, 'gene_symbol']
    )
    adata_visium_rbp = adata_visium[:, adata_visium.var['is_rbp']].copy()

    return adata_visium_rbp

def _run_hsic_visium(adata_visium):
    # calculate spatial kernel
    coords = adata_visium.obs.loc[:, ['array_row', 'array_col']] # n_spots x 2
    corr_sp = get_cov_sp(coords, k = 4, rho=0.99) # n_spots x n_spots

    # eigenvalues of the spatial kernel
    n_spots = coords.shape[0]
    H = torch.eye(n_spots) - 1/n_spots
    K_sp = H @ corr_sp @ H # centered spatial kernel
    lambda_sp = torch.linalg.eigvalsh(K_sp) # eigenvalues of length n_spots
    lambda_sp = lambda_sp[lambda_sp > 1e-5] # remove small eigenvalues

    # extract gene counts
    counts = adata_visium.layers['counts'].toarray() # n_spots x n_genes
    y = counts - counts.mean(axis = 0, keepdims = True) # centering per gene
    y = torch.from_numpy(y).T.unsqueeze(1).float() # n_genes x 1 x n_spots
    lambda_y = (y @ y.transpose(1, 2)).squeeze() # n_genes

    # calculate HSIC
    hsic_scaled = (y @ K_sp @ y.transpose(1, 2)).squeeze() # n_genes
    hsic = hsic_scaled / (n_spots - 1) ** 2 # n_genes

    # calculate p-values
    pval_list = []
    for _hsic_scaled, _lambda_y in zip(hsic_scaled, lambda_y):
        pval_list.append(liu_sf((_hsic_scaled * n_spots).numpy(), (lambda_sp * _lambda_y).numpy()))

    # save SV test results
    df_sv_genes = pd.DataFrame({
        'statistic': hsic.numpy(),
        'pvalue': pval_list,
        'pvalue_adj': false_discovery_control(np.array(pval_list)),
        'method': 'hsic-gc'
    })
    df_sv_genes.index = adata_visium.var.index

    return df_sv_genes

def _run_sparkx_visium(adata_visium):
    # for running SPARK-X
    import rpy2
    import rpy2.robjects as ro
    from rpy2.robjects import r
    from rpy2.robjects import numpy2ri
    from rpy2.robjects import pandas2ri
    from rpy2.robjects.packages import importr
    numpy2ri.activate()
    pandas2ri.activate()
    spark = importr('SPARK')

    # prepare robject inputs
    coords_r = ro.conversion.py2rpy(
        adata_visium.obs.loc[:, ['array_row', 'array_col']]
    ) # n_spots x 2

    counts_r = ro.conversion.py2rpy(adata_visium.layers['counts'].toarray().T) # n_genes x n_spots
    counts_r.colnames = ro.vectors.StrVector(r['rownames'](coords_r))

    # run SPARK-X and extract outputs
    sparkx_res = spark.sparkx(counts_r, coords_r)
    df_sv_genes = pd.DataFrame({
        'statistic': ro.conversion.rpy2py(sparkx_res.rx['stats'][0]).mean(1),
        'pvalue': ro.conversion.rpy2py(sparkx_res.rx['res_mtest'][0])['combinedPval'].values,
        'pvalue_adj': ro.conversion.rpy2py(sparkx_res.rx['res_mtest'][0])['adjustedPval'].values,
        'method': 'spark-x'
    })
    df_sv_genes.index = adata_visium.var.index

    return df_sv_genes


# set up the data and results directories
data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/visium_mouse_cbs/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/visium_mouse_cbs/"
date = '1119'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/anndata", f"{res_dir}/du_results"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    ## (0) prepare RBP expression data and run spatial variability tests
    # load visium expression data for RBP genes
    adata_visium_rbp = _prepare_visium_rbp_anndata(data_dir)

    # run spatial variability tests on RBP expression
    df_sve_hsic = _run_hsic_visium(adata_visium_rbp)
    df_sve_sparkx = _run_sparkx_visium(adata_visium_rbp)

    # compare the results
    print(f"Number of significant genes (HSIC): {(df_sve_hsic['pvalue_adj'] <= 0.01).sum()}")
    print(f"Number of significant genes (SPARK-X): {(df_sve_sparkx['pvalue_adj'] <= 0.01).sum()}")
    print(
        f"Number of overlapping significant genes: {
        df_sve_hsic[df_sve_hsic['pvalue_adj'] <= 0.01].index.intersection(
            df_sve_sparkx[df_sve_sparkx['pvalue_adj'] <= 0.01].index
        ).shape[0]}"
    )

    # save SV results
    adata_visium_rbp.var['pvalue_adj_hsic'] = df_sve_hsic['pvalue_adj']
    adata_visium_rbp.var['pvalue_adj_sparkx'] = df_sve_sparkx['pvalue_adj']
    adata_visium_rbp.var['is_visium_sve'] = adata_visium_rbp.var['pvalue_adj_sparkx'] < 0.01

    ## (1) load the filtered peak counts anndata and SV test results
    adata_peak = sc.read(f"{res_dir}/anndata/cbs_peak_filtered_{date}.h5ad")
    df_sv_pval = pd.read_csv(f"{res_dir}/sv_results/cbs_peak_sv_combined_{date}.csv", index_col=0)

    # reorder rbp adata spot order to match peak adata
    adata_visium_rbp = adata_visium_rbp[adata_peak.obs.index, :].copy()
    adata_visium_rbp.write_h5ad(f"{res_dir}/anndata/cbs_visium_rbp_{date}.h5ad")

    # focus on SVS genes only
    _is_sig_svs = (df_sv_pval['padj_hsic-ir'] < 0.01) & (df_sv_pval['count_avg'] > 0.5) & \
        (~ df_sv_pval.index.str.startswith('Gm')) & \
        (~ df_sv_pval.index.str.endswith('Rik')) # exclude unannotated genes
    _iso_ind = adata_peak.var['gene_symbol'].isin(df_sv_pval.index[_is_sig_svs])
    adata_peak_sv = adata_peak[:, _iso_ind].copy()

    # extract lists of isoform counts and ratios
    counts_list, _, gene_name_list, _ = extract_counts_n_ratios(adata_peak_sv, layer = 'counts', group_iso_by = 'gene_symbol')
    coords = adata_peak_sv.obs.loc[:, ['array_row', 'array_col']]

    ## (2) run the DU tests for RBP expression
    # prepare covariates inputs
    covariates = adata_visium_rbp[adata_peak_sv.obs.index, adata_visium_rbp.var['is_visium_sve']].layers['log1p'].toarray()
    covariate_names = adata_visium_rbp.var.loc[adata_visium_rbp.var['is_visium_sve'], 'gene_symbols']
    design_mtx = torch.from_numpy(covariates).float()

    print(f"Testing {adata_peak_sv.var['gene_symbol'].nunique()} SVS genes with {design_mtx.shape[1]} RBPs")
    print(f"Average number of isoforms per gene: {adata_peak_sv.shape[1] / adata_peak_sv.var['gene_symbol'].nunique()}")

    # set up the model
    model_np = SplisosmNP()
    model_np.setup_data(
        counts_list, coords, design_mtx = design_mtx, 
        gene_names = gene_name_list, covariate_names = covariate_names
    )

    # run all DU tests (unconditional HSIC or HSIC-GP)
    df_du_res = {}
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
    df_rbp_pval.to_csv(f"{res_dir}/du_results/cbs_rbp_du_np_combined_{date}.csv", index=False)


