# Test for spatially variable expression in RBP genes, 12 DLPFC samples
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
    # df_rbp = pd.read_table("/Users/jysumac/Projects/SPLISOSM_paper/data/Homo_sapiens.RBP.txt", header=None, sep='\t')
    # focus on RBPs with CISBP-RNA motif or POSTAR CLIP
    rbp_with_motif = pd.read_table(
        "~/reference/cisbp-rna/human_pm/rbp_with_known_cleaned_motif.txt", header = None
    )[0].to_list()
    rbp_with_clip = pd.read_table(
        '/Users/jysumac/reference/POSTAR3/human.rbp_with_clip.txt', header=None
    )[0].to_list()

    # add additional RBPs to test
    rbp_custom = [f'CELF{i+1}' for i in range(6)] + [f'RBFOX{i+1}' for i in range(3)] + \
        ['QKI', 'ADAR1', 'ADAR2', 'ADAR3', 'ARPP21', 'SON']

    # merge the list
    df_rbp = pd.DataFrame({
        'name': list(set(rbp_with_motif + rbp_with_clip + rbp_custom))
    })
    df_rbp['has_motif'] = df_rbp['name'].isin(rbp_with_motif)
    df_rbp['has_clip'] = df_rbp['name'].isin(rbp_with_clip)

    # sort dataframe by name
    df_rbp = df_rbp.sort_values('name')

    return df_rbp

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

    # loop over all samples
    for sample_id in sample_list:
        # load visium expression data
        adata_visium = sc.read_h5ad(f"{res_dir}/anndata/visium/{sample_id}.visium_{date}.h5ad")

        # load the RBP list and keep only non-ribosomal canonical RBPs
        df_rbp = _load_rbp_list()

        # subset to RBP genes
        adata_visium.var['is_rbp'] = adata_visium.var_names.isin(df_rbp['name'])
        adata_visium_rbp = adata_visium[:, adata_visium.var['is_rbp']].copy()
        adata_visium_rbp.var['has_motif'] = adata_visium_rbp.var_names.isin(df_rbp.query('has_motif')['name'])
        adata_visium_rbp.var['has_clip'] = adata_visium_rbp.var_names.isin(df_rbp.query('has_clip')['name'])

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
        adata_visium_rbp.var['is_visium_sve'] = (adata_visium_rbp.var['pvalue_adj_hsic'] < 0.01) & \
            (adata_visium_rbp.var['pvalue_adj_sparkx'] < 0.01)

        # save the anndata
        adata_visium_rbp.write_h5ad(f"{res_dir}/anndata/visium/{sample_id}.rbp_{date}.h5ad")

