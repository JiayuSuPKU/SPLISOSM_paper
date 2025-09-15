### Run spatial variability tests on the Xenium 5K mouse brain CBS sample
import os
import torch
import numpy as np
import pandas as pd
import itertools
import pickle
from splisosm.utils import load_visium_sp_meta, extract_counts_n_ratios, extract_gene_level_statistics
from splisosm.hyptest_np import SplisosmNP

def _prepare_input_data(
    pickle_file_path: str,
):
    # Load the dictionary of density maps from the pickle file
    with open(pickle_file_path, "rb") as f:
        cw_density_dict = pickle.load(f)

    # Process the data to map genes to their corresponding codewords
    gene_cwid_dict = {}
    if cw_density_dict:
        codeword_lists = list(cw_density_dict.keys())
        for codeword in codeword_lists:
            # Assumes format 'GeneName-CodewordID'
            gene = codeword.split('-')[0]
            if gene not in gene_cwid_dict:
                gene_cwid_dict[gene] = []
            gene_cwid_dict[gene].append(codeword)

    # Filter for genes that have more than one codeword, as these are the
    # only ones for which a log ratio plot is meaningful.
    genes_with_multiple_codewords = sorted([
        gene for gene, codewords in gene_cwid_dict.items() if len(codewords) > 1
    ])

    print(f"Found {len(genes_with_multiple_codewords)} genes with multiple codewords available for comparison.")
    
    # Convert the per-codeword dictionary into a per-gene list
    gene_matrices_list = []
    ordered_gene_names = [] # This list will hold the gene names in the same order as the matrices.

    # Compute the spot coordinates
    # Get the shape of the density maps from the first entry in the dictionary
    # We assume all density maps have the same dimensions.
    first_codeword = next(iter(cw_density_dict))
    map_height, map_width = cw_density_dict[first_codeword].shape
    num_spots = map_height * map_width

    print(f"\nDensity map dimensions: {map_height} (height) x {map_width} (width)")
    print(f"Total number of spatial bins (spots): {num_spots}")

    # Create a meshgrid to represent the X and Y coordinates of each bin
    # np.arange(map_width) -> [0, 1, ..., width-1]
    # np.arange(map_height) -> [0, 1, ..., height-1]
    x_coords, y_coords = np.meshgrid(np.arange(map_width), np.arange(map_height))

    # Flatten the coordinate grids and stack them into a single (N_spots, 2) matrix
    # The .flatten() method unravels the arrays in 'C' (row-major) order,
    # which is consistent with how we will flatten the density maps.
    spot_coordinates = np.c_[x_coords.flatten(), y_coords.flatten()]

    print(f"\nGenerated spot coordinate matrix with shape: {spot_coordinates.shape}")

    # Compute the per-gene spot-by-codeword density matrix
    print("\nProcessing genes to create density matrices...")

    for gene in genes_with_multiple_codewords:
        ordered_gene_names.append(gene)
        
        # Get the list of codewords for the current gene
        codewords_for_gene = gene_cwid_dict[gene]
        
        # This list will temporarily hold the flattened density vector for each codeword
        codeword_density_vectors = []
        
        for codeword in codewords_for_gene:
            # Retrieve the sparse matrix, convert to dense, and flatten to a 1D vector
            sparse_map = cw_density_dict[codeword]
            dense_vector = sparse_map.toarray().flatten()
            codeword_density_vectors.append(dense_vector)
            
        # Stack the list of 1D vectors into a single 2D matrix.
        # Each vector becomes a column in the final matrix.
        # Shape will be (#spots, #codewords_for_this_gene)
        gene_matrix = np.stack(codeword_density_vectors, axis=1)
        gene_matrices_list.append(gene_matrix)
        
        # Print progress for a few genes
        if len(ordered_gene_names) <= 3:
            print(f"  - Gene '{gene}': created matrix of shape {gene_matrix.shape}")

    print(f"\nCompleted. Generated {len(gene_matrices_list)} per-gene spot-by-codeword count matrices.")

    return gene_matrices_list, ordered_gene_names, spot_coordinates

data_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/xenium_5k_mouse_brain/"
res_dir = "/Users/jysumac/Projects/SPLISOSM_paper/results/xenium_5k_mouse_brain/"
date = '0808'

if __name__ == '__main__':
    # create the results directory if it does not exist
    for _dir in [res_dir, f"{res_dir}/sv_results", f"{res_dir}/figures"]:
        if not os.path.exists(_dir):
            os.makedirs(_dir)

    # load Xenium 5K probe gene annotations
    df_gene_meta = pd.read_csv(
        f"{data_dir}/XeniumPrimeMouse5Kpan_tissue_pathways_metadata.csv"
    ).set_index('gene_name')

    # extract lists of isoform counts and ratios
    density_map_path = f"{data_dir}/density_maps"
    pickle_file_path = f"{density_map_path}/all_multi_codeword_genes.pkl"
    counts_list, gene_name_list, coords = _prepare_input_data(pickle_file_path)

    # non-parametric testings
    model_np = SplisosmNP()
    model_np.setup_data(
        data = counts_list, 
        coordinates = coords, 
        gene_names = gene_name_list, 
        approx_rank=100
    )

    # run all SV tests
    df_sv_res = {}
    for _test_method in ['hsic-ir', 'hsic-ic', 'hsic-gc']:
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
        }, index=gene_name_list), how = 'right')
    df_sv_pval = df_sv_pval.sort_values('pvalue_hsic-ir', ascending=True)
    df_sv_pval.index.name = 'gene'

    # save results
    df_sv_pval.to_csv(f"{res_dir}/sv_results/all_multi_codeword_genes_{date}.csv", index=True)
