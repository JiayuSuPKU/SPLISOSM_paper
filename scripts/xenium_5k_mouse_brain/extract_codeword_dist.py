import zarr
import numpy as np
import pandas as pd
import os
import shutil
import scipy
from scipy.sparse import save_npz
from tqdm import tqdm
import warnings
import pickle

def extract_density_per_codeword(
    zarr_path: str,
    target_codewords = None,
    target_genes = None,
    quality_threshold: float = 20.0,
    bin_size_um: float = 10.0
):
    """
    Parses a transcripts.zarr.zip file to extract transcript locations for
    specified codewords and calculates their 2D density.

    Args:
        zarr_path (str): Path to the transcripts.zarr.zip file.
        target_codewords (list[str]): A list of codeword names to extract.
        target_genes (list[str]): A list of gene names to extract (all codewords of these genes will be included).
        quality_threshold (float): The minimum Q-Score to include a transcript.
        bin_size_um (float): The size of the bins in microns for density calculation.

    Returns:
        dict: A dictionary mapping codeword names to their 2D density matrices.
    """
    if not os.path.exists(zarr_path):
        raise FileNotFoundError(f"The file was not found at the specified path: {zarr_path}")

    print(f"Opening Zarr archive: {zarr_path}")
    # Open the zipped Zarr archive in read-only mode
    root = zarr.open(zarr_path, mode='r')

    # --- Step 1: Sanity check for codewords ---
    codeword2gene = root.density.codeword.attrs.get('codeword_names', [])
    all_codeword_ids = [i for i, _ in enumerate(codeword2gene)]
    
    if target_codewords is None:
        # use gene names to find all codewords
        if target_genes is None:
            raise ValueError("Either target_codewords or target_genes must be specified.")

        # Find all codewords for the specified genes
        target_codewords = []
        for gene in target_genes:
            if gene not in codeword2gene:
                warnings.warn(f"Gene '{gene}' not found in the Zarr file. Skipping.")
                continue
                # raise ValueError(f"Gene '{gene}' not found in the Zarr file.")
            # Get all codewords for this gene
            target_codewords.extend([i for i, g in enumerate(codeword2gene) if g == gene])

        target_codewords = list(set(target_codewords))  # Remove duplicates
        if (len(target_codewords) == 0):
            raise ValueError(f"No codewords found for the specified gene set.")
    else:
        # Make sure all target codewords are present in the Zarr file
        if not all(name in all_codeword_ids for name in target_codewords):
            raise ValueError("One or more target codewords are not present in the Zarr file.")

    # --- Step 2: Iterate through chunks and extract transcript data ---
    locations_by_codeword = {name: [] for name in target_codewords}
    
    # The grid keys for level 0 list all chunks to be processed.
    level0_keys = root['grids'].attrs['grid_keys'][0]
    print(f"Processing {len(level0_keys)} chunks from grid level 0...")

    for key in tqdm(level0_keys):
        chunk_path = f'grids/0/{key}'
        if chunk_path not in root:
            continue
        chunk = root[chunk_path]

        # Load necessary data arrays from the chunk
        q_scores = chunk['quality_score'][:].flatten()
        codeword_ids = chunk['codeword_identity'][:, 0] # First column is the index [1]
        locations = chunk['location'][:]

        # Create a boolean mask for quality filtering
        quality_mask = q_scores >= quality_threshold

        # Find transcripts for each target codeword
        for name in target_codewords:
            codeword_mask = codeword_ids == name
            # Combine masks to find final indices
            final_mask = quality_mask & codeword_mask
            if np.any(final_mask):
                locations_by_codeword[name].append(locations[final_mask])

    # --- Step 3: Calculate density ---
    # Create a common grid for all codewords
    grid_size = root.density.gene.attrs['grid_size'] # [10, 10]
    grid_origin = root.density.gene.attrs['origin'] # {'x': 0.0, 'y': 0.0}
    nrows = root.density.gene.attrs['rows']
    ncols = root.density.gene.attrs['cols']
    x_bins = np.linspace(
        grid_origin['x'], 
        grid_origin['x'] + ncols * grid_size[0], 
        ncols + 1
    )
    y_bins = np.linspace(
        grid_origin['y'], 
        grid_origin['y'] + nrows * grid_size[1], 
        nrows + 1
    )    

    # Adjust resize the grid if the bin size is different from the grid size
    if bin_size_um != grid_size[0]:
        # Calculate new grid size based on the bin size
        new_ncols = int((x_bins[-1] - x_bins[0]) / bin_size_um)
        new_nrows = int((y_bins[-1] - y_bins[0]) / bin_size_um)
        x_bins = np.linspace(x_bins[0], x_bins[-1], new_ncols + 1)
        y_bins = np.linspace(y_bins[0], y_bins[-1], new_nrows + 1)
        ncols, nrows = new_ncols, new_nrows
    
    print(f"Grid size: {bin_size_um} µm, Origin: {grid_origin}, Rows: {nrows}, Cols: {ncols}")

    # Initialize a dictionary to hold the sparse density maps for each codeword
    density_maps = {}
    for name, loc_list in tqdm(locations_by_codeword.items()):
        gene_name = codeword2gene[name]
        key = f"{gene_name}-{name}"

        if loc_list:
            locs = np.vstack(loc_list)  # Combine all locations for this codeword
            # Calculate 2D histogram for this specific codeword
            hist, _, _ = np.histogram2d(
                locs[:, 0], locs[:, 1], bins=[x_bins, y_bins]
            )
            # Convert the dense histogram to a sparse matrix (CSC format is efficient for column slicing)
            # Transpose for intuitive plotting (X on columns, Y on rows)
            sparse_hist = scipy.sparse.csc_matrix(hist.T)
            density_maps[key] = sparse_hist
            print(f"Created sparse density map for '{key}' with {sparse_hist.nnz} non-zero bins.")
        else:
            # Store an empty sparse matrix for codewords with no transcripts
            shape = (len(y_bins) - 1, len(x_bins) - 1)
            density_maps[key] = scipy.sparse.csc_matrix(shape)

    return density_maps


zarr_path = '/gpfs/commons/home/jsu/data/xenium_5k_mouse_brain/codeword-update/outs/transcripts.zarr.zip'
output_dir = '/gpfs/commons/home/jsu/data/xenium_5k_mouse_brain/density_maps/'

if __name__ == '__main__':
    # target_codewords = [18067, 18587]
    # target_genes = ['Gnao1', 'Gnal', 'Map4']
    # target_genes = ['Sept5', 'Dtnbp1', 'Aldoa', 'Sept8', 'Bin1', 'Rtn4', 'Ly6h']

    # Get the root data of transcripts.zarr.zip
    root = zarr.open(zarr_path, mode='r')

    # Create output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Extract density maps for all genes with more than one codeword
    all_codeword_genes = root.density.codeword.attrs.get('codeword_names', [])
    avail_genes = list(set(all_codeword_genes))
    # Remove genes starting with 'UnassignedCodeword'
    avail_genes = [gene for gene in avail_genes if not gene.startswith('UnassignedCodeword')]
    multi_codeword_genes = [gene for gene in avail_genes if all_codeword_genes.count(gene) > 1]
    print(f"Number of genes with multiple codewords: {len(multi_codeword_genes)} / {len(avail_genes)}")

    # Extract density maps for the target multi-codeword genes
    density_maps = extract_density_per_codeword(
        zarr_path=zarr_path,
        target_genes=multi_codeword_genes,
        bin_size_um=20.0,
    )

    # Save the dictionary of density maps as a single file using pickle
    with open(os.path.join(output_dir, 'all_multi_codeword_genes.pkl'), 'wb') as f:
        pickle.dump(density_maps, f)

    print(f"Saved all density maps to '{output_dir}/all_multi_codeword_genes.pkl'")
