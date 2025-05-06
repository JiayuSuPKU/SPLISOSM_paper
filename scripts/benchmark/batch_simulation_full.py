"""Script for simulating benchmark isoform expression datasets with varing spatial structures
"""

import os
import itertools
import pickle
import numpy as np
import torch
from isosde.simulation import _sample_multinom_sp_single_gene
from isosde.utils import get_cov_sp
from pyro.distributions import Normal, Multinomial, MultivariateNormal, Poisson

### simulation functions
### =============================
def _simulate_coords_n_design(grid_size, rho, beta_scale):
	torch.manual_seed(42) # fix random seed for reproducibility

	# simulate spatial covariates
	coords = np.array(list(itertools.product(
		np.linspace(0, 1, grid_size[0]),
		np.linspace(0, 1, grid_size[1])
	)))
	n_spots = coords.shape[0]

	# extract the spatial covariance matrix
	cov_sp = get_cov_sp(coords, rho=rho)

	# simulate donut-shaped covariates
	dist_center = np.linalg.norm(coords - coords.mean(0), axis=1)

	# two continuous covariates (need to add white noise to avoid perfect spatial correlation)
	# c1 = np.exp(-2 * np.sin(np.pi * dist_center / (0.5 * dist_center.max()))**2) + 0.1 * np.random.randn(n_spots)
	c1 = 1 - dist_center / dist_center.max() + 0.1 * np.random.randn(n_spots)
	c2 = MultivariateNormal(torch.zeros(n_spots), covariance_matrix=cov_sp + 1e-3 * torch.eye(n_spots)).sample()

	# two binary covariates
	b1 = (dist_center >= 0.25) * (dist_center < 0.5) * 1.0
	b2 = (dist_center >= 0.5) * 1.0

	# merge all covariates
	X_design_full = torch.tensor(np.stack(
		[c1, c2, b1, b2],
		axis = 1
	)).float() # n_spots x 3

	# normalize covariates to have zero mean and unit variance
	X_design_full = (X_design_full - X_design_full.mean(0, keepdim=True)) / X_design_full.std(0, keepdim=True)

	# specify true effect sizes
	beta_true = torch.tensor([
		[1, 0],
		[1, 1],
		[-1, 0],
		[-1, -1]
	]) * beta_scale # 4 x 2

	return coords, cov_sp, X_design_full, beta_true

def _simulate_data(
	gene_variability = 'none',
	iso_usage_variability = 'none',
	n_genes = 1000, # number of genes to simulate
	grid_size = (30, 30), # spatial grid size
	n_isos = 3, # number of isoforms
	total_counts_expected = 5, # mean total counts per spot
	var_sp = 0.1, # variance of spatial component
	var_nsp = 0.1, # variance of non-spatial component
	rho = 0.99, # spatial autocorrelation coefficient
	beta_scale = 0.5, # effect size of covariates
	return_params=True
):
	torch.manual_seed(42) # fix random seed for reproducibility

	assert gene_variability in ['none', 'donut']
	assert iso_usage_variability in ['none', 'mvn', 'donut']

	n_x, n_y = grid_size
	n_spots = n_x * n_y
	n_isos = int(n_isos)

	## simulate spatial covariates and design matrix
	coords, cov_sp, X_design_full, beta_true = _simulate_coords_n_design(grid_size, rho, beta_scale)
	cov_sp_full = var_sp * cov_sp + var_nsp * torch.eye(n_spots)
	prop_var_sp = var_sp / (var_sp + var_nsp) # proportion of re variance explained by spatial component

	## specify the distribution for expected total gene counts (Poisson mean)
	if gene_variability == 'none': # no spatial variance
		tc_mean = torch.ones(n_spots) * total_counts_expected
		tc_dist = Normal(torch.zeros(n_spots), torch.ones(n_spots) * 5)
	else: # Gaussian-process-like spatial variance with donut-shaped covariates
		tc_mean = X_design_full @ torch.tensor([0.0, 0, -1, 1]) * 2.0
		tc_mean = tc_mean - tc_mean.mean() + total_counts_expected
		tc_dist = MultivariateNormal(torch.zeros(n_spots), covariance_matrix=cov_sp * (var_sp + var_nsp) / 2)

	## specify the distribution for isoform usage ratios
	if iso_usage_variability == 'none':
		eta_mean = torch.zeros(n_spots, n_isos - 1)
		eta_dist = Normal(torch.zeros(n_isos - 1, n_spots), torch.ones(n_isos - 1, n_spots) * (var_sp + var_nsp))
		prop_var_sp = 0 # no spatial variance
	elif iso_usage_variability == 'mvn':
		eta_mean = torch.zeros(n_spots, n_isos - 1)
		eta_dist = MultivariateNormal(torch.zeros(n_isos - 1, n_spots), covariance_matrix=cov_sp_full.unsqueeze(0))
	else:
		eta_mean = X_design_full @ beta_true # n_spots x (n_isos - 1)
		eta_dist = MultivariateNormal(torch.zeros(n_isos - 1, n_spots), covariance_matrix=cov_sp_full.unsqueeze(0))

	## sampling for each gene and isoform
	tc_expected_list = []
	props_list = []
	counts_list = []
	for g in range(n_genes):
		# sample total gene counts
		tc_expected = torch.abs(tc_mean + tc_dist.sample()) # n_spots
		tc_expected_list.append(tc_expected)

		# sample isoform usage ratio
		eta = eta_mean + eta_dist.sample().T # n_spots x (n_isos - 1)
		# convert linear predictor to proportions
		props = torch.concat([eta, torch.zeros(n_spots, 1)], dim=-1)
		props = torch.softmax(props, dim=-1) # n_spots x n_isos
		props_list.append(props)

		# sample from multinomial distributions
		counts = _sample_multinom_sp_single_gene(props, tc_expected)
		counts_list.append(counts)

	tc_expected = torch.stack(tc_expected_list, dim = 0) # n_genes x n_spots
	props = torch.stack(props_list, dim = 0) # n_genes x n_spots x n_isos
	counts = torch.stack(counts_list, dim = 0) # n_genes x n_spots x n_isos

	## store simulated data and ground truth
	data = {"counts": counts, "coords": coords, "cov_sp": cov_sp, "design_mtx": X_design_full}

	## return simulation parameters
	if return_params:
		params = {
			"grid_size": grid_size,
			"n_spots": n_spots,
			"n_isos": n_isos,
			"total_counts_expected": tc_expected,
			"rho": rho,
			"var_sp": var_sp,
			"var_nsp": var_nsp,
			"prop_var_sp": prop_var_sp,
			"design_mtx": X_design_full,
			"beta_true": beta_true,
			"iso_ratio_expected": props,
		}
		return data, params
	else:
		return data


if __name__ == "__main__":

	### (1) simulate the general benchmark dataset
	### ==========================================
	# specify simulation parameters
	n_genes = 1000 # number of genes to simulate per scenario
	grid_size = (30, 30) # spatial grid size
	n_isos = 3 # number of isoforms
	total_counts_expected = 5 # mean total counts per spot
	var_sp = 0.1 # variance of spatial component
	var_nsp = 0.1 # variance of non-spatial component
	rho = 0.99 # spatial autocorrelation coefficient
	beta_scale = 0.5 # effect size of covariates

	# specify output directory
	output_dir = "/Users/jysumac/Projects/SPLISOSM_paper/data/simulation_data/general_six_scenarios"
	os.makedirs(output_dir, exist_ok=True)

	# run simulation for all six scenarios
	for gene_variability, iso_usage_variability in itertools.product(
		['none', 'donut'],
		['none', 'mvn', 'donut']
	):
		print(f"=== Simulating the general dataset (gene={gene_variability}, isoform={iso_usage_variability}) ===")

		# simulate genes with spatial variance
		data, params = _simulate_data(
			gene_variability=gene_variability,
			iso_usage_variability=iso_usage_variability,
			n_genes=n_genes,
			grid_size=grid_size,
			n_isos=n_isos,
			total_counts_expected=total_counts_expected,
			var_sp=var_sp,
			var_nsp=var_nsp,
			rho=rho,
			beta_scale=beta_scale
		)

		# save data and params
		with open(f"{output_dir}/data_gene-{gene_variability}_iso-{iso_usage_variability}.pkl", "wb") as f:
			pickle.dump(data, f)

		with open(f"{output_dir}/param_gene-{gene_variability}_iso-{iso_usage_variability}.pkl", "wb") as f:
			pickle.dump(params, f)

	### (2) simulate the SV-specific benchmark dataset (G: donut, I: none -> mvn)
	### =========================================================================
	# specify simulation parameters
	n_genes = 1000 # number of genes to simulate per scenario
	grid_size = (30, 30) # spatial grid size
	n_isos = 3 # number of isoforms
	total_counts_expected = 5 # mean total counts per spot
	rho = 0.99 # spatial autocorrelation coefficient
	beta_scale = 0.5 # effect size of covariates

	# specify output directory
	output_dir = f"/Users/jysumac/Projects/SPLISOSM_paper/data/simulation_data/sv_donut_gene"
	os.makedirs(output_dir, exist_ok=True)

	# run simulation for all six scenarios
	for prop_var_sp in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
		# prop_var_sp = 0.0 <=> general G-donut, I-none
		# prop_var_sp = 0.5 <=> general G-donut, I-mvn
		print(f"=== Simulating the SV dataset (prop_var_sp={prop_var_sp}) ===")

		# specify spatial and non-spatial variance
		var_sp = prop_var_sp * 0.2
		var_nsp = (1 - prop_var_sp) * 0.2

		# simulate genes with spatial variance
		data, params = _simulate_data(
			gene_variability='donut',
			iso_usage_variability='mvn',
			n_genes=n_genes,
			grid_size=grid_size,
			n_isos=n_isos,
			total_counts_expected=total_counts_expected,
			var_sp=var_sp,
			var_nsp=var_nsp,
			rho=rho,
			beta_scale=beta_scale
		)

		# save data and params
		with open(f"{output_dir}/data_prop_var_sp-{prop_var_sp}.pkl", "wb") as f:
			pickle.dump(data, f)

		with open(f"{output_dir}/param_prop_var_sp-{prop_var_sp}.pkl", "wb") as f:
			pickle.dump(params, f)


	### (3) simulate the DU-specific benchmark dataset (G: donut, I: mvn -> donut)
	### =========================================================================
	# specify simulation parameters
	n_genes = 1000 # number of genes to simulate per scenario
	grid_size = (30, 30) # spatial grid size
	n_isos = 3 # number of isoforms
	total_counts_expected = 5 # mean total counts per spot
	rho = 0.99 # spatial autocorrelation coefficient
	var_sp = 0.1 # variance of spatial component
	var_nsp = 0.1 # variance of non-spatial component

	# specify output directory
	output_dir = f"/Users/jysumac/Projects/SPLISOSM_paper/data/simulation_data/du_donut_gene"
	os.makedirs(output_dir, exist_ok=True)

	# run simulation for all six scenarios
	for beta_scale in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]:
		# beta_scale = 0.0 <=> general G-donut, I-mvn
		# beta_scale = 0.5 <=> general G-donut, I-donut
		print(f"=== Simulating the DU dataset (beta_scale={beta_scale}) ===")

		# simulate genes with spatial variance
		data, params = _simulate_data(
			gene_variability='donut',
			iso_usage_variability='donut',
			n_genes=n_genes,
			grid_size=grid_size,
			n_isos=n_isos,
			total_counts_expected=total_counts_expected,
			var_sp=var_sp,
			var_nsp=var_nsp,
			rho=rho,
			beta_scale=beta_scale
		)

		# save data and params
		with open(f"{output_dir}/data_beta-{beta_scale}.pkl", "wb") as f:
			pickle.dump(data, f)

		with open(f"{output_dir}/param_beta-{beta_scale}.pkl", "wb") as f:
			pickle.dump(params, f)
