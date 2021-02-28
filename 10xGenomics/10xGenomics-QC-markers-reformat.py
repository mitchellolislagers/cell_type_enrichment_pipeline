#!/usr/bin/python

# ABOUT: 	Quality control, marker gene detection and data formatting (cell-type specificity of genes) for
#		randomized representative subset of ~108,000 cells of the 10x Genomics mouse brain single-cell
#		RNA-sequencing dataset
# REQUIRED: 	- Input gene / cell matrix in H5 format (filename_input_h5)
#		- Clustering identities of cells (filename_cluster_ids)
#		- Annotated gene / cell matrix with all cells in H5AD format (filename_h5ad)
#		- Annotated gene / cell matrix with QC and marker gene detection applied in H5AD format (filename_h5ad_qc_markers)
#		- Ensembl mouse gene symbol > mouse gene identifier mapping (filename_mouse_mapping)
#		- Ensembl mouse gene identifier > human gene identifier mapping (filename_mouse_to_human)
#		- Output of reformatted data (filename_reformatted_out)
#		- Output logging file (filename_log_out)
# AUTHOR:	Koen Rademaker, GitHub repository 'perslab-sc-library' (https://github.com/perslab/perslab-sc-library, customized code for own purposes)
# DATE:	5 June 2019, updated 28 February 2021 for publication of paper


########## Import packages ##########
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib
import matplotlib.pyplot as plt
import os
import sys


########## Set variables ##########
filename_input_h5 = sys.argv[1]
filename_cluster_ids = sys.argv[2]
filename_h5ad = sys.argv[3]
filename_h5ad_qc_markers = sys.argv[4]
filename_mouse_mapping = sys.argv[5]
filename_mouse_to_human = sys.argv[6]
filename_reformatted_out = sys.argv[7]
filename_log_out = sys.argv[8]
genome = 'mm10'


########## Set Matplotlib & Scanpy settings ##########
matplotlib.use('Agg')
sc.logging.print_versions()
sc.settings.verbosity = 5
sc.settings.logfile = filename_log_out
sc.settings.autosave = True
sc.settings.max_memory = 128
sc.settings.n_jobs = 8
sc.settings.set_figure_params(dpi=600)


########## Function declaration ##########
# Function to map from Mm symbol to Hs Ensembl gene identifiers (from https://github.com/perslab/perslab-sc-library/blob/master/gene_mapping.py)
def get_mapping(df_mm2mm,df_mm2hs,mm_symbol):
	if mm_symbol in df_mm2mm.index:							# (Koen: Removed call to PDB debugging)
		# Discard many-to-many mappings (e.g. Pou6f1, which maps to ENSMUSG00000009739 and ENSMUSG00000098598)
		matches_mm = df_mm2mm.ix[mm_symbol,'Ensembl Gene ID']
		if isinstance(matches_mm, str) and matches_mm in df_mm2hs.index: # isinstance check fails if several matches for mm_symbol
			# Keep one-to-many mappings, by selecting most idential human gene
			matches_hs = df_mm2hs.ix[matches_mm,:]
			if len(matches_hs.shape) > 1:
				row_index = np.argmax(matches_hs.ix[:,'% Identity with respect to Human gene'].tolist())
				return matches_hs.ix[row_index,'Human Ensembl Gene ID']
			else:
				return matches_hs['Human Ensembl Gene ID']
	return 'not_found'

# Function to convert mouse gene identifiers to Ensembl human identifiers (from https://github.com/perslab/perslab-sc-library/blob/master/gene_mapping.py)
def to_ensembl(df_mm2mm,df_mm2hs,df):
	rows2drop = []
	mapping = []
	for ix in df.index.tolist():
		hs_ensembl_id = str(get_mapping(df_mm2mm,df_mm2hs,ix))
		if hs_ensembl_id == 'not_found' or hs_ensembl_id == 'nan':		        # (Koen: Updated to exclude both not_found and nan values)
			rows2drop.append(ix)
		else:
			mapping.append(hs_ensembl_id)
	df.drop(rows2drop,axis=0,inplace=True)
	df['Ensembl Gene ID'] = pd.Series(mapping, index=df.index)
	return df.groupby('Ensembl Gene ID',sort=False).mean(),rows2drop

# Function to normalize to 10k UMI and take log (from https://github.com/perslab/perslab-sc-library/blob/master/dropseq.py)
def normalize(df):
	dge = df.values								        # (Koen: Replaced .as_matrix() with .values as the former is deprecated)
	col_sums = np.apply_along_axis(sum,0,dge)
	mat_dge_norm =  np.log( dge/[float(x) for x in col_sums] * 10000 + 1 )
	df_dge_norm = pd.DataFrame(mat_dge_norm,index=df.index,columns=df.columns)
	df_dge_norm.drop(df_dge_norm.index[df_dge_norm.sum(axis=1) == 0],axis=0,inplace=True)
	return df_dge_norm

# Function to average cells by cluster (from https://github.com/perslab/perslab-sc-library/blob/master/dropseq.py)
def get_average_by_celltype(df_dge,df_cluster):
	df_cluster = df_cluster.merge(df_dge.transpose(),left_index=True,right_index=True,how='inner').groupby('cluster_id',sort=False).mean().transpose()
	df_cluster.columns = df_cluster.columns.astype(int)			        	# (Koen: Added code to sort columns in ascending order)
	df_cluster.columns.sort_index(axis=1, inplace=True)
	return df_cluster

# Function to standardize genes' expresion across cell types (from https://github.com/perslab/perslab-sc-library/blob/master/dropseq.py)
def standardize(df):
        return df.sub(df.mean(axis=1),axis=0).div(df.std(axis=1),axis=0)

# Function to remove cells of specific clusters from the AnnData object
def remove_cluster_cells(sc_data_obj, clusters_to_remove):
	cells_to_remove = []
	for cluster in clusters_to_remove:
		cells = sc_data_obj[sc_data_obj.obs['cell_labels']==str(cluster)].obs_names.tolist()
		for cell in cells:
			cells_to_remove.append(cell)
	cells_to_keep = [cell for cell in sc_data_obj.obs_names if (cell not in cells_to_remove)]
	return sc_data_obj[cells_to_keep, :]




# Module to prepare all the data
def prepare_data():
	########## Load mouse and human gene mapping data ##########
	mouse_mapping = pd.read_csv(filename_mouse_mapping, compression='gzip', index_col=1, sep='\t')
	mouse_mapping = mouse_mapping[[isinstance(x, str) for x in mouse_mapping.index]]
	mouse_to_human = pd.read_csv(filename_mouse_to_human, compression='gzip', index_col=0, sep='\t')


	########## Load data from H5 ##########
	sc_data=sc.read_10x_h5(filename_input_h5, genome)
	sc_data.var_names_make_unique()
	sc_data.obs['cell_labels'] = pd.read_csv(filename_cluster_ids, header=None, skiprows=1, dtype='category')[1].values


	########## Save data to H5AD for later use (plotting) ##########
	sc_data.write_h5ad(filename_h5ad)


	########## Calculate cell/gene metrics ##########
	mito_genes=sc_data.var_names.str.startswith('mt-')
	sc_data.obs['percent_mito'] = np.sum(sc_data[:, mito_genes].X, axis=1).A1 / np.sum(sc_data.X, axis=1).A1
	sc_data.obs['n_counts'] = sc_data.X.sum(axis=1).A1
	sc.pp.filter_cells(sc_data, min_genes=0)
	sc.pp.filter_genes(sc_data, min_cells=0)

	mito_mean=np.mean(sc_data.obs['percent_mito'])
	mito_sd=np.std(sc_data.obs['percent_mito'])

	umi_mean=np.mean(sc_data.obs['n_counts'])
	umi_sd=np.std(sc_data.obs['n_counts'])

	gene_mean=np.mean(sc_data.obs['n_genes'])
	gene_sd=np.std(sc_data.obs['n_genes'])

	cell_mean=np.mean(sc_data.var['n_cells'])
	cell_sd=np.std(sc_data.var['n_cells'])


	########## Plot cell quality metrics ##########
	sc.pl.violin(sc_data, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='_pre_qc.png')
	sc.pl.scatter(sc_data, x='n_counts', y='percent_mito', save='_counts_mito_pre_qc.png')
	sc.pl.scatter(sc_data, x='n_counts', y='n_genes', save='_counts_genes_pre_qc.png')


	########## Apply QC ##########
	sc.pp.filter_cells(sc_data, min_counts=1) 						        # Remove cells with no UMI counts
	sc.pp.filter_cells(sc_data, min_counts=umi_mean-3*umi_sd)				        # Remove cells with a UMI count 3 SDs below the mean
	sc.pp.filter_cells(sc_data, max_counts=umi_mean+3*umi_sd)				        # Remove cells with a UMI count 3 SDs above the mean
	sc.pp.filter_cells(sc_data, min_genes=1000)						        # Remove cells with less than 1000 genes
	sc.pp.filter_cells(sc_data, max_genes=gene_mean+3*gene_sd)				        # Remove cells with a gene count 3 SDs above the mean
	sc_data=sc_data[sc_data.obs['percent_mito'] < mito_mean+3*mito_sd,:]			        # Remove cells with a % mtDNA 3 SDs above the mean
	sc.pp.filter_genes(sc_data, min_counts=1)						        # Remove genes with no expression


	########## Remove cells from specific clusters ##########
	sc_data = remove_cluster_cells(sc_data, [12, 18, 19, 20])					# Remove cells from clusters 12, 18, 19 and 20


	########## Plot cell quality metrics ##########
	sc.pl.violin(sc_data, ['n_genes', 'n_counts', 'percent_mito'], jitter=0.4, multi_panel=True, save='_post_qc.png')
	sc.pl.scatter(sc_data, x='n_counts', y='percent_mito', save='_counts_mito_post_qc.png')
	sc.pl.scatter(sc_data, x='n_counts', y='n_genes', save='_counts_genes_post_qc.png')


	########## Detect marker genes ##########
	sc.tl.rank_genes_groups(sc_data, 'cell_labels', method='wilcoxon', corr_method='benjamini-hochberg')


	########## Save data to H5AD for later use ##########
	sc_data.write_h5ad(filename_h5ad_qc_markers)

#############################################################



# Module to reformat data into specificity
def reformat_data():
	########## Reformat data ##########
	cluster_id_cells = sc_data.obs['cell_labels'].to_frame()
	cluster_id_cells.columns = ['cluster_id']
	normalized = normalize(sc_data.to_df().T)							# Normalize cells
	cluster_averaged = get_average_by_celltype(normalized, cluster_id_cells)			# Average gene expression per cluster
	mapped_to_human, not_mapped = to_ensembl(mouse_mapping, mouse_to_human, cluster_averaged) 	# Map to human genes
	standardized = standardize(mapped_to_human)			 				# Standardize gene expression across cell types
	standardized.to_csv(filename_reformatted_out, sep='\t', header=True, index=True)		# Export final matrix
	

def main():
	# prepare_data()
	sc_data=sc.read_h5ad(filename_h5ad_qc_markers)
	reformat_data()

if __name__ == '__main__':
	main()
