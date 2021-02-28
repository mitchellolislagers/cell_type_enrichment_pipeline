## Running the 10x Genomics cell type annotation pipeline

### (1) Extract compressed data to .h5ad object
Explain

cat data/10x_Genomics_data_blocks.*.tar.gz | tar xzvf -

### (2) 
Explain

Note on RAM usage

python 10xGenomics-QC-markers-reformat.py {1} {2} {3} {108k_kmeans_20_qc_clusters_markers.h5ad} {mouse_mapping?} {mouse_to_human?} {specificity} {log} 

### (3) 
