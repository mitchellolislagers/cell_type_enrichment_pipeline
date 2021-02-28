## Running the 10x Genomics cell type annotation pipeline

### (0) Requirements
Code was run using Python version 3.7 and requires the complete installation of the [SCANPY](https://scanpy.readthedocs.io/en/stable/) library to function properly. A cluster computing environment is highly recommended since the large transformations of data matrices requires a lot of memory.

### (1) Extract compressed data to .h5ad object
H5AD data objects are a SCANPY implementation of the [hdf5 file type](https://www.h5py.org/), which enables the storage of large transcriptomics datasets together with all forms of QC, annotation and results.

Given the size of 1.7GB for the uncompressed file for 100436 cells (see manuscript for details on quality control, marker gene identification) and GitHub file size limits, we compressed the file into 50Mb chunks which you can extract into one single .h5ad file using:

```cat cell_type_enrichment_pipeline/10xGenomics/data/10x_Genomics_data_blocks.*.tar.gz | tar xzvf -```

### (2) Reformat data for cell type enrichment analysis
The script _10xGenomics-QC-markers-reformat.py_ contains functionality for quality control and marker gene detection, however specific functions have been disabled since the resulting .h5ad file is already provided. Hence, not all script parameters have to be furfilled for proper functioning (marked with "OBSOLETE").

Note that sufficient RAM is required to load and process the .h5ad data object, attempt at least 8GB first and otherwise scale up for what your computing cluster allows for.

**See also commentary in the script for details of input/output files!**

```python 10xGenomics-QC-markers-reformat.py OBSOLETE_H5 OBSOLETE_CLUSTERING OBSOLETE_H5AD data/108k_kmeans_20_qc_clusters_markers.h5ad data/ensembl_v96_ensembl_genename_Mm.txt.gz data/ensembl_v96_Mm_Hs_GRCh37.txt.gz {give output specificity file} {give output log file}```

### (3) Questions
For further details about the code / approach please ask the corresponding author.
