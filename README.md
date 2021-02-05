## Integrating scRNA-seq with GWAS using LDSC, MAGMA, DEPICT and FUMA

Last update: 05 Feb 2020

### Requirements
- [LDSC](https://github.com/bulik/ldsc)
- [MAGMA](https://ctg.cncr.nl/software/magma)
- [DEPICT](https://data.broadinstitute.org/mpg/depict/depict_download/bundles/)
- [PLINK](https://www.cog-genomics.org/plink/)
- [Specificity metrices of scRNA-seq data](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/scRNA_datasets).

### The pipeline
All of the code neccesary to reproduce our findings are found within this repository.

### LDSC
Follow the [steps](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/LDSC/README.md), including:
- Install [LDSC](https://github.com/bulik/ldsc).
- [Munge sumstats](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/LDSC/munge_sumstats/munge_sumstats.sh).
- Estimate [heritability](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/LDSC/heritability/calculate_and_collect_h2.sh) and [genetic correlations](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/LDSC/bivariate_correlations/calculate_and_collect_rg.sh).
- Run the [cell type enrichment analyses](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/LDSC/celltype_enrichment).

### MAGMA
Follow the [steps](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/MAGMA/README.md), including:
- [Preparing summary statistics](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/MAGMA/sum_stats/prepare_MAGMA_sumstats.sh).
- A MAGMA [sensitivity analysis](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/MAGMA/sum_stats/QC_MAGMA.sh).
- Running MAGMA in Top 10% and Linear mode.

### DEPICT
- [Configuration files](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/DEPICT/config_files) of all summary statistics and the [mapping files](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/DEPICT/mapping) used for DEPICT are included in this repository.
- For DEPICT analyses, summary statistics from the [MAGMA summary statistics preparation](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/MAGMA/sum_stats/prepare_MAGMA_sumstats.sh) can be used.

### FUMA
FUMA was run on the [web-based platform](https://fuma.ctglab.nl/). 
- Summary statistics from the [MAGMA summary statistics preparation](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/MAGMA/sum_stats/prepare_MAGMA_sumstats.sh) can be used as input for the [SNP2GENE](https://fuma.ctglab.nl/snp2gene) process.
- The gene analysis result can then be used as input for the [Cell type](https://fuma.ctglab.nl/celltype) process.