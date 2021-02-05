## Integrating scRNA-seq with GWAS using LDSC, MAGMA, DEPICT and FUMA

Last update: 05 Feb 2020

### Requirements
- [LDSC](https://github.com/bulik/ldsc)
- [MAGMA](https://ctg.cncr.nl/software/magma)
- [DEPICT](https://data.broadinstitute.org/mpg/depict/depict_download/bundles/)

### The pipeline
All of the code neccesary to reproduce our findings are found withing this repository.

### LDSC
Follow the [steps](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/LDSC/README.md), including:
- Install [LDSC](https://github.com/bulik/ldsc)
- [Munge sumstats](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/LDSC/munge_sumstats/munge_sumstats.sh).
- Estimate [heritability](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/LDSC/heritability/calculate_and_collect_h2.sh) and [genetic correlations](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/blob/master/LDSC/bivariate_correlations/calculate_and_collect_rg.sh).
- Run the [cell type enrichment analyses](https://github.com/mitchellolislagers/cell_type_enrichment_pipeline/tree/master/LDSC/celltype_enrichment).
