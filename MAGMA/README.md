## How to run MAGMA

### (1) Prepare summary statistics
Because summary statistics are not uniformly released, preparation differs per summary statistics.
```prepare_MAGMA_sumstats.sh```

### (2) (optional) Sensitivity analyses (HapMap3 restriction and MHC exlcusion)
```QC_MAGMA.sh```

### (3) Run MAGMA Top 10% mode and Linear mode using 10x Genomics and KI level 1 scRNA data
```MAGMA_v.1.08_celltype_analysis_hm3_mhc_restricted.R```

### (4) Run MAGMA Top 10% mode and Linear mode using KI level 2 scRNA data
```MAGMA_v.1.08_celltype_analysis_hm3_mhc_restricted_KI_level2_v33.R```

### (5) (optional) Run sensitivity analyses in schizophrenia
```MAGMA_v108_celltype_analysis_SCZ_all_QC.R```

### Reference files
To map rsids for SNPs withour CHRPOS, 1kg .bim files were used.

Reference data used for MAGMA analyses can be found [here](https://ctg.cncr.nl/software/MAGMA/ref_data/).

Also see the [original repository](https://github.com/NathanSkene/MAGMA_Celltyping)

[Skene, et al. Genetic identification of brain cell types underlying schizophrenia. Nature Genetics, 2018.](https://www.nature.com/articles/s41588-018-0129-5)