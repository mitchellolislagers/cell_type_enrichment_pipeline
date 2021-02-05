## Setting up LDSC environment

### (1) Clone LDSC repository from GitHub
```git clone https://github.com/bulik/ldsc.git```
```cd ldsc```

### (2) Set up designated environment
```wget https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh```
```bash Miniconda2-latest-Linux-x86_64.sh```
```conda env create --file environment.yml```
```source activate ldsc```

### (3) Munge sum stats
- ```munge_sumstats.sh```

### (4) Heritability
- ```calculate_and_collect_h2_phase3.sh```

### (5) Bivariate correlations
- ```calculate_and_collect_rg_phase3.sh```

### (6) order of scripts (found in respective folders for KI level 1, KI level 2 and 10X genomics
- ```make_gene_set_files```
- ```make_annotation_files```
- ```calculate_ld_scores```
- ```partition_h2```
