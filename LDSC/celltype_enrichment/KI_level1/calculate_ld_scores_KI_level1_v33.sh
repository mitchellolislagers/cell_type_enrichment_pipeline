#!/bin/bash

# TITLE:	calculate_ld_scores.sh
# GOAL:		Calculate LD scores for specificity deciles of cell types
# INPUT:	LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types
# INPUT:	1000 Genomes Project phase 3 plink files, HapMap3 SNPs (will automatically be downloaded and organized)
# OUTPUT:	LD score output files; .M (total number of SNPs), .M_5_50 (number of SNPs with MAF > 5%), .ldscore.gz (SNP LD scores)
# AUTHOR:	Koen Rademaker, adapted by Mitchell Olislagers
# DATE:		20 February 2020


########## Initialize script ##########
# Set folder paths
ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment/
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level1/annotations_v33
ld_scores_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level1/ld_scores_v33
tmp_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level1/temp_v33
miniconda_dir=~/miniconda2

# Create folder paths
mkdir -p ${tmp_dir}

# Declare variables
declare -a cell_types=("astrocytes_ependymal" "Dopaminergic_Adult" "Dopaminergic_Neuroblast" "Embryonic_Dopaminergic_Neuron" "Embryonic_GABAergic_Neuron" "Embryonic_midbrain_nucleus_neurons" "endothelial-mural" "Hypothalamic_Dopaminergic_Neurons" "Hypothalamic_GABAergic_Neurons" "Hypothalamic_Glutamatergic_Neurons" "interneurons" "Medium_Spiny_Neuron" "microglia" "Neural_Progenitors" "Neuroblasts" "Oligodendrocyte_Precursor" "Oligodendrocytes" "Oxytocin_and_Vasopressin_Expressing_Neurons" "pyramidal_CA1" "pyramidal_SS" "Radial_glia_like_cells" "Serotonergic_Neuron" "Striatal_Interneuron" "Vascular_Leptomeningeal_Cells")

########## Organize LDSC ##########
# Get SNPs from hapmap3
tail -n+2 /hpc/hers_en/molislagers/LDSC/ref_data/regression/w_hm3.snplist | cut -f1 > ${ref_dir}/w_hm3.snplist.snps

# Declare variables
phase3_1000G_dir=${ref_dir}/1000G_EUR_Phase3_plink
ld_window=1
dataset="KI"
# Activate LDSC environment
source ${miniconda_dir}/bin/activate ldsc

########## Calculate LD scores with annotation files ##########
# (1) Iterate over cell types
for ct in ${cell_types[@]}; do
	# (2) Iterate over chromosomes
	for chr in {1..22}; do
		# (3) Calculate LD scores
		python ${ldsc_dir}/ldsc.py \
			--l2 \
			--bfile ${phase3_1000G_dir}/1000G.EUR.QC.${chr} \
			--ld-wind-cm ${ld_window} \
			--annot ${annotation_dir}/${dataset}_${ct}.${chr}.annot.gz \
			--thin-annot \
			--print-snps ${ref_dir}/w_hm3.snplist.snps \
			--out ${ld_scores_dir}/${dataset}_${ct}.${chr}
	done
done
