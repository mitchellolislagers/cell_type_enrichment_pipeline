#!/bin/bash

# TITLE:	make_annotation_files.sh
# GOAL:		Make LDSC annotation files for cell types per autosomal chromosome, with specificity deciles constituting as sub-annotations
# INPUT:	Gene sets of Ensembl IDs for human genes in a given sub-annotation, as generated by the 'make_gene_set_files.R' script
# INPUT:	Ensembl gene coordinates, 1000 Genomes Project Phase 3 plink files
# OUTPUT:	LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types
# AUTHOR:	Koen Rademaker, adapted by Mitchell Olislagers
# DATE:		20 February 2020


########## Initialize script ##########
# Set folder paths
ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment
genesets_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/10X_genomics/geneset
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/10X_genomics/annotations
tmp_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/10X_genomics/temp
miniconda_dir=~/miniconda2
# Create folders
mkdir -p ${tmp_dir}
# Declare variables
declare -a cell_types=("Glutamatergic_neurons" "Neuroblasts_1" "Astrocytes_1" "Neuroblasts_2" "Intermediate_progenitors" "Enteric_glial_cells" "Interneurons_1" "Neurons" 
"Interneurons_2" "Vascular_endothelial_cells" "Enteric_neurons" "Cajal_Retzius_cells" "Oligodendrocytes" "Astrocytes_2" "Endothelial_cells" "Microglia")
declare -a specificity_deciles=(N 1 2 3 4 5 6 7 8 9 10)

########## Organize LDSC ##########
# Download background data
#echo "-- Downloading Ensembl gene coordinates files --"
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/make_annot_sample_files/ENSG_coord.txt -P ${ref_dir}
#echo "-- Downloading 1000 Genomes Project phase 3 plink files--"
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_plinkfiles.tgz -P ${ref_dir}
#tar -xvzf ${ref_dir}/1000G_Phase3_plinkfiles.tgz
# Declare variables
gene_coord_file=${ref_dir}/ENSG_coord.txt
phase3_1000G_dir=${ref_dir}/1000G_EUR_Phase3_plink
window_size=100000
dataset="10x_Genomics"
# Activate LDSC environment
source ${miniconda_dir}/bin/activate ldsc

########## Make LDSC annotation files ##########
# (1) Iterate over cell types
for ct in ${cell_types[@]}; do
	# (2) Iterate over chromosomes	
	for chr in {1..22}; do
		# (3) Copy gene set files for specificity deciles to temporary folder
		for gene_set_file in ${genesets_dir}/${dataset}.${ct}.${chr}.*.GeneSet; do cp "$gene_set_file" "${tmp_dir}"; done
		# (4) Iterate over specificity deciles
		for decile in "${specificity_deciles[@]}"; do
			# (5) Create sub-annotation file
			echo "-- Creating sub-annotation file for ${ct} (cell type), chromosome ${chr}, specificity decile ${decile} --"
			gene_set_file=${tmp_dir}/${dataset}.${ct}.${chr}.${decile}.GeneSet
			python ${ldsc_dir}/make_annot.py \
				--gene-set-file ${gene_set_file} \
				--gene-coord-file ${gene_coord_file} \
				--windowsize ${window_size} \
				--bimfile ${phase3_1000G_dir}/1000G.EUR.QC.${chr}.bim \
				--annot-file ${tmp_dir}/${dataset}_${ct}_${chr}_${decile}.annot
			# (6) Update sub-annotation file header
			echo "-- Updating sub-annotation file header --"
			sed -i "1 s/^.*$/D_${decile}/" ${tmp_dir}/${dataset}_${ct}_${chr}_${decile}.annot
		done
		# (7) Merge sub-annotation files into single annotation file
		echo "-- Merging sub-annotation files to single annotation file for ${ct} (cell type) on chromosome ${chr} --"
		paste -d "\t" ${tmp_dir}/${dataset}_${ct}_${chr}_{N,1,2,3,4,5,6,7,8,9,10}.annot > ${annotation_dir}/${dataset}_${ct}.${chr}.annot
		gzip ${annotation_dir}/${dataset}_${ct}.${chr}.annot
		# (8) Clear temporary folder
		rm ${tmp_dir}/*
	done	
done

