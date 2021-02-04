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
genesets_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/geneset_v33
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/annotations_v33
tmp_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/temp_v33
miniconda_dir=~/miniconda2
# Create folders
mkdir -p ${tmp_dir}
# Declare variables
declare -a cell_types=("Astro1" "Astro2" "CA1Pyr1" "CA1Pyr2" "CA1PyrInt" "CA2Pyr2" "Choroid" "ClauPyr" "Committed_oligodendrocyte_precursors" "Dopaminergic_Adult-Substantia_nigra" 
"Dopaminergic_Adult-Ventral_tegmental_area1" "Dopaminergic_Adult-Ventral_tegmental_area2" "Dopaminergic_Adult-Ventral_tegmental_area3" "Dopaminergic_Adult-Ventral_tegmental_area4" 
"Dopaminergic_Neuroblast" "Embryonic_Dopaminergic_Neuron_0" "Embryonic_Dopaminergic_Neuron_1" "Embryonic_Dopaminergic_Neuron_2" "Embryonic_GABAergic_Neuron_1a" 
"Embryonic_GABAergic_Neuron_1b" "Embryonic_GABAergic_Neuron_2" "Epend" "Hypothalamic_Adcyap1;_Tac1_(VMH)_Neuron" "Hypothalamic_Adcyap1(VMH)_Neuron" "Hypothalamic_Avp-high;Gal;Oxt-low_Neuron" 
"Hypothalamic_Avp-high;Oxt-low_Neuron" "Hypothalamic_Avp-medium_Neuron" "Hypothalamic_Calcr-high;VMAT+and-;Lhx1_(Sch)_Neuron" "Hypothalamic_Crh+and-;_Lhx6;GABA_(BST;_MPO)_Neuron" 
"Hypothalamic_Crh+and-;GABA;Pgr15l_(PVH-PVH-galo)_Neuron" "Hypothalamic_Crh+and-;Gal-low_(PVH)_Neuron" "Hypothalamic_Dopamine;_Dat;Nmur2;GABA_Neuron" 
"Hypothalamic_Dopamine;Tac+and-_Gad1;_GABA_Neuron" "Hypothalamic_Dopamine;Tac1;Ghrh;Pnoc;_Dat+and-;_GABA_Neuron" "Hypothalamic_GABA_1_Neuron" "Hypothalamic_GABA_2_Neuron" 
"Hypothalamic_GABA;_Gucy1a3_Neuron" "Hypothalamic_GABA;Pnoc;Tac2+and-_Neuron" "Hypothalamic_Gad-low;Gnrh-and+_Neuron" "Hypothalamic_Gad1;Gad2;VGAT+and-;Pnoc_Neuron" 
"Hypothalamic_Galanin_Neuron" "Hypothalamic_Ghrh-high;Th;Vglut2_Neuron" "Hypothalamic_Hcrt_Neuron" "Hypothalamic_Hmit+and-_Neuron" 
"Hypothalamic_Nms;VIP+and-;_Avp-low+and-;circadian_(SCH)_Neuron" "Hypothalamic_Npr2;Gm5595;4930422G04Rik;Tnr_Neuron" 
"Hypothalamic_Npvf_Neuron" "Hypothalamic_Npy-medium;Gad1;Gad2_Neuron" "Hypothalamic_Npy;Agrp_(ARH)_Neuron" "Hypothalamic_Nts;Gal;Pnoc;Tac2+and-;GABA_Neuron" 
"Hypothalamic_Nts;Pnoc;Tac2+and-;GABA_Neuron" "Hypothalamic_Otof;Lhx1+and-_Neuron" "Hypothalamic_Oxt;Avp-low_Neuron" "Hypothalamic_Oxt;Avp-low;Th;Cacna1H_Neuron" 
"Hypothalamic_Oxt;Avp-medium;Gad2-low_Neuron" "Hypothalamic_Oxt;Avp-medium;Th+and-_PDYN+and_Neuron" "Hypothalamic_Per2;circadian_(SCH)_Neuron" "Hypothalamic_Pmch_Neuron" 
"Hypothalamic_Pomc+and-_(ARH)_Neuron" "Hypothalamic_Qrfp_Neuron" "Hypothalamic_Sst_low-medium_Neuron" "Hypothalamic_Sst-high;Cartpt;Galr1_Neuron" "Hypothalamic_Sst-medium;Cartpt_Neuron" 
"Hypothalamic_Th;Tac1;Ghrh;Vmat2-lowand-_GABA_Neuron" "Hypothalamic_Trh-high;Adcyap1;Cartpt_Neuron" "Hypothalamic_Trh-low_Neuron" "Hypothalamic_Trh-medium_Neuron" 
"Hypothalamic_Vglut2_;Col9a2_(PVH)_Neuron" "Hypothalamic_Vglut2_A_Neuron" "Hypothalamic_Vglut2;A930013F10Rik;Pou2f2_Neuron" "Hypothalamic_Vglut2;Bdnf;Gad1_Neuron" 
"Hypothalamic_Vglut2;Cck;2AG;AEA;Npr2;Npy1r_Neuron" "Hypothalamic_Vglut2;Cnr1;Ninl;Rfx5;Zfp346_Neuron" "Hypothalamic_Vglut2;Cnr1;Npr2;Cacna1h_Neuron" 
"Hypothalamic_Vglut2;Gad1-low;Crh-lowand-_(VMH)_Neuron" "Hypothalamic_Vglut2;Gpr149_(PVH;VMH)_Neuron" "Hypothalamic_Vglut2;Hcn1;6430411K18Rik_Neuron" 
"Hypothalamic_Vglut2;Morn4;Prrc2a_Neuron" "Hypothalamic_Vglut2;Myt1;Lhx9_Neuron" "Hypothalamic_Vglut2;Penk;Oprk1_Neuron" "Hypothalamic_Vglut2;Pgam;Snx12_Neuron" 
"Hypothalamic_Vglut2;Prmt8;Ugdh_(VMH)_Neuron" "Hypothalamic_Vglut2;Zfp458;Ppp1r12b;Cacna1e-highest_Neuron" "Hypothalamic_Vip;Grp+and-;circadian_(SCH)_Neuron" "Int_Pvalb" "
Int1" "Int10" "Int11" "Int12" "Int13" "Int14" "Int15" "Int16" "Int2" "Int4" "Int5" "Int6" "Int7" "Int8" "Int9" "Lateral_Neuroblasts_1" "Lateral_Neuroblasts_2" 
"Mature_oligodendrocytes_1" "Mature_oligodendrocytes_2" "Mature_oligodendrocytes_3" "Mature_oligodendrocytes_4" "Mature_oligodendrocytes_5" "Mature_oligodendrocytes_6" 
"Medial_Neuroblasts_" "Mediolateral_Neuroblasts_1" "Mediolateral_Neuroblasts_2" "Mediolateral_Neuroblasts_3" "Mediolateral_Neuroblasts_4" "Mediolateral_Neuroblasts_5" 
"Medium_Spiny_Neuron_D1R" "Medium_Spiny_Neuron_D2R" "Mgl1" "Mgl2" "Myelin-forming_oligodendrocytes_1" "Myelin-forming_oligodendrocytes_2" "Neural_Progenitors" 
"Newly_formed_oligodendrocytes_1" "Newly_formed_oligodendrocytes_2" "Oculomotor_and_Trochlear_nucleus_embryonic_neurons" "Oligodendrocyte_Precursor" "Peric" "Pvm1" "Pvm2" 
"Radial_glia_like_cells_1" "Radial_glia_like_cells_2" "Radial_glia_like_cells_3" "Red_nucleus_embryonic_neurons" "S1PyrDL" "S1PyrL23" "S1PyrL4" "S1PyrL5" "S1PyrL5a" "S1PyrL6" 
"S1PyrL6b" "Serotonergic_Neuron" "Striatal_CHAT_Interneuron" "Striatal_Interneurons_(other)" "Striatal_Pvalb_Interneuron" "Striatal_Sst_Interneuron" "SubPyr" 
"Vascular_Leptomeningeal_Cells" "Vend1" "Vend2" "Vsmc")
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
dataset="KI"
# Activate LDSC environment
source ${miniconda_dir}/bin/activate ldsc

########## Replace all starting positions of genes lower than 100001 to 100001 to allow a window size of 100000
awk 'BEGIN{FS="\t"}{OFS="\t"} $3<100001 {$3=100001}1' ${gene_coord_file} > ${ref_dir}/ENSG_coord_windowsize.txt

gene_coord_file=${ref_dir}/ENSG_coord_windowsize.txt

########## Replace whitespace with underscore in filenames ###########
for file in ${genesets_dir}/*; do
	if [[ "${file}" == *" "* ]]; then
		mv "${file}" `echo ${file} | tr ' ' '_'`
	fi
done


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
