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
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/annotations_v33
ld_scores_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/ld_scores_v33
tmp_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/temp_v33
miniconda_dir=~/miniconda2

# Create folder paths
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
