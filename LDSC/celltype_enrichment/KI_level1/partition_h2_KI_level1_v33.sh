#!/bin/bash
  
# TITLE:        partition_h2.sh
# GOAL:         Partition heritability (h2) to specificity deciles of cell types.
# INPUT:        LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types, LDSC LD score files 
# INPUT:	1000 Genomes Project phase 3 baseline model LD scores, regression weights, allele frequencies (see README)
# OUTPUT:       Partitioned h2 output files
# AUTHOR:       Koen Rademaker
# DATE:         26 March 2020


########## Initialize script ##########
# Set paths
ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment
sum_stats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics/munged_sumstats
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level1/annotations_v33
ld_scores_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level1/ld_scores_v33
out_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level1/partitioned_h2_v33
miniconda_dir=~/miniconda2
# Declare variables
declare -a sum_stats=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance" "neuroticism")
declare -a cell_types=("astrocytes_ependymal" "Dopaminergic_Adult" "Dopaminergic_Neuroblast" "Embryonic_Dopaminergic_Neuron" "Embryonic_GABAergic_Neuron" "Embryonic_midbrain_nucleus_neurons" "endothelial-mural" "Hypothalamic_Dopaminergic_Neurons" "Hypothalamic_GABAergic_Neurons" "Hypothalamic_Glutamatergic_Neurons" "interneurons" "Medium_Spiny_Neuron" "microglia" "Neural_Progenitors" "Neuroblasts" "Oligodendrocyte_Precursor" "Oligodendrocytes" "Oxytocin_and_Vasopressin_Expressing_Neurons" "pyramidal_CA1" "pyramidal_SS" "Radial_glia_like_cells" "Serotonergic_Neuron" "Striatal_Interneuron" "Vascular_Leptomeningeal_Cells")
dataset="KI"


########## Organize LDSC ##########
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz -P ${ref_dir}
#tar -xvzf 1000G_Phase3_frq.tgz
# Copy annotation files and LD score files to reference folder
ref_ld_dir=${ref_dir}/ref_ld_v33
cp ${ld_scores_dir}/*l2* ${ref_ld_dir}
cp ${annotation_dir}/*.annot.gz ${ref_ld_dir}
# Declare variables
weights_dir=/hpc/hers_en/molislagers/LDSC/ref_data/regression/1000G_Phase3_weights_hm3_no_MHC
baseline_dir=/hpc/hers_en/molislagers/LDSC/ref_data/regression/1000G_EUR_Phase3_baseline
frq_dir=${ref_dir}/1000G_Phase3_frq
# Activate LDSC environment
source ${miniconda_dir}/bin/activate ldsc

###Renaming all files which are generated in the calculate_ld_score.sh script from ${ref_ld_dir}/${dataset}_${ct}_${chr} to ${ref_ld_dir}/${dataset}_${ct}.${chr} (_ to .)
###This is unnecessary to do if files from the calculate_ld_score.sh script are already in the ${ref_ld_dir}/${dataset}_${ct}.${chr} format
###This was fixed in the script for the level 2 KI data.
for ct in ${cell_types[@]}; do
	for chr in {1..22}; do
		mv ${ref_ld_dir}/${dataset}_${ct}_${chr}.l2.ldscore.gz ${ref_ld_dir}/${dataset}_${ct}.${chr}.l2.ldscore.gz
		mv ${ref_ld_dir}/${dataset}_${ct}_${chr}.l2.M ${ref_ld_dir}/${dataset}_${ct}.${chr}.l2.M
		mv ${ref_ld_dir}/${dataset}_${ct}_${chr}.l2.M_5_50 ${ref_ld_dir}/${dataset}_${ct}.${chr}.l2.M_5_50
	done
done

########## Partition heritability for combined cell types across specificity deciles ##########
for gwas in ${sum_stats[@]}; do
        python ${ldsc_dir}/ldsc.py \
                --h2 ${sum_stats_dir}/${gwas}.sumstats.gz \
                --w-ld-chr ${weights_dir}/weights.hm3_noMHC. \
                --ref-ld-chr ${ref_ld_dir}/${dataset}_${cell_types[0]}.,${ref_ld_dir}/${dataset}_${cell_types[1]}.,${ref_ld_dir}/${dataset}_${cell_types[2]}.,${ref_ld_dir}/${dataset}_${cell_types[3]}.,${ref_ld_dir}/${dataset}_${cell_types[4]}.,${ref_ld_dir}/${dataset}_${cell_types[5]}.,${ref_ld_dir}/${dataset}_${cell_types[6]}.,${ref_ld_dir}/${dataset}_${cell_types[7]}.,${ref_ld_dir}/${dataset}_${cell_types[8]}.,${ref_ld_dir}/${dataset}_${cell_types[9]}.,${ref_ld_dir}/${dataset}_${cell_types[10]}.,${ref_ld_dir}/${dataset}_${cell_types[11]}.,${ref_ld_dir}/${dataset}_${cell_types[12]}.,${ref_ld_dir}/${dataset}_${cell_types[13]}.,${ref_ld_dir}/${dataset}_${cell_types[14]}.,${ref_ld_dir}/${dataset}_${cell_types[15]}.,${ref_ld_dir}/${dataset}_${cell_types[16]}.,${ref_ld_dir}/${dataset}_${cell_types[17]}.,${ref_ld_dir}/${dataset}_${cell_types[18]}.,${ref_ld_dir}/${dataset}_${cell_types[19]}.,${ref_ld_dir}/${dataset}_${cell_types[20]}.,${ref_ld_dir}/${dataset}_${cell_types[21]}.,${ref_ld_dir}/${dataset}_${cell_types[22]}.,${ref_ld_dir}/${dataset}_${cell_types[23]}.,${baseline_dir}/baseline. \
                --thin-annot \
                --overlap-annot \
                --frqfile-chr ${frq_dir}/1000G.EUR.QC. \
                --out ${out_dir}/${gwas}_combined_model \
                --print-coefficients
done

