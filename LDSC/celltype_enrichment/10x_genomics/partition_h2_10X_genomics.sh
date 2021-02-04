#!/bin/bash
  
# TITLE:        partition_h2.sh
# GOAL:         Partition heritability (h2) to specificity deciles of cell types.
# INPUT:        LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types, LDSC LD score files 
# INPUT:	1000 Genomes Project phase 3 baseline model LD scores, regression weights, allele frequencies (see README)
# OUTPUT:       Partitioned h2 output files
# AUTHOR:       Koen Rademaker
# DATE:         21 February 2020


########## Initialize script ##########
# Set paths
ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment
sum_stats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics/munged_sumstats
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/10X_genomics/annotations
ld_scores_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/10X_genomics/ld_scores
out_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/10X_genomics/partitioned_h2
miniconda_dir=~/miniconda2
# Declare variables
declare -a sum_stats=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance" "neuroticism")

declare -a cell_types=("Glutamatergic_neurons" "Neuroblasts_1" "Astrocytes_1" "Neuroblasts_2" "Intermediate_progenitors" "Enteric_glial_cells" "Interneurons_1" "Neurons" 
"Interneurons_2" "Vascular_endothelial_cells" "Enteric_neurons" "Cajal_Retzius_cells" "Oligodendrocytes" "Astrocytes_2" "Endothelial_cells" "Microglia")
dataset="10x_Genomics"


########## Organize LDSC ##########
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz -P ${ref_dir}
#tar -xvzf 1000G_Phase3_frq.tgz
# Copy annotation files and LD score files to reference folder
ref_ld_dir=${ref_dir}/ref_ld_10x
cp ${ld_scores_dir}/*l2* ${ref_ld_dir}
cp ${annotation_dir}/*.annot.gz ${ref_ld_dir}
# Declare variables
weights_dir=/hpc/hers_en/molislagers/LDSC/ref_data/regression/1000G_Phase3_weights_hm3_no_MHC
baseline_dir=/hpc/hers_en/molislagers/LDSC/ref_data/regression/1000G_EUR_Phase3_baseline
frq_dir=${ref_dir}/1000G_Phase3_frq
# Activate LDSC environment
source ${miniconda_dir}/bin/activate ldsc


########## Partition heritability for combined cell types across specificity deciles ##########
for gwas in ${sum_stats[@]}; do
        python ${ldsc_dir}/ldsc.py \
                --h2 ${sum_stats_dir}/${gwas}.sumstats.gz \
                --w-ld-chr ${weights_dir}/weights.hm3_noMHC. \
                --ref-ld-chr ${ref_ld_dir}/${dataset}_${cell_types[0]}.,${ref_ld_dir}/${dataset}_${cell_types[1]}.,${ref_ld_dir}/${dataset}_${cell_types[2]}.,${ref_ld_dir}/${dataset}_${cell_types[3]}.,${ref_ld_dir}/${dataset}_${cell_types[4]}.,${ref_ld_dir}/${dataset}_${cell_types[5]}.,${ref_ld_dir}/${dataset}_${cell_types[6]}.,${ref_ld_dir}/${dataset}_${cell_types[7]}.,${ref_ld_dir}/${dataset}_${cell_types[8]}.,${ref_ld_dir}/${dataset}_${cell_types[9]}.,${ref_ld_dir}/${dataset}_${cell_types[10]}.,${ref_ld_dir}/${dataset}_${cell_types[11]}.,${ref_ld_dir}/${dataset}_${cell_types[12]}.,${ref_ld_dir}/${dataset}_${cell_types[13]}.,${ref_ld_dir}/${dataset}_${cell_types[14]}.,${ref_ld_dir}/${dataset}_${cell_types[15]}.,${baseline_dir}/baseline. \
                --thin-annot \
                --overlap-annot \
                --frqfile-chr ${frq_dir}/1000G.EUR.QC. \
                --out ${out_dir}/${gwas}_combined_model \
                --print-coefficients
done

