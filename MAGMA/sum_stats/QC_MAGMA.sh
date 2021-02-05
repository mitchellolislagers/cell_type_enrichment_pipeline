!#bin/bash

# Goal: Restrict summary statistics to hapmap3 SNPs and exclude the MHC region (chr6:25-34 Mb)
# Required: MAGMA summary statistics, as generated in the prepare_MAGMA_sumstats.sh script
# By: Mitchell Olislagers
# Last updated: 14 May 2020

input_sums=/hpc/hers_en/molislagers/MAGMA/sumstats/MAGMA_sumstats
hm3_dir=/hpc/hers_en/molislagers/MAGMA/sumstats/QC_MAGMA/hm3_restricted
hm3_mhc_dir=/hpc/hers_en/molislagers/MAGMA/sumstats/QC_MAGMA/hm3_restricted_MHC_excluded
mhc_dir=/hpc/hers_en/molislagers/MAGMA/sumstats/QC_MAGMA/MHC_excluded
ref_dir=/hpc/hers_en/molislagers/MAGMA/ref_data
temp_dir=/hpc/hers_en/molislagers/MAGMA/sumstats/QC_MAGMA/temp
mhc_bp_start=25000000
mhc_bp_end=34000000

phenotypes=("ADHD" "age_initiation" "alcohol_dependence" "alcohol_use" "all_epilepsy" "all_stroke" "ALS" "alzheimers" "AN" "anxiety" "ASD" "BIP" "BMI" "cannabis" "cardioembolic_stroke" "chronotype" "cigarettes_pd" "cognitive_performance" "cross" "daytime_sleepiness" "drinks_pw" "educational_attainment" "ever_smoked" "focal_epilepsy" "generalised_epilepsy" "height" "insomnia" "intelligence" "ischemic_stroke" "large_artery_stroke" "long_sleep_duration" "MDD" "neuroticism" "OCD" "overall_sleep_duration" "parkinson" "PTSD" "SCZ" "short_sleep_duration" "small_vessel_stroke" "smoking_cessation" "TS")

#Remove log files if they exist
if [ -f ${hm3_dir}/hm3_QC.log ]; then
	rm ${hm3_dir}/hm3_QC.log
	rm ${hm3_mhc_dir}/hm3_mhc_QC.log
	rm ${mhc_dir}/mhc_QC.log
fi


## QC TO RESTRICT TO HM3 SNPS AND EXCLUDE MHC REGION (CHR6:25-34 MB) ##
#Loop over phenotypes
for phenotype in ${phenotypes[@]}; do
#Restrict to Hapmap3 SNPs 
awk -v OFS="\t" 'FNR==NR{a[$1]=$1;next} ($1 in a) {print $0}' ${ref_dir}/w_hm3.snplist ${input_sums}/${phenotype}_MAGMA.txt > ${hm3_dir}/${phenotype}_hm3_restricted_MAGMA.txt
#Create file that contains all SNPs within MHC region from Hapmap3 restricted GWAS
awk -v OFS="\t" '$2==6' ${hm3_dir}/${phenotype}_hm3_restricted_MAGMA.txt | awk '$3 >= 25000000 && $3 <= 34000000' > ${temp_dir}/${phenotype}_rows_to_remove.txt
#Remove SNPs within the MHC region from Hapmap3 restricted GWAS
awk -v OFS="\t" 'FNR==NR{a[$1]=$1;next} !($1 in a) {print $0}' ${temp_dir}/${phenotype}_rows_to_remove.txt ${hm3_dir}/${phenotype}_hm3_restricted_MAGMA.txt > ${hm3_mhc_dir}/${phenotype}_hm3_mhc_restricted_MAGMA.txt

##Create file that contains all SNPs within MHC region from raw GWAS
awk -v OFS="\t" '$2==6' ${input_sums}/${phenotype}_MAGMA.txt | awk '$3 >= 25000000 && $3 <= 34000000' > ${temp_dir}/${phenotype}_rows_to_remove.txt
#Remove SNPs within the MHC region from raw GWAS
awk -v OFS="\t" 'FNR==NR{a[$1]=$1;next} !($1 in a) {print $0}' ${temp_dir}/${phenotype}_rows_to_remove.txt ${input_sums}/${phenotype}_MAGMA.txt > ${mhc_dir}/${phenotype}_mhc_restricted_MAGMA.txt
#Clean up rows to remove file
rm ${temp_dir}/${phenotype}_rows_to_remove.txt

#Create files to make a log file
echo ${phenotype} >> ${temp_dir}/GWAS_name.txt
sed -e '1d' ${input_sums}/${phenotype}_MAGMA.txt | wc -l >> ${temp_dir}/snps_before_QC.txt
sed -e '1d' ${hm3_dir}/${phenotype}_hm3_restricted_MAGMA.txt | wc -l  >> ${temp_dir}/snps_after_hm3_QC.txt
sed -e '1d' ${hm3_mhc_dir}/${phenotype}_hm3_mhc_restricted_MAGMA.txt | wc -l >> ${temp_dir}/snps_after_hm3_mhc_QC.txt
wc -l < ${mhc_dir}/${phenotype}_mhc_restricted_MAGMA.txt >> ${temp_dir}/snps_after_mhc_QC.txt
done

#Paste all log files together to create final hm3 log file
paste ${temp_dir}/GWAS_name.txt ${temp_dir}/snps_before_QC.txt ${temp_dir}/snps_after_hm3_QC.txt > ${temp_dir}/no_headers_hm3_QC.txt
#Paste all log files together to create final hm3 and MHC log file
paste ${temp_dir}/GWAS_name.txt ${temp_dir}/snps_before_QC.txt ${temp_dir}/snps_after_hm3_mhc_QC.txt > ${temp_dir}/no_headers_hm3_mhc_QC.txt
#Paste all log files together to create final MHC log file
paste ${temp_dir}/GWAS_name.txt ${temp_dir}/snps_before_QC.txt ${temp_dir}/snps_after_mhc_QC.txt > ${temp_dir}/no_headers_mhc_QC.txt

#Create header for hm3 log file
echo -e "phenotype\tsnps_before_qc\tsnps_after_hm3_restriction" | cat - ${temp_dir}/no_headers_hm3_QC.txt > ${temp_dir}/headers_hm3_QC.txt
#Create header for hm3 and MHC log file
echo -e "phenotype\tsnps_before_qc\tsnps_after_hm3_and_mhc_restriction" | cat - ${temp_dir}/no_headers_hm3_mhc_QC.txt > ${temp_dir}/headers_hm3_mhc_QC.txt
#Create header for MHC log file
echo -e "phenotype\tsnps_before_qc\tsnps_after_mhc_restriction" | cat - ${temp_dir}/no_headers_mhc_QC.txt > ${temp_dir}/headers_mhc_QC.txt
#Add 'snps_removed' column in hm3 log file
awk 'NR == 1 { $4 = "snps_removed" } NR >= 2 { $5 = $2 - $3 } 1' < ${temp_dir}/headers_hm3_QC.txt >> ${hm3_dir}/hm3_QC.log
#Add 'snps_removed' column in hm3 and MHC log file
awk 'NR == 1 { $4 = "snps_removed" } NR >= 2 { $5 = $2 - $3 } 1' < ${temp_dir}/headers_hm3_mhc_QC.txt >> ${hm3_mhc_dir}/hm3_mhc_QC.log
#Add 'snps_removed' column in mhc log file
awk 'NR == 1 { $4 = "snps_removed" } NR >= 2 { $5 = $2 - $3 } 1' < ${temp_dir}/headers_mhc_QC.txt >> ${mhc_dir}/mhc_QC.log

#Clean up files
rm ${temp_dir}/*
