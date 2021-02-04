#!bin/bash

##Script to estimate genetic correlations among phenotypes
##Requirements: Munged sumstats
##By: Mitchell Olislagers
##Last updated: 10 Feb 2020

ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
sumstats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/regression
munged_dir=${sumstats_dir}/munged_sumstats
rg_dir=/hpc/hers_en/molislagers/LDSC/bivariate_correlations/analysis_phase3
rg_output_dir=/hpc/hers_en/molislagers/LDSC/bivariate_correlations/output_phase3
conda activate ldsc

cd ${munged_dir}
phenotypes=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance" "neuroticism")
phenotypes_munged=( "${phenotypes[@]/%/.sumstats.gz}" )

for phenotype in "${phenotypes_munged[@]}"; do
	#Break the loop if only 1 variable in array
	if [ "${#phenotypes_munged[@]}" -eq 1 ]; then
		break
	fi
	#Join array by comma
	rg_sumstats=$(printf ",%s" "${phenotypes_munged[@]}")
	rg_sumstats=(${rg_sumstats:1})
	#Run bivariate correlations
	python ${ldsc_dir}/ldsc.py --rg $rg_sumstats --ref-ld-chr ${ref_dir}/1000G_EUR_Phase3_baseline/baseline. --w-ld-chr ${ref_dir}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --out ${rg_dir}/$phenotypes
	#Remove first variable of array
	phenotypes_munged=(${phenotypes_munged[@]:1})
	phenotypes=(${phenotypes[@]:1})
done

python ${ldsc_dir}/ldsc.py --rg smoking_initiation.sumstats.gz,ever_smoked.sumstats.gz --ref-ld-chr ${ref_dir}/1000G_EUR_Phase3_baseline/baseline. --w-ld-chr ${ref_dir}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --out ${rg_dir}/test.log


cd ${rg_dir}
phenotypes=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance")
#Create header file 
awk '/Summary/{y=1;next}y' ${phenotypes}.log > ${phenotypes}_no_top.log
head -n 1 ${phenotypes}_no_top.log > header_file.txt
rm ${phenotypes}_no_top.log

for phenotype in "${phenotypes[@]}"; do
	#Only select summary of results
	awk '/Summary/{y=1;next}y' ${phenotype}.log > ${phenotype}_no_top.log
	sed -i -e "1d" ${phenotype}_no_top.log
	head -n -3 ${phenotype}_no_top.log > ${phenotype}_no_top_no_bottom.log
	#Append to header file
	cat ${phenotype}_no_top_no_bottom.log >> header_file.txt
done
cp ${rg_dir}/header_file.txt ${rg_output_dir}/all_rg_phase3.log

