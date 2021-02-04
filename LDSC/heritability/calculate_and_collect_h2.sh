#!bin/bash

##Script to estimate SNP-heritability of phenotypes
##Requirements: Munged sumstats
##By: Mitchell Olislagers
##Last updated: 10 Feb 2020

ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
sumstats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data
munged_dir=${sumstats_dir}/munged_sumstats
h2_dir=/hpc/hers_en/molislagers/LDSC/heritability/analysis_phase3
h2_output_dir=/hpc/hers_en/molislagers/LDSC/heritability/output_phase3
conda activate ldsc

cd ${munged_dir}

phenotypes=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance" "neuroticism")

phenotypes_munged=( "${phenotypes[@]/%/.sumstats.gz}" )

for phenotype in "${phenotypes_munged[@]}"; do
	python ${ldsc_dir}/ldsc.py --h2 ${phenotype} --ref-ld-chr ${ref_dir}/1000G_EUR_Phase3_baseline/baseline. --w-ld-chr ${ref_dir}/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. --out ${h2_dir}/$phenotypes
	phenotypes=(${phenotypes[@]:1})
done

cd ${h2_dir}
#Create files with headers
echo "phenotype" > phen_header.txt
echo "h2" > h2_header.txt
echo "se" > se_header.txt
echo "lambda_gc" > lambda_gc_header.txt
echo "chi2" > chi2_header.txt
echo "intercept" > intercept_header.txt
echo "intercept_se" > intercept_se_header.txt

#Create log file of all phenotypes
phenotypes=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance" "neuroticism")
for phenotype in "${phenotypes[@]}"; do
	grep 'Total Observed' ${phenotype}.log > ${phenotype}_temp.log
	cut -d " " -f 5 ${phenotype}_temp.log >> h2_header.txt
	cut -d " " -f 6 ${phenotype}_temp.log >${phenotype}_se.log
	sed 's|[(),]||g' ${phenotype}_se.log >> se_header.txt
	echo ${phenotype} >> phen_header.txt
	rm ${phenotype}_se.log
	rm ${phenotype}_temp.log
	grep 'Lambda' ${phenotype}.log > ${phenotype}_lambda_line.txt
	cut -d " " -f 3 ${phenotype}_lambda_line.txt >> lambda_gc_header.txt
	rm ${phenotype}_lambda_line.txt
	grep 'Chi' ${phenotype}.log > ${phenotype}_chi_line.txt
	cut -d " " -f 3 ${phenotype}_chi_line.txt >> chi2_header.txt
	rm ${phenotype}_chi_line.txt
	grep 'Intercept' ${phenotype}.log > ${phenotype}_intercept_line.txt
	cut -d " " -f 2 ${phenotype}_intercept_line.txt >> intercept_header.txt
	cut -d " " -f 3 ${phenotype}_intercept_line.txt > ${phenotype}_intercept_se_tmp.txt
	sed 's|[(),]||g' ${phenotype}_intercept_se_tmp.txt >> intercept_se_header.txt
	rm ${phenotype}_intercept_se_tmp.txt
	rm ${phenotype}_intercept_line.txt
done

paste phen_header.txt h2_header.txt se_header.txt lambda_gc_header.txt chi2_header.txt intercept_header.txt intercept_se_header.txt | column -s $'\t' -t > ${h2_output_dir}/all_h2_phase3.log
