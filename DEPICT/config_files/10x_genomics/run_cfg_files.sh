#!/bin/bash

phenotypes=("ADHD" "age_initiation" "alcohol_dependence" "alcohol_use" "all_epilepsy" "all_stroke" "ALS" "alzheimers" "AN" "anxiety" "ASD" "BIP" "BMI" "cannabis" "cardioembolic_stroke" "chronotype" "cigarettes_pd" "cognitive_performance" "cross" "daytime_sleepiness" "drinks_pw" "educational_attainment" "ever_smoked" "focal_epilepsy" "generalised_epilepsy" "height" "insomnia" "intelligence" "ischemic_stroke" "large_artery_stroke" "long_sleep_duration" "MDD" "neuroticism" "OCD" "overall_sleep_duration" "parkinson" "PTSD" "SCZ" "short_sleep_duration" "small_vessel_stroke" "smoking_cessation" "TS")

# load modules
module load pre2019
module load Python/2.7.12-foss-2016b
module load plink/1.90b6.9

# execute DEPICT using custom .cfg file
for phenotype in ${phenotypes[@]}; do python $HOME/DEPICT/src/python/depict.py $HOME/Mitchell/10x_genomics/config_scripts/10x_genomics_${phenotype}.cfg; done
