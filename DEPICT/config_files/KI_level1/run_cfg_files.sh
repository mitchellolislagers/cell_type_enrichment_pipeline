#!/bin/bash

phenotypes=("ADHD" "age_initiation" "alcohol_dependence" "alcohol_use" "all_epilepsy" "all_stroke" "ALS" "alzheimers" "AN" "anxiety" "ASD" "BIP" "BMI" "cannabis" "cardioembolic_stroke" "chronotype" "cigarettes_pd" "cognitive_performance" "cross" "daytime_sleepiness" "drinks_pw" "educational_attainment" "ever_smoked" "focal_epilepsy" "generalised_epilepsy" "height" "insomnia" "intelligence" "ischemic_stroke" "large_artery_stroke" "long_sleep_duration" "MDD" "neuroticism" "OCD" "overall_sleep_duration" "parkinson" "PTSD" "SCZ" "short_sleep_duration" "small_vessel_stroke" "smoking_cessation" "TS")

# execute DEPICT using custom .cfg file
for phenotype in ${phenotypes[@]}; do python $HOME/DEPICT/src/python/depict.py $HOME/Mitchell/KI_level1/config_scripts/KI_${phenotype}_lvl1.cfg; done
