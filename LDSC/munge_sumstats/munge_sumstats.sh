!#bin/bash

# Goal: Munge GWAs summary statistics for LDSC analyses
# Required:	GWAS summary statistics files, working LDSC environment
# By: Mitchell Olislagers
# Last updated: Feb 10 2020

qlogin -l h_rt=08:00:00 
## Set up environment and variables
ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
sumstats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/regression
munged_dir=${sumstats_dir}/munged_sumstats
h2_dir=/hpc/hers_en/molislagers/LDSC/heritability
rg_dir=/hpc/hers_en/molislagers/LDSC/bivariate_correlations

conda activate ldsc

##Download HapMap3, regression weights and baseline LD scores
cd ${ref_dir}
##Download SNPLIST
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2
bunzip2 w_hm3.snplist.bz2

##Download regression weights and baseline LD scores
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_weights_hm3_no_MHC.tgz
tar -xvzf 1000G_Phase3_weights_hm3_no_MHC.tgz
wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_baseline_ldscores.tgz
tar -xvzf 1000G_Phase3_baseline_ldscores.tgz

## Schizophrenia
cd ${sumstats_dir}/psychiatric_diseases/SCZ/
gunzip ${sumstats_dir}/psychiatric_diseases/SCZ/clozuk_pgc2.meta.sumstats.txt.gz
##Preprocessing file
#Make file with header
head -n 1 ${sumstats_dir}/psychiatric_diseases/SCZ/raw_GWAS/clozuk_pgc2.meta.sumstats.txt > ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/tmp_clozuk_header
#Select all SNPs with rsnumber
grep 'rs' ${sumstats_dir}/psychiatric_diseases/SCZ/raw_GWAS/clozuk_pgc2.meta.sumstats.txt > ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/clozuk_pgc2.meta.sumstats.only_rs.txt
#Create extra column with only rs#
awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print $0, s}' ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/clozuk_pgc2.meta.sumstats.only_rs.txt > ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/clozuk_pgc2.meta.sumstats.only_rs.rs_col.txt
sed -i 's/\r//' ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/clozuk_pgc2.meta.sumstats.only_rs.rs_col.txt
sed -i 's/\r//' ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/tmp_clozuk_header
#Create new header with SNP_RS columnname
sed 's/$/\tSNP_RS/' ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/tmp_clozuk_header > ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/tmp_clozuk_header_2
#Merge new header with file containing rs# column
cat ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/tmp_clozuk_header_2 ${sumstats_dir}/psychiatric_diseases/SCZ/processing_sumstats/clozuk_pgc2.meta.sumstats.only_rs.rs_col.txt > ${sumstats_dir}/psychiatric_diseases/SCZ/final_sumstats/clozuk_pgc2.meta.sumstats.only_rs.rs_col.header.txt
#Munge SCZ
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/SCZ/final_sumstats/clozuk_pgc2.meta.sumstats.only_rs.rs_col.header.txt --out ${munged_dir}/SCZ --merge-alleles ${ref_dir}/w_hm3.snplist --N 105317 --snp SNP_RS --ignore SNP --p P --frq Freq.A1

##Anorexia
cd ${sumstats_dir}/psychiatric_diseases/AN/
gunzip pgcAN2.2019-07.vcf.tsv.gz
##Preprocess file
#Only select summary statistics
tail -n +71 ${sumstats_dir}/psychiatric_diseases/AN/raw_GWAS/pgcAN2.2019-07.vcf.tsv > ${sumstats_dir}/psychiatric_diseases/AN/processing_sumstats/only_sumstats_pgcAN2.2019-07.vcf.tsv
awk '$0=$0"\t"(NR==1?"N":$12+$13)' ${sumstats_dir}/psychiatric_diseases/AN/processing_sumstats/only_sumstats_pgcAN2.2019-07.vcf.tsv > ${sumstats_dir}/psychiatric_diseases/AN/final_sumstats/only_sumstats_ntot_pgcAN2.2019-07.vcf.tsv
#Munge AN
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/AN/final_sumstats/only_sumstats_ntot_pgcAN2.2019-07.vcf.tsv --out ${munged_dir}/AN --merge-alleles ${ref_dir}/w_hm3.snplist --snp ID --a1 ALT --a2 REF --p PVAL --info IMPINFO

##ADHD
cd ${sumstats_dir}/psychiatric_diseases/ADHD/
gunzip adhd_eur_jun2017.gz
#Munge ADHD
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/ADHD/adhd_eur_jun2017 --out ${munged_dir}/ADHD --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --info INFO --p P --N-cas 19099 --N-con 34194

##Anxiety disorder
cd ${sumstats_dir}/psychiatric_diseases/anxiety_disorder/angst.study.results
gunzip anxiety.meta.full.cc.tbl.gz
#Change column name of TotalN to N
sed -e '1s/TotalN/N/' anxiety.meta.full.cc.tbl > n_anxiety.meta.full.cc.tbl
#Munge anxiety disorder
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/anxiety_disorder/angst.study.results/n_anxiety.meta.full.cc.tbl --out ${munged_dir}/anxiety --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNPID --a1 Allele1 --a2 Allele2 --frq Freq1

##ASD
cd ${sumstats_dir}/psychiatric_diseases/ASD
gunzip iPSYCH-PGC_ASD_Nov2017.gz
#Munge ASD
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/ASD/iPSYCH-PGC_ASD_Nov2017 --out ${munged_dir}/ASD --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --info INFO --N-cas 18382 --N-con 27969

##Bipolar disorder
cd ${sumstats_dir}/psychiatric_diseases/BIP
gunzip daner_PGC_BIP32b_mds7a_0416a.gz
#Create column with N
awk '$0=$0"\t"(NR==1?"N":$17+$18)' ${sumstats_dir}/psychiatric_diseases/BIP/daner_PGC_BIP32b_mds7a_0416a > ${sumstats_dir}/psychiatric_diseases/BIP/ntot_daner_PGC_BIP32b_mds7a_0416a
#Munge bipolar disorder
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/BIP/ntot_daner_PGC_BIP32b_mds7a_0416a --out ${munged_dir}/BIP --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --info INFO

##Cross-disorders
cd ${sumstats_dir}/psychiatric_diseases/cross_disorders
gunzip pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt.gz
#Make column with total N
awk '$0=$0"\t"(NR==1?"N":$14+$15)' ${sumstats_dir}/psychiatric_diseases/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt > ${sumstats_dir}/psychiatric_diseases/cross_disorders/ntot_pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt
#Munge cross disorders
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/cross_disorders/ntot_pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt --out ${munged_dir}/cross --merge-alleles ${ref_dir}/w_hm3.snplist --snp ID --a1 ALT --a2 REF

##Major depression disorder
cd ${sumstats_dir}/psychiatric_diseases/MDD/DS_10283_3203
#Change column name of log odds ratio to LOG_ODDS
sed -e '1s/LogOR/LOG_ODDS/' ${sumstats_dir}/psychiatric_diseases/MDD/DS_10283_3203/PGC_UKB_depression_genome-wide.txt > ${sumstats_dir}/psychiatric_diseases/MDD/DS_10283_3203/logOR_PGC_UKB_depression_genome-wide.txt
#Munge major depression disorder
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/MDD/DS_10283_3203/logOR_PGC_UKB_depression_genome-wide.txt --out ${munged_dir}/MDD --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 A1 --a2 A2 --frq Freq --N-cas 170756 --N-con 329443

##Obsessive compulsive disorder
cd ${sumstats_dir}/psychiatric_diseases/OCD/PGC_OCD_Aug2017
gunzip ocd_aug2017.gz
#Munge obsessive compulsive disorder
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/OCD/PGC_OCD_Aug2017/ocd_aug2017 --out ${munged_dir}/OCD --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --info INFO --N-cas 2688 --N-con 7037

##Post-traumatic stress disorder
cd ${sumstats_dir}/psychiatric_diseases/PTSD
gunzip pts_eur_freeze2_overall.results.gz
#Create total N column
awk '$0=$0"\t"(NR==1?"N":$17+$18)' ${sumstats_dir}/psychiatric_diseases/PTSD/pts_eur_freeze2_overall.results > ${sumstats_dir}/psychiatric_diseases/PTSD/ntot_pts_eur_freeze2_overall.results
#Munge post-traumatic stress disorder
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/PTSD/ntot_pts_eur_freeze2_overall.results --out ${munged_dir}/PTSD --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --info INFO

##Tourette
cd ${sumstats_dir}/psychiatric_diseases/TS
gunzip TS_Oct2018.gz
#Munge Tourette
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/psychiatric_diseases/TS/TS_Oct2018 --out ${munged_dir}/TS --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --info INFO --N-cas 4819 --N-con 9488

##Drinks per week
cd ${sumstats_dir}/substance_use/drinks_pw
gunzip DrinksPerWeek.txt.gz
#Munge drinks per week
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/drinks_pw/DrinksPerWeek.txt --out ${munged_dir}/drinks_pw --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 ALT --a2 REF

##Cannabis use
cd ${sumstats_dir}/substance_use/cannabis
gunzip Cannabis_ICC_UKB_het.txt.gz
#Create total N column
awk '$0=$0"\t"(NR==1?"N":$15+$16)' ${sumstats_dir}/substance_use/cannabis/Cannabis_ICC_UKB_het.txt > ${sumstats_dir}/substance_use/cannabis/ntot_Cannabis_ICC_UKB_het.txt
#Munge cannabis use (Ignore Z column as there is a BETA column)
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/cannabis/ntot_Cannabis_ICC_UKB_het.txt --out ${munged_dir}/cannabis --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --frq FRQ --ignore Z

##Alcohol use disorder
cd ${sumstats_dir}/substance_use/alcohol_use
gunzip AUDIT_UKB_2018_AJP.txt.gz
##Preprocessing file
#Select all SNPs with rsnumber
grep 'rs' ${sumstats_dir}/substance_use/alcohol_use/AUDIT_UKB_2018_AJP.txt > ${sumstats_dir}/substance_use/alcohol_use/AUDIT_UKB_2018_AJP_only_rs.txt
#Change column name of beta_T to BETA
sed -e '1s/beta_T/BETA/' ${sumstats_dir}/substance_use/alcohol_use/AUDIT_UKB_2018_AJP_only_rs.txt > ${sumstats_dir}/substance_use/alcohol_use/beta_AUDIT_UKB_2018_AJP_only_rs.txt
#Munge alcohol use disorder
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/alcohol_use/beta_AUDIT_UKB_2018_AJP_only_rs.txt --out ${munged_dir}/alcohol_use --merge-alleles ${ref_dir}/w_hm3.snplist --snp rsid --a1 a_1 --a2 a_0 --info info --p P_T

##Alcohol dependence
cd ${sumstats_dir}/substance_use/alcohol_dependence
tar -xzvf pgc_alcdep.aug2018_release.tar.gz
gunzip pgc_alcdep.eur_discovery.aug2018_release.txt.gz
##Preprocessing file
#Make file with header
head -n 1 ${sumstats_dir}/substance_use/alcohol_dependence/input_GWAS/pgc_alcdep.eur_discovery.aug2018_release.txt > ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/header_pgc_alcdep.eur_discovery.aug2018_release.txt
#Select all SNPs with rsnumber
grep 'rs' ${sumstats_dir}/substance_use/alcohol_dependence/input_GWAS/pgc_alcdep.eur_discovery.aug2018_release.txt > ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs.txt
#Create extra column with only rs#
awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print $0, s}' ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs.txt > ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt
sed -i 's/\r//' ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt
sed -i 's/\r//' ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/header_pgc_alcdep.eur_discovery.aug2018_release.txt
#Create new header with SNP_RS columnname
sed -i 's/$/ wrong_col/' ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/header_pgc_alcdep.eur_discovery.aug2018_release.txt
sed 's/$/ SNP_RS/' ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/header_pgc_alcdep.eur_discovery.aug2018_release.txt > ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/header2_pgc_alcdep.eur_discovery.aug2018_release.txt
#Merge new header with file containing rs# column
cat ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/header2_pgc_alcdep.eur_discovery.aug2018_release.txt ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt > ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/wrong_col_pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt
#Remove wrong column
awk '{$9=""; print $0}' ${sumstats_dir}/substance_use/alcohol_dependence/processing_sumstats/wrong_col_pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt > ${sumstats_dir}/substance_use/alcohol_dependence/final_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt
#Munge alcohol dependence
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/alcohol_dependence/final_sumstats/pgc_alcdep.eur_discovery.aug2018_release_only_rs_rs_col.txt --out ${munged_dir}/alcohol_dependence --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP_RS --N-cas 11569 --N-con 34999 --ignore SNP,Weight --a1 A1 --a2 A2 --p P

##Age smoking intitiation
cd ${sumstats_dir}/substance_use/smoking_initiation
gunzip AgeofInitiation.txt.gz
#Munge age smoking initiation
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/smoking_initiation/AgeofInitiation.txt --out ${munged_dir}/smoking_initiation --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 ALT --a2 REF

##Ever smoked regularly
cd ${sumstats_dir}/substance_use/ever_smoked
gunzip SmokingInitiation.txt.gz
#Munge ever smoked regularly
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/ever_smoked/SmokingInitiation.txt --out ${munged_dir}/ever_smoked --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 ALT --a2 REF


##Cigarettes per day
cd ${sumstats_dir}/substance_use/cigarettes_pd
gunzip CigarettesPerDay.txt.gz
#Munge cigarettes per day
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/cigarettes_pd/CigarettesPerDay.txt --out ${munged_dir}/cigarettes_pd --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 ALT --a2 REF


##Smoking cessation
cd ${sumstats_dir}/substance_use/smoking_cessation
gunzip SmokingCessation.txt.gz
#Munge smoking cessation
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/substance_use/smoking_cessation/SmokingCessation.txt --out ${munged_dir}/smoking_cessation --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 ALT --a2 REF

##Reference data
#Prepare reference data
awk '{print $1":"$4 "\t" $2}' ${ref_dir}/SNPref/1kg_rs.bim > ${ref_dir}/SNPref/new_1kg_rs.bim
awk '{print $1":"$2 "\t" $3}' ${ref_dir}/SNPref/rsids > ${ref_dir}/SNPref/new_rsids

##Parkinson
#Remove 'chr' from positions
sed 's/chr//g' ${sumstats_dir}/neurological_disorders/parkinson/nallsEtAl2019_excluding23andMe_allVariants.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
#Replace SNP with position
sed -e '1s/SNP/position/' ${sumstats_dir}/neurological_disorders/parkinson/processing_data/nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
#Grep header
head -n 1 ${sumstats_dir}/neurological_disorders/parkinson/processing_data/pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/header_pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
#Add column names position and RSID in header
sed 's/$/\tpos\tRSID/' ${sumstats_dir}/neurological_disorders/parkinson/processing_data/header_pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_header_pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
#Merge chr:pos column to RSIDs in reference file new_1_kg_rs.bim
awk 'FNR==NR{a[$1]=$0;next}{ print $0 "\t" a[$1]}' ${ref_dir}/SNPref/new_1kg_rs.bim ${sumstats_dir}/neurological_disorders/parkinson/processing_data/pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab 
awk 'FNR==NR{a[$1]=$0;next}{ print $0 "\t" a[$1]}' ${ref_dir}/SNPref/new_rsids ${sumstats_dir}/neurological_disorders/parkinson/processing_data/nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab 
#Remove header and replace with new header
sed -e "1d" ${sumstats_dir}/neurological_disorders/parkinson/processing_data/1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/no_header_1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
sed -e "1d" ${sumstats_dir}/neurological_disorders/parkinson/processing_data/rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/no_header_rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab 
cat ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_header_pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab ${sumstats_dir}/neurological_disorders/parkinson/processing_data/no_header_1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
cat ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_header_pos_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab ${sumstats_dir}/neurological_disorders/parkinson/processing_data/no_header_rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab  > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab 
#Change column b to BETA
sed -e '1s/b/BETA/' ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/beta_final_1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
sed -e '1s/b/BETA/' ${sumstats_dir}/neurological_disorders/parkinson/processing_data/final_rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/processing_data/beta_final_rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab 
#Create total N column
awk '$0=$0"\t"(NR==1?"N":$8+$9)' ${sumstats_dir}/neurological_disorders/parkinson/processing_data/beta_final_1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/final_data/1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab
awk '$0=$0"\t"(NR==1?"N":$8+$9)' ${sumstats_dir}/neurological_disorders/parkinson/processing_data/beta_final_rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab > ${sumstats_dir}/neurological_disorders/parkinson/final_data/rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab 

#Munge parkinson; 1kg: 1137731 SNPs remain; rsids: 1136406 SNPs remain --> 1kg used for analysis
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/parkinson/final_data/1kg_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab --out ${munged_dir}/parkinson --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 A1 --a2 A2 --frq freq 
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/parkinson/final_data/rsids_nallsEtAl2019_excluding23andMe_allVariants_nochrom.tab --out ${munged_dir}/parkinson --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 A1 --a2 A2 --frq freq

##ALS
#Grep header 
head -n 1 ${sumstats_dir}/neurological_disorders/ALS/alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/processing_data/header_alsMetaSummaryStats_march21st2018.tab
#Add column names position and RSID in header
sed 's/$/\tpos\tRSID/' ${sumstats_dir}/neurological_disorders/ALS/processing_data/header_alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/processing_data/final_header_alsMetaSummaryStats_march21st2018.tab
#Merge chr:pos column to RSIDs in reference file new_1_kg_rs.bim
awk 'FNR==NR{a[$1]=$0;next}{ print $0 "\t" a[$1]}' ${ref_dir}/SNPref/new_1kg_rs.bim ${sumstats_dir}/neurological_disorders/ALS/alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/processing_data/merged_1kg_alsMetaSummaryStats_march21st2018.tab
awk 'FNR==NR{a[$1]=$0;next}{ print $0 "\t" a[$1]}' ${ref_dir}/SNPref/new_rsids ${sumstats_dir}/neurological_disorders/ALS/alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/processing_data/merged_rsids_alsMetaSummaryStats_march21st2018.tab
#Remove header and replace with new header
sed -e "1d" ${sumstats_dir}/neurological_disorders/ALS/processing_data/merged_1kg_alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/processing_data/no_header_merged_1kg_alsMetaSummaryStats_march21st2018.tab
sed -e "1d" ${sumstats_dir}/neurological_disorders/ALS/processing_data/merged_rsids_alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/processing_data/no_header_merged_rsids_alsMetaSummaryStats_march21st2018.tab
cat ${sumstats_dir}/neurological_disorders/ALS/processing_data/final_header_alsMetaSummaryStats_march21st2018.tab ${sumstats_dir}/neurological_disorders/ALS/processing_data/no_header_merged_1kg_alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/final_data/1kg_alsMetaSummaryStats_march21st2018.tab
cat ${sumstats_dir}/neurological_disorders/ALS/processing_data/final_header_alsMetaSummaryStats_march21st2018.tab ${sumstats_dir}/neurological_disorders/ALS/processing_data/no_header_merged_rsids_alsMetaSummaryStats_march21st2018.tab > ${sumstats_dir}/neurological_disorders/ALS/final_data/rsids_alsMetaSummaryStats_march21st2018.tab
#Munge ALS; 1kg: 1191988 SNPs remain; rsids: 1190217 SNPs remain --> 1kg used for analysis
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/ALS/final_data/1kg_alsMetaSummaryStats_march21st2018.tab --out ${munged_dir}/ALS --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --ignore SNP --a1 Allele1 --a2 Allele2 --frq effectAlleleFreq --N-cas 20806 --N-con 59804
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/ALS/final_data/rsids_alsMetaSummaryStats_march21st2018.tab --out ${munged_dir}/ALS --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --ignore SNP --a1 Allele1 --a2 Allele2 --frq effectAlleleFreq --N-cas 20806 --N-con 59804

##Alzheimers
cd ${sumstats_dir}/neurological_disorders/alzheimers
gunzip AD_sumstats_Jansenetal_2019sept.txt.gz
#Change column name of Nsum to N
sed -e '1s/Nsum/N/' ${sumstats_dir}/neurological_disorders/alzheimers/AD_sumstats_Jansenetal_2019sept.txt > ${sumstats_dir}/neurological_disorders/alzheimers/N_AD_sumstats_Jansenetal_2019sept.txt
#Munge Alzheimers (ignore Z as there is already a BETA column)
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/alzheimers/N_AD_sumstats_Jansenetal_2019sept.txt --out ${munged_dir}/alzheimers --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --frq EAF --ignore Z

##all eilepsy
cd ${sumstats_dir}/neurological_disorders/epilepsy/all_epilepsy
gunzip all_epilepsy_METAL.gz
#Munge all epilepsy
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/epilepsy/all_epilepsy/all_epilepsy_METAL --out ${munged_dir}/all_epilepsy --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 15212 --N-con 29677 --ignore Weight

##Generalized epilepsy
cd ${sumstats_dir}/neurological_disorders/epilepsy/generalized
gunzip generalised_epilepsy_METAL.gz
#Munge generalized epilepsy
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/epilepsy/generalized/generalised_epilepsy_METAL --out ${munged_dir}/generalized --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 3769 --N-con 29677 --ignore Weight

##Focal epilepsy
cd ${sumstats_dir}/neurological_disorders/epilepsy/focal
gunzip focal_epilepsy_METAL.gz
#Munge focal epilepsy
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/epilepsy/focal/focal_epilepsy_METAL --out ${munged_dir}/focal --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 9671 --N-con 29677 --ignore Weight

## All stroke
cd ${sumstats_dir}/neurological_disorders/stroke/all_stroke
#Munge all stroke
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/stroke/all_stroke/MEGASTROKE.1.AS.EUR.out --out ${munged_dir}/all_stroke --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 40585 --N-con 406111

##Cardioembolic stroke
cd ${sumstats_dir}/neurological_disorders/stroke/cardioembolic
#Munge carioembolic stroke
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/stroke/cardioembolic/MEGASTROKE.4.CES.EUR.out --out ${munged_dir}/cardioembolic --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 34217 --N-con 406111

##Ischemic stroke
cd ${sumstats_dir}/neurological_disorders/stroke/ischemic
#Munge ischemic stroke
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/stroke/ischemic/MEGASTROKE.2.AIS.EUR.out --out ${munged_dir}/ischemic --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 4373 --N-con 406111

##Large artery stroke
cd ${sumstats_dir}/neurological_disorders/stroke/large_artery
#Munge large artery stroke
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/stroke/large_artery/MEGASTROKE.3.LAS.EUR.out --out ${munged_dir}/large_artery --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 7193 --N-con 406111

##Small vessel stroke
cd ${sumstats_dir}/neurological_disorders/stroke/small_vessel
#Munge small vessel stroke
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/neurological_disorders/stroke/small_vessel/MEGASTROKE.5.SVS.EUR.out --out ${munged_dir}/small_vessel --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 Allele1 --a2 Allele2 --frq Freq1 --N-cas 5386 --N-con 406111

##Height
cd ${sumstats_dir}/behaviour/height
gunzip Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt.gz
#munge height
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/height/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt --out ${munged_dir}/height --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 Tested_Allele --a2 Other_Allele --frq Freq_Tested_Allele_in_HRS

##BMI
cd ${sumstats_dir}/behaviour/BMI
gunzip Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt.gz
#Munge BMI
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/BMI/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt --out ${munged_dir}/BMI --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 Tested_Allele --a2 Other_Allele --frq Freq_Tested_Allele_in_HRS

##Chronotype
cd ${sumstats_dir}/behaviour/chronotype
gunzip chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt.gz
#Munge chronotype
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/chronotype/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt --out ${munged_dir}/chronotype --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 ALLELE1 --a2 ALLELE0 --info INFO --frq A1FREQ --p P_BOLT_LMM --N 449734

##Excessive daytime sleepiness
cd ${sumstats_dir}/behaviour/daytime_sleepiness
unzip Saxena.fullUKBB.DaytimeSleepiness.sumstats.zip
#Munge excessive daytime sleepiness
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/daytime_sleepiness/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt --out ${munged_dir}/daytime_sleepiness --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --N 452071 --info INFO

##Overall sleep duration
cd ${sumstats_dir}/behaviour/sleep_duration/overall_sleep_duration
unzip sleepdurationsumstats.txt.zip
#Change BETA_SLEEPDURATION and P_SLEEPDURATION to BETA and P, respectively
sed -e '1s/BETA_SLEEPDURATION/BETA/' ${sumstats_dir}/behaviour/sleep_duration/overall_sleep_duration/sleepdurationsumstats.txt > ${sumstats_dir}/behaviour/sleep_duration/overall_sleep_duration/P_BETA_sleepdurationsumstats.txt
sed -i -e '1s/P_SLEEPDURATION/P/' ${sumstats_dir}/behaviour/sleep_duration/overall_sleep_duration/P_BETA_sleepdurationsumstats.txt
#Munge overall sleep duration
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/sleep_duration/overall_sleep_duration/P_BETA_sleepdurationsumstats.txt --out ${munged_dir}/overall_sleep_duration --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --info INFO --N 446118

##Short sleep duration
cd ${sumstats_dir}/behaviour/sleep_duration/short_sleep_duration
unzip shortsumstats.txt.zip
#Change BETA_SHORTSLEEP and P_SHORTSLEEP to BETA and P, respectively
sed -e '1s/BETA_SHORTSLEEP/BETA/' ${sumstats_dir}/behaviour/sleep_duration/short_sleep_duration/shortsumstats.txt > ${sumstats_dir}/behaviour/sleep_duration/short_sleep_duration/P_BETA_shortsumstats.txt
sed -i -e '1s/P_SHORTSLEEP/P/' ${sumstats_dir}/behaviour/sleep_duration/short_sleep_duration/P_BETA_shortsumstats.txt
#Munge short sleep duration
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/sleep_duration/short_sleep_duration/P_BETA_shortsumstats.txt --out ${munged_dir}/short_sleep_duration --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --info INFO --N-cas 106129 --N-con 305742

##Long sleep duration
cd ${sumstats_dir}/behaviour/sleep_duration/long_sleep_duration
unzip longsumstats.txt.zip
#Change BETA_LONGSLEEP and P_LONGSLEEP to BETA and P, respectively
sed -e '1s/BETA_LONGSLEEP/BETA/' ${sumstats_dir}/behaviour/sleep_duration/long_sleep_duration/longsumstats.txt > ${sumstats_dir}/behaviour/sleep_duration/long_sleep_duration/P_BETA_longsumstats.txt
sed -i -e '1s/P_LONGSLEEP/P/' ${sumstats_dir}/behaviour/sleep_duration/long_sleep_duration/P_BETA_longsumstats.txt
#Munge long sleep duration
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/sleep_duration/long_sleep_duration/P_BETA_longsumstats.txt --out ${munged_dir}/long_sleep_duration --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --info INFO --N-cas 34184 --N-con 305742

##Insomnia
cd ${sumstats_dir}/behaviour/insomnia
unzip Saxena_fullUKBB_Insomnia_summary_stats.zip
#Change BETA_INSOMNIA and P_INSOMNIA to BETA and P, respectively
sed -e '1s/BETA_INSOMNIA/BETA/' ${sumstats_dir}/behaviour/insomnia/Saxena_fullUKBB_Insomnia_summary_stats.txt > ${sumstats_dir}/behaviour/insomnia/P_BETA_Saxena_fullUKBB_Insomnia_summary_stats.txt
sed -i -e '1s/P_INSOMNIA/P/' ${sumstats_dir}/behaviour/insomnia/P_BETA_Saxena_fullUKBB_Insomnia_summary_stats.txt
#Munge insomnia
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/insomnia/P_BETA_Saxena_fullUKBB_Insomnia_summary_stats.txt --out ${munged_dir}/insomnia --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 ALLELE1 --a2 ALLELE0 --frq A1FREQ --info INFO --N 453379

##Intelligence
cd ${sumstats_dir}/behaviour/intelligence
unzip SavageJansen_IntMeta_sumstats.zip
#Change N_analyzed, stdBeta, minINFO to N, BETA and INFO, respectively
sed -e '1s/N_analyzed/N/' ${sumstats_dir}/behaviour/intelligence/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt > ${sumstats_dir}/behaviour/intelligence/sumstats/N_BETA_INFO_SavageJansen_2018_intelligence_metaanalysis.txt
sed -i -e '1s/minINFO/INFO/' ${sumstats_dir}/behaviour/intelligence/sumstats/N_BETA_INFO_SavageJansen_2018_intelligence_metaanalysis.txt
sed -i -e '1s/stdBeta/BETA/' ${sumstats_dir}/behaviour/intelligence/sumstats/N_BETA_INFO_SavageJansen_2018_intelligence_metaanalysis.txt
#Munge intelligence (ignore Zscore as there is BETA)
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/intelligence/sumstats/N_BETA_INFO_SavageJansen_2018_intelligence_metaanalysis.txt --out ${munged_dir}/intelligence --merge-alleles ${ref_dir}/w_hm3.snplist --snp SNP --a1 A1 --a2 A2 --frq EAF_HRC --info INFO --ignore Zscore

##Educational attainment
cd ${sumstats_dir}/behaviour/educational_attainment
#Munge educational attainment
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/educational_attainment/GWAS_EA_excl23andMe.txt --out ${munged_dir}/educational_attainment --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 A1 --a2 A2 --frq EAF --N 766345

##Cognitive performance
cd ${sumstats_dir}/behaviour/cognitive_performance
#Munge cognitive performance
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/cognitive_performance/GWAS_CP_all.txt --out ${munged_dir}/cognitive_performance --merge-alleles ${ref_dir}/w_hm3.snplist --snp MarkerName --a1 A1 --a2 A2 --frq EAF --N 257828

##Neuroticism
cd ${sumstats_dir}/behaviour/neuroticism
gunzip sumstats_neuroticism_ctg_format.txt.gz
#Munge neuroticism
python ${ldsc_dir}/munge_sumstats.py --sumstats ${sumstats_dir}/behaviour/neuroticism/sumstats_neuroticism_ctg_format.txt --out ${munged_dir}/neuroticism --merge-alleles ${ref_dir}/w_hm3.snplist --snp RSID --a1 A1 --a2 A2 --frq EAF_UKB --info INFO_UKB --ignore SNP


#H2
python ${ldsc_dir}/ldsc.py --h2 ${munged_dir}/SCZ.sumstats.gz --ref-ld-chr ${ref_dir}/eur_ref_ld_chr/ --w-ld-chr ${ref_dir}/eur_w_ld_chr/ --out ${h2_dir}/scz_h2
python ${ldsc_dir}/ldsc.py --h2 ${munged_dir}/AN.sumstats.gz --ref-ld-chr ${ref_dir}/eur_ref_ld_chr/ --w-ld-chr ${ref_dir}/eur_w_ld_chr/ --out ${h2_dir}/an_h2
python ${ldsc_dir}/ldsc.py --h2 ${munged_dir}/ADHD.sumstats.gz --ref-ld-chr ${ref_dir}/eur_ref_ld_chr/ --w-ld-chr ${ref_dir}/eur_w_ld_chr/ --out ${h2_dir}/adhd_h2


#Genetic correlations
python ${ldsc_dir}/ldsc.py --rg ${munged_dir}/ADHD.sumstats.gz,${munged_dir}/AN.sumstats.gz,${munged_dir}/anxiety.sumstats.gz,${munged_dir}/ASD.sumstats.gz --ref-ld-chr ${ref_dir}/eur_w_ld_chr/ --w-ld-chr ${ref_dir}/eur_w_ld_chr/ --out ${rg_dir}/an_scz_adhd
