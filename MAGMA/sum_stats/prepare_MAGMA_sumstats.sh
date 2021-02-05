!#bin/bash

# Goal: Prepare summary statistics for MAGMA analyses
# Processed summary statistics files should:
# have SNP, CHR, BP as first three columns;
# have at least one of these columns: (“Z”,“OR”,“BETA”,“LOG_ODDS”,“SIGNED_SUMSTAT”);
# have all of these columns: (“SNP”,“CHR”,“BP”,“P”,“A1”,“A2”).
# Required:	GWAS summary statistics files
# By: Mitchell Olislagers
# Last updated: April 15 2020

sumstats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics
magma_sums_dir=/hpc/hers_en/molislagers/MAGMA/sumstats
magma_ref=/hpc/hers_en/molislagers/MAGMA/ref_data

## SCHIZOPHRENIA ##
## Input columns: SNP	Freq.A1	CHR	BP	A1	A2	OR	SE	P	SNP_RS
## Output columns: SNP	CHR	BP	Freq.A1	A1	A2	OR	SE	P

cp ${sumstats_dir}/psychiatric_diseases/SCZ/final_sumstats/clozuk_pgc2.meta.sumstats.only_rs.rs_col.header.txt ${magma_sums_dir}/SCZ

#Select columns of interest in the right order
awk -v OFS="\t" '{print $10,$3,$4,$2,$5,$6,$7,$8,$9}' ${magma_sums_dir}/SCZ/clozuk_pgc2.meta.sumstats.only_rs.rs_col.header.txt > ${magma_sums_dir}/SCZ/clozuk_pgc2.meta.sumstats_columns_of_interest.txt
#Rename SNP_RS to SNP
sed -e '1s/SNP_RS/SNP/' ${magma_sums_dir}/SCZ/clozuk_pgc2.meta.sumstats_columns_of_interest.txt > ${magma_sums_dir}/SCZ/clozuk_pgc2.meta.sumstats_MAGMA.txt
#Copy to MAGMA sumstats dir
cp ${magma_sums_dir}/SCZ/clozuk_pgc2.meta.sumstats_MAGMA.txt ${magma_sums_dir}/MAGMA_sumstats/SCZ_MAGMA.txt

## ANOREXIA NERVOSA ##
## Input columns: CHROM	POS	ID	REF	ALT	BETA	SE	PVAL	NGT	IMPINFO	NEFFDIV2	NCAS	NCON	DIRE
## Output columns: SNP	CHR	BP	A1	A2	BETA	SE	P	NGT	IMPINFO	NEFFDIV2	NCAS	NCON	DIRE
cp ${sumstats_dir}/psychiatric_diseases/AN/raw_GWAS/pgcAN2.2019-07.vcf.tsv ${magma_sums_dir}/AN
#Strip first lines to only retain the summary statistics
tail -n +71 ${magma_sums_dir}/AN/pgcAN2.2019-07.vcf.tsv > ${magma_sums_dir}/AN/only_sumstats_pgcAN2.2019-07.vcf.tsv
#Reorder columns
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14}' ${magma_sums_dir}/AN/only_sumstats_pgcAN2.2019-07.vcf.tsv > ${magma_sums_dir}/AN/only_sumstats_pgcAN2.2019-07_reordered.vcf.tsv
#Rename column names
sed -e '1s/ID/SNP/' -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/ALT/A1/' -e '1s/REF/A2/' -e '1s/PVAL/P/' ${magma_sums_dir}/AN/only_sumstats_pgcAN2.2019-07_reordered.vcf.tsv > ${magma_sums_dir}/AN/only_sumstats_pgcAN2.2019-07_reordered_MAGMA.vcf.tsv
#Copy to MAGMA sumstats dir
cp ${magma_sums_dir}/AN/only_sumstats_pgcAN2.2019-07_reordered_MAGMA.vcf.tsv ${magma_sums_dir}/MAGMA_sumstats/AN_MAGMA.txt

## ANXIETY ##
## Input columns: SNPID	CHR	BP	Allele1	Allele2	Freq1	Effect	StdErr	P.value	TotalN
## Output columns: SNP	CHR	BP	A1	A2	Freq1	BETA	StdErr	P.value	TotalN
cp ${sumstats_dir}/psychiatric_diseases/anxiety_disorder/angst.study.results/anxiety.meta.full.cc.tbl ${magma_sums_dir}/anxiety
#Rename column names
sed -e '1s/SNPID/SNP/' -e '1s/Allele1/A1/' -e '1s/Allele2/A2/' -e '1s/Effect/BETA/' -e '1s/P.value/P/' ${magma_sums_dir}/anxiety/anxiety.meta.full.cc.tbl > ${magma_sums_dir}/anxiety/anxiety.meta.full.cc_MAGMA.tbl
#Copy to MAGMA sumstats dir
cp ${magma_sums_dir}/anxiety/anxiety.meta.full.cc_MAGMA.tbl ${magma_sums_dir}/MAGMA_sumstats/anxiety_MAGMA.txt

## AUTISM SPECTRUM DISORDER ##
## Input columns: CHR	SNP	BP	A1	A2	INFO	OR	SE	P
## Output columns: SNP	CHR	BP	A1	A2	INFO	OR	SE	P
cp ${sumstats_dir}/psychiatric_diseases/ASD/iPSYCH-PGC_ASD_Nov2017 ${magma_sums_dir}/ASD
#Reorder columns
awk -v OFS="\t" '{print $2,$1,$3,$4,$5,$6,$7,$8,$9}' ${magma_sums_dir}/ASD/iPSYCH-PGC_ASD_Nov2017 > ${magma_sums_dir}/ASD/iPSYCH-PGC_ASD_Nov2017_MAGMA.tsv
#Copy to MAGMA sumstats dir
cp ${magma_sums_dir}/ASD/iPSYCH-PGC_ASD_Nov2017_MAGMA.tsv ${magma_sums_dir}/MAGMA_sumstats/ASD_MAGMA.txt

## ADHD ##
## Input columns: CHR	SNP	BP	A1	A2	INFO	OR	SE	P
## Output columns: SNP	CHR	BP	A1	A2	INFO	OR	SE	P
cp ${sumstats_dir}/psychiatric_diseases/ADHD/adhd_eur_jun2017 ${magma_sums_dir}/ADHD
#Reorder columns
awk -v OFS="\t" '{print $2,$1,$3,$4,$5,$6,$7,$8,$9}' ${magma_sums_dir}/ADHD/adhd_eur_jun2017 > ${magma_sums_dir}/ADHD/adhd_eur_jun2017_MAGMA.tsv
#Copy to MAGMA sumstats dir
cp ${magma_sums_dir}/ADHD/adhd_eur_jun2017_MAGMA.tsv ${magma_sums_dir}/MAGMA_sumstats/ADHD_MAGMA.txt

## BIPOLAR DISORDER ##
## Input columns: CHR	SNP	BP	A1	A2	FRQ_A_20352	FRQ_U_31358	INFO	OR	SE	P	ngt	Direction	HetISqt	HetDf	HetPVa	Nca	Nco	Neff
## Output columns: SNP	CHR	BP	A1	A2	FRQ_A_20352	FRQ_U_31358	INFO	OR	SE	P	ngt	Direction	HetISqt	HetDf	HetPVa	Nca	Nco	Neff
cp ${sumstats_dir}/psychiatric_diseases/BIP/daner_PGC_BIP32b_mds7a_0416a ${magma_sums_dir}/BIP
# Reorder columns
awk -v OFS="\t" '{ t = $1; $1 = $2; $2 = t; print; }' ${magma_sums_dir}/BIP/daner_PGC_BIP32b_mds7a_0416a > ${magma_sums_dir}/BIP/daner_PGC_BIP32b_mds7a_0416a_MAGMA.tsv
#Copy to MAGMA sumstats dir
cp ${magma_sums_dir}/BIP/daner_PGC_BIP32b_mds7a_0416a_MAGMA.tsv ${magma_sums_dir}/MAGMA_sumstats/BIP_MAGMA.txt

## CROSS DISORDERS ##
## Input columns: CHROM	POS	ID	REF	ALT	BETA	SE	PVAL	NGT	FCAS	FCON	IMPINFO	NEFFDIV2	NCAS	NCON	DIRE	ASSET	m.SCZ	m.BIP	m.MD	m.ASD	m.ADHD	m.TS	m.AN	m.OCD
## Output columns SNP	CHR	BP	A1	A2	BETA	SE	P	NGT	FCAS	FCON	IMPINFO	NEFFDIV2	NCAS	NCON	DIRE	ASSET	m.SCZ	m.BIP	m.MD	m.ASD	m.ADHD	m.TS	m.AN	m.OCD
cp ${sumstats_dir}/psychiatric_diseases/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt ${magma_sums_dir}/cross_disorders
#Reorder columns
awk -v OFS="\t" '{snp=$3; chr=$1; bp=$2; a1=$5; a2=$4; $1=snp; $2=chr; $3=bp; $4=a1; $5=a2; print;}' ${magma_sums_dir}/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner.txt > ${magma_sums_dir}/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner_reordered.txt
#Rename column names
sed -e '1s/ID/SNP/' -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/ALT/A1/' -e '1s/REF/A2/' -e '1s/PVAL/P/' ${magma_sums_dir}/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner_reordered.txt > ${magma_sums_dir}/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner_reordered_MAGMA.txt
cp ${magma_sums_dir}/cross_disorders/pgc_cdg2_meta_no23andMe_oct2019_v2.txt.daner_reordered_MAGMA.txt ${magma_sums_dir}/MAGMA_sumstats/cross_MAGMA.txt

## MAJOR DEPRESSIVE DISORDER
## Input columns: MarkerName A1 A2 Freq LogOR StdErrLogOR P
## Output columns: SNP	CHR	BP	A1	A2	Freq	LogOR	StdErrLogOR	P
cp ${sumstats_dir}/psychiatric_diseases/MDD/DS_10283_3203/PGC_UKB_depression_genome-wide.txt ${magma_sums_dir}/MDD
#Summary statistics only contain rsno, and no chromosomal position. 1kg_rs.bim is used to map rsno's to chromosomal position.
# 1kg_rs.bim is the original file, chrpos_1k_rs.bim is a new file with the chr and pos column concatenated by a ":"
cp /hpc/hers_en/molislagers/LDSC/ref_data/regression/SNPref/1kg_rs.bim ${magma_ref}
cp /hpc/hers_en/molislagers/LDSC/ref_data/regression/SNPref/new_1kg_rs.bim ${magma_ref}/chrpos_1kg_rs.bim
#Select SNP, CHR, and BP from reference file
awk -v OFS="\t" '{print $2,$1,$4}' ${magma_ref}/colnames_1kg_rs.bim > ${magma_ref}/snpchrbp_1kg_rs.bim
#Add CHR to summary stats
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/MDD/PGC_UKB_depression_genome-wide.txt > ${magma_sums_dir}/MDD/chr_PGC_UKB_depression_genome-wide.txt
#Add BP to summary stats
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/MDD/chr_PGC_UKB_depression_genome-wide.txt > ${magma_sums_dir}/MDD/chr_bp_PGC_UKB_depression_genome-wide.txt
#Create columns names
echo -e "SNP\tCHR\tBP\tA1\tA2\tFreq\tLOG_ODDS\tLOG_SE\tP" | cat - ${magma_sums_dir}/MDD/chr_bp_PGC_UKB_depression_genome-wide.txt > ${magma_sums_dir}/MDD/PGC_UKB_depression_genome-wide_MAGMA.txt
cp ${magma_sums_dir}/MDD/PGC_UKB_depression_genome-wide_MAGMA.txt ${magma_sums_dir}/MAGMA_sumstats/MDD_MAGMA.txt

## OBSESSIVE COMPULSIVE DISORDER ##
##Input columns: CHR	SNP	BP	A1	A2	INFO	OR	SE	P
##Output columns: SNP	CHR	BP	A1	A2	INFO	OR	SE	P
cp ${sumstats_dir}/psychiatric_diseases/OCD/PGC_OCD_Aug2017/ocd_aug2017 ${magma_sums_dir}/OCD
#Reorder columns
awk -v OFS="\t" '{snp=$2; chr=$1; $1=snp; $2=chr; print;}' ${magma_sums_dir}/OCD/ocd_aug2017 > ${magma_sums_dir}/OCD/ordered_ocd_aug2017.txt
cp ${magma_sums_dir}/OCD/ordered_ocd_aug2017.txt ${magma_sums_dir}/MAGMA_sumstats/OCD_MAGMA.txt

## POST-TRAUMATIC STRESS DISORDER ##
##Input columns: CHR	SNP	BP	A1	A2	FRQ_A_23212	FRQ_U_151447	INFO	OR	SE	P	ngt	Direction	HetISqt	HetDf	HetPVa	Nca	Nco	Neff
##Output columns: SNP	CHR	BP	A1	A2	FRQ_A_23212	FRQ_U_151447	INFO	OR	SE	P	ngt	Direction	HetISqt	HetDf	HetPVa	Nca	Nco	Neff
cp ${sumstats_dir}/psychiatric_diseases/PTSD/pts_eur_freeze2_overall.results ${magma_sums_dir}/PTSD
#Reorder columns
awk -v OFS="\t" '{snp=$2; chr=$1; $1=snp; $2=chr; print;}' ${magma_sums_dir}/PTSD/pts_eur_freeze2_overall.results > ${magma_sums_dir}/PTSD/ordered_pts_eur_freeze2_overall.results
cp ${magma_sums_dir}/PTSD/ordered_pts_eur_freeze2_overall.results ${magma_sums_dir}/MAGMA_sumstats/PTSD_MAGMA.txt

## TOURETTE SYNDROME ##
##Input columns: SNP CHR BP A1 A2 INFO OR SE P
##Output columns: SNP	CHR	BP	A1	A2	INFO	OR	SE	P (tab-delimited)

cp ${sumstats_dir}/psychiatric_diseases/TS/TS_Oct2018 ${magma_sums_dir}/TS
awk -v OFS="\t" '{print $1,$2,$3,$4,$5,$6,$7,$8,$9}' ${magma_sums_dir}/TS/TS_Oct2018 > ${magma_sums_dir}/TS/TS_MAGMA.txt
cp ${magma_sums_dir}/TS/TS_MAGMA.txt ${magma_sums_dir}/MAGMA_sumstats/TS_MAGMA.txt

## ALS ##
## Input columns: SNP	Allele1	Allele2	Effect	StdErr	P	Direction	HetISq	HetChiSq	HetDf	HetPVal	SNPforpos	CHR	BP	effectAlleleFreq	effectAlleleMinFreq	effectAlleleMaxFreq	effectAlleleFreqStdErr
## Ouput columns: SNP	CHR	BP	A1	A2	BETA	SE	P	CHRPOS	Direction	HetISq	HetChiSq	HetDf	HetPVal	SNPforpos	effectAlleleFreq	effectAlleleMinFreq	effectAlleleMaxFreq	effectAlleleFreqStdErr

cp ${sumstats_dir}/neurological_disorders/ALS/alsMetaSummaryStats_march21st2018.tab ${magma_sums_dir}/ALS
#grep header
head -n 1 ${magma_sums_dir}/ALS/alsMetaSummaryStats_march21st2018.tab > ${magma_sums_dir}/ALS/header_alsMetaSummaryStats_march21st2018.tab
#Change columns names (SNP = CHRPOS)
sed -e '1s/SNP/CHRPOS/' -e '1s/Effect/BETA/' -e '1s/Allele1/A1/' -e '1s/Allele2/A2/' -e '1s/StdErr/SE/' ${magma_sums_dir}/ALS/header_alsMetaSummaryStats_march21st2018.tab > ${magma_sums_dir}/ALS/renamed_header_alsMetaSummaryStats_march21st2018.tab
#Add SNP columns name for rsids
sed -e 's/^/SNP\t/' ${magma_sums_dir}/ALS/renamed_header_alsMetaSummaryStats_march21st2018.tab > ${magma_sums_dir}/ALS/all_renamed_header_alsMetaSummaryStats_march21st2018.tab
#Merge chr:pos column to RSIDs in reference file new_1_kg_rs.bim
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print a[$1],$1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16,$17,$18}' ${magma_ref}/chrpos_1kg_rs.bim ${magma_sums_dir}/ALS/alsMetaSummaryStats_march21st2018.tab >  ${magma_sums_dir}/ALS/snp_alsMetaSummaryStats_march21st2018.tab
#Add header
cat ${magma_sums_dir}/ALS/all_renamed_header_alsMetaSummaryStats_march21st2018.tab  ${magma_sums_dir}/ALS/snp_alsMetaSummaryStats_march21st2018.tab  > ${magma_sums_dir}/ALS/all_columns_and_header_alsMetaSummaryStats_march21st2018
#Reorder columns
awk -v OFS="\t" '{print $1,$14,$15,$3,$4,$5,$6,$7,$2,$8,$9,$10,$11,$12,$13,$16,$17,$18,$19}' ${magma_sums_dir}/ALS/all_columns_and_header_alsMetaSummaryStats_march21st2018 > ${magma_sums_dir}/ALS/reordered_all_columns_and_header_alsMetaSummaryStats_march21st2018
cp ${magma_sums_dir}/ALS/reordered_all_columns_and_header_alsMetaSummaryStats_march21st2018 ${magma_sums_dir}/MAGMA_sumstats/ALS_MAGMA.txt

## ALZHEIMERS ##
##Input columns: uniqID.a1a2	CHR	BP	A1	A2	SNP	Z	P	Nsum	Neff	dir	EAF	BETA	SE
##Output columns: SNP	CHR	BP	A1	A2	BETA	SE	P	uniqID.a1a2	Z	Nsum	Neff	dir	EAF
cp ${sumstats_dir}/neurological_disorders/alzheimers/AD_sumstats_Jansenetal_2019sept.txt ${magma_sums_dir}/alzheimers
#Reorder columns
awk -v OFS="\t" '{print $6,$2,$3,$4,$5,$13,$14,$8,$1,$7,$9,$10,$11,$12}' ${magma_sums_dir}/alzheimers/AD_sumstats_Jansenetal_2019sept.txt > ${magma_sums_dir}/alzheimers/ordered_AD_sumstats_Jansenetal_2019sept.txt
cp ${magma_sums_dir}/alzheimers/ordered_AD_sumstats_Jansenetal_2019sept.txt ${magma_sums_dir}/MAGMA_sumstats/alzheimers_MAGMA.txt

## EPILEPSY ##
#ALL EPILEPSY, GENERALIZED EPILEPSY AND FOCAL EPILEPSY
#Input columns: CHR	BP	MarkerName	Allele1	Allele2	Freq1	FreqSE	Weight	Zscore	P-value	Direction	HetISq	HetChiSq	HetDf	HetPVal
#Output columns: SNP	CHR	BP	A1	A2	Z	P	Freq1	FreqSE	Weight	Direction	HetISq	HetChiSq	HetDf	HetPVal
cp ${sumstats_dir}/neurological_disorders/epilepsy/all_epilepsy/all_epilepsy_METAL ${magma_sums_dir}/epilepsy
cp ${sumstats_dir}/neurological_disorders/epilepsy/focal/focal_epilepsy_METAL ${magma_sums_dir}/epilepsy
cp ${sumstats_dir}/neurological_disorders/epilepsy/generalized/generalised_epilepsy_METAL ${magma_sums_dir}/epilepsy
#Change columns names
sed -e '1s/MarkerName/SNP/' -e'1s/Allele1/A1/' -e'1s/Allele2/A2/' -e '1s/Zscore/Z/' -e'1s/P-value/P/' ${magma_sums_dir}/epilepsy/all_epilepsy_METAL > ${magma_sums_dir}/epilepsy/renamed_all_epilepsy_METAL
sed -e '1s/MarkerName/SNP/' -e'1s/Allele1/A1/' -e'1s/Allele2/A2/' -e '1s/Zscore/Z/' -e'1s/P-value/P/' ${magma_sums_dir}/epilepsy/focal_epilepsy_METAL > ${magma_sums_dir}/epilepsy/renamed_focal_epilepsy_METAL
sed -e '1s/MarkerName/SNP/' -e'1s/Allele1/A1/' -e'1s/Allele2/A2/' -e '1s/Zscore/Z/' -e'1s/P-value/P/' ${magma_sums_dir}/epilepsy/generalised_epilepsy_METAL > ${magma_sums_dir}/epilepsy/renamed_generalised_epilepsy_METAL
#Reorder columns
awk -v OFS="\t" '{print $3,$1,$2,$4,$5,$9,$10,$6,$7,$8,$11,$12,$13,$14,$15}' ${magma_sums_dir}/epilepsy/renamed_all_epilepsy_METAL > ${magma_sums_dir}/epilepsy/ordered_renamed_all_epilepsy_METAL
awk -v OFS="\t" '{print $3,$1,$2,$4,$5,$9,$10,$6,$7,$8,$11,$12,$13,$14,$15}' ${magma_sums_dir}/epilepsy/renamed_focal_epilepsy_METAL > ${magma_sums_dir}/epilepsy/ordered_renamed_focal_epilepsy_METAL
awk -v OFS="\t" '{print $3,$1,$2,$4,$5,$9,$10,$6,$7,$8,$11,$12,$13,$14,$15}' ${magma_sums_dir}/epilepsy/renamed_generalised_epilepsy_METAL > ${magma_sums_dir}/epilepsy/ordered_renamed_generalised_epilepsy_METAL
cp ${magma_sums_dir}/epilepsy/ordered_renamed_all_epilepsy_METAL ${magma_sums_dir}/MAGMA_sumstats/all_epilepsy_MAGMA.txt
cp ${magma_sums_dir}/epilepsy/ordered_renamed_focal_epilepsy_METAL ${magma_sums_dir}/MAGMA_sumstats/focal_epilepsy_MAGMA.txt
cp ${magma_sums_dir}/epilepsy/ordered_renamed_generalised_epilepsy_METAL ${magma_sums_dir}/MAGMA_sumstats/generalised_epilepsy_MAGMA.txt

## PARKINSON ##
#Input columns: SNP	A1	A2	freq	b	se	p	N_cases	N_controls
## Input SNP columns is CHR:POS format, so will need to be converted to RSID
#Output columns: SNP	CHR	BP	A1	A2	BETA	P	SE	CHRPOS(old SNP colum)	freq	N_cases	N_controls
cp ${sumstats_dir}/neurological_disorders/parkinson/nallsEtAl2019_excluding23andMe_allVariants.tab ${magma_sums_dir}/parkinson
#Grep header
head -n 1 ${magma_sums_dir}/parkinson/nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/header_nallsEtAl2019_excluding23andMe_allVariants.tab
#Rename columns headers
sed -e '1s/SNP/CHRPOS/' -e '1s/b/BETA/' -e '1s/se/SE/' -e '1s/p/P/' ${magma_sums_dir}/parkinson/header_nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/renamed_header_nallsEtAl2019_excluding23andMe_allVariants.tab
#Create new header
sed -e 's/^/SNP\tCHR\tBP\t/' ${magma_sums_dir}/parkinson/renamed_header_nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/finalheader_nallsEtAl2019_excluding23andMe_allVariants.tab
#Remove 'chr' from positions
sed 's/chr//g' ${magma_sums_dir}/parkinson/nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/nochr_nallsEtAl2019_excluding23andMe_allVariants.tab
#Map chr:pos to rsIDs
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print a[$1],$1,$2,$3,$4,$5,$6,$7,$8,$9}' ${magma_ref}/chrpos_1kg_rs.bim ${magma_sums_dir}/parkinson/nochr_nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab
#Add CHR column
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7,$8,$9, $10}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/parkinson/rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab
#Add BP column
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8,$9, $10,$11}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/parkinson/chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/bpcol_chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab
#Add header
cat ${magma_sums_dir}/parkinson/finalheader_nallsEtAl2019_excluding23andMe_allVariants.tab ${magma_sums_dir}/parkinson/bpcol_chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/final_bpcol_chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab
#Reorder columns
awk -v OFS="\t" '{print $1,$2,$3,$5,$6,$8,$10,$9,$4,$7,$11,$12}' ${magma_sums_dir}/parkinson/final_bpcol_chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab > ${magma_sums_dir}/parkinson/ordered_final_bpcol_chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab
cp ${magma_sums_dir}/parkinson/ordered_final_bpcol_chrcol_rsids_nochr__nallsEtAl2019_excluding23andMe_allVariants.tab ${magma_sums_dir}/MAGMA_sumstats/parkinson_MAGMA.txt

##STROKE
#ALL STROKE, CARDIOEMBOLIC STROKE, ISCHEMIC STROKE, LARGE ARTERY STROKE, SMALL VESSEL STROKE
#Input columns: MarkerName Allele1 Allele2 Freq1 Effect StdErr P-value
#Output columns: SNP	CHR	BP	A1	A2	P	SE	BETA	Freq1
cp ${sumstats_dir}/neurological_disorders/stroke/all_stroke/MEGASTROKE.1.AS.EUR.out ${magma_sums_dir}/stroke
cp ${sumstats_dir}/neurological_disorders/stroke/ischemic/MEGASTROKE.2.AIS.EUR.out ${magma_sums_dir}/stroke
cp ${sumstats_dir}/neurological_disorders/stroke/small_vessel/MEGASTROKE.5.SVS.EUR.out ${magma_sums_dir}/stroke
cp ${sumstats_dir}/neurological_disorders/stroke/cardioembolic/MEGASTROKE.4.CES.EUR.out ${magma_sums_dir}/stroke
cp ${sumstats_dir}/neurological_disorders/stroke/large_artery/MEGASTROKE.3.LAS.EUR.out ${magma_sums_dir}/stroke
#Grep header
head -n 1 ${magma_sums_dir}/stroke/MEGASTROKE.1.AS.EUR.out > ${magma_sums_dir}/stroke/header_MEGASTROKE.out
#Replace whitespaces with tabs in header
awk -v OFS="\t" '$1=$1' ${magma_sums_dir}/stroke/header_MEGASTROKE.out > ${magma_sums_dir}/stroke/tab_header_MEGASTROKE.out
#Rename header
sed -e '1s/MarkerName/SNP/' -e '1s/Allele1/A1/' -e '1s/Allele2/A2/' -e '1s/Effect/BETA/' -e '1s/StdErr/SE/' -e '1s/P-value/P/' ${magma_sums_dir}/stroke/tab_header_MEGASTROKE.out > ${magma_sums_dir}/stroke/renamed_tab_header_MEGASTROKE.out
sed -e '1s/SNP/SNP\tCHR\tBP/' ${magma_sums_dir}/stroke/renamed_tab_header_MEGASTROKE.out > ${magma_sums_dir}/stroke/final_renamed_tab_header_MEGASTROKE.out
#Add CHR column
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/MEGASTROKE.1.AS.EUR.out > ${magma_sums_dir}/stroke/chr_MEGASTROKE.1.AS.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/MEGASTROKE.2.AIS.EUR.out > ${magma_sums_dir}/stroke/chr_MEGASTROKE.2.AIS.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/MEGASTROKE.5.SVS.EUR.out > ${magma_sums_dir}/stroke/chr_MEGASTROKE.5.SVS.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/MEGASTROKE.4.CES.EUR.out > ${magma_sums_dir}/stroke/chr_MEGASTROKE.4.CES.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$2;next} ($1 in a) {print $1,a[$1],$2,$3,$4,$5,$6,$7}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/MEGASTROKE.3.LAS.EUR.out > ${magma_sums_dir}/stroke/chr_MEGASTROKE.3.LAS.EUR.out 
#Add BP column
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/chr_MEGASTROKE.1.AS.EUR.out  > ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.1.AS.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/chr_MEGASTROKE.2.AIS.EUR.out  > ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.2.AIS.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/chr_MEGASTROKE.5.SVS.EUR.out  > ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.5.SVS.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/chr_MEGASTROKE.4.CES.EUR.out > ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.4.CES.EUR.out
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/stroke/chr_MEGASTROKE.3.LAS.EUR.out  > ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.3.LAS.EUR.out
#Add header
cat ${magma_sums_dir}/stroke/final_renamed_tab_header_MEGASTROKE.out ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.1.AS.EUR.out > ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.1.AS.EUR.out
cat ${magma_sums_dir}/stroke/final_renamed_tab_header_MEGASTROKE.out ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.2.AIS.EUR.out > ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.2.AIS.EUR.out
cat ${magma_sums_dir}/stroke/final_renamed_tab_header_MEGASTROKE.out ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.5.SVS.EUR.out > ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.5.SVS.EUR.out
cat ${magma_sums_dir}/stroke/final_renamed_tab_header_MEGASTROKE.out ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.4.CES.EUR.out > ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.4.CES.EUR.out
cat ${magma_sums_dir}/stroke/final_renamed_tab_header_MEGASTROKE.out ${magma_sums_dir}/stroke/chr_bp_MEGASTROKE.3.LAS.EUR.out > ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.3.LAS.EUR.out
cp ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.1.AS.EUR.out ${magma_sums_dir}/MAGMA_sumstats/all_stroke_MAGMA.txt
cp ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.2.AIS.EUR.out ${magma_sums_dir}/MAGMA_sumstats/ischemic_stroke_MAGMA.txt
cp ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.5.SVS.EUR.out ${magma_sums_dir}/MAGMA_sumstats/small_vessel_epilepsy_MAGMA.txt
cp ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.4.CES.EUR.out ${magma_sums_dir}/MAGMA_sumstats/cardioembolic_stroke_MAGMA.txt
cp ${magma_sums_dir}/stroke/final_chr_bp_MEGASTROKE.3.LAS.EUR.out ${magma_sums_dir}/MAGMA_sumstats/large_artery_stroke_MAGMA.txt

## ALCOHOL USE ##
#Input columns: chr rsid a_0 a_1 info beta_T se_T p_T beta_C se_C p_C beta_P se_P p_P N
#Output columns: SNP	CHR	BP	A1	A2	info	BETA	SE	P	beta_C	se_C	p_C	beta_P	se_P	p_P	N
cp ${sumstats_dir}/substance_use/alcohol_use/AUDIT_UKB_2018_AJP.txt ${magma_sums_dir}/alcohol_use
#Replace whitespaces with tabs
awk -v OFS="\t" '$1=$1' ${magma_sums_dir}/alcohol_use/AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/tab_AUDIT_UKB_2018_AJP.txt
#Rename column names
sed -e '1s/rsid/SNP/' -e '1s/chr/CHR\tBP/' -e '1s/a_1/A1/' -e '1s/a_0/A2/' -e '1s/beta_T/BETA/' -e '1s/se_T/SE/' -e '1s/p_T/P/' ${magma_sums_dir}/alcohol_use/tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/renamed_tab_AUDIT_UKB_2018_AJP.txt
#Grep header
head -n 1 ${magma_sums_dir}/alcohol_use/renamed_tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/header_renamed_tab_AUDIT_UKB_2018_AJP.txt
#Grep rsids
grep 'rs' ${magma_sums_dir}/alcohol_use/renamed_tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt
#Change column order
awk -v OFS="\t" '{chr=$1; snp=$2; a1=$4; a2=$3; $1=snp; $2=chr; $3=a1; $4=a2; print;}' ${magma_sums_dir}/alcohol_use/rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/ordered_rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt
#Add BP column
awk -v OFS="\t" 'FNR==NR{a[$1]=$3;next} ($1 in a) {print $1,$2,a[$1],$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15}' ${magma_ref}/snpchrbp_1kg_rs.bim ${magma_sums_dir}/alcohol_use/ordered_rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/bp_ordered_rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt
#Change header order
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$6,$7,$8,$9,$10,$11,$12,$13,$14,$15,$16}' ${magma_sums_dir}/alcohol_use/header_renamed_tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/finalheader_renamed_tab_AUDIT_UKB_2018_AJP.txt
#Add header
cat ${magma_sums_dir}/alcohol_use/finalheader_renamed_tab_AUDIT_UKB_2018_AJP.txt ${magma_sums_dir}/alcohol_use/bp_ordered_rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt > ${magma_sums_dir}/alcohol_use/final_bp_ordered_rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt
cp ${magma_sums_dir}/alcohol_use/final_bp_ordered_rsids_renamed_tab_AUDIT_UKB_2018_AJP.txt ${magma_sums_dir}/MAGMA_sumstats/alcohol_use_MAGMA.txt

## CANNABIS ##
#Input columns: CHR	SNP	BP	A1	A2	FRQ	BETA	SE	Z	P	Direction	HetISq	HetDf	HetPVa	Nca	Nco	Neff
#Output columns: SNP	CHR	BP	A1	A2	FRQ	BETA	SE	Z	P	Direction	HetISq	HetDf	HetPVa	Nca	Nco	Neff
cp ${sumstats_dir}/substance_use/cannabis/Cannabis_ICC_UKB_het.txt ${magma_sums_dir}/cannabis
#Reorder columns
awk -v OFS="\t" '{snp=$2; chr=$1; $1=snp; $2=chr; print;}' ${magma_sums_dir}/cannabis/Cannabis_ICC_UKB_het.txt > ${magma_sums_dir}/cannabis/ordered_Cannabis_ICC_UKB_het.txt
cp ${magma_sums_dir}/cannabis/ordered_Cannabis_ICC_UKB_het.txt ${magma_sums_dir}/MAGMA_sumstats/cannabis_MAGMA.txt

## ALCOHOL DEPENDENCE ##
#Input columns: CHR SNP BP A1 A2 Z P Weight
#Output columns: SNP	CHR	BP	A1	A2	Z	P	Weight
cp ${sumstats_dir}/substance_use/alcohol_dependence/input_GWAS/pgc_alcdep.eur_discovery.aug2018_release.txt ${magma_sums_dir}/alcohol_dependence
#Create extra column with only rs#
awk -v n=1 '{s = gensub("^(([^:]*:){"n"}).*$", "\\1", 1); sub(".$","",s); print $0, s}' ${magma_sums_dir}/alcohol_dependence/pgc_alcdep.eur_discovery.aug2018_release.txt > ${magma_sums_dir}/alcohol_dependence/rscol_pgc_alcdep.eur_discovery.aug2018_release.txt
sed -i 's/\r//' ${magma_sums_dir}/alcohol_dependence/rscol_pgc_alcdep.eur_discovery.aug2018_release.txt
#Remove header
sed '1d' ${magma_sums_dir}/alcohol_dependence/rscol_pgc_alcdep.eur_discovery.aug2018_release.txt > ${magma_sums_dir}/alcohol_dependence/nohead_rscol_pgc_alcdep.eur_discovery.aug2018_release.txt
#Order columns
awk -v OFS="\t" '{print $10,$1,$3,$4,$5,$6,$7,$8}' ${magma_sums_dir}/alcohol_dependence/nohead_rscol_pgc_alcdep.eur_discovery.aug2018_release.txt > ${magma_sums_dir}/alcohol_dependence/ordered_nohead_rscol_pgc_alcdep.eur_discovery.aug2018_release.txt
#Create header
echo -e "SNP\tCHR\tBP\tA1\tA2\tZ\tP\tWeight" | cat - ${magma_sums_dir}/alcohol_dependence/ordered_nohead_rscol_pgc_alcdep.eur_discovery.aug2018_release.txt > ${magma_sums_dir}/alcohol_dependence/final_ordered_rscol_pgc_alcdep.eur_discovery.aug2018_release.txt
cp ${magma_sums_dir}/alcohol_dependence/final_ordered_rscol_pgc_alcdep.eur_discovery.aug2018_release.txt ${magma_sums_dir}/MAGMA_sumstats/alcohol_dependence_MAGMA.txt

## SMOKING PHENOTYPES ##
#Input columns: CHROM	POS	RSID	REF	ALT	AF	STAT	PVALUE	BETA	SE	N	EFFECTIVE_N	Number_of_Studies	ANNO	ANNOFULL
#Output columns: SNP	CHR	BP	A1	A2	P	BETA	SE	AF	STAT	N	EFFECTIVE_N	Number_of_Studies	Anno	ANNOFULL
# AGE OF INITIATION, CIGARETTES PER DAY, DRINKS PER WEEK, SMOKING CESSATION, SMOKING INITIATION
cp ${sumstats_dir}/substance_use/ever_smoked/SmokingInitiation.txt ${magma_sums_dir}/ever_smoked
cp ${sumstats_dir}/substance_use/cigarettes_pd/CigarettesPerDay.txt ${magma_sums_dir}/cigarettes_pd
cp ${sumstats_dir}/substance_use/drinks_pw/DrinksPerWeek.txt ${magma_sums_dir}/drinks_pw
cp ${sumstats_dir}/substance_use/smoking_cessation/SmokingCessation.txt ${magma_sums_dir}/smoking_cessation
cp ${sumstats_dir}/substance_use/smoking_initiation/AgeofInitiation.txt ${magma_sums_dir}/age_initiation
#Change columns names
sed -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/RSID/SNP/' -e '1s/REF/A2/' -e '1s/ALT/A1/' -e '1s/PVALUE/P/' ${magma_sums_dir}/ever_smoked/SmokingInitiation.txt > ${magma_sums_dir}/ever_smoked/renamed_SmokingInitiation.txt
sed -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/RSID/SNP/' -e '1s/REF/A2/' -e '1s/ALT/A1/' -e '1s/PVALUE/P/' ${magma_sums_dir}/cigarettes_pd/CigarettesPerDay.txt > ${magma_sums_dir}/cigarettes_pd/renamed_CigarettesPerDay.txt
sed -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/RSID/SNP/' -e '1s/REF/A2/' -e '1s/ALT/A1/' -e '1s/PVALUE/P/' ${magma_sums_dir}/drinks_pw/DrinksPerWeek.txt > ${magma_sums_dir}/drinks_pw/renamed_DrinksPerWeek.txt
sed -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/RSID/SNP/' -e '1s/REF/A2/' -e '1s/ALT/A1/' -e '1s/PVALUE/P/' ${magma_sums_dir}/smoking_cessation/SmokingCessation.txt > ${magma_sums_dir}/smoking_cessation/renamed_SmokingCessation.txt
sed -e '1s/CHROM/CHR/' -e '1s/POS/BP/' -e '1s/RSID/SNP/' -e '1s/REF/A2/' -e '1s/ALT/A1/' -e '1s/PVALUE/P/' ${magma_sums_dir}/age_initiation/AgeofInitiation.txt > ${magma_sums_dir}/age_initiation/renamed_AgeofInitiation.txt
#Order columns
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$8,$9,$10,$6,$7,$11,$12,$13,$14,$15}' ${magma_sums_dir}/ever_smoked/renamed_SmokingInitiation.txt > ${magma_sums_dir}/ever_smoked/ordered_renamed_SmokingInitiation.txt
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$8,$9,$10,$6,$7,$11,$12,$13,$14,$15}' ${magma_sums_dir}/cigarettes_pd/renamed_CigarettesPerDay.txt > ${magma_sums_dir}/cigarettes_pd/ordered_renamed_CigarettesPerDay.txt
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$8,$9,$10,$6,$7,$11,$12,$13,$14,$15}' ${magma_sums_dir}/drinks_pw/renamed_DrinksPerWeek.txt > ${magma_sums_dir}/drinks_pw/ordered_renamed_DrinksPerWeek.txt
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$8,$9,$10,$6,$7,$11,$12,$13,$14,$15}' ${magma_sums_dir}/smoking_cessation/renamed_SmokingCessation.txt > ${magma_sums_dir}/smoking_cessation/ordered_renamed_SmokingCessation.txt
awk -v OFS="\t" '{print $3,$1,$2,$5,$4,$8,$9,$10,$6,$7,$11,$12,$13,$14,$15}' ${magma_sums_dir}/age_initiation/renamed_AgeofInitiation.txt > ${magma_sums_dir}/age_initiation/ordered_renamed_AgeofInitiation.txt
cp ${magma_sums_dir}/ever_smoked/ordered_renamed_SmokingInitiation.txt ${magma_sums_dir}/MAGMA_sumstats/ever_smoked_MAGMA.txt
cp ${magma_sums_dir}/cigarettes_pd/ordered_renamed_CigarettesPerDay.txt ${magma_sums_dir}/MAGMA_sumstats/cigarettes_pd_MAGMA.txt
cp ${magma_sums_dir}/drinks_pw/ordered_renamed_DrinksPerWeek.txt ${magma_sums_dir}/MAGMA_sumstats/drinks_pw_MAGMA.txt
cp ${magma_sums_dir}/smoking_cessation/ordered_renamed_SmokingCessation.txt ${magma_sums_dir}/MAGMA_sumstats/smoking_cessation_MAGMA.txt
cp ${magma_sums_dir}/age_initiation/ordered_renamed_AgeofInitiation.txt ${magma_sums_dir}/MAGMA_sumstats/age_initiation_MAGMA.txt

## BODY MASS INDEX ##
#Input columns: CHR	POS	SNP	Tested_Allele	Other_Allele	Freq_Tested_Allele_in_HRS	BETA	SE	P	N
#Output columns: SNP	CHR	BP	A1	A2	BETA	P	Freq_Tested_Allele_in_HRS	SE	N
cp ${sumstats_dir}/behaviour/BMI/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt ${magma_sums_dir}/BMI
#Rename columns
sed -e '1s/POS/BP/' -e '1s/Tested_Allele/A1/' -e '1s/Other_Allele/A2/' ${magma_sums_dir}/BMI/Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt > ${magma_sums_dir}/BMI/renamed_Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt
#Order columns
awk -v OFS="\t" '{print $3,$1,$2,$4,$5,$7,$9,$6,$8,$10}' ${magma_sums_dir}/BMI/renamed_Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt > ${magma_sums_dir}/BMI/ordered_renamed_Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt
cp ${magma_sums_dir}/BMI/ordered_renamed_Meta-analysis_Locke_et_al%2BUKBiobank_2018_UPDATED.txt ${magma_sums_dir}/MAGMA_sumstats/BMI_MAGMA.txt

## HEIGHT ##
#Input columns: CHR	POS	SNP	Tested_Allele	Other_Allele	Freq_Tested_Allele_in_HRS	BETA	SE	P	N
#Output columns: SNP	CHR	BP	A1	A2	BETA	P	Freq_Tested_Allele_in_HRS	SE	N
cp ${sumstats_dir}/behaviour/height/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt ${magma_sums_dir}/height
#Rename columns
sed -e '1s/POS/BP/' -e '1s/Tested_Allele/A1/' -e '1s/Other_Allele/A2/' ${magma_sums_dir}/height/Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt > ${magma_sums_dir}/height/renamed_Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt
#Order columns
awk -v OFS="\t" '{print $3,$1,$2,$4,$5,$7,$9,$6,$8,$10}' ${magma_sums_dir}/height/renamed_Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt > ${magma_sums_dir}/height/ordered_renamed_Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt
cp ${magma_sums_dir}/height/ordered_renamed_Meta-analysis_Wood_et_al%2BUKBiobank_2018.txt ${magma_sums_dir}/MAGMA_sumstats/height_MAGMA.txt

## CHRONOTYPE ##
#Input columns: SNP	CHR	BP	ALLELE1	ALLELE0	A1FREQ	INFO	BETA	SE	P_BOLT_LMM	HWE_P
#Output columns: SNP	CHR	BP	A1	A2	A1FREQ	INFO	BETA	SE	P	HWE_P
cp ${sumstats_dir}/behaviour/chronotype/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt ${magma_sums_dir}/chronotype
#Rename columns
sed -e '1s/ALLELE1/A1/' -e '1s/ALLELE0/A2/' -e '1s/P_BOLT_LMM/P/' ${magma_sums_dir}/chronotype/chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt > ${magma_sums_dir}/chronotype/renamed_chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt
cp ${magma_sums_dir}/chronotype/renamed_chronotype_raw_BOLT.output_HRC.only_plus.metrics_maf0.001_hwep1em12_info0.3.txt ${magma_sums_dir}/MAGMA_sumstats/chronotype_MAGMA.txt

## EXCESSIVE DAYTIME SLEEPINESS ##
#Input columns: SNP	CHR	BP	ALLELE1	ALLELE0	A1FREQ	INFO	BETA	SE	P
#Output columns: SNP	CHR	BP	A1	A2	A1FREQ	INFO	BETA	SE	P
cp ${sumstats_dir}/behaviour/daytime_sleepiness/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt ${magma_sums_dir}/daytime_sleepiness
#Rename columns
sed -e '1s/ALLELE1/A1/' -e '1s/ALLELE0/A2/' -e '1s/P_BOLT_LMM/P/' ${magma_sums_dir}/daytime_sleepiness/Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt > ${magma_sums_dir}/daytime_sleepiness/renamed_Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt
cp ${magma_sums_dir}/daytime_sleepiness/renamed_Saxena.fullUKBB.DaytimeSleepiness.sumstats.txt ${magma_sums_dir}/MAGMA_sumstats/daytime_sleepiness_MAGMA.txt

## SLEEP DURATION ##
## OVERALL SLEEP DURATION, SHORT SLEEP DURATION, LONG SLEEP DURATION
#Input columns: SNP	CHR	BP	ALLELE1	ALLELE0	A1FREQ	INFO	BETA_SLEEPDURATION	SE_SLEEPDURATION	P_SLEEPDURATION
#Output columns: SNP	CHR	BP	A1	A2A1FREQ	INFO	BETA	SE	P
cp ${sumstats_dir}/behaviour/sleep_duration/overall_sleep_duration/sleepdurationsumstats.txt ${magma_sums_dir}/sleep_duration
cp ${sumstats_dir}/behaviour/sleep_duration/long_sleep_duration/longsumstats.txt ${magma_sums_dir}/sleep_duration
cp ${sumstats_dir}/behaviour/sleep_duration/short_sleep_duration/shortsumstats.txt ${magma_sums_dir}/sleep_duration
#Rename columns
sed -e '1s/ALLELE1/A1/' -e '1s/ALLELE0/A2/' -e '1s/BETA_SLEEPDURATION/BETA/' -e '1s/SE_SLEEPDURATION/SE/' -e '1s/P_SLEEPDURATION/P/' ${magma_sums_dir}/sleep_duration/sleepdurationsumstats.txt > ${magma_sums_dir}/sleep_duration/renamed_sleepdurationsumstats.txt
sed -e '1s/ALLELE1/A1/' -e '1s/ALLELE0/A2/' -e '1s/BETA_LONGSLEEP/BETA/' -e '1s/SE_LONGSLEEP/SE/' -e '1s/P_SLONGSLEEP/P/' ${magma_sums_dir}/sleep_duration/longsumstats.txt > ${magma_sums_dir}/sleep_duration/renamed_longsumstats.txt
sed -e '1s/ALLELE1/A1/' -e '1s/ALLELE0/A2/' -e '1s/BETA_SHORTSLEEP/BETA/' -e '1s/SE_SHORTSLEEP/SE/' -e '1s/P_SHORTSLEEP/P/' ${magma_sums_dir}/sleep_duration/shortsumstats.txt > ${magma_sums_dir}/sleep_duration/renamed_shortsumstats.txt
cp ${magma_sums_dir}/sleep_duration/renamed_sleepdurationsumstats.txt ${magma_sums_dir}/MAGMA_sumstats/overall_sleep_duration_MAGMA.txt
cp ${magma_sums_dir}/sleep_duration/renamed_longsumstats.txt ${magma_sums_dir}/MAGMA_sumstats/long_sleep_duration_MAGMA.txt
cp ${magma_sums_dir}/sleep_duration/renamed_shortsumstats.txt ${magma_sums_dir}/MAGMA_sumstats/short_sleep_duration_MAGMA.txt

## INSOMNIA ##
#Input columns: SNP CHR BP ALLELE1 ALLELE0 A1FREQ INFO BETA_INSOMNIA SE_INSOMNIA P_INSOMNIA
#Output columns: SNP CHR BP A1 A2 A1FREQ INFO BETA SE P
cp ${sumstats_dir}/behaviour/insomnia/Saxena_fullUKBB_Insomnia_summary_stats.txt ${magma_sums_dir}/insomnia
#Replace whitespace with tab
awk -v OFS="\t" '{$1=$1; print;}' ${magma_sums_dir}/insomnia/Saxena_fullUKBB_Insomnia_summary_stats.txt > ${magma_sums_dir}/insomnia/tab_Saxena_fullUKBB_Insomnia_summary_stats.txt
#Rename columns
sed -e '1s/ALLELE1/A1/' -e '1s/ALLELE0/A2/' -e '1s/BETA_INSOMNIA/BETA/' -e '1s/SE_INSOMNIA/SE/' -e '1s/P_INSOMNIA/P/' ${magma_sums_dir}/insomnia/tab_Saxena_fullUKBB_Insomnia_summary_stats.txt > ${magma_sums_dir}/insomnia/renamed_tab_Saxena_fullUKBB_Insomnia_summary_stats.txt
cp ${magma_sums_dir}/insomnia/renamed_tab_Saxena_fullUKBB_Insomnia_summary_stats.txt ${magma_sums_dir}/MAGMA_sumstats/insomnia_MAGMA.txt

## EDUCATIONAL ATTAINMENT ##
#Input columns: MarkerName	CHR	POS	A1	A2	EAF	Beta	SE	Pval
#Output columns: SNP	CHR	POS	A1	A2	EAF	BETA	SE	P
cp ${sumstats_dir}/behaviour/educational_attainment/GWAS_EA_excl23andMe.txt ${magma_sums_dir}/educational_attainment
#Rename columns
sed -e '1s/MarkerName/SNP/' -e '1s/POS/BP/' -e '1s/Beta/BETA/' -e '1s/Pval/P/' ${magma_sums_dir}/educational_attainment/GWAS_EA_excl23andMe.txt > ${magma_sums_dir}/educational_attainment/renamed_GWAS_EA_excl23andMe.txt
cp ${magma_sums_dir}/educational_attainment/renamed_GWAS_EA_excl23andMe.txt ${magma_sums_dir}/MAGMA_sumstats/educational_attainment_MAGMA.txt

## COGNITIVE PERFORMANCE ##
#Input columns: MarkerName	CHR	POS	A1	A2	EAF	Beta	SE	Pval
#Output columns: SNP	CHR	POS	A1	A2	EAF	BETA	SE	P
cp ${sumstats_dir}/behaviour/cognitive_performance/GWAS_CP_all.txt ${magma_sums_dir}/cognitive_performance
#Rename columns
sed -e '1s/MarkerName/SNP/' -e '1s/POS/BP/' -e '1s/Beta/BETA/' -e '1s/Pval/P/' ${magma_sums_dir}/cognitive_performance/GWAS_CP_all.txt > ${magma_sums_dir}/cognitive_performance/renamed_GWAS_CP_all.txt
cp ${magma_sums_dir}/cognitive_performance/renamed_GWAS_CP_all.txt ${magma_sums_dir}/MAGMA_sumstats/cognitive_performance_MAGMA.txt

## INTELLIGENCE ##
#Input columns: SNP	UNIQUE_ID	CHR	POS	A1	A2	EAF_HRC	Zscore	stdBeta	SE	P	N_analyzed	minINFO	EffectDirection
#Output columns: SNP	CHR	BP	A1	A2	UNIQUE_ID	EAF_HRC	Z	stdBeta	SE	P	N_analyzed	minINFO	EffectDirection
cp ${sumstats_dir}/behaviour/intelligence/sumstats/SavageJansen_2018_intelligence_metaanalysis.txt ${magma_sums_dir}/intelligence
#Rename columns
sed -e '1s/POS/BP/' -e '1s/Zscore/Z/' ${magma_sums_dir}/intelligence/SavageJansen_2018_intelligence_metaanalysis.txt > ${magma_sums_dir}/intelligence/renamed_SavageJansen_2018_intelligence_metaanalysis.txt
#Reorder columns
awk -v OFS="\t" '{print $1,$3,$4,$5,$6,$2,$7,$8,$9,$10,$11,$12,$13,$14}' ${magma_sums_dir}/intelligence/renamed_SavageJansen_2018_intelligence_metaanalysis.txt > ${magma_sums_dir}/intelligence/ordered_renamed_SavageJansen_2018_intelligence_metaanalysis.txt
cp ${magma_sums_dir}/intelligence/ordered_renamed_SavageJansen_2018_intelligence_metaanalysis.txt ${magma_sums_dir}/MAGMA_sumstats/intelligence_MAGMA.txt

## NEUROTICISM ##
#Input columns: RSID	SNP	CHR	POS	A1	A2	EAF_UKB	MAF_UKB	Z	P	N	INFO_UKB
#Output columns: SNP	CHR	BP	A1	A2 UNIQUE_ID	EAF_UKB	MAF_UKB	Z	P	N	INFO_UKB
cp ${sumstats_dir}/behaviour/neuroticism/sumstats_neuroticism_ctg_format.txt ${magma_sums_dir}/neuroticism
#Rename columns
sed -e '1s/SNP/UNIQUE_ID/' -e '1s/RSID/SNP/' -e '1s/POS/BP/' ${magma_sums_dir}/neuroticism/sumstats_neuroticism_ctg_format.txt >  ${magma_sums_dir}/neuroticism/renamed_sumstats_neuroticism_ctg_format.txt 
#Reorder columns
awk -v OFS="\t" '{print $1,$3,$4,$5,$6,$2,$7,$8,$9,$10,$11,$12}' ${magma_sums_dir}/neuroticism/renamed_sumstats_neuroticism_ctg_format.txt > ${magma_sums_dir}/neuroticism/ordered_renamed_sumstats_neuroticism_ctg_format.txt 
cp ${magma_sums_dir}/neuroticism/ordered_renamed_sumstats_neuroticism_ctg_format.txt ${magma_sums_dir}/MAGMA_sumstats/neuroticism_MAGMA.txt

## ORGANIZE ALL SUMSTAST ##
for gwas in ${magma_sums_dir}/*
do
	if [ $gwas == "${magma_sums_dir}/MAGMA_sumstats" ] || [ $gwas == "${magma_sums_dir}/GWAS" ] ; then
		continue;
	fi
	mv $gwas ${magma_sums_dir}/GWAS
done
