#!/bin/bash
  
# TITLE:        partition_h2.sh
# GOAL:         Partition heritability (h2) to specificity deciles of cell types.
# INPUT:        LDSC thin-annotation files (annot.gz) per autosomal chromosome for all cell types, LDSC LD score files 
# INPUT:	1000 Genomes Project phase 3 baseline model LD scores, regression weights, allele frequencies (see README)
# OUTPUT:       Partitioned h2 output files
# AUTHOR:       Koen Rademaker
# DATE:         20 February 2020


########## Initialize script ##########
# Set paths
ldsc_dir=/hpc/hers_en/molislagers/LDSC/ldsc
ref_dir=/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment
sum_stats_dir=/hpc/hers_en/molislagers/LDSC/summary_statistics/munged_sumstats
annotation_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/annotations_v33
ld_scores_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/ld_scores_v33
out_dir=/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/partitioned_h2_v33
miniconda_dir=~/miniconda2

# Declare variables
declare -a sum_stats=("ADHD" "AN" "anxiety" "ASD" "BIP" "cross" "MDD" "OCD" "PTSD" "SCZ" "TS" "alcohol_use" "alcohol_dependence" "drinks_pw" "cannabis" "smoking_initiation" "ever_smoked" "cigarettes_pd" "smoking_cessation" "ALS" "alzheimers" "all_epilepsy" "generalized" "focal" "all_stroke" "cardioembolic" "ischemic" "large_artery" "small_vessel" "parkinson" "height" "BMI" "chronotype" "daytime_sleepiness" "overall_sleep_duration" "short_sleep_duration" "long_sleep_duration" "insomnia" "intelligence" "educational_attainment" "cognitive_performance" "neuroticism")

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
"Hypothalamic_Vglut2;Prmt8;Ugdh_(VMH)_Neuron" "Hypothalamic_Vglut2;Zfp458;Ppp1r12b;Cacna1e-highest_Neuron" "Hypothalamic_Vip;Grp+and-;circadian_(SCH)_Neuron" "Int_Pvalb" 
"Int1" "Int10" "Int11" "Int12" "Int13" "Int14" "Int15" "Int16" "Int2" "Int4" "Int5" "Int6" "Int7" "Int8" "Int9" "Lateral_Neuroblasts_1" "Lateral_Neuroblasts_2" 
"Mature_oligodendrocytes_1" "Mature_oligodendrocytes_2" "Mature_oligodendrocytes_3" "Mature_oligodendrocytes_4" "Mature_oligodendrocytes_5" "Mature_oligodendrocytes_6" 
"Medial_Neuroblasts_" "Mediolateral_Neuroblasts_1" "Mediolateral_Neuroblasts_2" "Mediolateral_Neuroblasts_3" "Mediolateral_Neuroblasts_4" "Mediolateral_Neuroblasts_5" 
"Medium_Spiny_Neuron_D1R" "Medium_Spiny_Neuron_D2R" "Mgl1" "Mgl2" "Myelin-forming_oligodendrocytes_1" "Myelin-forming_oligodendrocytes_2" "Neural_Progenitors" 
"Newly_formed_oligodendrocytes_1" "Newly_formed_oligodendrocytes_2" "Oculomotor_and_Trochlear_nucleus_embryonic_neurons" "Oligodendrocyte_Precursor" "Peric" "Pvm1" "Pvm2" 
"Radial_glia_like_cells_1" "Radial_glia_like_cells_2" "Radial_glia_like_cells_3" "Red_nucleus_embryonic_neurons" "S1PyrDL" "S1PyrL23" "S1PyrL4" "S1PyrL5" "S1PyrL5a" "S1PyrL6" 
"S1PyrL6b" "Serotonergic_Neuron" "Striatal_CHAT_Interneuron" "Striatal_Interneurons_(other)" "Striatal_Pvalb_Interneuron" "Striatal_Sst_Interneuron" "SubPyr" 
"Vascular_Leptomeningeal_Cells" "Vend1" "Vend2" "Vsmc")

dataset="KI"


########## Organize LDSC ##########
#wget https://data.broadinstitute.org/alkesgroup/LDSCORE/1000G_Phase3_frq.tgz -P ${ref_dir}
#tar -xvzf 1000G_Phase3_frq.tgz
# Copy annotation files and LD score files to reference folder
ref_ld_dir=${ref_dir}/ref_ld_level2
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
				--ref-ld-chr ${ref_ld_dir}/${dataset}_${cell_types[0]}.,${ref_ld_dir}/${dataset}_${cell_types[1]}.,${ref_ld_dir}/${dataset}_${cell_types[2]}.,${ref_ld_dir}/${dataset}_${cell_types[3]}.,${ref_ld_dir}/${dataset}_${cell_types[4]}.,${ref_ld_dir}/${dataset}_${cell_types[5]}.,${ref_ld_dir}/${dataset}_${cell_types[6]}.,${ref_ld_dir}/${dataset}_${cell_types[7]}.,${ref_ld_dir}/${dataset}_${cell_types[8]}.,${ref_ld_dir}/${dataset}_${cell_types[9]}.,${ref_ld_dir}/${dataset}_${cell_types[10]}.,${ref_ld_dir}/${dataset}_${cell_types[11]}.,${ref_ld_dir}/${dataset}_${cell_types[12]}.,${ref_ld_dir}/${dataset}_${cell_types[13]}.,${ref_ld_dir}/${dataset}_${cell_types[14]}.,${ref_ld_dir}/${dataset}_${cell_types[15]}.,${ref_ld_dir}/${dataset}_${cell_types[16]}.,${ref_ld_dir}/${dataset}_${cell_types[17]}.,${ref_ld_dir}/${dataset}_${cell_types[18]}.,${ref_ld_dir}/${dataset}_${cell_types[19]}.,${ref_ld_dir}/${dataset}_${cell_types[20]}.,${ref_ld_dir}/${dataset}_${cell_types[21]}.,${ref_ld_dir}/${dataset}_${cell_types[22]}.,${ref_ld_dir}/${dataset}_${cell_types[23]}.,${ref_ld_dir}/${dataset}_${cell_types[24]}.,${ref_ld_dir}/${dataset}_${cell_types[25]}.,${ref_ld_dir}/${dataset}_${cell_types[26]}.,${ref_ld_dir}/${dataset}_${cell_types[27]}.,${ref_ld_dir}/${dataset}_${cell_types[28]}.,${ref_ld_dir}/${dataset}_${cell_types[29]}.,${ref_ld_dir}/${dataset}_${cell_types[30]}.,${ref_ld_dir}/${dataset}_${cell_types[31]}.,${ref_ld_dir}/${dataset}_${cell_types[32]}.,${ref_ld_dir}/${dataset}_${cell_types[33]}.,${ref_ld_dir}/${dataset}_${cell_types[34]}.,${ref_ld_dir}/${dataset}_${cell_types[35]}.,${ref_ld_dir}/${dataset}_${cell_types[36]}.,${ref_ld_dir}/${dataset}_${cell_types[37]}.,${ref_ld_dir}/${dataset}_${cell_types[38]}.,${ref_ld_dir}/${dataset}_${cell_types[39]}.,${ref_ld_dir}/${dataset}_${cell_types[40]}.,${ref_ld_dir}/${dataset}_${cell_types[41]}.,${ref_ld_dir}/${dataset}_${cell_types[42]}.,${ref_ld_dir}/${dataset}_${cell_types[43]}.,${ref_ld_dir}/${dataset}_${cell_types[44]}.,${ref_ld_dir}/${dataset}_${cell_types[45]}.,${ref_ld_dir}/${dataset}_${cell_types[46]}.,${ref_ld_dir}/${dataset}_${cell_types[47]}.,${ref_ld_dir}/${dataset}_${cell_types[48]}.,${ref_ld_dir}/${dataset}_${cell_types[49]}.,${ref_ld_dir}/${dataset}_${cell_types[50]}.,${ref_ld_dir}/${dataset}_${cell_types[51]}.,${ref_ld_dir}/${dataset}_${cell_types[52]}.,${ref_ld_dir}/${dataset}_${cell_types[53]}.,${ref_ld_dir}/${dataset}_${cell_types[54]}.,${ref_ld_dir}/${dataset}_${cell_types[55]}.,${ref_ld_dir}/${dataset}_${cell_types[56]}.,${ref_ld_dir}/${dataset}_${cell_types[57]}.,${ref_ld_dir}/${dataset}_${cell_types[58]}.,${ref_ld_dir}/${dataset}_${cell_types[59]}.,${ref_ld_dir}/${dataset}_${cell_types[60]}.,${ref_ld_dir}/${dataset}_${cell_types[61]}.,${ref_ld_dir}/${dataset}_${cell_types[62]}.,${ref_ld_dir}/${dataset}_${cell_types[63]}.,${ref_ld_dir}/${dataset}_${cell_types[64]}.,${ref_ld_dir}/${dataset}_${cell_types[65]}.,${ref_ld_dir}/${dataset}_${cell_types[66]}.,${ref_ld_dir}/${dataset}_${cell_types[67]}.,${ref_ld_dir}/${dataset}_${cell_types[68]}.,${ref_ld_dir}/${dataset}_${cell_types[69]}.,${ref_ld_dir}/${dataset}_${cell_types[70]}.,${ref_ld_dir}/${dataset}_${cell_types[71]}.,${ref_ld_dir}/${dataset}_${cell_types[72]}.,${ref_ld_dir}/${dataset}_${cell_types[73]}.,${ref_ld_dir}/${dataset}_${cell_types[74]}.,${ref_ld_dir}/${dataset}_${cell_types[75]}.,${ref_ld_dir}/${dataset}_${cell_types[76]}.,${ref_ld_dir}/${dataset}_${cell_types[77]}.,${ref_ld_dir}/${dataset}_${cell_types[78]}.,${ref_ld_dir}/${dataset}_${cell_types[79]}.,${ref_ld_dir}/${dataset}_${cell_types[80]}.,${ref_ld_dir}/${dataset}_${cell_types[81]}.,${ref_ld_dir}/${dataset}_${cell_types[82]}.,${ref_ld_dir}/${dataset}_${cell_types[83]}.,${ref_ld_dir}/${dataset}_${cell_types[84]}.,${ref_ld_dir}/${dataset}_${cell_types[85]}.,${ref_ld_dir}/${dataset}_${cell_types[86]}.,${ref_ld_dir}/${dataset}_${cell_types[87]}.,${ref_ld_dir}/${dataset}_${cell_types[88]}.,${ref_ld_dir}/${dataset}_${cell_types[89]}.,${ref_ld_dir}/${dataset}_${cell_types[90]}.,${ref_ld_dir}/${dataset}_${cell_types[91]}.,${ref_ld_dir}/${dataset}_${cell_types[92]}.,${ref_ld_dir}/${dataset}_${cell_types[93]}.,${ref_ld_dir}/${dataset}_${cell_types[94]}.,${ref_ld_dir}/${dataset}_${cell_types[95]}.,${ref_ld_dir}/${dataset}_${cell_types[96]}.,${ref_ld_dir}/${dataset}_${cell_types[97]}.,${ref_ld_dir}/${dataset}_${cell_types[98]}.,${ref_ld_dir}/${dataset}_${cell_types[99]}.,${ref_ld_dir}/${dataset}_${cell_types[100]}.,${ref_ld_dir}/${dataset}_${cell_types[101]}.,${ref_ld_dir}/${dataset}_${cell_types[102]}.,${ref_ld_dir}/${dataset}_${cell_types[103]}.,${ref_ld_dir}/${dataset}_${cell_types[104]}.,${ref_ld_dir}/${dataset}_${cell_types[105]}.,${ref_ld_dir}/${dataset}_${cell_types[106]}.,${ref_ld_dir}/${dataset}_${cell_types[107]}.,${ref_ld_dir}/${dataset}_${cell_types[108]}.,${ref_ld_dir}/${dataset}_${cell_types[109]}.,${ref_ld_dir}/${dataset}_${cell_types[110]}.,${ref_ld_dir}/${dataset}_${cell_types[111]}.,${ref_ld_dir}/${dataset}_${cell_types[112]}.,${ref_ld_dir}/${dataset}_${cell_types[113]}.,${ref_ld_dir}/${dataset}_${cell_types[114]}.,${ref_ld_dir}/${dataset}_${cell_types[115]}.,${ref_ld_dir}/${dataset}_${cell_types[116]}.,${ref_ld_dir}/${dataset}_${cell_types[117]}.,${ref_ld_dir}/${dataset}_${cell_types[118]}.,${ref_ld_dir}/${dataset}_${cell_types[119]}.,${ref_ld_dir}/${dataset}_${cell_types[120]}.,${ref_ld_dir}/${dataset}_${cell_types[121]}.,${ref_ld_dir}/${dataset}_${cell_types[122]}.,${ref_ld_dir}/${dataset}_${cell_types[123]}.,${ref_ld_dir}/${dataset}_${cell_types[124]}.,${ref_ld_dir}/${dataset}_${cell_types[125]}.,${ref_ld_dir}/${dataset}_${cell_types[126]}.,${ref_ld_dir}/${dataset}_${cell_types[127]}.,${ref_ld_dir}/${dataset}_${cell_types[128]}.,${ref_ld_dir}/${dataset}_${cell_types[129]}.,${ref_ld_dir}/${dataset}_${cell_types[130]}.,${ref_ld_dir}/${dataset}_${cell_types[131]}.,${ref_ld_dir}/${dataset}_${cell_types[132]}.,${ref_ld_dir}/${dataset}_${cell_types[133]}.,${ref_ld_dir}/${dataset}_${cell_types[134]}.,${ref_ld_dir}/${dataset}_${cell_types[135]}.,${ref_ld_dir}/${dataset}_${cell_types[136]}.,${ref_ld_dir}/${dataset}_${cell_types[137]}.,${ref_ld_dir}/${dataset}_${cell_types[138]}.,${ref_ld_dir}/${dataset}_${cell_types[139]}.,${ref_ld_dir}/${dataset}_${cell_types[140]}.,${ref_ld_dir}/${dataset}_${cell_types[141]}.,${ref_ld_dir}/${dataset}_${cell_types[142]}.,${ref_ld_dir}/${dataset}_${cell_types[143]}.,${ref_ld_dir}/${dataset}_${cell_types[144]}.,${ref_ld_dir}/${dataset}_${cell_types[145]}.,${ref_ld_dir}/${dataset}_${cell_types[146]}.,${ref_ld_dir}/${dataset}_${cell_types[147]}.,${ref_ld_dir}/${dataset}_${cell_types[148]}.,${baseline_dir}/baseline. \
				--thin-annot \
                --overlap-annot \
                --frqfile-chr ${frq_dir}/1000G.EUR.QC. \
                --out ${out_dir}/${gwas}_combined_model \
                --print-coefficients
done

