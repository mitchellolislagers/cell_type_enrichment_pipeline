!#/usr/bin/bash
#	JOB DEFINITIONS
#########################################
#$ -N run_MAGMA
#$ -S /bin/bash
#$ -cwd
#$ -e /home/hers_en/molislagers/errors/
#$ -o /home/hers_en/molislagers/output/
#$ -M mitchellolislagers@gmail.com
#$ -m beas
#$ -l h_rt=40:00:00
#$ -l h_vmem=30G
#########################################

module load R/3.6.3
cd /hpc/hers_en/molislagers/MAGMA/sumstats/analysis/v1.08/scripts
Rscript ./MAGMA_v108_celltype_analysis_SCZ_all_QC.R
