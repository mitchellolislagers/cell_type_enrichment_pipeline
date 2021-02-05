#!/usr/bin/Rscript

# TITLE:    MAGMA_celltype_analysis_hm3_mhc_restricted_KI_level2_v33.R
# ABOUT:    Script to run MAGMA Celltyping analysis 42 phenotypes using KI level 2 data in both linear and top 10 mode using MAGMA v1.08
# INPUT:    10X Genomics mean expression data
#           42 phenotypes as generated using the prepare_MAGMA_sumstats.sh and QC_MAGMA.sh scripts
# AUTHOR:   Mitchell Olislagers
# DATE:     10 Nov 2020

########## Load required libraries ##########
library(One2One)
library(MAGMA.Celltyping)
library(ggplot2)
library(readxl)

#Creating a costum function to map SNPs to genes
map.snps.to.genes_costum <- function(path_formatted,upstream_kb=10,downstream_kb=1.5,N=NULL,genome_ref_path){
    path_formatted = path.expand(path_formatted)
    magmaPaths = get.magma.paths(path_formatted,upstream_kb,downstream_kb)
    
    # Check whether there is an N column in the sumstats file (if it wasn't provided as an argument)
    if(is.null(N)){
        con <- file(path_formatted,"r") ; first_line <- readLines(con,n=1) ; close(con)
        column_headers = strsplit(first_line,"\t")[[1]]
        if("N" %in% column_headers){n_arg = "ncol=N"}else{
            nval <- as.numeric(readline("There is no N column within the sumstats file. What is the N value for this GWAS?"))
            
            if(is.na(nval)){stop(sprintf("%s provided but value of N for the GWAS must be numeric",nval))}
            if(nval<1000){stop("Value of N provided is less than 1000. This seems unlikely.")}
            if(nval>100000000){stop("Value of N provided is over than 100000000. In 2018 this seems unlikely.")}
            n_arg = sprintf("N=%s",nval)
        }
    }else{
        n_arg = sprintf("N=%s",N)
    }
    
    # Determine which genome build it uses & get path to gene loc file
    genome_build = get_genomebuild_for_sumstats(path_formatted)
    gene_loc_dir = sprintf("%s/extdata",system.file(package="MAGMA.Celltyping"))
    if(genome_build == "GRCh37"){genomeLocFile=sprintf("%s/NCBI37.3.gene.loc",gene_loc_dir)}
    if(genome_build == "GRCh38"){genomeLocFile=sprintf("%s/NCBI38.gene.loc",gene_loc_dir)}
    print(sprintf("GWAS Sumstats appear to come from genome build: %s",genome_build))
    
    #sumstatsPrefix = sprintf("%s.%sUP.%sDOWN",path_formatted,upstream_kb,downstream_kb)
    #sumstatsPrefix = 
    outPath = gsub("\\/$","",magmaPaths$filePathPrefix) # Remove a trailing slash to avoid errors on windows
    magma_cmd = sprintf("magma --annotate window=%s,%s --snp-loc '%s' --gene-loc '%s' --out '%s'",upstream_kb,downstream_kb,path_formatted,genomeLocFile,outPath)
    system(magma_cmd)
    
    # SCHIZ CLOZUK N=35802
    magma_cmd = sprintf("magma --bfile '%s' --pval '%s' %s --gene-annot '%s.genes.annot' --out '%s'",path.expand(genome_ref_path),path_formatted,n_arg,magmaPaths$filePathPrefix,outPath)
    #magma_cmd
    system(magma_cmd)
    
    # Return path to genes.out file
    return(sprintf("%s.genes.out",path_formatted))
}

########## Set variables ##########
genome_ref_dir = '/hpc/hers_en/molislagers/MAGMA/ref_data/g1000_eur'
if(!file.exists(sprintf('%s/g1000_eur.bed', genome_ref_dir))){
    download.file('https://ctg.cncr.nl/software/MAGMA/ref_data/g1000_eur.zip',destfile=sprintf('%s.zip', genome_ref_dir))
    unzip(sprintf('%s.zip',genome_ref_dir), exdir=genome_ref_dir)
}

genome_ref_path = sprintf('%s/g1000_eur', genome_ref_dir)
sum_stats_dir = '/hpc/hers_en/molislagers/MAGMA/sumstats'
output_dir = '/hpc/hers_en/molislagers/MAGMA/sumstats/analysis/v1.08'
formal_cell_type_names_KI_level2 <- c("Astro1", "Astro2", "CA1Pyr1", "CA1Pyr2", "CA1PyrInt", "CA2Pyr2", "Choroid", "ClauPyr", "Committed oligodendrocyte precursors", "Dopaminergic Adult-Substantia nigra", 
                                    "Dopaminergic Adult-Ventral tegmental area1", "Dopaminergic Adult-Ventral tegmental area2", "Dopaminergic Adult-Ventral tegmental area3", "Dopaminergic Adult-Ventral tegmental area4", 
                                    "Dopaminergic Neuroblast", "Embryonic Dopaminergic Neuron 0", "Embryonic Dopaminergic Neuron 1", "Embryonic Dopaminergic Neuron 2", "Embryonic GABAergic Neuron 1a", 
                                    "Embryonic GABAergic Neuron 1b", "Embryonic GABAergic Neuron 2", "Epend", "Hypothalamic Adcyap1; Tac1 (VMH) Neuron", "Hypothalamic Adcyap1(VMH) Neuron", "Hypothalamic Avp-high;Gal;Oxt-low Neuron", 
                                    "Hypothalamic Avp-high;Oxt-low Neuron", "Hypothalamic Avp-medium Neuron", "Hypothalamic Calcr-high;VMAT+and-;Lhx1 (Sch) Neuron", "Hypothalamic Crh+and-; Lhx6;GABA (BST; MPO) Neuron", 
                                    "Hypothalamic Crh+and-;GABA;Pgr15l (PVH-PVH-galo) Neuron", "Hypothalamic Crh+and-;Gal-low (PVH) Neuron", "Hypothalamic Dopamine; Dat;Nmur2;GABA Neuron", 
                                    "Hypothalamic Dopamine;Tac+and- Gad1; GABA Neuron", "Hypothalamic Dopamine;Tac1;Ghrh;Pnoc; Dat+and-; GABA Neuron", "Hypothalamic GABA 1 Neuron", "Hypothalamic GABA 2 Neuron", 
                                    "Hypothalamic GABA; Gucy1a3 Neuron", "Hypothalamic GABA;Pnoc;Tac2+and- Neuron", "Hypothalamic Gad-low;Gnrh-and+ Neuron", "Hypothalamic Gad1;Gad2;VGAT+and-;Pnoc Neuron", 
                                    "Hypothalamic Galanin Neuron", "Hypothalamic Ghrh-high;Th;Vglut2 Neuron", "Hypothalamic Hcrt Neuron", "Hypothalamic Hmit+and- Neuron", 
                                    "Hypothalamic Nms;VIP+and-; Avp-low+and-;circadian (SCH) Neuron", "Hypothalamic Npr2;Gm5595;4930422G04Rik;Tnr Neuron", 
                                    "Hypothalamic Npvf Neuron", "Hypothalamic Npy-medium;Gad1;Gad2 Neuron", "Hypothalamic Npy;Agrp (ARH) Neuron", "Hypothalamic Nts;Gal;Pnoc;Tac2+and-;GABA Neuron", 
                                    "Hypothalamic Nts;Pnoc;Tac2+and-;GABA Neuron", "Hypothalamic Otof;Lhx1+and- Neuron", "Hypothalamic Oxt;Avp-low Neuron", "Hypothalamic Oxt;Avp-low;Th;Cacna1H Neuron", 
                                    "Hypothalamic Oxt;Avp-medium;Gad2-low Neuron", "Hypothalamic Oxt;Avp-medium;Th+and- PDYN+and Neuron", "Hypothalamic Per2;circadian (SCH) Neuron", "Hypothalamic Pmch Neuron", 
                                    "Hypothalamic Pomc+and- (ARH) Neuron", "Hypothalamic Qrfp Neuron", "Hypothalamic Sst low-medium Neuron", "Hypothalamic Sst-high;Cartpt;Galr1 Neuron", "Hypothalamic Sst-medium;Cartpt Neuron", 
                                    "Hypothalamic Th;Tac1;Ghrh;Vmat2-lowand- GABA Neuron", "Hypothalamic Trh-high;Adcyap1;Cartpt Neuron", "Hypothalamic Trh-low Neuron", "Hypothalamic Trh-medium Neuron", 
                                    "Hypothalamic Vglut2 ;Col9a2 (PVH) Neuron", "Hypothalamic Vglut2 A Neuron", "Hypothalamic Vglut2;A930013F10Rik;Pou2f2 Neuron", "Hypothalamic Vglut2;Bdnf;Gad1 Neuron", 
                                    "Hypothalamic Vglut2;Cck;2AG;AEA;Npr2;Npy1r Neuron", "Hypothalamic Vglut2;Cnr1;Ninl;Rfx5;Zfp346 Neuron", "Hypothalamic Vglut2;Cnr1;Npr2;Cacna1h Neuron", 
                                    "Hypothalamic Vglut2;Gad1-low;Crh-lowand- (VMH) Neuron", "Hypothalamic Vglut2;Gpr149 (PVH;VMH) Neuron", "Hypothalamic Vglut2;Hcn1;6430411K18Rik Neuron", 
                                    "Hypothalamic Vglut2;Morn4;Prrc2a Neuron", "Hypothalamic Vglut2;Myt1;Lhx9 Neuron", "Hypothalamic Vglut2;Penk;Oprk1 Neuron", "Hypothalamic Vglut2;Pgam;Snx12 Neuron", 
                                    "Hypothalamic Vglut2;Prmt8;Ugdh (VMH) Neuron", "Hypothalamic Vglut2;Zfp458;Ppp1r12b;Cacna1e-highest Neuron", "Hypothalamic Vip;Grp+and-;circadian (SCH) Neuron", "Int Pvalb", 
                                    "Int1", "Int10", "Int11", "Int12", "Int13", "Int14", "Int15", "Int16", "Int2", "Int4", "Int5", "Int6", "Int7", "Int8", "Int9", "Lateral Neuroblasts 1", "Lateral Neuroblasts 2", 
                                    "Mature oligodendrocytes 1", "Mature oligodendrocytes 2", "Mature oligodendrocytes 3", "Mature oligodendrocytes 4", "Mature oligodendrocytes 5", "Mature oligodendrocytes 6", 
                                    "Medial Neuroblasts ", "Mediolateral Neuroblasts 1", "Mediolateral Neuroblasts 2", "Mediolateral Neuroblasts 3", "Mediolateral Neuroblasts 4", "Mediolateral Neuroblasts 5", 
                                    "Medium Spiny Neuron D1R", "Medium Spiny Neuron D2R", "Mgl1", "Mgl2", "Myelin-forming oligodendrocytes 1", "Myelin-forming oligodendrocytes 2", "Neural Progenitors", 
                                    "Newly formed oligodendrocytes 1", "Newly formed oligodendrocytes 2", "Oculomotor and Trochlear nucleus embryonic neurons", "Oligodendrocyte Precursor", "Peric", "Pvm1", "Pvm2", 
                                    "Radial glia like cells 1", "Radial glia like cells 2", "Radial glia like cells 3", "Red nucleus embryonic neurons", "S1PyrDL", "S1PyrL23", "S1PyrL4", "S1PyrL5", "S1PyrL5a", "S1PyrL6", 
                                    "S1PyrL6b", "Serotonergic Neuron", "Striatal CHAT Interneuron", "Striatal Interneurons (other)", "Striatal Pvalb Interneuron", "Striatal Sst Interneuron", "SubPyr", 
                                    "Vascular Leptomeningeal Cells", "Vend1", "Vend2", "Vsmc")

phenotypes <- read_excel("/hpc/hers_en/molislagers/MAGMA/sumstats/analysis/phenotypes_totalN.xlsx")
phenotype_names <- phenotypes$phenotype
phenotype_totalN <- phenotypes$totalN

########## Remove level 1 data from KI gene expression ##########
ctd_KI = prepare.quantile.groups(ctd_allKI, specificity_species = 'mouse', numberOfBins = 40)
ctd_KI[[1]] <- NULL

#Loop over phenotypes
for (i in 1:length(phenotype_names)) {
    hm3_mhc_filename <- paste(phenotype_names[i], "_hm3_mhc_restricted_MAGMA.txt", sep = "")
    #TEMPORARY SUMSTATS FOR PARRALELL RUNNING
    hm3_mhc_sumstats <- paste(sum_stats_dir, "/QC_MAGMA/hm3_restricted_MHC_excluded_level2/", hm3_mhc_filename, sep="")
    ## Map SNPs to genes
    hm3_mhc_mapped_genes = map.snps.to.genes_costum(hm3_mhc_sumstats, genome_ref_path = genome_ref_path, N = phenotype_totalN[i])
    #Perform linear cell type association analysis using KI level 2 data
    hm3_mhc_linear_assoc_KI_level2 = calculate_celltype_associations(ctd_KI, hm3_mhc_sumstats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
    #Perform top 10 cell type association analysis using KI data
    hm3_mhc_top_10_assoc_KI_level2 = calculate_celltype_associations(ctd_KI, hm3_mhc_sumstats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
    #Update cell types
    hm3_mhc_linear_assoc_KI_level2[[1]]$results$Celltype = formal_cell_type_names_KI_level2
    hm3_mhc_top_10_assoc_KI_level2[[1]]$results$Celltype = formal_cell_type_names_KI_level2
    #Save results
    write.table(hm3_mhc_linear_assoc_KI_level2[[1]]$results, file = paste(output_dir, "/KI_level2/hm3_mhc_restricted/linear/", phenotype_names[i], ".tsv", sep = ""), sep = '\t', row.names = FALSE)
    write.table(hm3_mhc_top_10_assoc_KI_level2[[1]]$results, file = paste(output_dir, "/KI_level2/hm3_mhc_restricted/top10/", phenotype_names[i], ".tsv", sep = ""), sep = '\t', row.names = FALSE)

}
