#!/usr/bin/Rscript

# TITLE:    MAGMA_celltype_analysis.R
# ABOUT:    Script to run MAGMA Celltyping analysis 42 phenotypes using KI and 10x data in both linear and top 10 mode using MAGMA v.1.08
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
formal_cell_type_names_KI <- c("Astrocytes/Ependymal", "Dopaminergic adult", "Dopaminergic neuroblast", 
                            "Embryonic dopaminergic neuron", "Embryonic GABAergic neuron", 
                            "Embryonic midbrain nucleus neurons", "Endothelial-mural", 
                            "Hypothalamic dopaminergic neurons", "Hypothalamic GABAergic neurons", 
                            "Hypothalamic glutamatergic neurons", "Interneurons", "Medium spiny neuron", 
                            "Microglia", "Neural progenitors", "Neuroblasts", "Oligodendrocyte precursor", 
                            "Oligodendrocytes", "Oxytocin/Vasopressin neurons", "Pyramidal (CA1)", 
                            "Pyramidal (SS)", "Radial glia like", "Serotonergic neuron", 
                            "Striatal interneuron", "Vascular leptomeningeal cells")
formal_cell_type_names_10X <- c("Glutamatergic neurons", "Neuroblasts 1", "Astrocytes 1", "Neuroblasts 2",
                                "Intermediate progenitors", "Enteric glial cells", "Interneurons 1", "Neurons",
                                "Interneurons 2", "Vascular endothelial cells", "Enteric neurons",
                                "Cajal Retzius cells", "Oligodendrocytes", "Astrocytes 2", "Endothelial cells",
                                "Microglia")

phenotypes <- read_excel("/hpc/hers_en/molislagers/MAGMA/sumstats/analysis/phenotypes_totalN.xlsx")
phenotype_names <- phenotypes$phenotype
phenotype_totalN <- phenotypes$totalN

########## Prepare quantile groups for cell types ##########
mean_expr_path_10xGenomics <- "/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment/10x_Genomics_16_cell_types_mean_expression.tsv"
mean_expr_data <- read.delim(mean_expr_path_10xGenomics, row.names = 1)
specificity_data <- mean_expr_data/rowSums(mean_expr_data)

########## Restructure ctd object for MAGMA cell type association analysis ##########
ctd_all10X <- MAGMA.Celltyping::ctd_allKI
# REMOVE IF DEEMED OBSOLETE
ctd_all10X[[1]]$mean_exp <- mean_expr_data
# REMOVE IF DEEMED OBSOLETE
ctd_all10X[[1]]$specificity <- specificity_data

ctd_KI = prepare.quantile.groups(ctd_allKI, specificity_species = 'mouse', numberOfBins = 40)
ctd_10X = prepare.quantile.groups(ctd_all10X, specificity_species = 'mouse', numberOfBins = 40)

#Loop over phenotypes
for (i in 1:length(phenotype_names)) {
    hm3_mhc_filename <- paste(phenotype_names[i], "_hm3_mhc_restricted_MAGMA.txt", sep = "")
    hm3_mhc_sumstats <- paste(sum_stats_dir, "/QC_MAGMA/hm3_restricted_MHC_excluded/", hm3_mhc_filename, sep="")
    ## Map SNPs to genes
    hm3_mhc_mapped_genes = map.snps.to.genes_costum(hm3_mhc_sumstats, genome_ref_path = genome_ref_path, N = phenotype_totalN[i])
    #Perform linear cell type association analysis using KI data
    hm3_mhc_linear_assoc_KI = calculate_celltype_associations(ctd_KI, hm3_mhc_sumstats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
    #Perform linear cell type association analysis using 10X data
    hm3_mhc_linear_assoc_10X = calculate_celltype_associations(ctd_10X, hm3_mhc_sumstats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Linear')
    #Perform top 10 cell type association analysis using KI data
    hm3_mhc_top_10_assoc_KI = calculate_celltype_associations(ctd_KI, hm3_mhc_sumstats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
    #Perform top 10 cell type association analysis using 10X data
    hm3_mhc_top_10_assoc_10X = calculate_celltype_associations(ctd_10X, hm3_mhc_sumstats, genome_ref_path = genome_ref_path, specificity_species = 'mouse', EnrichmentMode = 'Top 10%')
    #Update cell types
    hm3_mhc_linear_assoc_KI[[1]]$results$Celltype = formal_cell_type_names_KI
    hm3_mhc_linear_assoc_10X[[1]]$results$Celltype = formal_cell_type_names_10X
    hm3_mhc_top_10_assoc_KI[[1]]$results$Celltype = formal_cell_type_names_KI
    hm3_mhc_top_10_assoc_10X[[1]]$results$Celltype = formal_cell_type_names_10X
    #Save results
    write.table(hm3_mhc_linear_assoc_KI[[1]]$results, file = paste(output_dir, "/KI/hm3_mhc_restricted/linear/", phenotype_names[i], ".tsv", sep = ""), sep = '\t', row.names = FALSE)
    write.table(hm3_mhc_linear_assoc_10X[[1]]$results, file = paste(output_dir, "/10x_genomics/hm3_mhc_restricted/linear/", phenotype_names[i], ".tsv", sep = ""), sep = '\t', row.names = FALSE)
    write.table(hm3_mhc_top_10_assoc_KI[[1]]$results, file = paste(output_dir, "/KI/hm3_mhc_restricted/top10/", phenotype_names[i], ".tsv", sep = ""), sep = '\t', row.names = FALSE)
    write.table(hm3_mhc_top_10_assoc_10X[[1]]$results, file = paste(output_dir, "/10x_genomics/hm3_mhc_restricted/top10/", phenotype_names[i], ".tsv", sep = ""), sep = '\t', row.names = FALSE)

}
