# TITLE:	make_gene_set_files.R
# ABOUT:	Script to prepare LDSC annotation data by creating gene set files for sub-annotations.
# INPUT:	mean_expr_path: Mean gene expression data for cell types.
# INPUT:    human_gene_reference_path: Reference file for human genes (hg19), including chromosomes, gene names and gene IDs.
# AUTHOR:	Koen Rademaker, adapted by Mitchell Olislagers
# DATE:		20 February 2020

########## Load required libraries ##########
library(MAGMA.Celltyping)
library(One2One)
library(ggplot2)
library(dplyr)
library(tibble)

########## Declare custom functions ##########
get_specificity_decile_genes <- function(ctd, cell_type, decile){
    genes <- as.data.frame(ctd[[1]]$specificity_quantiles[, cell_type]) %>%
        rownames_to_column('gene') %>%
        filter(ctd[[1]]$specificity_quantiles[, cell_type] == decile) %>%
        pull(gene)
    return(genes)
}

get_homolog_translation <- function(){
    # Load orthology data
    all_homologs = load.homologs()
    ortholog_data = analyse.orthology('mouse', 'human', all_homologs)$orthologs_one2one
    # Select human and mouse entries
    selected_homologs <- all_homologs %>%
        filter(all_homologs$'HomoloGene ID' %in% ortholog_data$'HomoloGene ID' & all_homologs$'Common Organism Name' %in% c('mouse, laboratory', 'human')) %>%
        select('HomoloGene ID', 'Common Organism Name', 'Symbol')
    # Join mouse-human homologs together
    selected_homologs <- full_join(selected_homologs %>% filter(selected_homologs$`Common Organism Name` == 'mouse, laboratory'),
                                   selected_homologs %>% filter(selected_homologs$`Common Organism Name` == 'human'),
                                   by = 'HomoloGene ID')
    # Rename columns and set mouse gene symbols as row names
    selected_homologs <- selected_homologs %>%
        select('HomoloGene ID', 'Symbol.x', 'Symbol.y') %>%
        rename('Symbol_Mouse' = 'Symbol.x', 'Symbol_Human' = 'Symbol.y')
    selected_homologs <- selected_homologs %>% distinct()
    row.names(selected_homologs) <- selected_homologs$'Symbol_Mouse'
    selected_homologs <- selected_homologs %>%
        select('HomoloGene ID', 'Symbol_Human')
    return(selected_homologs)
}

translate_mouse_symbol_to_ensembl <- function(genes_mouse_symbol){
    # Convert mouse gene symbols to human gene symbols
    genes_mouse_symbol_df <- as.data.frame(genes_mouse_symbol)
    genes_mouse_symbol_df[] <- setNames(homolog_translation$Symbol_Human, row.names(homolog_translation))[as.character(unlist(genes_mouse_symbol_df))]
    genes_human_symbol <- as.vector(unlist(genes_mouse_symbol_df))
    # Convert human gene symbols to human gene Ensembl IDs
    genes_human_symbol_df <- as.data.frame(genes_human_symbol)
    genes_human_symbol_df[] <- setNames(human_gene_reference$GENE_ID, human_gene_reference$GENE_NAME)[as.character(unlist(genes_human_symbol_df))]
    genes_human_ensembl <- as.vector(unlist(genes_human_symbol_df))
    return(genes_human_ensembl)
}

partition_genes_by_chromosome <- function(genes, dataset, cell_type, decile, out_dir=getwd()){
    if(!dir.exists(out_dir)){
        dir.create(out_dir)
    }
    for (chr in 1:22){
        # Partition genes per chromosome
        chr_partitioned_genes <- human_gene_reference %>%
            filter(human_gene_reference$CHROMOSOME == chr & human_gene_reference$GENE_ID %in% genes) %>%
            select('GENE_ID')
        # Write file
        f <- file(paste(out_dir, '/', dataset, '.', cell_type, '.', chr, '.', decile, '.GeneSet', sep = ''), open = 'w')
        write.table(chr_partitioned_genes, file = f, sep = '\t', row.names = F, col.names = F, quote = F)
        close(f)
    }
}

########## Set variables ##########
human_gene_reference_path <- '/hpc/hers_en/molislagers/LDSC/ref_data/celltype_enrichment/hg19_reference/gencode.v33lift37_genes_only.tab'
geneset_out_path = '/hpc/hers_en/molislagers/LDSC/celltype_enrichment/KI_level2/geneset_v33'

########## Load human gene reference data and remove Ensembl ID version numbers ##########
human_gene_reference <- read.delim(human_gene_reference_path)
human_gene_reference$GENE_ID <- sapply(strsplit(as.character(human_gene_reference$GENE_ID), '\\.'), `[`, 1)

########## Load homolog translation table ##########
homolog_translation <- get_homolog_translation()

#ANALYSIS FOR KAROLINSKA INSTITUTE DATASET

########## Load KI data ########## 
data(ctd_allKI)
mean_expr_data <- ctd_allKI[[2]]$mean_exp
specificity <- ctd_allKI[[2]]$specificity

########## Calculate specificity deciles (least-to-most specific expression) ##########
ctd_allKI = prepare.quantile.groups(ctd_allKI, specificity_species = 'mouse', numberOfBins = 10)
ctd_allKI[[1]] <- NULL

########## Compose gene sets for specificity deciles of cell types per chromosome and store to file ##########
for (cell_type in colnames(ctd_allKI[[1]]$specificity)){
    # Genes not expressed in a cell type
    cell_type_non_expressed <- get_specificity_decile_genes(ctd_allKI, cell_type, 0)
    ensembl_cell_type_non_expressed <- translate_mouse_symbol_to_ensembl(cell_type_non_expressed)
    partition_genes_by_chromosome(genes = ensembl_cell_type_non_expressed,
                                  dataset = 'KI',
                                  cell_type = cell_type,
                                  decile = 'N',
                                  out_dir = geneset_out_path)
    # Genes in specificity decile of cell type 
    for (decile_n in 1:10){
        cell_type_decile_expressed <- get_specificity_decile_genes(ctd_allKI, cell_type, decile_n)
        ensembl_cell_type_decile_expressed <- translate_mouse_symbol_to_ensembl(cell_type_decile_expressed)
        partition_genes_by_chromosome(genes = ensembl_cell_type_decile_expressed,
                                      dataset = 'KI',
                                      cell_type = cell_type,
                                      decile = toString(decile_n),
                                      out_dir = geneset_out_path)
    }
}
