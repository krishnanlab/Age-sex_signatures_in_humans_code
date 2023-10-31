# findEnrichedGOBP --------------------------------------------------------
#' @name findEnrichedGOBP
#' @description uses topGO to find enriched GOBP terms for a list of genes.
#' ------------ Pvalues from weight01 fisher test
#' @param goi character vector of genes of interest
#' @param background character vector of all background genes
#' @param org.db the org.db for the organism of interest from bioconductor. 
#' ------------- Must already be installed. (human = org.Hs.eg.db)
#' @param id_type gene identifier to use. 
#' ------------- c("entrez","genbank","alias","ensembl","symbol","genename")
#' @param min_size smallest number of genes for a GO term to be used
#' @param max_size largest number of genes for a GO term to be used
#' @return a data frame with the pvalue, bh FDR, and overlapping genes 
#' ------- between the genes of interest and GOBPs. 

findEnrichedGOBP = function(goi, 
                            background, 
                            org.db, 
                            id_type,
                            min_size,
                            max_size){
  
  require(paste(org.db), character.only = TRUE)
  require(tidyverse)
  require(topGO)
  
  submit = factor(ifelse(background %in% goi, 1, 0))
  names(submit) = background
  
  GO = new("topGOdata",
           ontology = "BP",
           allGenes = submit,
           nodeSize = min_size,
           annot = annFUN.org,
           mapping = org.db,
           ID = id_type)
  
  nTerms = length(usedGO(GO))
  
  Fisher <- runTest(GO,
                    algorithm = "weight01", 
                    statistic = "fisher")
  
  results = GenTable(GO, 
                     pval = Fisher, 
                     topNodes = nTerms,
                     numChar=1000)
  
  results = results %>%
    filter(Annotated <=  max_size)
  
  results$FDR = p.adjust(results$pval, method = "BH")
  
  ann.genes = genesInTerm(GO, results$GO.ID)
  symbol = names(submit)[submit == 1]
  
  results$OverlappingGenes = sapply(results$GO.ID, 
                                    function(x){
                                      sig = ann.genes[[x]][ann.genes[[x]] %in% symbol]
                                      sig = paste(sig, collapse = ", ")
                                      return(sig)}
  )
  return(results)
}