# code from Stephanie
# overlapSets ------------------------------------------------------------
#' @name overlap_sets
#' @description performs pairwise hypergeometric tests between multiple sets 
#' of genes i.e. cluster marker genes from two different single cell
#' different single cell data sets
#' @param table1 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 1      
#' @param table2 a dataframe with two columns, 
#' ie gene and cluster assignment, from dataset 2
#' format example
#' --------------  
#' Gene  Cluster
#' Gene1 cluster1
#' Gene2 cluster1
#' Gene3 cluster2
#' Gene4 cluster3
#' @param background the number of genes in the gene universe 
#' ie. all expressed genes
#' @return a list with two elements: 
#' [1]intersecting.genes: a tidy table where each row shows two groups 
#' (one from table1 and one from table2) and an intersecting gene or "none"
#' [2] pval.scores: a tidy table where each row is a pair of groups 
#' (one from table1 and one from table2) with the number of genes in each group,
#' the nunber of intersecting genes, the uncorrected pval, and the enrichment score

overlapSets <- function(table1, table2, background){
  require(dplyr)
  
  # make table2 into a list divided by group
  colnames(table1) = c("Gene","Group")
  tab1_groups = unique(table1$Group)
  tab1_list = list()
  
  for(group in tab1_groups){
    tab1_list[[group]] = table1 %>%
      filter(Group == group)
  }
  
  # make table2 into a list divided by group
  colnames(table2) = c("Gene","Group")
  tab2_groups = unique(table2$Group)
  tab2_list = list()
  
  for(group in tab2_groups){
    tab2_list[[group]] = table2 %>%
      filter(Group == group)
  }
  
  # make a list of dfs with the names of genes shared between
  # every group in table1 and every group in table2
  INT = list()
  # make a list of dfs with the number of genes shared between
  # every group in table1 and every group in table2
  nINT = list()
  
  for(group in tab1_groups){
    genes = tab1_list[[group]]$Gene
    INT_list = lapply(tab2_list, 
                      function(x){
                        y = intersect(x$Gene, genes);
                        if(length(y) > 0){
                          z = data.frame(SharedGene = y)}
                        else{z = data.frame(SharedGene = "none")};
                        z
                      })
    
    for(i in 1:length(INT_list)){
      INT_list[[i]]$Group2 = names(INT_list)[i]
    }
    
    INT[[group]] = do.call(rbind,INT_list)
    INT[[group]]$Group1 = group
    
    INT[[group]] = INT[[group]] %>%
      select(Group1, Group2, SharedGene)
    
    nINT_list = lapply(tab2_list, 
                       function(x){
                         y = intersect(x$Gene, genes);
                         z = length(y);
                         z
                       })
    
    nINT_vec = unlist(nINT_list)
    nGroup2 = unlist(lapply(tab2_list, nrow))
    
    nINT[[group]] = data.frame(Group1 = group,
                               nGroup1 = length(genes),
                               Group2 = names(nINT_vec),
                               nGroup2  = nGroup2,
                               nSharedGenes = nINT_vec)
    
    
  }
  
  INTdf = do.call(rbind, INT)
  nINTdf = do.call(rbind, nINT)
  
  nINTdf = nINTdf %>%
    mutate(Pval = phyper(nSharedGenes-1,
                         nGroup2,
                         background-nGroup2,
                         nGroup1,
                         lower.tail= FALSE)
    )
  
  nINTdf = nINTdf %>%
    mutate(Enrichment = 
             log2(
               (nSharedGenes/nGroup2)/
                 (nGroup1/background))
    )
  
  return(list(intersecting.genes = INTdf, 
              pval.scores = nINTdf))
}  