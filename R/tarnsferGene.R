#' @import org.Hs.eg.db
#' @import tidyverse
to.id <- function(gene.table){
  eg <- clusterProfiler::bitr(
    geneID = gene.table$gene,
    fromType = "SYMBOL",
    toType=c("ENTREZID","ENSEMBL"),
    OrgDb="org.Hs.eg.db")
  return(eg)
}

add.id <- function(gene.table){
  eg <- clusterProfiler::bitr(
    geneID = gene.table$gene,
    fromType = "SYMBOL",
    toType=c("ENTREZID"),
    OrgDb="org.Hs.eg.db")
  colnames(eg) <- c("gene", "ENTREZID")
  gene.table <- gene.table %>%
    dplyr::left_join(eg,by = "gene")
  return(gene.table)
}
