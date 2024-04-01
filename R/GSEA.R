#' @import tidyverse

kegg.gsea <- function(gene.table,
                      pvalueCutoff,
                      pAdjustMethod,
                      minGSSize,
                      maxGSSize) {
  gene.table <- add.id(gene.table)
  genes.gsea <- gene.table %>% tidyr::drop_na() %>%
    with(setNames(.$log2FoldChange, .$ENTREZID))
  genes.gsea <- sort(genes.gsea, decreasing = T)

  kegg <- clusterProfiler::gseKEGG(
    geneList = genes.gsea,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    use_internal_data = FALSE
  )
  return(kegg)

}

go.gsea <- function(gene.table,
                    pvalueCutoff,
                    pAdjustMethod,
                    minGSSize,
                    maxGSSize,
                    ont) {
  gene.table <- add.id(gene.table)
  genes.gsea <- gene.table %>% tidyr::drop_na() %>%
    with(setNames(.$log2FoldChange, .$ENTREZID))
  genes.gsea <- sort(genes.gsea, decreasing = T)

  go <- clusterProfiler::gseGO(
    geneList = genes.gsea,
    OrgDb = org.Hs.eg.db,
    ont= ont,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff = pvalueCutoff,
    keyType = 'ENTREZID')
  return(go)
}



plot.gsea <- function(enrichobj, pathway_name) {
  return(enrichplot::gseaplot2(enrichobj,pathway_name, pvalue_table = TRUE))
}
