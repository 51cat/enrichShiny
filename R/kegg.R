#' @import tidyverse
kegg <- function(gene.table,
                 pvalueCutoff,
                 pAdjustMethod,
                 minGSSize,
                 maxGSSize,
                 qvalueCutoff) {
  eg <- to.id(gene.table)
  kegg <- clusterProfiler::enrichKEGG(
    gene = eg$ENTREZID,
    organism = "hsa",
    keyType = "kegg",
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = pAdjustMethod,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    qvalueCutoff = qvalueCutoff,
    use_internal_data = FALSE)
  return(kegg)
}
