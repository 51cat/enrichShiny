#' @import org.Hs.eg.db
#' @import tidyverse
go <- function(gene.table,
               pvalueCutoff,
               pAdjustMethod,
               minGSSize,
               maxGSSize,
               qvalueCutoff,
               ont){
  eg <- to.id(gene.table)
  go <- clusterProfiler::enrichGO(
    eg$ENTREZID,
    OrgDb = org.Hs.eg.db,
    ont=ont,
    minGSSize = minGSSize,
    maxGSSize = maxGSSize,
    pAdjustMethod = pAdjustMethod,
    pvalueCutoff = pvalueCutoff,
    qvalueCutoff = qvalueCutoff,
    keyType = 'ENTREZID')
  return(go)
}
