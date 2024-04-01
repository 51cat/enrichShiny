# analysis ----------------------------------------------------------------

enrich.Select <- function(df, enrich.method, pvalueCutoff, pAdjustMethod, minGSSize, maxGSSize, qvalueCutoff, ont) {
  if (enrich.method == "KEGG"){
    return(kegg(
      gene.table = df,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      qvalueCutoff = qvalueCutoff
    ))
  }else if (enrich.method == "GO") {
    return(go(
      gene.table = df,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      qvalueCutoff = qvalueCutoff,
      ont = ont
    ))
  }else if (enrich.method == "gseKEGG"){
    return(kegg.gsea(
      gene.table = df,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize
    ))
  }else if (enrich.method == "gseGO"){
    return(go.gsea(
      gene.table = df,
      pvalueCutoff = pvalueCutoff,
      pAdjustMethod = pAdjustMethod,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      ont = ont
    ))
  }
}

enrich.p <- function(enrichobj, show = 10, type = "barplot") {
  if( type == "barplot") {
    return(barplot(enrichobj,showCategory = show))
  }else{
    if ("ONTOLOGY" %in% colnames(enrichobj@result)) {
      return(clusterProfiler::dotplot(enrichobj,showCategory = show,split = "ONTOLOGY") + facet_wrap(.~ONTOLOGY,scales = "free_y"))
    } else {
      return(clusterProfiler::dotplot(enrichobj,showCategory = show))
    }
  }
}

enrich.df <- function(enrichobj) {
  return(enrichobj@result)
}
