Genes <- "Cdh1"

symbol_to_entrez_mouse <- function(x){
  require(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  # query biomart
  results <- getBM(attributes = c("entrezgene_id", "mgi_symbol"),
                   filters = "mgi_symbol", values = x,
                   mart = mart)
  results
}
common_mouse_entrez <- symbol_to_entrez_mouse(Genes)
head(common_mouse_entrez)