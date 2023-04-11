# Description: Convert gene symbols to entrez gene IDs

Genes <- "Cdh1"


symbol_to_entrez_mouse <- function(x){
  require(biomaRt) # load biomaRt package
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") # select mouse mart
  # query biomart
  results <- getBM(attributes = c("entrezgene_id", "mgi_symbol"), # select entrez gene ID and gene symbol
                   filters = "mgi_symbol", # filter by gene symbol
                   values = x, # gene symbols
                   mart = mart) # select mouse mart
  results # return
}
# convert gene symbols to entrez gene IDs
common_mouse_entrez <- symbol_to_entrez_mouse(Genes)
# print top 6 rows
head(common_mouse_entrez)