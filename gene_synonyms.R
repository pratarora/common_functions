Genes_id="ENSG00000115414"
gene_synonyms <- function(x){
  require(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol",
                                   "entrezgene_id", "description", "external_synonym"),
                   filters = "ensembl_gene_id", values = x,
                   mart = mart)
  results
}
syn_genes <- gene_synonyms(Genes_id)
head(syn_genes)
