# Description: Get gene synonyms from Ensembl Biomart
# Example: Genes_id="ENSG00000115414"
Genes_id="ENSG00000115414"

# function to get gene synonyms
gene_synonyms <- function(x){
  require(biomaRt) # load biomaRt package
  hs_mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl") # select human mart
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol", # select entrez gene ID and gene symbol
                                   "entrezgene_id", "description", "external_synonym"), # select entrez gene ID and gene symbol
                   filters = "ensembl_gene_id", # filter by gene symbol
                   values = x, # gene symbols
                   mart = hs_mart) # select human mart
  results # return
}
# get gene synonyms
syn_genes <- gene_synonyms(Genes_id)
# print top 6 rows
head(syn_genes)
