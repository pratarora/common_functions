zgGenes <- rownames(SCT_normalizeddata)
# Basic function to convert zebrafish to human gene names


convertDanioGeneList_Human <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl")
  genesV2 = getLDS(attributes = c("ensembl_gene_id", "zfin_id_symbol"), filters = "ensembl_gene_id", 
                                values = x , mart = danio, attributesL = c("hgnc_symbol", "ensembl_gene_id", "description"), 
                                martL = human, uniqueRows=T)
  
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish"
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID.1"] <- "EnsmblID_Human"
  
  # Print the first 6 genes found to the screen
  return(genesV2)
}

Human_Genes <- convertDanioGeneList_Human(zgGenes)
