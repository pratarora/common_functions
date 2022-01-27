zgGenes <- "cdh1"
# Basic function to convert zebrafish to human gene names


convertDanioGeneList_Mouse <- function(x){
  require("biomaRt")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl")
  genesV2 = getLDS(attributes = c("ensembl_gene_id", "zfin_id_symbol"), filters = "ensembl_gene_id", 
                                values = x , mart = danio, attributesL = c("mgi_symbol", "ensembl_gene_id", "description"), 
                                martL = mouse, uniqueRows=T)
  
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish"
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID.1"] <- "EnsmblID_Mouse"
  
  # Print the first 6 genes found to the screen
  return(genesV2)
}

Mouse_Genes <- convertDanioGeneList_Mouse(zgGenes)
Mouse_Genes
