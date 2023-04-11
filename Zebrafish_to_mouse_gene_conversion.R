# Basic function to convert zebrafish to human gene names

zgGenes <- "cdh1"

# Function to convert zebrafish gene names to human gene names
convertDanioGeneList_Mouse <- function(x){ 
  require("biomaRt") # load biomaRt package
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") # use human mart
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl") # use zebrafish mart
  genesV2 = getLDS(attributes = c("ensembl_gene_id", "zfin_id_symbol"), 
                  filters = "ensembl_gene_id", # get zebrafish gene names
                  values = x , # use the zebrafish gene names
                  mart = danio,  # use the zebrafish mart
                  attributesL = c("mgi_symbol", "ensembl_gene_id", "description"), # get human gene names
                  martL = mouse, uniqueRows=T) # use the human mart
  
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish" # rename columns
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID.1"] <- "EnsmblID_Mouse" # rename columns
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2) # return the genes
}
# Run the function
Mouse_Genes <- convertDanioGeneList_Mouse(zgGenes)
# print the first 6 genes
print(head(Mouse_Genes))