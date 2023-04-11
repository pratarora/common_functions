
# Basic function to convert zebrafish to human gene names
# This function takes a list of zebrafish gene names and converts them to human gene names
# It uses the biomaRt package to do this
# It requires the user to have the biomaRt package installed
# It requires the user to have an internet connection
# Example: 
zgGenes <- "ENSDARG00000102750"

# Function to convert zebrafish gene names to human gene names
convertDanioGeneList_Human <- function(x){
  require("biomaRt") # load biomaRt package
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl") # use human mart
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl") # use zebrafish mart
  genesV2 = getLDS(attributes = c("ensembl_gene_id", "zfin_id_symbol"), filters = "ensembl_gene_id",  # get zebrafish gene names
                                values = x ,  # use the zebrafish gene names
                                mart = danio, # use the zebrafish mart
                                attributesL = c("hgnc_symbol", "ensembl_gene_id", "description"), # get human gene names
                                martL = human, # use the human mart
                                uniqueRows=T)# remove duplicates
  
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish"# rename columns
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID.1"] <- "EnsmblID_Human"# rename columns
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  return(genesV2)
}
# Run the function
Human_Genes <- convertDanioGeneList_Human(zgGenes)
print(head(Human_Genes))

