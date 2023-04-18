# Basic function to convert zebrafish to mouse gene names

zgGenes <- "cdh1"

# This function uses biomaRt to convert zebrafish gene names to mouse gene names. 
# This is done by using the zebrafish and mouse ensembl mart datasets. 
# The function takes a vector of zebrafish gene names as input and returns 
# a data frame of zebrafish and mouse gene names.

convertDanioGeneList_Mouse <- function(x){ 
  require("biomaRt") # load biomaRt package
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") # use mouse mart
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl") # use zebrafish mart
  
  genesV2 = getLDS(attributes = c("ensembl_gene_id", "zfin_id_symbol"), 
                   filters = "ensembl_gene_id", # get zebrafish gene names
                   values = x , # use the zebrafish gene names
                   mart = danio,  # use the zebrafish mart
                   attributesL = c("mgi_symbol", "ensembl_gene_id", "description"), # get mouse gene names
                   martL = mouse, uniqueRows=T) # use the mouse mart
  
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish" # rename columns
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID.1"] <- "EnsmblID_Mouse" # rename columns
  
  # Check if the gene is not found
  if (length(genesV2) == 0) {
    print("No gene found for this input")
  } else { 
    return(genesV2) # return the genes
  }
}

# Run the function
Mouse_Genes <- convertDanioGeneList_Mouse(zgGenes)
# print the first 6 genes
print(head(Mouse_Genes))