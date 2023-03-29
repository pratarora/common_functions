# convert gene from zebrafish to mouse using biomart
# https://www.bioconductor.org/packages/release/bioc/vignettes/biomart/inst/doc/biomart.html
zgGenes <- "ENSDARG00000102750"
# Basic function to convert zebrafish to human gene names

# Function to convert zebrafish gene names to mouse gene names
convertDanioGeneList_Mouse <- function(x){
  require("biomaRt")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl") # mouse mart
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl") # zebrafish mart
  # Get the gene names for the genes in the list
  genesV2 = getLDS(attributes = c("ensembl_gene_id", "zfin_id_symbol"), # zebrafish gene name
  filters = "ensembl_gene_id", # zebrafish gene id
    values = x, # zebrafish gene id
    mart = danio, # zebrafish mart
    attributesL = c("mgi_symbol", "ensembl_gene_id", "description"), # mouse gene name
    martL = mouse, # mouse mart
    uniqueRows = T) # only return unique rows

  # Rename the columns
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish"
  colnames(genesV2)[colnames(genesV2)== "Gene.stable.ID.1"] <- "EnsmblID_Mouse"
  
  # Print the first 6 genes found to the screen
  print(head(genesV2))
  # Return the list of genes
  return(genesV2)
}

# Run the function
Mouse_Genes <- convertDanioGeneList_Mouse(zgGenes)

# Print the first 6 genes found to the screen
Mouse_Genes
