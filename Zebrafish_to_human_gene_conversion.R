zgGenes <- "ENSDARG00000102750"
# Basic function to convert zebrafish to human gene names

# Function to convert zebrafish to human gene names
convertDanioGeneList_Human <- function(x) {
    # Load the biomaRt package
    require("biomaRt")
    # Load the ensembl mart for human and zebrafish
    human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl") # Human
    danio <- useMart("ensembl", dataset = "drerio_gene_ensembl") # Zebrafish
    # Get the human gene names for the zebrafish genes
    genesV2 <- getLDS(
        filters = "ensembl_gene_id", # Filter by ensembl gene id
        attributes = c("ensembl_gene_id", "zfin_id_symbol"), # Get the ensembl gene id and zebrafish gene name
        attributesL = c("hgnc_symbol", "ensembl_gene_id", "description"), # Get the human gene name, ensembl gene id, and description
        values = x, # The zebrafish gene names
        mart = danio, # The zebrafish mart
        martL = human, # The human mart
        uniqueRows = T
    ) # Only return unique rows
    # Rename the columns
    colnames(genesV2)[colnames(genesV2) == "Gene.stable.ID"] <- "EnsmblID_Zebrafish"
    colnames(genesV2)[colnames(genesV2) == "Gene.stable.ID.1"] <- "EnsmblID_Human"

    # Print the first 6 genes found to the screen
    print(head(genesV2))
    # Return the data frame
    return(genesV2)
}
# Use the function to convert the zebrafish gene names to human gene names
Human_Genes <- convertDanioGeneList_Human(zgGenes)
# print the first 6 genes to the screen
head(Human_Genes)
