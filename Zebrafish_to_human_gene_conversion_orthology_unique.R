zgGenes <- c("ENSDARG00000000966", "ENSDARG00000012314")

# Function to convert zebrafish to human gene names with maximum one-one orthology 

convertDanioGeneList_Mouse <- function(x){
  require("biomaRt")
  danio = useMart("ensembl", dataset = "drerio_gene_ensembl",  host="useast.ensembl.org") # query organism (from this organism)
  
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl",  host="useast.ensembl.org") # target organism (convert genes to this organism)
  
  
  danio_ensembl_id_attributes <-c("ensembl_gene_id", "mmusculus_homolog_perc_id", # attributes to query from danio database (query organism) 
                                  "mmusculus_homolog_goc_score") #orthology score to mouse from zebrfaish side
                                                                
  
  
  #(you can't get both ensembl_id zfin_id and mouse orthology scores together, so make two dfs and merge)
  danio_zfin_id_attributes <-c("ensembl_gene_id", "zfin_id_symbol") # attributes to convert ensembl id to zfin id 
  
  
  mmus_attributes <- c("mgi_symbol", "ensembl_gene_id","entrezgene_id") # attributes that you want from mouse database
  
  #make a converted dataframe with ensemble ids and mouse orthology scores 
  genes_conv_ensembl = getLDS(values = x , filters = "ensembl_gene_id",  # query using ensemble ids
                              attributes = danio_ensembl_id_attributes, # query these attributes from danio mart
                              mart = danio, # use danio mart
                              attributesL = mmus_attributes, # convert to attributes of mouse mart 
                              martL = mouse, uniqueRows=T) # use mouse mart
  
  # make dataframe for zf ensembl id and zfin
  genes_conv_zfin = getBM(attributes = c("ensembl_gene_id", "zfin_id_symbol"),
                          filters = "ensembl_gene_id", values = x,
                          mart = danio)
  
  # merge conversion df with zfin id df
  gene_merge = merge(genes_conv_ensembl, genes_conv_zfin, by.x= "Gene.stable.ID",by.y="ensembl_gene_id", all=TRUE)
  
  # rename column names to readable format
  colnames(gene_merge)[colnames(gene_merge)== "Gene.stable.ID"] <- "EnsmblID_Zebrafish"
  colnames(gene_merge)[colnames(gene_merge)== "X.id..target.Mouse.gene.identical.to.query.gene"] <- "Mouse_homology_percentage"
  colnames(gene_merge)[colnames(gene_merge)== "Mouse.Gene.order.conservation.score"] <- "Mouse_Gene_order_conservation_score"
  colnames(gene_merge)[colnames(gene_merge)== "NCBI.gene..formerly.Entrezgene..ID"] <- "EntrezID_Mouse"
  colnames(gene_merge)[colnames(gene_merge)== "Gene.stable.ID.1"] <- "EnsmblID_Mouse"
  
  #select the top ortholog with max homology percentage 
  
  gene_merge <- gene_merge %>% group_by(EnsmblID_Zebrafish) %>% # group by ensembl id (since ensemble ids are the one used for query)
    arrange(desc(Mouse_homology_percentage),desc(Mouse_Gene_order_conservation_score),.by_group = TRUE) %>% #arrange each of them by homology score then by GOC score
    mutate(num_orthologues=n()) %>% #write in a column how many orthologs are possible
    dplyr::slice(n= 1) %>%  # if two orthologues have same homology percentage and goc then choose the top one
    arrange(desc(num_orthologues)) # arrange by the genes which have maximum orthogues
  
  
  # Print the first 6 genes found to the screen
  return(gene_merge)
}


df <- convertDanioGeneList_Mouse(zgGenes)
head(df)