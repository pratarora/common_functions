



# This function converts gene symbol to ensembl id for human, mouse and zebrafish
Genes = "CDH1"

# Human
symbol_to_ensmbl_human <- function(x){
  require(biomaRt) # load biomaRt package
  ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org") # use the ensembl mart
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"), # get the ensembl gene id and the hgnc symbol
                 filters = "hgnc_symbol", # filter by the hgnc symbol
                 values = x, # the gene symbol
                 mart = ensembl)# the ensembl mart
  results # return the results
                                    }
symbol_to_ensmbl_human(Genes) # run the function

# Zebrafish
Genes="cdh1"


symbol_to_ensmbl_danio <- function(x){ 
  require(biomaRt) # load biomaRt package
  zf_mart <- useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl") # use the ensembl mart
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "zfin_id_symbol"), # get the ensembl gene id and the zfin symbol
                                  filters = "zfin_id_symbol", # filter by the zfin symbol
                                  values = x, # the gene symbol
                                  mart = zf_mart) # the ensembl mart
  results # return the results
                                      }
symbol_to_ensmbl_danio(Genes) # run the function

# Mouse
Genes="Cdh1"

symbol_to_ensmbl_mouse <- function(x){
  require(biomaRt) # load biomaRt package
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl") # use the ensembl mart
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"), # get the ensembl gene id and the mgi symbol
                                  filters = "mgi_symbol", # filter by the mgi symbol
                                  values = x, # the gene symbol
                                  mart = mart) # the ensembl mart
  results # return the results
                                    }
symbol_to_ensmbl_mouse(Genes) # run the function

