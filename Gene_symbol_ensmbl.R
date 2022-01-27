




Genes = "CDH1"

symbol_to_ensmbl_human <- function(x){
  require(biomaRt)
  ensembl = useMart(biomart="ensembl", dataset="hsapiens_gene_ensembl", host="useast.ensembl.org")
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                 filters = "hgnc_symbol", values = x,
                 mart = ensembl)
  results
                                    }
symbol_to_ensmbl_human(Genes)

?useMart


Genes="Cdh1"

symbol_to_ensmbl_danio <- function(x){
  require(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = "drerio_gene_ensembl")
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "zfin_id_symbol"),
                                  filters = "zfin_id_symbol", values = x,
                                  mart = mart)
  results
                                      }
symbol_to_ensmbl_danio(Genes)

Genes="cdh1"

symbol_to_ensmbl_mouse <- function(x){
  require(biomaRt)
  mart <- useMart(biomart = "ensembl", dataset = "mmusculus_gene_ensembl")
  # query biomart
  results <- getBM(attributes = c("ensembl_gene_id", "mgi_symbol"),
                                  filters = "mgi_symbol", values = x,
                                  mart = mart)
  results
                                    }
symbol_to_ensmbl_mouse(Genes)

