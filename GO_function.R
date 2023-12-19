GO_function <- function(gene_list, pval  = 0.05, onto= "MF", prefix="", org="mouse"){

if(org=="mouse"){
    orgdb <- "org.Mm.eg.db"
    org_reactome <- "mouse"
}else if(org=="human"){
    orgdb <- "org.Hs.eg.db"
    org_reactome <- "human"
}else{
    message("Please enter a valid organism (mouse or human)")
}

    if(onto %in% c("MF", "CC", "BP")){
        
    
    compGO <- enrichGO(gene = gene_list, pvalueCutoff  = pval,keyType = "SYMBOL",
                             pAdjustMethod = "BH",OrgDb = orgdb, ont = onto)
    }else if(onto=="reactome"){
        gene_list <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
        gene_list <- gene_list$ENTREZID
        compGO <- enrichPathway(gene = gene_list, pvalueCutoff  = 0.05, organism= org_reactome,      readable = TRUE)
    
    }else{
        message("Please enter a valid GO term")
    }

      
    
    if(is.null(compGO)){
    message(paste0("No GO:",onto, " obtained"))
    message(paste0("****************************************************************************************"))
    message(paste0("\n"))
   
    }else {

        compGO_df <- as.data.frame(compGO)
        compGO_df$GeneRatio_decimal <- compGO_df$GeneRatio
        compGO_df$GeneRatio_decimal <- sapply(compGO_df$GeneRatio_decimal, 
                                                      function(x) (eval(parse(text = as.character(x)))))
        compGO_df$BgRatio_decimal <- compGO_df$BgRatio
        compGO_df$BgRatio_decimal <- sapply(compGO_df$BgRatio_decimal, 
                                                    function(x) (eval(parse(text = as.character(x)))))
        compGO_df <- compGO_df %>% tidyr::separate_rows(geneID, sep = "/", convert = FALSE) %>%
          arrange(desc(GeneRatio_decimal))
        compGO_df %>% head 

            if(nrow(compGO_df)==0){
            message(paste0("No GO:",onto, " obtained"))
            message(paste0("****************************************************************************************"))
            message(paste0("\n"))

            } else{

            write.csv(compGO_df, paste0(prefix,"_GO_",onto, "_pathways.csv"))


            full_name= switch(onto,
                    MF= "Moleuclar Function",
                    CC= "Cellular Components",
                    BP= "Biological Pathways",
                    reactome= "Reactome Pathways"
                             )                                           

            print(dotplot(compGO, showCategory = 15, title = paste0("GO Pathway Enrichment Analysis \n",full_name), 
                    font.size = 12))
            dev.copy(
            pdf,
            file = paste0(prefix,"_GO_",onto, "_pathways.pdf"),
            width = 10,
            height = 12
            )
            dev.off ()




        message(paste0("Pathway analysis GO:",onto, " done"))
        message(paste0("****************************************************************************************"))
        message(paste0("\n"))

        }
    }
}
