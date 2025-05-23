# function to convert gene symbol to ensembl id
GO_function <- function(gene_list, pval  = 0.05, onto= "MF", prefix="", org="mouse", gmt_file = NULL){

if(org=="mouse"){
    orgdb <- "org.Mm.eg.db"
    org_reactome <- "mouse"
}else if(org=="human"){
    orgdb <- "org.Hs.eg.db"
    org_reactome <- "human"
}else{
    message("Please enter a valid organism (mouse or human)")
    return(NULL) # Return NULL to stop execution if organism is invalid
}

if(!is.null(gmt_file)){
    # Use custom GMT file for enrichment analysis
    message(paste0("Using custom GMT file: ", gmt_file))
    gmt <- clusterProfiler::read.gmt(gmt_file)
    compGO <- enricher(gene = gene_list,
                         pvalueCutoff  = pval,
                         pAdjustMethod = "BH",
                         TERM2GENE = gmt
                         )
    full_name = "TH2_vs_TH1" # Set a name for custom GMT
    onto_for_filename = "TH2_vs_TH1" # For filename purposes

} else {
    # Use default GO or Reactome databases
    if(onto %in% c("MF", "CC", "BP")){

        compGO <- enrichGO(gene = gene_list, pvalueCutoff  = pval,keyType = "SYMBOL",
                                 pAdjustMethod = "BH",OrgDb = orgdb, ont = onto)
        full_name= switch(onto,
                        MF= "Moleuclar Function",
                        CC= "Cellular Components",
                        BP= "Biological Process"
                                 )
        onto_for_filename = onto


    }else if(onto=="reactome"){
        gene_list_entrez <- bitr(gene_list, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = orgdb)
        if(nrow(gene_list_entrez) == 0){
            message("No ENTREZID found for the input gene list. Reactome analysis skipped.")
            return(NULL) # Return NULL to stop execution if no ENTREZID found
        }
        gene_list_entrez <- gene_list_entrez$ENTREZID
        compGO <- enrichPathway(gene = gene_list_entrez, pvalueCutoff  = pval, organism= org_reactome,      readable = TRUE)
        print("reactome")
        full_name= "Reactome Pathways"
        onto_for_filename = "reactome"

    }else{
        message("Please enter a valid GO term (MF, CC, BP) or 'reactome' or provide a GMT file.")
        return(NULL) # Return NULL to stop execution if onto is invalid and no GMT
    }
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

        compGO_df %>% head  %>% print()


            if(nrow(compGO_df)==0){
            message(paste0("No GO:",onto, " obtained"))
            message(paste0("****************************************************************************************"))
            message(paste0("\n"))

            } else{
            Genes=compGO_df$`geneID`
            go_annotations <- symbol_to_ensmbl_mouse(Genes) # run the function
            go_annotations_tf <- go_annotations[grepl("DNA-binding transcription factor activity, RNA polymerase II-specific", go_annotations$name_1006),] %>% pull(mgi_symbol)
            # go_annotations_tf %>% print()
            compGO_df$Transcription_Factor <- ifelse((compGO_df$`geneID` %in% go_annotations_tf), "Yes", "No")



            compGO_df %>% arrange(desc(GeneRatio_decimal)) %>% distinct(Description) %>% head(20) %>% print()
            compGO_df %>% arrange(desc(BgRatio_decimal)) %>% distinct(Description) %>% head(20) %>% print()
            compGO_df %>% arrange(p.adjust) %>% distinct(Description) %>% head(20) %>% print()
            compGO_df %>% arrange(pvalue) %>% distinct(Description) %>% head(20) %>% print()
            compGO_df %>% arrange(qvalue) %>% distinct(Description) %>% head(20) %>% print()
            compGO_df %>% arrange(desc(Count)) %>% distinct(Description) %>% head(20) %>% print()
            write.csv(compGO_df, paste0(prefix,"_GO_",onto, "_pathways.csv"))


            print(dotplot(compGO, showCategory = 10, title = paste0("Pathway Enrichment Analysis \n",full_name),
                    font.size = 12))
            dev.copy(
            pdf,
            file = paste0(prefix,"_GO_",onto_for_filename, "_pathways.pdf"),
            width = 7,
            height = 7
            )
            dev.off ()




        message(paste0("Pathway analysis GO:",onto, " done"))
        message(paste0("****************************************************************************************"))
        message(paste0("\n"))
        return(compGO)
        }
    }
}
