GSEA_function <- function(df,cell_type, pval_deg  = 0.05,pval_enrich=0.2, onto= "MF", prefix=""){
    
    gene_list_df <- df[df$cell_type==cell_type & df$p_val<=pval_deg,]
    gene_list_df <- gene_list_df %>% arrange(desc(avg_logFC))
    gene_list_df <- gene_list_df[!is.na(gene_list_df$gene),]
    gene_list_df <- gene_list_df[!duplicated(gene_list_df$gene),]
    gene_list <- gene_list_df %>% pull(avg_logFC)
    names(gene_list) <- gene_list_df %>% pull(gene)

    gene_list <- gene_list[!duplicated(gene_list)]
    print(head(gene_list)) 



    compGO <- gseGO(gene = gene_list, pvalueCutoff  = pval_enrich,,keyType = "SYMBOL",
                             pAdjustMethod = "BH",OrgDb = "org.Mm.eg.db", ont = onto)
    # print(compGO)
    if(is.null(compGO)|nrow(compGO@result)==0){
    message(paste0("No GP:",onto, " obtained for Cell Type : ", cell_type, " done"))
    message(paste0("****************************************************************************************"))
    message(paste0("\n"))
   
    }else {

        compGO_df <- as.data.frame(compGO)
          compGO_df <- compGO_df %>% tidyr::separate_rows(core_enrichment, sep = "/", convert = FALSE) %>%  arrange((p.adjust))

            if(nrow(compGO_df)==0){
            message(paste0("No GP:",onto, " obtained for Cell Type : ", cell_type, " done"))
            message(paste0("****************************************************************************************"))
            message(paste0("\n"))

            } else{

            write.csv(compGO_df, paste0(prefix,cell_type,"_GO_",onto, "_pathways.csv"))


            full_name= switch(onto,
                    MF= "Moleuclar Function",
                    CC= "Cellular Components",
                    BP= "Biological Pathways"
                             )                                           

            print(dotplot(compGO, showCategory = 15, title = paste0("DEG GO Pathway Enrichment Analysis \n",full_name, " for ",cell_type, " Cells"), 
                    font.size = 12) + facet_grid(.~.sign))
            dev.copy(
            svg,
            file = paste0(prefix,cell_type,"_GO_",onto, "_pathways.svg"),
            width = 10,
            height = 8
            )
            dev.off ()




        message(paste0("DEG Pathway analysis GO:",onto, " for Cell Type : ", cell_type, " done"))
        message(paste0("****************************************************************************************"))
        message(paste0("\n"))

        }
    }
}



dir.create("GSEA_subcluster", showWarnings = FALSE)
dir.create("GSEA_subcluster/Molecular_function", showWarnings = FALSE)
dir.create("GSEA_subcluster/Cellular_component", showWarnings = FALSE)
dir.create("GSEA_subcluster/Biological_pathways", showWarnings = FALSE)

cell_type <- DEG_Libra_edger_lrt$cell_type %>% unique %>% sort  %>% as.character
cell_type
lapply(X = cell_type,  function(x) {GSEA_function(df = DEG_Libra_edger_lrt, x, onto = "MF", prefix="GSEA_subcluster/Molecular_function/")})
