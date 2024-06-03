GO_function <- function(df, pval_DEG  = 0.05, pval_enrichment  = 0.05, onto= "MF", prefix="", org="mouse", add_to_graph_title="", num_pathways=15){

gene_list <- df %>% dplyr::filter(p_Value<=pval_DEG) %>% pull(Gene)
message(paste0("Number of DEG: ",length(gene_list)))

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
        
    
    compGO <- enrichGO(gene = gene_list, pvalueCutoff  = pval_enrichment,keyType = "SYMBOL",
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
        compGO_df <- compGO_df %>%  tidyr::separate_rows(geneID, sep = "/", convert = FALSE) %>%
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

            # print(dotplot(compGO, showCategory = 15, title = paste0("GO Pathway Enrichment Analysis \n",full_name, "\nfor DEGs in ", add_to_graph_title), 
            #         font.size = 12))



        pathways <- unique(compGO_df$Description)
        pathways_selected <- pathways[1:num_pathways] #select top 15 pathways
        
        compGO_df_subset <- compGO_df[compGO_df$Description %in% pathways_selected,] #subset to only include top 15 pathways
        compGO_df_subset <- compGO_df_subset %>% arrange((GeneRatio_decimal)) #arrange by GeneRatio_decimal
        
        compGO_df_subset$Description <- factor(compGO_df_subset$Description, levels=unique(compGO_df_subset$Description)) 
        # print(head(compGO_df_subset))

        p <- ggplot2::ggplot(data = (compGO_df_subset), # make ggplot object 
                            aes(x = GeneRatio_decimal, y = Description,  # x axis is GeneRatio_decimal, y axis is Description
                            color = p.adjust, size = Count)) + # color by p.adjust, size of dot by Count
            geom_point(alpha = 0.5) + # add points with alpha 0.5
            scale_color_gradient(low = "red", high ="blue" ,guide=guide_colourbar(reverse = TRUE)) + # set color gradient to blue to red
            # facet_grid(vars(GO_class), scales = "free",  space = "free_y") + # facet by GO_class with free scales and free y axis
            theme_bw() + # set theme to black and white
            labs(title = paste0("GO Pathway Enrichment Analysis \n", full_name, "\nfor DEGs in ", add_to_graph_title), x = "GeneRatio", y = "") + # set title and axis labels
            theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 10), # set theme for plot title and axis text
                axis.title.x = element_text(size = 10), axis.title.y = element_text(size = 10),
                axis.text.x=element_text(colour="black"),
                axis.text.y=element_text(colour="black")) +
            # change color of Counts in legend to black
            guides(size = guide_legend(override.aes = list(color = "black", alpha=1))) 

        print(p)

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

