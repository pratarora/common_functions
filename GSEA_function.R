GSEA_function_bulk <- function(df, pval_deg = 0.05, pval_enrich = 0.2, onto = "MF", prefix = "", gmt = NULL, organism = "mouse", msig_category = "C7") {
    # Filter DEGs based on p-value and arrange by log2FoldChange
    gene_list_df <- df[df$padj <= pval_deg, ]
    gene_list_df <- gene_list_df %>% arrange(desc(log2FoldChange))

    # Use 'gene' column consistently, handle NA and duplicates in gene symbols
    if (!"gene" %in% colnames(gene_list_df)) {
        if ("Symbol" %in% colnames(gene_list_df)) {
            colnames(gene_list_df)[colnames(gene_list_df) == "Symbol"] <- "gene"
        } else {
            stop("Error: 'gene' or 'Symbol' column not found in input dataframe.")
        }
    }
    gene_list_df <- gene_list_df[!is.na(gene_list_df$gene), ]
    gene_list_df <- gene_list_df[!duplicated(gene_list_df$gene), ]

    # Create gene list ranked by log2FoldChange
    gene_list <- gene_list_df %>% pull(log2FoldChange)
    names(gene_list) <- gene_list_df %>% pull(gene)
    gene_list <- gene_list[!duplicated(gene_list)]
    print(head(gene_list))

    # Select OrgDb based on organism
    OrgDb <- if (organism == "human") {
        "org.Hs.eg.db"
    } else if (organism == "mouse") {
        "org.Mm.eg.db"
    } else {
        stop("Unsupported organism. Please use 'human' or 'mouse'.")
    }

    # Perform GSEA based on 'onto' parameter
    if (onto %in% c("MF", "CC", "BP")) {
        compGO <- gseGO(gene = gene_list, pvalueCutoff = pval_enrich, keyType = "SYMBOL",
                        pAdjustMethod = "BH", OrgDb = OrgDb, ont = onto)
    } else if (onto == "reactome") {
        if (!requireNamespace("ReactomePA", quietly = TRUE)) {
            stop("Package 'ReactomePA' is needed for reactome analysis. Please install it.")
        }
        gene_entrez <- bitr(names(gene_list), fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb)
        gene_list_entrez <- gene_list[names(gene_list) %in% gene_entrez$SYMBOL] # Use gene_list_entrez to avoid modifying original gene_list names before gseGO/GSEA if other ontos are run later
        names(gene_list_entrez) <- gene_entrez$ENTREZID[match(names(gene_list_entrez), gene_entrez$SYMBOL)]
        compGO <- gsePathway(geneList = gene_list_entrez, pvalueCutoff = pval_enrich,
                             pAdjustMethod = "BH", organism = organism, verbose = FALSE)
    } else if (onto == "msigdb") {
        if (!requireNamespace("msigdbr", quietly = TRUE)) {
            stop("Package 'msigdbr' is needed for MSigDB analysis. Please install it.")
        }
        msig_org <- if (organism == "mouse") "Mus musculus" else "Homo sapiens"
        msigdb_set <- msigdbr(species = msig_org, category = msig_category)
        gmt <- msigdb_set %>% dplyr::select(gs_name, gene_symbol)
        compGO <- GSEA(geneList = gene_list, TERM2GENE = gmt, pvalueCutoff = pval_enrich, pAdjustMethod = "BH")
    } else if (onto == "gmt") {
        if (is.null(gmt)) {
            stop("Error: GMT file must be provided when onto='gmt'.")
        }
        compGO <- GSEA(geneList = gene_list, TERM2GENE = gmt, pvalueCutoff = pval_enrich, pAdjustMethod = "BH")
    }
    else {
        stop("Unsupported ontology type. Use MF, CC, BP, reactome, msigdb, or gmt.")
    }

    # Handle cases with no significant pathways
    if (is.null(compGO) || nrow(compGO@result) == 0) {
        message(paste0("No pathways obtained for: ", onto))
        message(paste0("****************************************************************************************"))
        message(paste0("\n"))
        return(NULL) # Return NULL explicitly when no results
    } else {
        # Add 'sign' column for activation/suppression status
        compGO@result$.sign <- ifelse(compGO@result$NES > 0, "Activated", "Suppressed")
        compGO@result$.sign <- factor(compGO@result$.sign, levels = c("Activated", "Suppressed"))
        print(table(compGO@result$.sign))

        # Prepare dataframe for output
        compGO_df <- as.data.frame(compGO)
        compGO_df <- compGO_df %>% tidyr::separate_rows(core_enrichment, sep = "/", convert = FALSE) %>% arrange((p.adjust))

        # Write results to CSV
        write.csv(compGO_df, paste0(prefix, "_GSEA_", onto, "_pathways.csv"))

        # Determine full ontology name for plot title
        full_name <- switch(onto,
            MF = "Molecular Function",
            CC = "Cellular Components",
            BP = "Biological Pathways",
            reactome = "Reactome Pathways",
            msigdb = paste0("MSigDB ", msig_category, " Signatures"),
            gmt = "Custom GMT Pathways"
        )

        # Generate and save dotplot
        p <- dotplot(compGO, showCategory = 15, title = paste0("GSEA Analysis\n", full_name),
                     split = ".sign", font.size = 12) + facet_grid(.~.sign)
        print(p)
        dev.copy(pdf, file = paste0(prefix, "_GSEA_", onto, "_pathways.pdf"), width = 15, height = 10)
        dev.off()

        # Messages for completion
        message(paste0("GSEA for ", full_name, " done"))
        message(paste0("****************************************************************************************"))
        message(paste0("\n"))

        return(compGO) # Return compGO object
    }
}
