# make bigger graph
options(repr.plot.width=6, repr.plot.height=8)
ha = HeatmapAnnotation(
    Genotype = rep(c("WT","Cox7a1 KO" ), each=3), 
    col = list(Genotype = c(
                # "WT"=paletteer::paletteer_d("calecopal::arbutus", n=2)[1],
                            "WT"="#226741",
                            # "Cox7a1_KO"=paletteer::paletteer_d("calecopal::arbutus", n=2)[2])
                            "Cox7a1 KO"="#f69431")
                
               ),
               # change font size of the annotation
                annotation_name_gp = grid::gpar(fontsize = 5),
                annotation_legend_param = list(
        # title = "Expression Levels",
        title_gp = gpar(fontsize = 5),
        labels_gp = gpar(fontsize = 5),
        grid_width	= unit(5,"point"),
        grid_height	= unit(5,"point")

        
))

Heatmap(heatmap_mat_v2_others, name="Expression Levels", 
 cluster_columns = FALSE,  
 cluster_rows = FALSE,  
row_dend_width = unit(0, "cm"), 
# col = viridis::viridis(100, direction = -1),
# col= paletteer::paletteer_c("grDevices::Rocket", n=100, direction=-1),
show_column_names = FALSE,
top_annotation = ha,
row_split = factor(heatmap_list_v2_others_unique$OXPHOS.Complex, levels = unique(heatmap_list_v2_others_unique$OXPHOS.Complex)),
gap=unit(2, "point"),
border= TRUE,
 border_gp = gpar(col = "grey40"),
row_title_gp = gpar(fontsize = 5), # change row split font size
row_title_rot = 0, 
bottom_annotation = ha,
heatmap_legend_param = list(
        title = "Expression Levels",
        title_gp = gpar(fontsize = 5),
        labels_gp = gpar(fontsize = 5),
        grid_width	= unit(5,"point"),
        grid_height	= unit(5,"point")
    ),
row_names_gp = grid::gpar(fontsize = 3))

dev.copy(pdf, file = "proteomics_n_3_muscle/muscle_heatmap_annotated_v5_others.pdf", width = 4, height = 8)
dev.off()
