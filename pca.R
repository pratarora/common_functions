PCA <- function(mat,color_pca="",shape_pca= "", label_pca= "",save_plot= "no", name_of_plot= "PCA", comp1=1, comp2=2){
  #Get the differential expressed values from the comparrison interested,
  #extract the normalized values from the assay of vsd and plot them.
  #Giving condition and group from your design table
  
  #1. Extract the counts.
  dt <- mat
  
  #2. Perform pca
  pca_dt <- prcomp(t(dt))
  cat("PCA running...\n")
  # Sys.sleep(0.2)
  
  #3. Extract percentVar data.
  percentVar_dt <- pca_dt$sdev^2/sum(pca_dt$sdev^2)
  cat("Percents calculated...\n")
  # Sys.sleep(0.2)
  
  #4. Create the new dataframe to plot.
  dt_f <- data.frame(PC1=pca_dt$x[,comp1],
                     PC2=pca_dt$x[,comp2],
                     color_pca=color_pca,
                     shape_pca=shape_pca,
                     label_pca= label_pca)
  cat("Data frame built...\n")
  # Sys.sleep(0.2)
  
  #5. Plot it 
  cat("Plotting...\n")
  # Sys.sleep(0.2)
  print(save_plot)
  require(ggplot2)
  require(ggrepel)
  if (save_plot== "no") {
    pca_p <- ggplot(data = dt_f, aes_string(x = paste0("PC1"),
                                          y = paste0("PC2"),
                                          color = "color_pca", 
                                          shape= "shape_pca", label="label_pca")) +
            geom_point(size = 3) +
            geom_text_repel(size= 3, max.overlaps = 50, 
                            box.padding   = 1.5,point.padding = 0.5,force = 50)+
            xlab(paste0("PC", comp1,": ", 
                        round(percentVar_dt[comp1] * 100), "% variance")) +
            ylab(paste0("PC",comp2,": ",
                        round(percentVar_dt[comp2] * 100), "% variance")) +
            # coord_fixed()+ 
      NULL
  }
  if (save_plot== "yes"){
  png(filename =paste0(name_of_plot,".png"),res = 300,width = 2560,height = 1440)
  cat("Saving plot as: ",paste0(name_of_plot,"...\n"))
  pca_p <- ggplot(data = dt_f, aes_string(x = paste0("PC",comp1),
                                          y = paste0("PC",comp2),
                                          color = "color_pca", 
                                          shape= "shape_pca", label="label_pca")) +
            geom_text_repel(size= 3, max.overlaps = 50, 
                            box.padding   = 1.5,
                            point.padding = 0.5,force = 50)+
            geom_point(size = 3) +
            xlab(paste0("PC", comp1,": ", round(percentVar_dt[comp1] * 100), "% variance")) +
            ylab(paste0("PC",comp2,": ", round(percentVar_dt[comp2] * 100), "% variance")) +
            # coord_fixed()+ 
    NULL
  print(pca_p)
  dev.off()
   }
  # Sys.sleep(0.2)
  cat("Done")
  print(pca_p)

  #return(pca_p)
}
