################################
#                              #
#       umap visualization     #
#                              #
#                              #
################################

umapTheMainData <- function(rawData=matrix_counts,normData=assay(vsd),dataSet,name,group=s2c$group) {
  
  #Plot the UMAP version of the data.
  #Normalized data
  #Raw data
  
  
  #1. Plot the raw counts umap 
  cat("Plotting raw counts UMAP...\n")
  Sys.sleep(0.2)
  #pdf(file = paste0("./03_Graphs/miRNAs/PCA/PCA_DExpressed_mirnas_",name,".pdf"),width = 10,height = 12)
  png(filename =paste0("./03_Graphs/miRNAs/UMAP/Origial_counts_umap_mirnas_",name,".png"),res = 300,width = 2560,height = 1440)
  cat("Saving plot as: ",paste0("./03_Graphs/miRNAs/UMAP/Original_counts_umap_mirnas_",name,"...\n"))
  umap_p <- umap(rawData,
                 labels=as.factor(group),
                 controlscale=TRUE,
                 scale=3)
  print(umap_p)
  rm(umap_p)
  dev.off()
  cat("Plotting raw counts UMAP...\nDone!\n\n")
  Sys.sleep(0.2)
  
  #2. Plot the normalized UMAP
  cat("Plotting normalized counts UMAP...\n")
  Sys.sleep(0.2)
  #png
  png(filename =paste0("./03_Graphs/miRNAs/UMAP/Normalized_counts_vsd_umap_mirnas_",name,".png"),res = 300,width = 2560,height = 1440)
  cat("Saving plot as: ",paste0("./03_Graphs/miRNAs/UMAP/Normalized_counts_vsd_umap_mirnas_",name,"...\n"))
  umap_n <- umap(normData,
                 labels=as.factor(group),
                 controlscale=TRUE,
                 scale=3)
  print(umap_n)
  umap_n
  dev.off()
  rm(umap_n)
  cat("Plotting VSt Normalize counts UMAP...\n\nDone!\n")
  Sys.sleep(0.2)
  
  
}

umapTheGroupData <- function(normData=assay(vsd),dataSet,name,group=s2c$group,condition=s2c$condition) {
  
  #Plot the comparisons data.  
  
  #1. Extract the counts from specific comparisons.
  dt <- normData
  dt_a <-dt[intersect(rownames(dt),dataSet),]
  
  cat("Counts values extracted...\n")
  Sys.sleep(0.2)
  
  #2. Plotting the data of the specific grouup comparisons.
  cat(paste0("Plotting the group:",name,"UMAP...\n"))
  Sys.sleep(0.2)
  #png
  png(filename =paste0("./03_Graphs/miRNAs/UMAP/Umap_mirnas_",name,".png"),res = 300,width = 2560,height = 1440)
  cat("Saving plot as: ",paste0("./03_Graphs/miRNAs/UMAP/Umap_mirnas_",name,"...\n"))
  umap_s <- umap(dt_a,
                 labels=as.factor(group),
                 controlscale=TRUE,
                 scale=3)
  print(umap_s)
  dev.off()
  cat("UMAP\n...\nDone!\n")
  
}

#All onrmalized and Raw
umapTheMainData(name = "all")


#aci vs caci
umapTheGroupData(name = "aci_vs_caci",
                 dataSet = row.names(mylist_res_cleaned_005_log2fc_1$ACI_vs_CACI))
#bci vs cbci
umapTheGroupData(name = "bci_vs_cbci",
                 dataSet = row.names(mylist_res_cleaned_005_log2fc_1$BCI_vs_CBCI))

#bci vs aci
umapTheGroupData(name = "bci_vs_cbci",
                 dataSet = row.names(mylist_res_cleaned_005_log2fc_1$BCI_vs_ACI))

#cbci vs caci
umapTheGroupData(name = "cbci_vs_caci",
                 dataSet = row.names(mylist_res_cleaned_005_log2fc_1$CBCI_vs_CACI))




#Specific genes
png(filename="./03_Graphs/miRNAs/UMAP/Umap_Normalized_mirnas_miR-107b.png",
    res = 300,
    width = 2560,height = 1440)

umap(mydata = assay(vsd),
     labels = scale(as.numeric(assay(vsd)[row.names(assay(vsd)) == "dre-miR-107b",])),
     controlscale = TRUE,scale = 2,legendtitle = "miR-107b")

dev.off()

