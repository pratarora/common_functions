setwd(dir = "C:/Users/marius/Documents/Phd/Projects/smallRNAs/dockerDoor/")
# --- Differentially miRNA expressed and get ready for annotation and downstream analysis.
files <- list.files(all.files = TRUE,path = "./02_Data/mirna/countsOutputs",
                    pattern = "*SM_all.tsv",
                    recursive = TRUE, 
                    full.names = TRUE)
#Sort the files by the number of output
files <- files[order(nchar(files))]
files


#Subsittute all the string for a shorter name.
#names_files <- gsub("/media/marius/Samsung_T5/MercaderLab/benne_andres/longRNAseq/sample_\\d*/","",files)
name_files <- gsub("02_Data/mirna/countsOutputs/counts_","",files)
#Subsittute all the string for a shorter name.
name_files <- gsub("_SM_all.tsv","",name_files)
name_files
#names_files
#Set the names
names(files) <- name_files
head(files)

lapply(names(files),function(x) head(x))
#lapply(files[,])

#STAR 2.7.3a
#column 1: gene ID
#column 2: counts for unstranded RNA-seq
#column 3: counts for the 1st read strand aligned with RNA (htseq-count option -s yes)
#column 4: counts for the 2nd read strand aligned with RNA (htseq-count option -s reverse)
#Get the gene names for the table of counts.

# --- read the genes from the fasta file.

comdwcFa <- 'zcat 02_Data/mirna/matureFasta/mature.fa.gz | grep -E -w "Danio rerio" | grep -Eo ".{0,0}dre-.{0,18}" | cut -f 1 -d" " | sort -V | uniq | wc -l'
system(command = comdwcFa)
# --- as this idk how to do it yet in windows use it with ubuntu and the saved file.
comdFa <- 'zcat   ./02_Data/mirna/matureFasta/mature.fa.gz | grep -E -w "Danio rerio" | grep -Eo ".{0,0}dre-.{0,18}" | cut -f 1 -d" " | sort -V | uniq'
genes_dre_fasta <- as.data.frame(system(command = comdFa, intern = TRUE))
#read it out of the table generated with the comman above
genes_dre_fasta <- read.table(file = "02_Data/mirna/matureFasta/miRBase_IDs.tsv")
colnames(genes_dre_fasta) <-c("genes_dre_fasta")
head(genes_dre_fasta)
dim(genes_dre_fasta)

#########################
#                       #
#   Matrix Counts       #
#                       #
#########################
#samples_vec <- c()
#for (i in length(list_counts)){
#  samples_vec[i] <- names(list_counts[i][[2]])
#}



list_counts <- lapply(files,function(fj) read.table(fj,header = 1))
lapply(list_counts,tail)
lapply(list_counts,dim)



#list_counts
samples_vec <- rep(paste0("X",seq(1:20),"_SM_all"))
outputs_vec <- rep(paste0("output_",seq(1:20)),1)
outputs_vec
samples_vec

#
#for (i in length(list_counts)){
#  outputs_[[i]] <- merge(x = genes_dre_fasta,y = paste("list_counts$output_",i), by.x = "genes_dre_fasta", by.y = paste0("X",i,"_SM_all"), all.x = TRUE)
#}



#Change names from files to be easier to merge and acces.
colnamesss <- c("Counts","SM_all")
list_counts <- lapply(list_counts,setNames, colnamesss)
lapply(list_counts,head)

head(genes_dre_fasta)

outputs_counts_fasta <- lapply(names(list_counts), function(x){
  merge(x = genes_dre_fasta, y = list_counts[[x]],
        by.x = "genes_dre_fasta",
        by.y = "SM_all",
        all.x = TRUE)})

#print(paste0("Completed with: ",names(list_counts[[x]][2])))})

names(outputs_counts_fasta) <- name_files
head(outputs_counts_fasta)

#Check IDs are the same for all the lists.
identical(x = outputs_counts_fasta$output_1$genes_dre_fasta,outputs_counts_fasta$output_19$genes_dre_fasta)


#outputs_counts_fasta

#Change NAs to 0.
outputs_counts_fasta <- lapply(outputs_counts_fasta, function(e) {
  e[is.na(e)] <- 0
  e
})
#Or with replace.
#outputs_counts_fasta <- lapply(outputs_counts_fasta,function(e) replace(e, is.na(e),0))


#Create the matrix counts now.

#CHange the names of the columns to recognize the sample outputs later.
#colnamesss <- c("IDs",outputs_vec)
outputs_counts_fasta <- lapply(outputs_counts_fasta,setNames, colnamesss)
lapply(outputs_counts_fasta,head)

#cbind them

#Now select only the columns we want to cbind.
outputs_counts_fasta_selected <- lapply(outputs_counts_fasta,function(x) x[,2])
outputs_counts_fasta_selected
matrix_counts <- do.call(cbind, outputs_counts_fasta_selected)


#Create the sample table for analyzing the conditions.
s2c<-data.frame(sample=paste0("output_",seq(1:length(files))),
                #state="injured",
                group=rep(c("BCI","CBCI","ACI","CACI"),each=5),
                condition=c(rep("uninjured",10),rep("injured",5),rep("uninjured",5)),
                parental=rep(c("sperm_before_injury","control_sample_sperm_before_injury","sperm_after_injury","control_sperm_after_injury"),each=5),
                group_easy=rep(c("1","2","3","4"),each=5),
                row.names = colnames(matrix_counts))
head(s2c)



# --- normalize the data with rlog or vsd  log parametric, normalized with respect to library size
vsd <- DESeq2::varianceStabilizingTransformation(dds, blind=FALSE,fitType = "local")
rld <- rlog(dds, blind=FALSE,fitType='local')
# this gives log2(n + 1)
ntd <- normTransform(dds)

png(filename = "./03_Graphs/miRNAs/Boxplot_Normalized_Counts_mirna_vsd.png",res = 300,width = 1440,height = 2560 )
boxplot(assay(vsd),main="Normalized")
dev.off()

png(filename = "./03_Graphs/miRNAs/Hist_Normalized_Counts_mirna_vsd.png",res = 300,width = 1440,height = 2560 )
hist(assay(vsd),main="Normalized")
dev.off()



mylist_res_cleaned_005_log2fc_1 <- lapply(mylist_res,function(x) subset(x,padj < 0.05,abs(log2FoldChange) > 1 ))



depiwiPCA <- function(dataSet,normalizeDataSet=assay(vsd),group,condition,name){
  #Get the differential expressed values from the comparrison interested,
  #extract the normalized values from the assay of vsd and plot them.
  #Giving condition and group from your design table
  
  #1. Extract the counts.
  dt <- normalizeDataSet
  dt_a <-dt[intersect(rownames(dt),dataSet),]
  cat("Counts values extracted...\n")
  Sys.sleep(0.2)
  
  #2. Perform pca
  pca_dt_a <- prcomp(t(dt_a))
  cat("PCA running...\n")
  Sys.sleep(0.2)
  
  #3. Extract percentVar data.
  percentVar_dt_a <- pca_dt_a$sdev^2/sum(pca_dt_a$sdev^2)
  cat("Percents calculated...\n")
  Sys.sleep(0.2)
  
  #4. Create the new dataframe to plot.
  dt_f <- data.frame(PC1=pca_dt_a$x[,1],
                     PC2=pca_dt_a$x[,2],
                     group=group,
                     condition=condition)
  cat("Data frame built...\n")
  Sys.sleep(0.2)
  
  #5. Plot it 
  cat("Plotting...\n")
  Sys.sleep(0.2)
  #pdf(file = paste0("./03_Graphs/piRNAs/PCA/PCA_DExpressed_pirnas_",name,".pdf"),width = 10,height = 12)
  png(filename =paste0("./03_Graphs/miRNAs/PCA/PCA_DExpressed_mirnas_",name,".png"),res = 300,width = 2560,height = 1440)
  cat("Saving plot as: ",paste0("./03_Graphs/miRNAs/PCA/PCA_DExpressed_mirnas_",name,"...\n"))
  pca_p <- ggplot(data = dt_f, aes_string(x = "PC1",
                                          y = "PC2",
                                          color = "group")) +
    geom_point(size = 3) +
    xlab(paste0("PC1: ", round(percentVar_dt_a[1] * 100), "% variance")) +
    ylab(paste0("PC2: ", round(percentVar_dt_a[2] * 100), "% variance")) + 
    coord_fixed()
  print(pca_p)
  dev.off()
  
  Sys.sleep(0.2)
  cat("Done")
  #return(pca_p)
  
  
  #plotly_pca_DE_pirna_acivscaci<- ggplotly(pca_depirnas_acivscaci)
  #htmlwidgets::saveWidget(plotly_pca_DE_pirna_acivscaci, "C:/Users/marius/Documents/Phd/Projects/smallRNAs/dockerDoor/03_Graphs/piRNAs/PCA/PCA_DE_pirna_acivscaci.html")
  #plotly_pca_DE_pirna_acivscaci
}


