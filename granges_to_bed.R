## Convert macs3 peaks to ".bed" for ChipSeeker


convert_peaks_to_granges <- function(df){ # df is a data frame with the macs3 peaks
  
  library(GenomicRanges) # load GenomicRanges package
  
  # convert peaks to GRanges
  peaks_GR <- GRanges(
    seqnames=df[,"chr"], # chromosome
    IRanges(df[,"start"], # start
            df[,"end"] # end
    ),
    strand= NULL, # strand 
    length= df$length, # length of peak
    abs_summit= df$abs_summit, # absolute summit
    pileup= df$pileup, # pileup of reads
    pval= df$`X.log10.pvalue.`, # p-value
    fdr= df$`X.log10.qvalue.`, # FDR
    fold_enrichment= df$fold_enrichment, # fold enrichment
    name = df$name # name of peak
  )
  
  # return peaks_GR
  peaks_GR
  
}

# convert peaks to GRanges
peakfile <-paste0("Test_1_peaks.xls") #this file is present in the git folder; it is an output file of macs3

# read peaks
peaks_DF <- read.delim2(peakfile, comment.char="#")

# convert to GRanges
peak_gr <- convert_peaks_to_granges(peaks_DF)

# add types
peak_gr <- add_types_to_peaks(peak_gr, "peaks")

peak_gr


# This code converts the bedgraph files into a GRanges object
# It does this by reading in the bedgraph file and using the data.frame function to
# convert the data into a GRanges object. It then saves the GRanges object as a
# variable with the same name as the sample name

# The variables that are used in this code are:
# peakfile = the name of the bedgraph file that is being read in
# peaks_DF = the data.frame that is created from the bedgraph file
# gr = the GRanges object that is created from the data.frame

# The purpose of this code is to convert the bedgraph files into GRanges objects
# so that they can be used in the next code chunk to create the UCSC tracks

# The context in which this code is used is in a for loop that iterates through
# all of the bedgraph files that are in the bedgraph_files directory

library(rtracklayer)
for (i in 1:length(peaks_files)){

  peakfile <-paste0("./macs3_peaks_03/",sample_names[i],"_peaks.xls") #this file is present in the git folder; it is an output file of macs3
  # print(peakfile)
  peaks_DF <- read.delim2(peakfile, comment.char="#") # read peaks file 
  assign(paste0(sample_names[i], "_granges"), convert_peaks_to_granges(peaks_DF)) # convert to GRanges and save as variable with same name as sample name
  gr= get(paste0(sample_names[i], "_granges")) # get the GRanges object
  gr %>% head() %>% print() # print the GRanges object
  df <- data.frame(seqnames=seqnames(gr), # chromosome 
         starts=start(gr)-1, # start
         ends=end(gr), # end
         names=gr$name, # name of peak 
         scores=gr$fold_enrichment, # fold enrichment
         length= elementMetadata(gr)$length, # length of peak 
         strands=strand(gr)# strand
          )
print("df created") # print that the df was created
write.table(df, file=paste0("./macs3_peaks_03/", sample_names[i], "_granges",".bed"), quote=F, sep="\t", row.names=F, col.names=F) # write the df to a bed file
  # print that the new variable is now stored in the global environment
  print(paste("new variable with peaks in granges called", paste0(sample_names[i], "_granges"), "is now stored in the global environment", sep = " "))
}
 # print that the new variables are now stored in the global environment
head(URV_1_granges)