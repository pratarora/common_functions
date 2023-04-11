## Function to convert macs peaks to granges
convert_peaks_to_granges <- function(df){
  
  library(GenomicRanges) # load GenomicRanges package
  
  peaks_GR <- GRanges( # convert peaks to GRanges
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
  
  
  peaks_GR # return peaks_GR
  
}

peakfile <-paste0("Test_1_peaks.xls") #this file is present in the git folder; it is an output file of macs3
peaks_DF <- read.delim2(peakfile, comment.char="#") # read peaks
peak_gr <- convert_peaks_to_granges(peaks_DF) # convert to GRanges
peak_gr # print peak_gr