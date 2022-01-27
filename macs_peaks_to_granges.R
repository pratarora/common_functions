## Function to convert macs peaks to granges
convert_peaks_to_granges <- function(df){
  
  library(GenomicRanges)
  
  peaks_GR <- GRanges(
    seqnames=df[,"chr"],
    IRanges(df[,"start"],
            df[,"end"]
    ),
    strand= NULL,
    length= df$length,
    abs_summit= df$abs_summit,
    pileup= df$pileup,
    pval= df$`X.log10.pvalue.`,
    fdr= df$`X.log10.qvalue.`,
    fold_enrichment= df$fold_enrichment,
    name = df$name
  )
  
  
  peaks_GR
  
}

peakfile <-paste0("Test_1_peaks.xls") #this file is present in the git folder; it is an output file of macs3
peaks_DF <- read.delim2(peakfile, comment.char="#")
peak_gr <- convert_peaks_to_granges(peaks_DF)
peak_gr