## Convert macs3 peaks to ".bed" for ChipSeeker

```{r}

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


library(rtracklayer)
for (i in 1:length(peaks_files)){

  peakfile <-paste0("./macs3_peaks_03/",sample_names[i],"_peaks.xls")
  # print(peakfile)
  peaks_DF <- read.delim2(peakfile, comment.char="#")
  assign(paste0(sample_names[i], "_granges"), convert_peaks_to_granges(peaks_DF))
  gr= get(paste0(sample_names[i], "_granges"))
  gr %>% head() %>% print()
  df <- data.frame(seqnames=seqnames(gr),
         starts=start(gr)-1,
         ends=end(gr),
         names=gr$name,
         scores=gr$fold_enrichment,
         length= elementMetadata(gr)$length,
         strands=strand(gr)
          )
print("df created")
write.table(df, file=paste0("./macs3_peaks_03/", sample_names[i], "_granges",".bed"), quote=F, sep="\t", row.names=F, col.names=F)

  print(paste("new variable with peaks in granges called", paste0(sample_names[i], "_granges"), "is now stored in the global environment", sep = " "))
}
head(URV_1_granges)