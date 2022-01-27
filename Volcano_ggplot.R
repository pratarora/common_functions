####################
#                  #
# VOLCANO FUNCTION #
#                  #
####################
require(ggplot2)
require(ggrepel)
require(clusterProfiler)
require(tidyverse)

draw_volcano<- function(fileinput, title) {
  # read input file
  # drawing plots
  ggplot(data =fileinput , aes(x = logFoldChange, y = -log10(padj))) +
    # draw lines
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", col = "chartreuse4") +
    geom_hline(yintercept = -log10(0.05), linetype = "dotted", col = "darkgoldenrod") +
    geom_vline(xintercept = 0, linetype = "dashed")+
    # draw points
    geom_point(x = fileinput$logFoldChange, y = -log10(fileinput$padj),alpha = 0.5, size = 2,color="grey51") +
    # draw coloured points
    geom_point(data = fileinput[which(fileinput$padj < 0.05 & fileinput$logFoldChange < -1),],
               aes(x=logFoldChange, y = -log10(padj)), shape = 21, color = "royalblue3", fill = "royalblue3",
               alpha = 0.5, size = 2) +
    geom_point(data = fileinput[which(fileinput$padj < 0.05 & fileinput$logFoldChange > 1),],
               aes(x=logFoldChange, y = -log10(padj)), shape = 21, color = "red3", fill = "red3",
               alpha = 0.5, size = 2) +
    # x axis scale
    scale_x_continuous(breaks = seq(round(min(fileinput$logFoldChange)- 0.5),round(max(fileinput$logFoldChange)+ 0.5),by = 1), limits = c(round(min(fileinput$logFoldChange)-1),round(max(fileinput$logFoldChange)+1))) + xlab("logFoldChange") + #ylab("-Log10(p.value)") +
    
    scale_y_continuous(breaks = seq(0,round(-log10(min(fileinput$padj))+1),by = 4), limits = c(0,round(-log10(min(fileinput$padj))+1))) + ylab("-Log10(pAdjusted)") + # ylab("-log10(p.value)")+
    
    # set title
    ggtitle(title)+
    # x and y axis limits
    # black and white theme
    theme_bw() +
    # center title
    theme(plot.title = element_text(hjust = 0.5), axis.text = element_text(size = 10), axis.title.x = element_text(size = 10),
          axis.title.y = element_text(size = 10))
}
###################

b = draw_volcano(so_uninjured_offspring_walt_table_volcano,"Uninjured hearts from\nnon-infarted parental vs infarted parentals")

#Set genes for marking 
wt_set<-c("WT1A","MPEG1.1","WT1B","TCF12","TNFA","IL1B","MMP14A","MMP14B","MYCA","MYCB")#"VCANB","MAFBA","TIMP2B","RLN1","TLR7","ARG2","C1QA","C1QB","MYCA","IRF8","SPI1A","LCP1","TLR9","MAFBB","RLN1") #tcf21
wt_set
wt_geneset<-as.data.frame(so_uninjured_offspring_walt_table_volcano[toupper(so_uninjured_offspring_walt_table_volcano$SYMBOL)%in%wt_set,])

require(ggrepel)
#Paint the genes in the plot
c= b + geom_point(data=so_uninjured_offspring_walt_table_volcano[toupper(so_uninjured_offspring_walt_table_volcano$SYMBOL)%in%wt_set,],color=ifelse(so_uninjured_offspring_walt_table_volcano[toupper(so_uninjured_offspring_walt_table_volcano$SYMBOL)%in%wt_set,]$logFoldChange > 0,"midnightblue","mediumvioletred"),size=2) +
  geom_text_repel(data = so_uninjured_offspring_walt_table_volcano[toupper(so_uninjured_offspring_walt_table_volcano$SYMBOL)%in%wt_set,],aes(label=so_uninjured_offspring_walt_table_volcano[toupper(so_uninjured_offspring_walt_table_volcano$SYMBOL)%in%wt_set,]$SYMBOL),
  )#nudge_x = 0,
#nudge_y = 2,segment.size = 0.1)
c


#Paint the top genes in the y Axis
d=c + geom_point(data=so_uninjured_offspring_walt_table_volcano[so_uninjured_offspring_walt_table_volcano$padj < min(so_uninjured_offspring_walt_table_volcano$padj*100),] ,color=ifelse(so_uninjured_offspring_walt_table_volcano[so_uninjured_offspring_walt_table_volcano$padj < min(so_uninjured_offspring_walt_table_volcano$padj*100),]$logFoldChange > 0,"midnightblue","mediumvioletred"),size=2) +
  geom_text_repel(data = so_uninjured_offspring_walt_table_volcano[so_uninjured_offspring_walt_table_volcano$padj < min(so_uninjured_offspring_walt_table_volcano$padj*100),],aes(label=so_uninjured_offspring_walt_table_volcano[so_uninjured_offspring_walt_table_volcano$padj < min(so_uninjured_offspring_walt_table_volcano$padj*100),]$SYMBOL),
  )#nudge_x = 0,
#nudge_y = 2,segment.size = 0.1)

d
#Save the plot
pdf("/home/marius/Desktop/epis_andres_old_rnaseq/R/uninjured_offspring/Volcano/so_uninjured_offspring_walt_table_volcano.pdf",height = 15,width = 15)
d
