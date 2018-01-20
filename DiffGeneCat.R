# Read arguments---------------------------------------------------------------
title <- 'siTIP60'
counttable.raw.file <-  'TIP60.raw.counttable.txt'
coordinate.file <-  'coordinates/Hela.transcripts.all.gtf'
groups <- c('siCTL','siCTL','siTIP60','siTIP60')
groupn <- length(groups)
min.rpkm <- 0.1
# Load packages and functions needed ---------------------------------------------------
library(GenomicRanges)
library(magrittr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stringr)
require(TxDb.Hsapiens.UCSC.hg19.knownGene)
require(org.Hs.eg.db)
for(i in list.files('r.func/')){
  load(paste0('r.func/',i))
}

# Change transcripts to GR -------------------------------------------------
transcripts.gr <- GTF2GR('coordinates/Hela.transcripts.all.gtf') %>% sort()

# Infomation from raw counttable -----------------------------------------------------
counttable <- fread(counttable.raw.file)
Info <- colnames(counttable)[1]

# combine coordinate info
counttable <- makeGRangesFromDataFrame(counttable[,-1],keep.extra.columns = T) %>%sort()
counttable <- cbind(as.data.table(transcripts.gr),mcols(counttable)[,tail(seq(ncol(mcols(counttable))),groupn)]) %>% as.data.table()

# select major gene_type.
total.reads <- tail(colnames(counttable),groupn) %>% str_extract('\\d+\\.0') %>% as.numeric() %>% sum()

temp <- counttable[,c(tail(colnames(counttable),groupn)),with = F] %>% rowSums()
genetype.info <- data.table(reads=temp,gene_type = counttable$gene_type)
genetype.info <- genetype.info[,.(readsperc= sum(reads)*100/total.reads,
                                  count=.N),by = gene_type]
genetype.info <- genetype.info[order(readsperc)]

# draw pie plot of gene_type.
genetype.major <- genetype.info[readsperc > 0.001 & count > 200,][order(-count)]$gene_type %>% as.character()

genetype.info <- genetype.info[gene_type %in% genetype.major,]
genetype.info <- rbind(list('other',100-sum(genetype.info$readsperc),0),genetype.info)
genetype.info$gene_type <- factor(genetype.info$gene_type,levels = genetype.info$gene_type)
genetype.info[gene_type == 'other',count:=NA]

library(plotly)
plot.table <- genetype.info
plot.table$gene_type <- paste0(plot.table$gene_type,'(',plot.table$count,')')
p <- plot_ly(plot.table,labels = ~gene_type, values = ~readsperc, type = 'pie', 
             textposition = 'outside',
             textinfo = 'percent',
             marker = list(colors = colors,
                           line = list(color = '#FFFFFF',
                                       width = 1)))
p <- layout(p,
       title = 'Distribution of GRO-seq reads on genomic elements',
       margin = list(l=100,t=80,b=50),
       showlegend = T)
p
plotly_IMAGE(p, width = 700, height = 700, out_file = 'figures/GRO-seq_reads_distribution.png')

#If error, need:
Sys.setenv("plotly_username"="Feng_Xiong") 
Sys.setenv("plotly_api_key"="cKUazDNFFbwzc7m8EGdj")


#Draw expressed percentage.
rpkm <- counttable[,tail(colnames(counttable),groupn),with=F] %>% rowSums()*10E9/(counttable$width*total.reads)
write.csv(counttable[rpkm >= 0.1,],file = 'expressed.csv',row.names = F)

#Plot the expressing file
plot.table <- counttable
plot.table$rpkm <- rpkm
plot.table[,`:=`(express = 'expressed')]
plot.table[rpkm < 0.1, express := 'non-expressed']
plot.table <- plot.table[,.N,by=.(gene_type,express)] %>% merge(genetype.info)
plot.table <- plot.table[order(-count)][,`:=`(experc=paste0(round(N/count*100),'%'),
                                              gene_type= factor(gene_type, 
                                                                 levels = unique(gene_type)))]
plot.table[express == 'non-expressed', experc := NA]

plot <- ggplot(plot.table[count > 200,],aes(x =gene_type, y= N, fill = express, label = experc))
plot + geom_bar(stat = 'identity') + theme_xf + 
  labs(title = 'expression of genomic elements (rpkm > 0.1)',x=NULL, y = 'Count') +
  geom_text(aes(y = count),vjust = -0.5, size =3.5)+
  scale_fill_manual(values = c('firebrick3','grey50'))+
  ylim(0,max(plot.table$count *1.1))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(.85,.75))
ggsave(filename = 'figures/genome_expression.jpg', width = 6, height = 6, dpi =150, units = 'in')


# DESeq2 analysis ---------------------------------------------------------
exptable <- counttable[rpkm >= min.rpkm,]
res.deseq2 <- DefExp(exptable,groups = c('siCTL','siCTL','siTIP60','siTIP60'),ref = 'siCTL')
#On expressed genes

plot.table <- res.deseq2[gene_type %in% genetype.major,][,gene_type:=factor(gene_type,
                                                                            levels = genetype.major)]
plot.table$regulation <- 'NC'
plot.table[pvalue < 0.05 & log2FoldChange > 0,regulation := 'Increase']
plot.table[pvalue < 0.05 & log2FoldChange < 0,regulation := 'Decrease']
plot.table[,regulation := factor(regulation, levels = c('Increase','Decrease','NC'))]
res.deseq2 <- plot.table

plot <- ggplot(plot.table,aes(x = gene_type, fill = regulation))
plot + geom_bar(width = 0.9) + theme_xf + theme(axis.text.x = element_text(angle = 90, hjust = 1))+
  labs(title = 'expression after siTIP60 (p < 0.05)',x =NULL)+
  scale_fill_manual(values = c('firebrick2','royalblue3','grey60')) +
    theme(legend.position = c(.85,.7),
          plot.title = element_text(hjust = 0.5))
ggsave('figures/dif_expressed_transcripts_siTIP60.jpg', width = 6, height = 6, dpi =150, units = 'in')




# Write Increased/decreased genes to bed generate matrix and plot-----------
GROPlotonTemp <- function(gr,inputfix,name,figurefile){
  Direction2Gene(paste0('matrix/Hela-siCTL-rpt1-GROpos_',inputfix,'_scaled.mx'),
                 paste0('matrix/Hela-siCTL-rpt1-GROneg_',inputfix,'_scaled.mx'),
                 temp,'siCTL-rpt1')
  
  Direction2Gene(paste0('matrix/Hela-siTIP60-rpt1-GROpos_',inputfix,'_scaled.mx'),
                 paste0('matrix/Hela-siTIP60-rpt1-GROneg_',inputfix,'_scaled.mx'),
                 temp,'siTIP60-rpt1')
  
  index <- as.data.table(temp)$seqnames != 'chrM'
  plot.dt <- data.table(siCTL = colSums(`siCTL-rpt1.sense.mx`[index,])/length(temp[index]),
                        siCTL_AS = -colSums(`siCTL-rpt1.antisense.mx`[index,])/length(temp[index]),
                        siTIP60 = colSums(`siTIP60-rpt1.sense.mx`[index,])/length(temp[index]),
                        siTIP_AS = -colSums(`siTIP60-rpt1.antisense.mx`[index,])/length(temp[index]))
  MxProfilePlot(plot.dt, breaks =c(1,61,261,320)) +
    scale_color_manual(values = c('firebrick2','firebrick3','royalblue2','royalblue3'))+
    labs(title = name) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(figurefile, height = 6, width = 8, dpi = 150, units = 'in')
}
ChIPPlotonTemp <- function(temp,antibody,inputfix,name,figurefile){
  siCTL.mx <- ReadMatrix(paste0('matrix/Hela-',antibody,'-SC-siCTL_',inputfix,'_scaled.mx'))
  siTIP60.mx <- ReadMatrix(paste0('matrix/Hela-',antibody,'-SC-siTIP60_',inputfix,'_scaled.mx'))
  index <- as.data.table(temp)$seqnames != 'chrM'
  plot.dt <- data.table(siCTL = colSums(siCTL.mx[index,])/length(temp[index]),
                        siTIP60 = colSums(siTIP60.mx[index,])/length(temp[index]))
  MxProfilePlot(plot.dt, breaks =c(1,61,261,320)) +
    scale_color_manual(values = c('firebrick2','royalblue2'))+
    labs(title = name) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(figurefile, height = 6, width = 8, dpi = 150, units = 'in')
}
#NC_coding_genes
temp <- res.deseq2[gene_type == 'protein_coding' & regulation == 'NC',] %>% makeGRangesFromDataFrame()
GR2BED(temp,'coordinates/NC_genes_siTIP60.bed')

GROPlotonTemp(temp,'NC_genes','NC_protein_coding_genes','figures/GRO-seq_on_NC_protein_coding.jpg')
ChIPPlotonTemp(temp,'H3K27ac','NC_genes','H3K27ac_NC_protein_coding_genes','figures/H3K27ac_on_NC_protein_coding.jpg')
ChIPPlotonTemp(temp,'H3K4me1','NC_genes','H3K4me1_NC_protein_coding_genes','figures/H3K4me1_on_NC_protein_coding.jpg')
ChIPPlotonTemp(temp,'Pol2','NC_genes','Pol2_NC_protein_coding_genes','figures/Pol2_on_NC_protein_coding.jpg')

#Increase_coding_genes
temp <- res.deseq2[gene_type == 'protein_coding' & regulation == 'Increase',] %>% makeGRangesFromDataFrame()
GR2BED(temp,'coordinates/Increase_genes_siTIP60.bed')

GROPlotonTemp(temp,'Increase_genes','Increased_protein_coding_genes','figures/GRO-seq_on_increase_protein_coding.jpg')
ChIPPlotonTemp(temp,'H3K27ac','Increase_genes','H3K27ac_Increase_protein_coding_genes','figures/H3K27ac_on_Increase_protein_coding.jpg')
ChIPPlotonTemp(temp,'H3K4me1','Increase_genes','H3K4me1_Increase_protein_coding_genes','figures/H3K4me1_on_Increase_protein_coding.jpg')
ChIPPlotonTemp(temp,'Pol2','Increase_genes','Pol2_Increase_protein_coding_genes','figures/Pol2_on_Increase_protein_coding.jpg')

temp <- res.deseq2[gene_type == 'protein_coding' & regulation == 'Decrease',] %>% makeGRangesFromDataFrame()
GR2BED(temp,'coordinates/Decrease_genes_siTIP60.bed')

GROPlotonTemp(temp,'Decrease_genes','Decreased_protein_coding_genes','figures/GRO-seq_on_decrease_protein_coding.jpg')
ChIPPlotonTemp(temp,'H3K27ac','Decrease_genes','H3K27ac_Decrease_protein_coding_genes','figures/H3K27ac_on_Decrease_protein_coding.jpg')
ChIPPlotonTemp(temp,'H3K4me1','Decrease_genes','H3K4me1_Decrease_protein_coding_genes','figures/H3K4me1_on_Decrease_protein_coding.jpg')
ChIPPlotonTemp(temp,'Pol2','Decrease_genes','Pol2_Decrease_protein_coding_genes','figures/Pol2_on_Decrease_protein_coding.jpg')



GROPlotonTemp <- function(gr,inputfix,name,figurefile){
  Direction2Gene(paste0('matrix/Hela-siCTL-rpt1-GROpos_',inputfix,'_scaled.mx'),
                 paste0('matrix/Hela-siCTL-rpt1-GROneg_',inputfix,'_scaled.mx'),
                 temp,'siCTL-rpt1')
  
  Direction2Gene(paste0('matrix/Hela-siTIP60-rpt1-GROpos_',inputfix,'_scaled.mx'),
                 paste0('matrix/Hela-siTIP60-rpt1-GROneg_',inputfix,'_scaled.mx'),
                 temp,'siTIP60-rpt1')
  
  plot.dt <- data.table(siCTL = colSums(`siCTL-rpt1.sense.mx`)/length(temp),
                        siCTL_AS = -colSums(`siCTL-rpt1.antisense.mx`)/length(temp),
                        siTIP60 = colSums(`siTIP60-rpt1.sense.mx`)/length(temp),
                        siTIP_AS = -colSums(`siTIP60-rpt1.antisense.mx`)/length(temp))
  MxProfilePlot(plot.dt, breaks =c(1,61,120),labels = c('-1.5kb','center','1.5kb')) +
    scale_color_manual(values = c('firebrick2','firebrick3','royalblue2','royalblue3'))+
    labs(title = name) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(figurefile, height = 6, width = 8, dpi = 150, units = 'in')
}
temp <- res.deseq2[gene_type == 'enhancer' & regulation == 'NC',] %>% makeGRangesFromDataFrame()
GR2BED(temp,'coordinates/NC_enhancer_siTIP60.bed')
GROPlotonTemp(temp,'NC_enhancer','NC_enhancer_eRNAs','figures/GRO-seq_on_NC_enhancers.jpg')

temp <- res.deseq2[gene_type == 'enhancer' & regulation == 'Increase',] %>% makeGRangesFromDataFrame()
GR2BED(temp,'coordinates/Increase_enhancer_siTIP60.bed')

GROPlotonTemp(temp,'Increase_enhancer','Increase_enhancer_eRNAs','figures/GRO-seq_on_increase_enhancers.jpg')


temp <- res.deseq2[gene_type == 'enhancer' & regulation == 'Decrease',] %>% makeGRangesFromDataFrame()
GR2BED(temp,'coordinates/Decrease_enhancer_siTIP60.bed')
GROPlotonTemp(temp,'Decrease_enhancer','Decrease_enhancer_eRNAs','figures/GRO-seq_on_decrease_enhancers.jpg')

# Differences in Gene fetures_length_exons --------------------------------
#RPKM
plot.table <- res.deseq2[gene_type == 'protein_coding' & seqnames != 'chrM',]
plot.table$rpkm <- select(plot.table,contains('Homer')) %>% rowSums()*10E9/(plot.table$width*total.reads)
plot.table[,regulation := factor(regulation, levels = c('Decrease','NC','Increase'))]
p <- ggplot(plot.table,aes(x= regulation, y = rpkm))
temp <- boxplot.stats(plot.table$rpkm)$stats[c(1,5)]
p + geom_boxplot(outlier.shape = NA, lwd =1) +
  ylim(temp*1.05) +
  theme_xf +
  labs(title = 'Expression of protein coding genes',x =NULL) +
  theme(plot.title = element_text(hjust =0.5))

ggsave(filename = 'figures/Expression_boxplot_coding.jpg', width = 4, height = 4, dpi =150, units = 'in')

#exons.
genes.region <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
exons.region <- exons(TxDb.Hsapiens.UCSC.hg19.knownGene)
test <- findOverlaps(genes.region,exons.region) %>% as.data.table()
test <- test[,.N, by = queryHits]
exon.count <- data.table(ENTREZID=genes.region$gene_id, exon.n = test$N)
plot.table <- res.deseq2[gene_type == 'protein_coding' & seqnames != 'chrM',]
temp <- AnnotationDbi::select(org.Hs.eg.db,plot.table$gene_name,'ENTREZID',keytype = 'SYMBOL')
plot.table <- merge(plot.table, temp,by.x = 'gene_name',by.y ='SYMBOL',all.x = F,all.y =F) 
plot.table <- merge(plot.table,exon.count,by = 'ENTREZID',all.x = T, all.y = F)


plot.table[exon.n > 30, exon.n := 30]
plot.table <- plot.table[!is.na(exon.n),]
plot.table[,exon.n := as.character(exon.n)][exon.n == '30', exon.n := '30+']
plot.table[,exon.n := factor(exon.n,levels = c(as.character(1:29),'30+'))]
p <- ggplot(plot.table,aes(x= exon.n, y = log2FoldChange))
p + geom_boxplot()

test <- plot.table[,.(.N,Mean = median(log2FoldChange),SD = sd(log2FoldChange)), by = exon.n]
test[,se:= SD/sqrt()]

# Differential gene analysis ----------------------------------------------
hk_gene <- fread('misc/HK_genes.txt',header = F) %>% `colnames<-`(c('gene_name','refseq'))
hk_gene$housekp <- 'housekeeping'
res.deseq2 <- merge(res.deseq2,hk_gene,all.x = T)
res.deseq2[is.na(refseq),housekp := 'regulated gene']

library(plotly)
plot.table <- res.deseq2[gene_type == 'protein_coding',][,.N,by=.(housekp)]
p <- plot_ly(plot.table,labels = ~housekp,values= ~N, type = 'pie', 
             textposition = 'inside',
             textinfo = 'value+label+percent',
             insidetextfont = list(color = '#FFFFFF',
                                   size = 20),
             marker = list(colors = colors,
                           line = list(color = '#FFFFFF',
                                       width = 1)))
p <- layout(p, title = 'All expressing genes',
       font = list(size =20),
       margin = list(t=100)
       )
plotly_IMAGE(p, width = 600, height = 800,format = 'png', out_file = 'figures/All_expressed_hkgenes.png')


plot.table <- res.deseq2[gene_type == 'protein_coding' & regulation == 'Decrease',][,.N,by=.(housekp)]
p <- plot_ly(plot.table,labels = ~housekp,values= ~N, type = 'pie', 
             textposition = 'inside',
             textinfo = 'value+label+percent',
             insidetextfont = list(color = '#FFFFFF',
                                   size = 20),
             marker = list(colors = colors,
                           line = list(color = '#FFFFFF',
                                       width = 1)))
p <- layout(p, title = 'Genes decreased by siTIP60',
            font = list(size =20),
            margin = list(t=100)
)
plotly_IMAGE(p, width = 600, height = 800,format = 'png', out_file = 'figures/Genes_decreased_hkgenes.png')



plot.table <- res.deseq2[gene_type == 'protein_coding' & regulation == 'Increase',][,.N,by=.(housekp)]
p <- plot_ly(plot.table,labels = ~housekp,values= ~N, type = 'pie', 
             textposition = 'inside',
             textinfo = 'value+label+percent',
             insidetextfont = list(color = '#FFFFFF',
                                   size = 20),
             marker = list(colors = colors,
                           line = list(color = '#FFFFFF',
                                       width = 1)))
p <- layout(p, title = 'Genes increased by siTIP60',
            font = list(size =20),
            margin = list(t=100)
)
plotly_IMAGE(p, width = 600, height = 800,format = 'png', out_file = 'figures/Genes_increased_hkgenes.png')


# GO analysis of increase/decresed genes ----------------------------------


