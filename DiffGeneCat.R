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

for(i in list.files('r.func/')){
  load(paste0('r.func/',i))
}

# Change transcripts to GR -------------------------------------------------
transcripts.gr <- GTF2GR('coordinates/Hela.transcripts.all.gtf')

# Infomation from raw counttable -----------------------------------------------------

counttable <- fread(counttable.raw.file)
Info <- colnames(counttable)[1]

# combine coordinate info
counttable <- makeGRangesFromDataFrame(counttable[,-1],keep.extra.columns = T) %>% sort()
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
#Sys.setenv("plotly_username"="Feng_Xiong") 
#Sys.setenv("plotly_api_key"="cKUazDNFFbwzc7m8EGdj")


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

# Differential gene analysis ----------------------------------------------
hk_gene <- fread('misc/HK_genes.txt',header = F) %>% `colnames<-`(c('gene_name','refseq'))
hk_gene$housekp <- 'housekeeping'
res.deseq2 <- merge(res.deseq2,hk_gene,all.x = T)
res.deseq2[is.na(refseq),housekp := 'other']

library(plotly)
plot.table <- res.deseq2[gene_type == 'protein_coding',][,.N,by=housekp]
p <- plot_ly(res.deseq2[gene_type == 'protein_coding',],labels = ~housekp, type = 'pie', 
             textposition = 'inside',
             textinfo = 'label',
             marker = list(colors = colors,
                           line = list(color = '#FFFFFF',
                                       width = 1)))
layout(p,
            margin = list(l=100,t=80,b=50),
            showlegend = T)
p
