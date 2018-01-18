# Read arguments---------------------------------------------------------------

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

# Read raw counttable -----------------------------------------------------

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
genetype.major <- genetype.info[readsperc > 0.001,]$gene_type
genetype.info <- genetype.info[gene_type %in% genetype.major,]
genetype.info <- rbind(list('other',100-sum(genetype.info$readsperc),0),genetype.info)
genetype.info$gene_type <- factor(genetype.info$gene_type,levels = genetype.info$gene_type)

p <- ggplot(genetype.info,aes(x='',y=readsperc, fill = gene_type))
p + geom_bar(stat = 'identity',width = 1) + 
  coord_polar(theta = 'y', start = 0) + 
  scale_fill_manual(values = rev(c(brewer.pal(9,'Set1'),brewer.pal(nrow(genetype.info)-10,'Set3'),'grey50')))+
  theme_xf +
  theme(axis.ticks = element_blank(),
        axis.title = element_blank(),
        axis.text = element_blank(),
        axis.line = element_blank())
library(plotly)

p <- plotly(genetype.info,labels = ~ gene_type, values = ~readsperc, type = 'pie', textposition = 'outside')
  

###################################################################stoped here

#Differential expressed
rpkm <- dif.tb[,tail(seq(ncol(dif.tb)),groupn),with=F] %>% rowSums()*1000*1000000/(dif.tb$width*total.reads)
write.csv(dif.tb[rpkm >= 0.1,],file = 'expressed.csv',row.names = F)

plot.table <- dif.tb
plot.table[,`:=`(rpkm=rpkm,express = 'expressed')]
plot.table[rpkm < 0.1, express := 'non-expressed']

plot <- ggplot(plot.table,aes(x =gene_type, fill = express))
plot + geom_bar() + theme_xf + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('expressed_transcripts_TIP60.jpg', width = 6, height = 6, dpi =150, units = 'in')

dif.deseq2 <- DefExp('expressed.csv',groups = c('siCTL','siCTL','siTIP60','siTIP60'),ref = 'siCTL')


#rpkm for expression in all samples threshold.

temp <- dif.deseq2
temp$regulation <- 'NC'
temp[padj < 0.5 & log2FoldChange > 0,]$regulation <- 'Increase'
temp[padj < 0.5 & log2FoldChange < 0,]$regulation <- 'Decrease'


test <- ggplot(temp,aes(x = gene_type, fill = regulation))
test + geom_bar(width = 0.9) + theme_xf + theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave('dif_expressed_transcripts_TIP60_loss.jpg', width = 6, height = 6, dpi =150, units = 'in')



