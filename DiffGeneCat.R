
# Library -----------------------------------------------------------------
library(GenomicRanges)
library(magrittr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stringr)

# Functions ---------------------------------------------------------------
load('r.func/ImportBiopacks.rdata')

diff.file = 'TIP60.diffoutput.txt'
coordinate.file = 'coordinates/gencode.v19.annotation.geneonly.gtf'


# Change coordinate to GR -------------------------------------------------
coordinate <- GTF2GR('coordinates/Hela.transcripts.all.gtf')
# Diff expression analysis (start from homer diff output) -----------------

groups <- c('siCTL','siCTL','siTIP60','siTIP60')
groupn <- length(groups)
min.exp <- 0.1
#Starting parameter

dif.exp <- fread('TIP60.raw.counttable.txt')
Info <- colnames(dif.exp)[1]

dif.exp <- makeGRangesFromDataFrame(dif.exp[,-1],keep.extra.columns = T) %>% sort()
dif.tb <- cbind(mcols(coordinate),as.data.table(dif.exp)) %>% as.data.table()



#Summary genetype.info
total.reads <- tail(colnames(dif.tb),groupn) %>% str_extract('\\d+\\.0') %>% as.numeric() %>% sum()
reads <- dif.tb[,c(tail(colnames(dif.tb),groupn)),with = F] %>% rowSums()
genetype.info <- data.table(reads=reads,gene_type = dif.tb$gene_type)
genetype.info <- genetype.info[,.(sumreads= sum(reads)*100/total.reads),by = gene_type]

temp <- table(dif.tb$gene_type)
temp <- data.table(count = as.numeric(temp), gene_type = names(temp))
genetype.info <- merge(genetype.info,temp)
genetype.major <- genetype.info[count > 200,]$gene_type

dif.tb <- dif.tb[gene_type %in% genetype.major,]

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



