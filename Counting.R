library(data.table)
library(ggplot2)
library(magrittr)
library(stringr)

# Read --------------------------------------------------------------------

genes <- fread('genecount.txt')
exons <- fread('exoncount.txt')
introns <- fread('introncount.txt')
utr3 <- fread('3utrcount.txt')
utr5 <- fread('5utrcount.txt')

total.reads <- str_extract(colnames(exons[,9:14]),"\\(\\d+") %>% str_replace("\\(","") %>% as.numeric()

reads.list <- list(genes,exons, introns, utr3, utr5)

count.list <- NULL
for (i in seq(length(reads.list))) {
  dt <- reads.list[[i]][,9:14] %>% colSums()
  dt <- dt/total.reads
  count.list[[i]] <- dt %>% as.data.table()
}
count.tb <- as.data.table(count.list)
colnames(count.tb) <- c('genes','exons','introns','3utr','5utr')
count.tb <- cbind(beads = c(rep('C1513',6)),count.tb)
count.tb <- cbind(samples = colnames(exons[,9:14]), count.tb )

write.csv(count.tb,file = './count.tb.csv')
