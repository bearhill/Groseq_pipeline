# libraries ---------------------------------------------------------------
library(GenomicRanges)
library(magrittr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stringr)

for(i in list.files('r.func/')){
  load(paste0('r.func/',i))
}

for(i in list.files('r.data/')){
  load(paste0('r.data/',i))
}

# Theme -------------------------------------------------------------------

save.image(file = ".RData")


SimpleVenn <- function(N1,N2,N12,title1,title2,labelposi= c(list(c(0.5,0),list(0.5,0)))){
  library(VennDiagram)
  temp <- draw.pairwise.venn(area1 = N1,
                             area2 = N2,
                             cross.area = N12, 
                             category = c(title1, title2),
                             col = c('black','black'),
                             fill = c('red','blue'),
                             cex = rep(1.5,3),
                             fontfamily = rep('sans',3),
                             cat.cex = rep(1.5,2),
                             cat.pos = c(-160,160),
                             cat.col = c('black','black'),
                             cat.just = labelposi,
                             cat.fontfamily = rep('sans',2),
                             margin = 0.05)
  grid.draw(temp)
}
save(SimpleVenn, file = 'r.func/SimpleVenn.rdata')

# Temp --------------------------------------------------------------------

MxPeakcenterPlot <- function(plot.dt,breaks = c(1,60,120), labels = c('-1.5kb','Peakcenter','1.5kb')){
  plot.dt <- cbind(plot.dt,bins = 1:nrow(plot.dt))
  plot.dt <- melt(plot.dt, id.vars = 'bins')
  p <- ggplot(plot.dt, aes(x=bins, y = value, group = variable, color = variable))
  p <- p + geom_line(size =1) +
    scale_x_continuous(breaks = breaks, labels = labels)+
    labs(x = NULL, y = 'Average Reads')+
    theme_xf+
    theme(legend.position = 'bottom')
  p  
}

save(MxPeakcenterPlot, file = 'r.func/MxPeakcenterPlot.rdata')


GRGeneInfo <- function(gr){
  require(GenomicRanges)
  require(TxDb.Hsapiens.UCSC.hg19.knownGene)
  require(org.Hs.eg.db)
  genes <- genes(TxDb.Hsapiens.UCSC.hg19.knownGene)
  symbols <- select(org.Hs.eg.db,genes$gene_id,'SYMBOL',keytype = 'ENTREZID')$SYMBOL
  temp <- nearest(gr,genes)
  gr$gene <- symbols[temp]
  gr$gene_id <- genes$gene_id[temp]
  gr
}

save(GRGeneInfo,file = 'r.func/GRGeneInfo.rdata')

Direction2Gene <- function(posmxfile, negmxfile, grmap, prefix){
  mx1 <- fread(posmxfile) %>% as.matrix()
  mx2 <- fread(negmxfile) %>% as.matrix()
  mx1[is.na(mx1)] <- 0
  mx1[is.nan(mx1)] <- 0
  mx2[is.na(mx2)] <- 0
  mx2[is.nan(mx2)] <- 0
  neg.index <- as.logical(strand(grmap) == '-')
  mx.swap <- mx1[neg.index,]
  mx1[neg.index,] <- mx2[neg.index,]
  mx2[neg.index,] <- mx.swap
  assign(paste0(prefix, '.sense.mx'),mx1, envir = parent.frame())
  assign(paste0(prefix, '.antisense.mx'),mx2, envir = parent.frame())
}
save(Direction2Gene, file = 'r.func/Direction2Gene.rdata')

MxProfilePlot <- function(plot.dt, breaks = c(1,61,260,320), labels = c('-1.5kb','TSS','TES','1.5kb')){
  plot.dt <- cbind(plot.dt,bins = 1:nrow(plot.dt))
  plot.dt <- melt(plot.dt, id.vars = 'bins')
  p <- ggplot(plot.dt, aes(x=bins, y = value, group = variable, color = variable))
  p <- p + geom_line(size =1) +
    scale_x_continuous(breaks = breaks, labels = labels)+
    labs(x = NULL, y = 'Average Reads')+
    theme_xf+
    theme(legend.position = 'bottom')
  p  
}

save(MxProfilePlot, file = 'r.func/MxProfilePlot.rdata')


RefseqHeatmap <- function(mx,tsscol=61,tescol=260,color = colorRampPalette(brewer.pal(n=7,'YlOrRd'))(100),labelsize = 20,sort = c('auto','index'),inputindex = NULL, returnindex = F,pictitle,picfile,...){
  library(RColorBrewer)
  library(pheatmap)
  mx[is.na(mx)] <- 0
  mx[is.nan(mx)] <- 0
  if(sort == 'auto'){
    rownames(mx) <- 1:nrow(mx)
    mx.sum <- rowSums(mx[1:nrow(mx), (tsscol + 1):tescol])
    mx <- mx[order(mx.sum, decreasing = T), 1:ncol(mx)]
  } else if(sort == 'index'){
    mx <- mx[inputindex,1:ncol(mx)]
  } else {
    stop('sort information needed!')
  }
  if(returnindex == T){
    cat('Object:heatmapindex generated.')
    assign('heatmapindex',as.integer(rownames(mx)),envir = parent.frame())
  }
  mx[mx>quantile(mx,0.95)] <- quantile(mx,0.95)
  colnames(mx) <- c('-1.5kb',rep('',tsscol-2),'TSS',rep('',tescol-tsscol-1),'TES',rep('',ncol(mx)-tescol-1),'1.5kb')
  jpeg(filename = picfile,width = 1000,height = 3000,units = 'px', res = 150)
  pheatmap(mx, cluster_rows = F, cluster_cols = F,
           annotation_names_col = F,
           annotation_names_row = F,
           show_colnames = T,
           show_rownames = F,
           fontsize = 20,
           main = pictitle,
           color = color)
  dev.off()
}

save(RefseqHeatmap, file = 'r.func/RefseqHeatmap.rdata')



GRPie <- function(GRlist,labels = names(GRlist),title,...){
  pie(elementNROWS(GRlist),
      clockwise = T,
      labels = paste(elementNROWS(GRlist),labels),
      main = paste(sum(elementNROWS(GRlist)), title),
      ...)
}

save(GRPie,file = 'r.func/GRPie.rdata')

SplitbyRef <- function(gr){
  gr <- unique(gr)
  if(!exists('refseq.hg19.gr')){
    if(!exists('hum.annotation')){
      require(AnnotationHub)
      hum.annotation <- subset(AnnotationHub(), species == 'Homo sapiens')
      assign('hum.annotation', hum.annotation,envir = parent.frame())
      cat('hum.annotation')
    }
    temp <- query(hum.annotation,'RefSeq','hg19')
    refseq.hg19.gr <- temp[[1]] %>% 
      granges() %>% 
      unique() %>% 
      keepStandardChromosomes(pruning.mode = 'coarse')
    assign('refseq.hg19.gr', refseq.hg19.gr, envir = parent.frame())
  }
  promoters.hg19.gr <- promoters(refseq.hg19.gr)
  gr$type <- 'intergenic'
  gr$type[from(findOverlaps(gr,refseq.hg19.gr))] <- 'intragenic'
  strand(gr)[from(findOverlaps(gr,refseq.hg19.gr))] <- strand(refseq.hg19.gr)[to(findOverlaps(gr,refseq.hg19.gr))]
  gr$type[from(findOverlaps(gr,promoters.hg19.gr))] <- 'promoter'
  strand(gr)[from(findOverlaps(gr,refseq.hg19.gr))] <- strand(refseq.hg19.gr)[to(findOverlaps(gr,refseq.hg19.gr))]
  grlist <- split(gr,gr$type)
  grlist
} 

save(SplitbyRef, file = 'r.func/SplitbyRef.rdata')



FindCommonPeaks <- function(..., mincommon = 2){
  peaks.gr.list <- GenomicRangesList(...)
  all.peaks <- peaks.gr.list[[1]]
  n <- length(peaks.gr.list)
  for(i in 2:n){
    all.peaks <- union(all.peaks,peaks.gr.list[[i]])
  }
  mx <- matrix(0,nrow = length(all.peaks), ncol = n)
  for(i in seq(n)){
    index <- findOverlaps(all.peaks,peaks.gr.list[[i]]) %>% from()
    mx[index,i] <- 1
  }
  index <- rowSums(mx) >= mincommon
  all.peaks[index]
}
save(FindCommonPeaks,file = 'FindCommonPeaks.rdata')

Bed2GR <- function(file){
  gr <- fread(file) %>% 
    `colnames<-`(c('chr','start','end','name','value')) %>% 
    makeGRangesFromDataFrame() %>% keepStandardChromosomes(pruning.mode = "coarse")
  gr
}
save(Bed2GR, file = 'r.func/BED2GR.rdata')

ImportBiopacks <- function(package){
  source("https://bioconductor.org/biocLite.R")
  biocLite(package)
}

save(ImportBiopacks, file = 'r.func/ImportBiopacks.rdata')

GR2GTF <- function(gr,file){
  dt <- mcols(gr)
  attribute <- NULL
  if(ncol(dt) > 1){
    attribute <- paste(colnames(dt)[1], dt[,1])
  }
  if(ncol(dt) >= 2){
    for (i in 2:ncol(dt)){
      temp <- paste(colnames(dt)[i], dt[,i])
      attribute <- paste(attribute,temp,sep = '; ')
    }}
  df <- data.frame(seqnames=seqnames(gr),
                   source='Custom',
                   feature = 'exon',
                   starts=start(gr)-1,
                   ends=end(gr),
                   scores=c(rep(".", length(gr))),
                   strands=strand(gr),
                   frame = '0',
                   attributes = attribute)
  write.table(df, file = file, quote=F, sep="\t", row.names=F, col.names=F)
}
# Expression on refgene ---------------------------------------------------
refseq.hg19.3000p.gr <- refseq.hg19.gr[width(refseq.hg19.gr) > 3000]
refseq.hg19.3000p.gr <- GRGeneInfo(refseq.hg19.3000p.gr)
GR2GTF(refseq.hg19.3000p.gr,file = 'coordinate/refseq.hg19.3000p.gtf')

tss.hg19.gr <- flank(refseq.hg19.gr,1000,both = T)
GR2BED(tss.hg19.gr,file = 'coordinates/tss1000.hg19.bed')


# Pausing ratio determine -------------------------------------------------
pausing.res.rpt1 <- fread('pausing.count.rpt1.txt')
pausing.res.rpt2 <- fread('pausing.count.rpt2.txt')
pausing.res.3000p.rpt1 <- fread('pausing.count.3000p.rpt1.txt')
pausing.res.3000p.rpt2 <- fread('pausing.count.3000p.rpt2.txt')

pausing.cols <-
  c('gene_id',
    'chr',
    'start',
    'end',
    'strand',
    'length',
    'copies',
    'gene_name',
    'siCTL_ratio',
    'siCTL_promoter_reads',
    'siCTL_genebody_reads',
    'siP400_ratio',
    'siP400_promoter_reads',
    'siP400_genebody_reads',
    'siTIP60_ratio',
    'siTIP60_promoter_reads',
    'siTIP60_genebody_reads'
  )
colnames(pausing.res.rpt1) <- pausing.cols
colnames(pausing.res.rpt2) <- pausing.cols
colnames(pausing.res.3000p.rpt1) <- pausing.cols
colnames(pausing.res.3000p.rpt2) <- pausing.cols

save(PausingAnalysis, file = 'r.func/PausingAnalysis.rdata')

PausingAnalysis(pausing.res.rpt1)
PausingAnalysis(pausing.res.rpt2)
PausingAnalysis(pausing.res.3000p.rpt1)
PausingAnalysis(pausing.res.3000p.rpt2)

# Define the enhancer region ----------------------------------------------
# Find P300 peaks:
require(AnnotationHub)
hum.annotation <- subset(AnnotationHub(), species == 'Homo sapiens')

temp <- query(hum.annotation,'P300','Hela','hg19') 
temp <- temp[str_detect(mcols(temp)$title,'Hela')]

hela.p300.awg <- temp[[1]] %>% sort() %>% keepStandardChromosomes(pruning.mode = "coarse")
hela.p300.Sydh <- temp[[2]] %>% sort() %>% keepStandardChromosomes(pruning.mode = "coarse")
hela.p300.Kuznetsova <- Bed2GR('coordinate/hela.p300.peaks.Kuznetsova.hg19.bed')
hela.p300.peaks.co <- FindCommonPeaks(hela.p300.awg,hela.p300.Sydh,hela.p300.Kuznetsova)
hela.p300.peaks.union <- union(hela.p300.awg,hela.p300.Sydh) %>% union(hela.p300.Kuznetsova)
rm(hela.p300.awg,hela.p300.Sydh,hela.p300.Kuznetsova)

#Find k27ac peaks
hela.k27ac.encode <- Bed2GR('coordinate//hela.k27ac.peaks.encode.hg19.bed')
hela.k27ac.Kuznetsova <- Bed2GR('coordinate//hela.k27ac.peaks.Kuznetsova.hg19.bed')
hela.k27ac.Lai <- Bed2GR('coordinate//hela.k27ac.peaks.Lai.rpt2.hg19.bed')

hela.k27ac.peaks.co <- FindCommonPeaks(hela.k27ac.encode,
                                       hela.k27ac.Kuznetsova,
                                       hela.k27ac.Lai)
hela.k27ac.peaks.union <- union(hela.k27ac.encode,hela.k27ac.Kuznetsova) %>% union(hela.k27ac.Lai)
rm(hela.k27ac.encode,hela.k27ac.Kuznetsova,hela.k27ac.Lai)

#Find Overlap between p300 and k27ac
hela.enhancers.p300k27 <- hela.p300.peaks.union[from(findOverlaps(hela.p300.peaks.union, hela.k27ac.peaks.union))] %>% 
  reduce()

#Download data from fantom and enhancer atlas, pick ther intersections
hela.enhancer.atlas <- fread('coordinate//hela.enhancers.atlas') %>%
  `colnames<-`(c('chr','start','end','target')) %>% 
  makeGRangesFromDataFrame(keep.extra.columns = T) %>% keepStandardChromosomes()
human.enhancer.fantom <- fread('coordinate//human.enhancers.fantom.bed') %>% 
  `colnames<-`(c('chr','start','end','','','')) %>%
  makeGRangesFromDataFrame() %>% keepStandardChromosomes()
hela.enhancer.public <- intersect(hela.enhancer.atlas,human.enhancer.fantom)

#Union public data.
hela.enhancer <- union(hela.enhancers.p300k27,hela.enhancer.public) %>% reduce()
hela.enhancer.list <- SplitbyRef(hela.enhancer)
GRPie(hela.enhancer.list,title = 'co-peaks of P300 and H3k27ac in Hela')

hela.enhancers.intergenic <- hela.enhancer.list[['intergenic']] %>% keepStandardChromosomes()
hela.enhancers.intragenic <- hela.enhancer.list[['intragenic']] %>% keepStandardChromosomes()

GR2BED(hela.enhancers.intergenic, file = 'coordinate/hela.enhancer.intergenic.bed')
GR2BED(hela.enhancers.intragenic, file = 'coordinate/hela.enhancer.intragenic.bed')
#Test
txdb <- TxDb.Hsapiens.UCSC.hg19.knownGene

# Matrix around genes --------------------------------------------


#siCTL GRO on refseq
Direction2Gene('matrix/Hela-siCTL-rpt1-GROpos.ref.scaled.mx',
               'matrix/Hela-siCTL-rpt1-GROneg.ref.scaled.mx',
               refseq.hg19.gr,
               "siCTL.GRO.rpt1")
Direction2Gene('matrix/Hela-siCTL-rpt2-GROpos.ref.scaled.mx',
               'matrix/Hela-siCTL-rpt2-GROneg.ref.scaled.mx',
               refseq.hg19.gr,
               "siCTL.GRO.rpt2")

#siP400 GRO on refseq
Direction2Gene('matrix/Hela-siP400-rpt1-GROpos.ref.scaled.mx',
               'matrix/Hela-siP400-rpt1-GROneg.ref.scaled.mx',
               refseq.hg19.gr,
               "siP400.GRO.rpt1")
Direction2Gene('matrix/Hela-siP400-rpt2-GROpos.ref.scaled.mx',
               'matrix/Hela-siP400-rpt2-GROneg.ref.scaled.mx',
               refseq.hg19.gr,
               "siP400.GRO.rpt2")

#siTIP60 GRO on refseq
Direction2Gene('matrix/Hela-siTIP60-rpt1-GROpos.ref.scaled.mx',
               'matrix/Hela-siTIP60-rpt1-GROneg.ref.scaled.mx',
               refseq.hg19.gr,
               "siTIP60.GRO.rpt1")
Direction2Gene('matrix/Hela-siTIP60-rpt2-GROpos.ref.scaled.mx',
               'matrix/Hela-siTIP60-rpt2-GROneg.ref.scaled.mx',
               refseq.hg19.gr,
               "siTIP60.GRO.rpt2")

## plot profile

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.sense.mx)/length(refseq.hg19.gr),
                      # siP400 = colSums(siP400.GRO.rpt1.sense.mx)/length(refseq.hg19.gr),
                      siTIP60 = colSums(siTIP60.GRO.rpt1.sense.mx)/length(refseq.hg19.gr))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))


plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.antisense.mx)/length(refseq.hg19.gr),
                      # siP400 = colSums(siP400.GRO.rpt1.antisense.mx)/length(refseq.hg19.gr),
                      siTIP60 = colSums(siTIP60.GRO.rpt1.antisense.mx)/length(refseq.hg19.gr))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))


#ChIPseq signal
siCTL.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siCTL.ref.scaled.mx') 
siP400.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siP400.ref.scaled.mx')
siTIP60.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siTIP60.ref.scaled.mx')

plot.dt <- data.table(siCTL = colSums(siCTL.H3k27ac.mx)/length(refseq.hg19.gr),
                      siP400 = colSums(siP400.H3k27ac.mx)/length(refseq.hg19.gr),
                      siTIP60 = colSums(siTIP60.H3k27ac.mx)/length(refseq.hg19.gr))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))


siCTL.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siCTL.ref.scaled.mx') 
siP400.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siP400.ref.scaled.mx')
siTIP60.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siTIP60.ref.scaled.mx')

plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx)/length(refseq.hg19.gr),
                      siP400 = colSums(siP400.Pol2.mx)/length(refseq.hg19.gr),
                      siTIP60 = colSums(siTIP60.Pol2.mx)/length(refseq.hg19.gr))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

siCTL.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siCTL.ref.scaled.mx') 
siP400.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siP400.ref.scaled.mx')
siTIP60.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siTIP60.ref.scaled.mx')

plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx)/length(refseq.hg19.gr),
                      siP400 = colSums(siP400.Pol2.mx)/length(refseq.hg19.gr),
                      siTIP60 = colSums(siTIP60.Pol2.mx)/length(refseq.hg19.gr))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))


# Identify diffexp genes by DESeq2 ----------------------------------------
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

def.exp.tip60 <- DefExp('TIP60.countTable.txt', groups = c('siCTL','siCTL','siTIP60','siTIP60'), ref = 'siCTL')
def.exp.p400 <- DefExp('P400.countTable.txt', groups = c('siCTL','siCTL','siP400','siP400'), ref = 'siCTL')

temp <- select(org.Hs.eg.db, def.exp.tip60[[1]],'SYMBOL','REFSEQ')$SYMBOL
def.exp.tip60[[1]] <- temp

temp <- select(org.Hs.eg.db, def.exp.p400[[1]],'SYMBOL','REFSEQ')$SYMBOL
def.exp.p400[[1]] <- temp

colnames(def.exp.tip60)[1] <- 'gene'
colnames(def.exp.p400)[1] <- 'gene'
##########################stop here
DefExpMA(def.exp.tip60, title = 'Gene Expression - siTIP60')
DefExpMA(def.exp.p400, title = 'Gene Expression - siP400')

upgene.tip60 <- def.exp.tip60[log2FoldChange > 1 & padj < 0.05, 1]
downgene.tip60 <- def.exp.tip60[log2FoldChange < -1 & padj < 0.05, 1]

upgene.p400 <- def.exp.p400[log2FoldChange > 1 & padj < 0.05, 1]
downgene.p400 <- def.exp.p400[log2FoldChange < -1 & padj < 0.05, 1]


N1 <- nrow(upgene.p400)
N2 <- nrow(upgene.tip60)
N12 <- sum(downgene.p400$Transcript %in% downgene.tip60$Transcript)
SimpleVenn(N1,N2,N12,'up in siP400','up in siTIP40')


N1 <- nrow(downgene.p400)
N2 <- nrow(downgene.tip60)
N12 <- sum(downgene.p400$Transcript %in% downgene.tip60$Transcript)
SimpleVenn(N1,N2,N12,'down in siP400','down in siTIP40')

# Replot expression on up or down genes ----------------------------------------------
def.exp.tip60.gr <- makeGRangesFromDataFrame(def.exp.tip60, keep.extra.columns = T)
up.tip60.gr <- def.exp.tip60.gr[def.exp.tip60.gr$gene %in% upgene.tip60$gene]
down.tip60.gr <- def.exp.tip60.gr[def.exp.tip60.gr$gene %in% downgene.tip60$gene]

refup.index.tip60 <- refseq.hg19.gr$gene %in% up.tip60.gr$gene
refdown.index.tip60 <- refseq.hg19.gr$gene %in% down.tip60.gr$gene
#up:
plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.sense.mx[refup.index.tip60,])/sum(refup.index.tip60),
                      siTIP60 = colSums(siTIP60.GRO.rpt1.sense.mx[refup.index.tip60,])/sum(refup.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.antisense.mx[refup.index.tip60,])/sum(refup.index.tip60),
                      siTIP60 = colSums(siTIP60.GRO.rpt1.antisense.mx[refup.index.tip60,])/sum(refup.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
#down:
plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.sense.mx[refdown.index.tip60,])/sum(refdown.index.tip60),
                      siTIP60 = colSums(siTIP60.GRO.rpt1.sense.mx[refdown.index.tip60,])/sum(refdown.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.antisense.mx[refdown.index.tip60,])/sum(refdown.index.tip60),
                      siTIP60 = colSums(siTIP60.GRO.rpt1.antisense.mx[refdown.index.tip60,])/sum(refdown.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

#p400
def.exp.p400.gr <- makeGRangesFromDataFrame(def.exp.p400, keep.extra.columns = T)
up.p400.gr <- def.exp.p400.gr[def.exp.p400.gr$gene %in% upgene.p400$gene]
down.p400.gr <- def.exp.p400.gr[def.exp.p400.gr$gene %in% downgene.p400$gene]

refup.index.p400 <- refseq.hg19.gr$gene %in% up.p400.gr$gene
refdown.index.p400 <- refseq.hg19.gr$gene %in% down.p400.gr$gene
#up:
plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.sense.mx[refup.index.p400,])/sum(refup.index.p400),
                      siP400 = colSums(siP400.GRO.rpt1.sense.mx[refup.index.p400,])/sum(refup.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.antisense.mx[refup.index.p400,])/sum(refup.index.p400),
                      siP400 = colSums(siP400.GRO.rpt1.antisense.mx[refup.index.p400,])/sum(refup.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

#down:
plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.sense.mx[refdown.index.p400,])/sum(refdown.index.p400),
                      siP400 = colSums(siP400.GRO.rpt1.sense.mx[refdown.index.p400,])/sum(refdown.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.rpt1.antisense.mx[refdown.index.p400,])/sum(refdown.index.p400),
                      siP400 = colSums(siP400.GRO.rpt1.antisense.mx[refdown.index.p400,])/sum(refdown.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

# Replot ChIP signals on up or down genes ----------------------------------------------
# H3K27ac
plot.dt <- data.table(siCTL = colSums(siCTL.H3k27ac.mx[refup.index.p400,])/sum(refup.index.p400),
                      siP400 = colSums(siP400.H3k27ac.mx[refup.index.p400,])/sum(refup.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.H3k27ac.mx[refdown.index.p400,])/sum(refdown.index.p400),
                      siP400 = colSums(siP400.H3k27ac.mx[refdown.index.p400,])/sum(refdown.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.H3k27ac.mx[refup.index.tip60,])/sum(refup.index.tip60),
                      siTIP60 = colSums(siTIP60.H3k27ac.mx[refup.index.tip60,])/sum(refup.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))

plot.dt <- data.table(siCTL = colSums(siCTL.H3k27ac.mx[refdown.index.tip60,])/sum(refdown.index.tip60),
                      siTIP60 = colSums(siTIP60.H3k27ac.mx[refdown.index.tip60,])/sum(refdown.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))


# Pol2
plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refup.index.p400,])/sum(refup.index.p400),
                      siP400 = colSums(siP400.Pol2.mx[refup.index.p400,])/sum(refup.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refup_Pol2_p400.jpeg',height = 5, width = 6, units = 'in', dpi = 100)

plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refdown.index.p400,])/sum(refdown.index.p400),
                      siP400 = colSums(siP400.Pol2.mx[refdown.index.p400,])/sum(refdown.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refdown_Pol2_p400.jpeg',height = 5, width = 6, units = 'in', dpi = 100)


plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refup.index.tip60,])/sum(refup.index.tip60),
                      siTIP60 = colSums(siTIP60.Pol2.mx[refup.index.tip60,])/sum(refup.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refup_Pol2_tip60.jpeg',height = 5, width = 6, units = 'in', dpi = 100)

plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refdown.index.tip60,])/sum(refdown.index.tip60),
                      siTIP60 = colSums(siTIP60.Pol2.mx[refdown.index.tip60,])/sum(refdown.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refdown_Pol2_tip60.jpeg',height = 5, width = 6, units = 'in', dpi = 100)


# Pol2
plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refup.index.p400,])/sum(refup.index.p400),
                      siP400 = colSums(siP400.Pol2.mx[refup.index.p400,])/sum(refup.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refup_Pol2_p400.jpeg',height = 5, width = 6, units = 'in', dpi = 100)

plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refdown.index.p400,])/sum(refdown.index.p400),
                      siP400 = colSums(siP400.Pol2.mx[refdown.index.p400,])/sum(refdown.index.p400))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refdown_Pol2_p400.jpeg',height = 5, width = 6, units = 'in', dpi = 100)


plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refup.index.tip60,])/sum(refup.index.tip60),
                      siTIP60 = colSums(siTIP60.Pol2.mx[refup.index.tip60,])/sum(refup.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refup_Pol2_tip60.jpeg',height = 5, width = 6, units = 'in', dpi = 100)

plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx[refdown.index.tip60,])/sum(refdown.index.tip60),
                      siTIP60 = colSums(siTIP60.Pol2.mx[refdown.index.tip60,])/sum(refdown.index.tip60))
MxProfilePlot(plot.dt, breaks = c(1,31,131,160))
ggsave('figure/refdown_Pol2_tip60.jpeg',height = 5, width = 6, units = 'in', dpi = 100)

# Matrix around enhancers ------------------------------------------------------------
# enhancer.intergenic.pos.expression
siCTL.enh.inter.pos.rpt1.mx <- ReadMatrix('matrix/Hela-siCTL-rpt1-GROpos.enhancer.intergenic.mx')
siP400.enh.inter.pos.rpt1.mx <- ReadMatrix('matrix/Hela-siP400-rpt1-GROpos.enhancer.intergenic.mx')
siTIP60.enh.inter.pos.rpt1.mx <- ReadMatrix('matrix/Hela-siTIP60-rpt1-GROpos.enhancer.intergenic.mx')


plot.dt <- data.table(siCTL = colSums(siCTL.enh.inter.pos.rpt1.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.enh.inter.pos.rpt1.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.enh.inter.pos.rpt1.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))

# enhancer.intergenic.neg.expression
siCTL.enh.inter.neg.rpt1.mx <- ReadMatrix('matrix/Hela-siCTL-rpt1-GROneg.enhancer.intergenic.mx')
siP400.enh.inter.neg.rpt1.mx <- ReadMatrix('matrix/Hela-siP400-rpt1-GROneg.enhancer.intergenic.mx')
siTIP60.enh.inter.neg.rpt1.mx <- ReadMatrix('matrix/Hela-siTIP60-rpt1-GROneg.enhancer.intergenic.mx')


plot.dt <- data.table(siCTL = colSums(siCTL.enh.inter.neg.rpt1.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.enh.inter.neg.rpt1.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.enh.inter.neg.rpt1.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))


# enhancer intergenic H3K27ac 
siCTL.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siCTL.enhancer.intergenic.mx') 
siP400.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siP400.enhancer.intergenic.mx')
siTIP60.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siTIP60.enhancer.intergenic.mx')


plot.dt <- data.table(siCTL = colSums(siCTL.H3k27ac.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.H3k27ac.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.H3k27ac.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))


# Test --------------------------------------------------------------------

test <- fread('coordinate/gencode.v19.annotation.gtf')
library(TxDb.Hsapiens.UCSC.hg19.knownGene)

temp <- makeTxDbFromGFF('coordinate/gencode.v19.annotation.gtf',format = 'gtf')

# enhancer intergenic H3K4me1 
siCTL.H3K4me1.mx <- ReadMatrix('matrix/Hela-H3K4me1-SC-siCTL.enhancer.intergenic.mx') 
siP400.H3K4me1.mx <- ReadMatrix('matrix/Hela-H3K4me1-SC-siP400.enhancer.intergenic.mx')
siTIP60.H3K4me1.mx <- ReadMatrix('matrix/Hela-H3K4me1-SC-siTIP60.enhancer.intergenic.mx')


plot.dt <- data.table(siCTL = colSums(siCTL.H3K4me1.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.H3K4me1.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.H3K4me1.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))



# enhancer intergenic Pol2 
siCTL.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siCTL.enhancer.intergenic.mx') 
siP400.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siP400.enhancer.intergenic.mx')
siTIP60.Pol2.mx <- ReadMatrix('matrix/Hela-Pol2-SC-siTIP60.enhancer.intergenic.mx')


plot.dt <- data.table(siCTL = colSums(siCTL.Pol2.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.Pol2.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.Pol2.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))

Direction2Gene('matrix/Hela-siCTL-rpt1-GROpos.enhancer.intergenic.mx',
               'matrix/Hela-siCTL-rpt1-GROneg.enhancer.intergenic.mx',
               hela.enhancers.intergenic,
               "siCTL.GRO.enh.inter.rpt1")
# siP400
Direction2Gene('matrix/Hela-siP400-rpt1-GROpos.enhancer.intergenic.mx',
               'matrix/Hela-siP400-rpt1-GROneg.enhancer.intergenic.mx',
               hela.enhancers.intergenic,
               "siP400.GRO.enh.inter.rpt1")

# siTIP60
Direction2Gene('matrix/Hela-siTIP60-rpt1-GROpos.enhancer.intergenic.mx',
               'matrix/Hela-siTIP60-rpt1-GROneg.enhancer.intergenic.mx',
               hela.enhancers.intergenic,
               "siTIP60.GRO.enh.inter.rpt1")

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.enh.inter.rpt1.sense.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.GRO.enh.inter.rpt1.sense.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.GRO.enh.inter.rpt1.sense.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))

plot.dt <- data.table(siCTL = colSums(siCTL.GRO.enh.inter.rpt1.antisense.mx)/length(hela.enhancers.intergenic),
                      siP400 = colSums(siP400.GRO.enh.inter.rpt1.antisense.mx)/length(hela.enhancers.intergenic),
                      siTIP60 = colSums(siTIP60.GRO.enh.inter.rpt1.antisense.mx)/length(hela.enhancers.intergenic))
MxPeakcenterPlot(plot.dt, breaks = c(1,41,80), labels = c('-2kb','enh.center','2kb'))


siCTL.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27a') 
siP400.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siP400.ref.scaled.mx')
siTIP60.H3k27ac.mx <- ReadMatrix('matrix/Hela-H3K27ac-SC-siTIP60.ref.scaled.mx')