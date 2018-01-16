library(GenomicRanges)
library(magrittr)
library(data.table)

# Get eRNA from defined enhancer region --------------------------------------------------
load('r.func/BED2GR.rdata')

enhancer.intergen.gr <- Bed2GR('coordinates/hela.enhancer.intergenic.bed')
enhancer.intergen.gr$gene_id <- paste0('intergenic',seq(length(enhancer.intergen.gr)))
enhancer.intergen.gr$gene_status <- 'intergenic'

enhancer.intragen.gr <- Bed2GR('coordinates/hela.enhancer.intragenic.bed')
enhancer.intragen.gr$gene_id <- paste0('intragenic',seq(length(enhancer.intragen.gr)))

enhancer.intragen.gr$gene_status <- 'intragenic'
enhancer.gr <- c(enhancer.intergen.gr,enhancer.intragen.gr)

enhancer.center.gr <- resize(enhancer.gr,width = 1, fix = 'center')
erna.plus.gr <- flank(enhancer.center.gr,width = 1500, start = F)
strand(erna.plus.gr) <- '+'
erna.minus.gr <- flank(enhancer.center.gr,width = 1500)
strand(erna.minus.gr) <- '-'
erna.gr <- c(erna.plus.gr,erna.minus.gr)
erna.gr$gene_type <- 'enhancer'

temp <- findOverlaps(erna.gr, refseq.hg19.gr)
erna.gr <- erna.gr[-from(temp)]
erna.gr$transcript_id <- paste0('eRNA',seq(length(erna.gr)))

GR2GTF(erna.gr,'coordinates/Hela.eRNA.gtf')



