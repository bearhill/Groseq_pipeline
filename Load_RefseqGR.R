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