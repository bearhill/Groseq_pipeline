GTF2GR <- function(file){
  require(data.table)
  require(GenomicRanges)
  coordinate <- fread(file) %>% `colnames<-`(c('seqname','source','feature','start','end','score','strand','frame','attribute'))
  attribute.type <- strsplit(coordinate$attribute[1],';') %>% unlist() %>% str_replace_all( '\\".*\\"','')
  attribute <- list(NULL)
  for (i in seq(length(attribute.type))) {
    name <- attribute.type[i]
    temp <- str_extract(coordinate$attribute, paste0(name, '\\".*?;')) %>%
      str_replace_all(name, '') %>%
      str_replace_all('\"', '') %>%
      str_replace_all(';', '') %>%
      str_replace_all(' ', '') 
    attribute[[i]] <- temp
  }
  attribute <- as.data.table(attribute) %>% `colnames<-`(str_replace_all(attribute.type,' ',''))
  coordinate$attribute <- NULL
  coordinate <- cbind(coordinate, attribute)
  coordinate <- makeGRangesFromDataFrame(coordinate,keep.extra.columns = T) %>% sort()
  coordinate
}
save(GTF2GR, file = 'r.func/GTF2GR.rdata')



