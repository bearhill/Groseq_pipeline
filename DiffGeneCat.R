
# Library -----------------------------------------------------------------
library(GenomicRanges)
library(magrittr)
library(ggplot2)
library(data.table)
library(RColorBrewer)
library(stringr)

# Functions ---------------------------------------------------------------
load('r.func/ImportBiopacks.rdata')

diff.file = 'P400.diffoutput.txt'
coordinate.file = 'coordinates/gencode.v19.annotation.geneonly.gtf'


# Change coordinate to GR -------------------------------------------------
coordinate <- fread(coordinate.file) %>% `colnames<-`(c('seqname','source','feature','start','end','score','strand','frame','attribute'))
coordinate <- makeGRangesFromDataFrame(coordinate,keep.extra.columns = T) %>% sort()
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
#Transform attribute type to data.table

coordinate$attribute <- NULL
mcols(coordinate) <- cbind(mcols(coordinate), attribute)
# Diff expression analysis (start from homer diff output) -----------------
dif.exp <- fread(diff.file)
Info <- colnames(dif.exp)[1]
dif.exp <- makeGRangesFromDataFrame(dif.exp[,-1],keep.extra.columns = T) %>% sort()
dif.tb <- cbind(as.data.table(dif.exp),mcols(coordinate))


test <- DefExp('P400.countTable.txt',groups = c('siCTL','siCTL','siP400','siP400'),ref = 'siCTL')
test <- makeGRangesFromDataFrame(test,keep.extra.columns = T) %>% sort()
