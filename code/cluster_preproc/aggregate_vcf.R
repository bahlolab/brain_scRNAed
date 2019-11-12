
# This file concatenates output for each sample, created in filter_vcf.R, into a single file.

library(tidyverse)

library(data.table)
library(tidyverse)
fileList <- list.files(pattern='DP.txt') 

l <- lapply(fileList, fread, sep="\t") #reads in 3K files. Takes time.
dt <- rbindlist( l , idcol = "listNum") %>% 
  mutate(chromosome=as.character(chromosome)) 

#join sampleNames
dt_label <- tibble(fileList) %>% rownames_to_column() %>% 
  mutate(sample = str_remove(fileList,'_processed.txt')) %>% 
  select(-fileList) %>% mutate(rowname=as.numeric(rowname)) %>% 
  left_join(dt,by=c("rowname"="listNum"))
  

saveRDS(dt_label, file='phs000834_aggregateSNPs_inclUncatalogd.Rds')



########### ########### ########### ########### ########### ###########
########### ########### ########### ########### ########### ###########
