

library(tidyverse)

library(here)
here()

#### READ DATA ####

BL_joinKey <- readRDS('data/phs000834/BL_metadata_highQualCells.Rds')

dt_filt <- readRDS("data/phs000834/dt_filt.Rds") 

dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")

dt_siteStats_TDjoin %>% count(siteID,sort=TRUE) #NB siteID in this table is redundant on Darmanis et al pxID

#### PLOT ####

dt_filt %>% filter(siteID %in% dt_siteStats_TDjoin$siteID) %>% 
  mutate(rdrType=ifelse(rdrType=='drnd_Only','uncatalog', rdrType)) %>% 
  mutate(rdrType=case_when(rdrType=='uncatalog' ~ 'Novel',
                           rdrType=='nonRep' ~ "Non-repetitive",
                           rdrType=='rep_nonAlu'~ 'Repetitive_non-Alu',
                           TRUE ~ rdrType)) %>% 
  mutate(rdrType=factor(rdrType,levels=c('Novel','Repetitive_non-Alu','Non-repetitive','Alu'))) %>% 
  rename(`RADAR class` = rdrType) %>% 
  ggplot(aes(x=altProp)) + geom_histogram(col='black',lwd=0.2, aes(fill=`RADAR class`), position='dodge') +
  guides(fill = guide_legend(reverse=T)) + 
  scale_fill_brewer()  + xlab('Frequency of inosine (FI)') + ylab("N. sites")

#### SAVE ####

ggsave('charts/FI_distribution.pdf',width=6,height=3.5)


