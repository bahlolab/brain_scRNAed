

library(tidyverse)

library(here)
here()



###### Read data ######

dt_filt     <-         readRDS("data/phs000834/dt_filt.Rds")
dt_siteStats_filt <-   readRDS("data/phs000834/dt_siteStats_filt.Rds")
dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")
##### ##### ##### ##### 



#Figure 1.

##  basic descriptors  ##
#1) sites per cell; 2) cells per site


#1) 
dt_filt %>% filter(siteID %in% (dt_siteStats_TDjoin$siteID)) %>% 
  count(sample) #3055 cells

F1a_table <- dt_filt %>% filter(siteID %in% (dt_siteStats_TDjoin$siteID)) %>% distinct() %>% 
  count(sample) 

F1b_table <- dt_siteStats_TDjoin %>% select(siteID, n_Cells_Lake) %>% distinct() 




#These charts are be limited to 83% of sites covered by a single non-OL gene region;

F1c_alt_table <- dt_siteStats_filt %>% 
  filter(str_detect(site_status,'single')) %>% 
  select(siteID,ENSGbioType,site_type,rdrType,genicFeature) %>% distinct() %>% 
  mutate(lineBound=case_when(site_type=='Alu' | site_type=='Alu_Novel' ~ 'A',
                             site_type=='nonRep' | site_type=='nonRep_Novel'~ 'B',
                             site_type=='rep_nonAlu' | site_type=='rep_nonAlu_Novel'~ 'C',
                             TRUE ~ 'NULL')) %>% 
  mutate(gene_type=ifelse(ENSGbioType=='protein_coding','Coding','Non-coding')) %>% 
  dplyr::rename(`RADAR class` = rdrType, `ENSG Bio-type` = ENSGbioType, `Genic feature` = genicFeature) 
  


F1d_table <- dt_siteStats_filt %>% filter(str_detect(site_status,'single')) %>%  
  select(siteID,ENSGbioType,site_type,rdrType,genicFeature) %>% distinct() %>% 
  mutate(rdrType=case_when(rdrType=='uncatalog' ~ 'Novel',
                           rdrType=='nonRep' ~ "Non-repetitive",
                           rdrType=='rep_nonAlu'~ 'Repetitive_non-Alu',
                           TRUE ~ rdrType)) %>% 
  mutate(rdrType=factor(rdrType,levels=c('Novel','Repetitive_non-Alu','Non-repetitive','Alu'))) %>% 
  mutate(tag=ifelse(str_detect(site_type,'Novel'),'Novel','Reported')) %>% 
  filter(ENSGbioType=="protein_coding") %>% #count(genicFeature)
  mutate(genicFeature=ifelse(genicFeature=='exon','exonic',genicFeature)) %>% 
  mutate(genicFeature=factor(genicFeature, levels = rev(c('five_prime_utr','exonic','intronic','stop_codon','three_prime_utr')))) %>% 
  filter(!is.na(genicFeature)) %>% 
  dplyr::rename(`RADAR class` = rdrType, `Genic feature` = genicFeature) 


###### ###### ###### ###### ######
###### ###### ###### ###### ######

#Write output

saveRDS(F1a_table, 'data/figure_data/F1a_table.Rds')
saveRDS(F1b_table, 'data/figure_data/F1b_table.Rds')
saveRDS(F1c_alt_table, 'data/figure_data/F1c_alt_table.Rds')
saveRDS(F1d_table, 'data/figure_data/F1d_table.Rds')
 
