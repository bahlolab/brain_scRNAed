#Test site-specific differential editing

library(tidyverse)

here::here()

BL_joinKey <- readRDS('data/phs000834/BL_metadata_highQualCells.Rds')


################################################################
############### Test individual editing sites ##################
################################################################

gte5_DPsites <- read_tsv("data/samtools_depth_output/samDepth_BQSR_sitesGTE5.out", col_names = FALSE,
                         col_types = cols('X2'=col_character())) %>%
  magrittr::set_colnames(c('sample','siteID','depth')) %>%
  filter(sample %in% BL_joinKey$SRA) %>% 
  filter(siteID %in% dt_siteStats_TDjoin$siteID)

### editing detection per cell per site covered ###

depth_edStatus <- gte5_DPsites %>% left_join(dt_filt, by=c('siteID','sample')) %>% 
  select(siteID,sample, rdrType, altProp, totalDP, depth) %>% 
  mutate(status=ifelse(is.na(altProp),'nonEdit','edit')) %>% 
  left_join(BL_joinKey %>% select(SRA,neuType,Group_ID,area),by=c('sample'='SRA')) %>% distinct()

saveRDS(depth_edStatus, "data/stat_tests/binom_testing/depth_edStatus.Rds")



#multidplyr

#Run Fisher test on gross neuronal phenotype within cortical regions.
# This requires an area ~ neuType ed_V_nonEd matrix.

library(multidplyr)

depth_edStatus <- readRDS('data/stat_tests/binom_testing/depth_edStatus.Rds') %>% arrange(siteID)

cluster <- new_cluster(32)

gp_3vars <- depth_edStatus %>%  group_by(siteID,Group_ID,area) %>% partition(cluster)

countAP_gp3 <- gp_3vars %>% count(AP_isNA = is.na(altProp)) %>% collect()
#runs in ~10-45 seconds

countAP_gp3 %>% ungroup() %>% count(siteID,sort=TRUE) #
countAP_gp3 %>% filter(siteID=="1_100003296")
countAP_gp3 <- countAP_gp3 %>% ungroup()


countAP_gp3_neuType <- countAP_gp3 %>% 
  mutate(neuType=case_when(str_detect(Group_ID, 'Ex') ~ 'Ex',
                           str_detect(Group_ID, 'In') ~ 'In',
                           str_detect(Group_ID, 'No') ~ 'NoN',
                           TRUE ~ NA_character_))


countAP_area_neuType <- countAP_gp3_neuType %>% group_by(siteID,area,neuType,AP_isNA) %>% summarize(nSummed=sum(n))



#################################
######## for chiSQ tables #######
#################################

#neuType (collapsed / 'summed' Ex vs In neurons); create edCount table

area_neuType_edCount <- countAP_area_neuType %>% 
  select(siteID,area,neuType,AP_isNA,nSummed) %>% 
  filter(neuType!="NoN") %>% 
  distinct() %>% 
  mutate(AP_isNA = ifelse(AP_isNA==FALSE,'edit','nonEdit')) %>%  
  mutate(forSpd=paste(area,neuType,AP_isNA,sep="_")) %>% 
  select(siteID,forSpd,nSummed) %>% distinct() %>% 
  spread(forSpd,nSummed,fill=0) 

saveRDS(area_neuType_edCount, file="data/stat_tests/binom_testing/area_neuType_edCount.Rds")


######## ###### ###### ###### ###### ###### ########
###### TEST cortical region  (area) X neuType ######
######## ###### ###### ###### ###### ###### ########

chiSQ_area_nType <- function(a,b,C,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x){
  data <- matrix(c(a,b,C,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x), nrow = 12)
  chisq.test(data, simulate.p.value = TRUE,B=999)$p.value
}


cluster_copy(cluster, "chiSQ_area_nType")

for_area_chisq <- area_neuType_edCount %>% 
  group_by(siteID) %>% partition(cluster)

as.tibble(
  paste(  
    rep(paste0('BA',c(8,10,17,21,22,41)),each=4),
    rep(c('Ex','In'),each=2 ),
    rep(c('edit','nonEdit'),12 ),sep="_")) %>% head()


area_nType_chisq <- for_area_chisq %>%  
  mutate(p=chiSQ_area_nType(BA8_Ex_edit, BA8_Ex_nonEdit, BA8_In_edit, BA8_In_nonEdit, BA10_Ex_edit,BA10_Ex_nonEdit,BA10_In_edit,BA10_In_nonEdit,
                            BA17_Ex_edit,BA17_Ex_nonEdit,BA17_In_edit,BA17_In_nonEdit,BA21_Ex_edit,BA21_Ex_nonEdit,BA21_In_edit,BA21_In_nonEdit,
                            BA22_Ex_edit,BA22_Ex_nonEdit,BA22_In_edit,BA22_In_nonEdit,BA41_Ex_edit,BA41_Ex_nonEdit,BA41_In_edit,BA41_In_nonEdit)[[1]]) %>% 
  collect() %>% ungroup() %>% mutate(FDR = p.adjust(p,method="BH"))

area_nType_chisq %>% count(p>0) #11273 sites are appropriate for testing.

saveRDS(area_nType_chisq, file='data/stat_tests/binom_testing/area_nType_chisq.Rds')


######## ########## ######### ##########
######## ########## ######### ##########

#postHoc testing of ChiSq results to identifiy significant pair-wise differences.

#if(!require(rcompanion)){install.packages("rcompanion")}
library(tidyverse)
library(rcompanion)

library(foreach)
library(doParallel)  

detectCores()

cluster <- new_cluster(detectCores())

#NB no significant sites = no post-hoc tests

area_nType_chisq_FDR0.1 <- area_nType_chisq %>% select(siteID,p,FDR,everything()) %>% filter(FDR < 0.1) 

area_neuType_forPostHoc <- area_nType_chisq_FDR0.1 %>% filter(FDR<0.05) %>% select(siteID, contains('edit'),contains('nonEdit')) %>% 
  gather(key,value,-siteID) %>% separate(key,into=c('area','neuType','status'), sep="_") %>%
  spread(status,value) 

registerDoParallel(32)


foreach(i = unique(area_neuType_forPostHoc$siteID), .combine=rbind) %dopar% {
  a <- area_neuType_forPostHoc %>% filter(siteID==i) 
  m <-a %>% select(edit,nonEdit) %>% as.matrix()
  rownames(m) <- paste(a$area,a$neuType,sep="_")
  o <- rcompanion::pairwiseNominalIndependence(m, fisher=TRUE, gtest=FALSE, chisq=FALSE, digits=5) %>% as_tibble()
  o %>% mutate(siteID=i)
} -> area_neuType_FDR0.05_postHoc

area_neuType_FDR0.05_postHoc %>% count(p.adj.Fisher < 0.05)

area_neuType_FDR0.05_postHoc %>% filter(p.adj.Fisher < 0.05)
area_neuType_FDR0.05_postHoc %>% filter(p.adj.Fisher < 0.05) %>% count(Comparison,sort=TRUE)


################## ################## ################## ################## ##################
################## ################## ################## ################## ##################

#Write output

saveRDS(area_neuType_FDR0.05_postHoc, file='data/stat_tests/binom_testing/area_neuType_FDR0.05_postHoc.Rds')

