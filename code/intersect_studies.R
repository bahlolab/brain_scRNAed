#p3. intersect putative editing sites with editing data from other studies

library(tidyverse)

library(here)
here()

source("code/dplyr_override.R")


##################
### Read files ###
##################


#Lake et al sites

dt_siteStats_filt <- readRDS("data/phs000834/dt_siteStats_filt.Rds")

#Darmanis et al GSE67835 gatk output

dt_siteStatsNeurons_GSE <- readRDS("data/GSE67835/dt_siteStatsNeurons_GSE.Rds")

dt_siteStatsNeurons_GSE %>% count(pxID,Group) %>% mutate(sum=sum(n))
dt_siteStatsNeurons_GSE %>% filter(n_Cells>=5) %>% count(pxID,Group) %>% mutate(sum=sum(n))

#Tan et al Brodmann Area 9 sites

Tan_FI_cortex <- readRDS("data/related_studies/Tan_etal_Nature/Tan_FI_cortexBA9_Hg38.Rds") 

Tan_FI_cortex %>% count(BA9>0)   #86K sites
Tan_FI_cortex %>% filter(BA9>0) %>% summary(.)
Tan_FI_cortex %>% filter(siteID=="4_157336723") 
Tan_FI_cortex %>% filter(region=="4", position==157336723) %>% head() 




################################################

#INTERSECT Darmanis etal and Lake etal neurons #

################################################


#with depth filter on Darmanis:
dt_siteStats_filt %>% select(siteID,basechange,rdrType) %>% distinct() %>% 
  left_join(dt_siteStatsNeurons_GSE %>% filter(n_Cells >= 3) %>% select(siteID,pxID), by='siteID', suffix=c('_Lake','_Darmanis')) %>% 
  count(pxID) 


#drop depth filter on Darmanis:
dt_siteStats_filt %>% select(siteID,basechange,rdrType) %>% distinct() %>% 
  left_join(dt_siteStatsNeurons_GSE %>% filter(n_Cells>=1) %>% select(siteID,pxID), by='siteID', suffix=c('_Lake','_Darmanis')) %>% 
  filter(!is.na(pxID)) %>% count(siteID)
#17667 phs000834 sites detected in at least one Darmanis etal individual.

dt_siteStats_filt %>% select(siteID,basechange,rdrType) %>% distinct() %>% 
  left_join(dt_siteStatsNeurons_GSE %>% select(siteID,pxID,Group,n_Cells) %>% distinct(), by=c('siteID')) %>% 
  filter(!is.na(Group)) %>% 
  count(Group,pxID,rdrType) %>% 
  rename(`RADAR class` = rdrType) %>% 
  ggplot(aes(x=pxID, y=n, fill=`RADAR class`)) + geom_col(col='black',lwd=0.2) + 
  geom_col(aes(x=pxID,y=n,fill=`RADAR class`), col='black',lwd=0.2) +
  xlab('Donor ID') + ylab('N. common sites') + coord_flip() + 
  scale_fill_brewer() +
  ggtitle('Darmanis_GSE67835_neurons', subtitle='Intersect with BlueLake neurons') 



#w depth filter on Darmanis:
dt_siteStats_filt %>% select(siteID,basechange,rdrType) %>% distinct() %>% 
  left_join(dt_siteStatsNeurons_GSE %>% select(siteID,pxID,Group,n_Cells) %>%
              filter(n_Cells >= 3) %>% distinct(), by=c('siteID')) %>% 
  filter(!is.na(Group)) %>% 
  count(Group,pxID,rdrType) %>% 
  rename(`RADAR class` = rdrType) %>% 
  ggplot(aes(x=pxID, y=n, fill=`RADAR class`)) + geom_col(col='black',lwd=0.2) + 
  geom_col(aes(x=pxID,y=n,fill=`RADAR class`), col='black',lwd=0.2) +
  xlab('Donor ID') + ylab('N. common sites') + coord_flip() + 
  scale_fill_brewer() +
  ggtitle('Darmanis_GSE67835_neurons', subtitle='Intersect with BlueLake neurons') 





##################################
##################################



### JOIN TAN ####

dt_siteStats_filt %>% select(siteID,rdrType,n_Cells) %>% distinct() %>% 
  left_join(Tan_FI_cortex %>% filter(BA9>0) %>% select(siteID,BA9), by='siteID') %>% 
  filter(!is.na(BA9))
#2060 sites of 40K
#NB distribution is totally different for bulk. majority of Tan sites have FI < 0.07.


Tan_FI_cortex %>% filter(BA9>0) %>% ggplot(aes(x=BA9)) + geom_histogram(fill='dodger blue', col='grey') + 
  ggtitle('Tan bulk brain RNAseq',subtitle = 'Brodmann area 9') +
  xlab('Alternate allele proportion')



#################################
######### UPSET PLOT ############
#################################



### Create object for UPSET plot ###

dt_siteStats_TDjoin <- dt_siteStats_filt %>% select(-c(source:ENSGbioType,ENSGID)) %>% distinct() %>% 
  left_join(dt_siteStatsNeurons_GSE, 
            by=c('siteID','basechange','rdrType'), suffix=c('_Lake','_Darmanis')) %>% 
  left_join(Tan_FI_cortex %>% select(siteID,BA9),by='siteID') %>% 
  rename(Tan_etal = BA9) 
#NB this is redundant on Darmanis pxIDs

dt_siteStats_TDjoin %>% select(siteID,basechange,site_type,pxID,Tan_etal) %>% distinct() %>% 
  count(site_type, Darm=!is.na(pxID), Tan=!is.na(Tan_etal), sort=TRUE) %>% mutate(sum=sum(n)) 


saveRDS(dt_siteStats_TDjoin, file="data/phs000834/dt_siteStats_TDjoin.Rds")


#### count site_type ####

dt_siteStats_TDjoin %>% 
  filter(str_detect(site_status,'single')) %>% 
  select(siteID,site_type,site_status) %>% distinct() %>% #count(siteID,sort=TRUE)
  count(site_type) %>% 
  mutate(sum=sum(n))



library(UpSetR)

forUPSET_p1 <- dt_siteStats_TDjoin %>% 
  mutate(dummy1=1) %>% spread(rdrType, dummy1) %>% 
  mutate(dummy2=1) %>% spread(pxID,dummy2) %>% 
  select(siteID,Alu,rep_nonAlu,nonRep,contains('AB_'),Tan_etal) %>% 
  distinct() %>% 
  mutate(Novel=ifelse(is.na(Alu) &is.na(rep_nonAlu) & is.na(nonRep) & 
                        is.na(AB_S11) & is.na(AB_S4) & is.na(AB_S5) & is.na(AB_S7) & is.na(Tan_etal), 1,NA)) 

forUPSET_p1 %>% dim() #46623 (redundant)
forUPSET_p1 %>% count(siteID, sort=TRUE)
forUPSET_p1 %>% filter(siteID=="1_100171601")

m <- as.matrix(forUPSET_p1 %>% select(-siteID))
rownames(m) <- forUPSET_p1$siteID
table(forUPSET_p1$Novel)
m[1:4, ]

m[!is.na(m)] <- 1
m[is.na(m)] <- 0


upset(as.data.frame(m), 
      sets=rev(c("Alu","rep_nonAlu","nonRep","Tan_etal","AB_S4", "AB_S5", "AB_S7","AB_S11","Novel")),
      keep.order = TRUE, nintersects = NA)

#collapse Darmanis etal

forUPSET_p1 %>% 
  mutate(Darmanis_etal=ifelse(!is.na(AB_S11)|!is.na(AB_S4)|!is.na(AB_S5)|!is.na(AB_S7),1,NA)) %>% 
  select(-contains('AB_')) %>% distinct() %>% dim() #40.9K

forUPSET_clpseDarm <- forUPSET_p1 %>% 
  mutate(Darmanis_etal=ifelse(!is.na(AB_S11)|!is.na(AB_S4)|!is.na(AB_S5)|!is.na(AB_S7),1,NA)) %>% 
  select(-contains('AB_')) %>% distinct()

#add intersect of novel sites in RepeatMasker database
forUPSET_clpseDarm_RM <- forUPSET_clpseDarm %>% 
  left_join(dt_siteStats_filt %>% 
              filter(str_detect(site_status,'single')) %>% 
              select(siteID,site_type) %>% mutate(dummy=1) %>% spread(site_type,dummy)) %>% 
  select(-Novel) 

#toggle off to revert to previous chart:
forUPSET_clpseDarm <- forUPSET_clpseDarm_RM


forUPSET_clpseDarm %>% group_by_at(vars(2:9,-Tan_etal)) %>% count(sort=TRUE) #
forUPSET_clpseDarm %>% count(Darmanis_etal) %>% mutate(prop=n/sum(n))
forUPSET_clpseDarm %>% mutate(Tan_etal_bin = ifelse(!is.na(Tan_etal), 1 , NA)) %>% 
   count(Tan_etal_bin) %>% mutate(prop=n/sum(n))


colSums(!is.na(forUPSET_clpseDarm[,-1] )) 

sum(colSums(!is.na(forUPSET_clpseDarm[,-1] )) )

m2 <- as.matrix(forUPSET_clpseDarm %>% select(-siteID))
m2[!is.na(m2)] <- 1
m2[ is.na(m2)] <- 0


m2_table <- as.data.frame(m2) %>% rename(`Repet_non-Alu` = rep_nonAlu, `Non-repet` = nonRep)
names(m2_table)


upset(m2_table, 
      sets=rev(c("Tan_etal","Darmanis_etal","Alu","Alu_Novel","Repet_non-Alu","rep_nonAlu_Novel","Non-repet","nonRep_Novel")),
      keep.order = TRUE, order.by = 'freq', nintersects = NA)
#save pdf charts/upset_studies.pdf 4x7 [--> Figure 2]



##### ##### ##### ##### ##### ##### ##### #####
##### ##### ##### ##### ##### ##### ##### #####


#How many sites are covered but not edited in these cells? [negative control for editing proportion]

#for samtools depth
write_tsv(dt_siteStats_TDjoin %>% select(siteID) %>% distinct() %>% 
  separate(siteID,into=c('region',"position"),sep="_", convert = TRUE) %>% 
  mutate(start=position-1) %>% select(1,3,2) ,
  col_names = FALSE,
  path='data/sam_depth_input/dt_siteStats_TDjoin_sites.bed')

