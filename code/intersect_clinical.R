#Intersect differentially-edited sites with clinically relevant variants

library(tidyverse)
here::here()

dt_siteStats_filt <-   readRDS("data/phs000834/dt_siteStats_filt.Rds")
dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")

BL_joinKey <- readRDS('data/phs000834/BL_metadata_highQualCells.Rds')

anno <- readRDS('data/expn_processing/GRCh38_genes.Rds')

VEP_out <- readRDS('data/VEP_output/VEP_out.Rds')


source('code/dplyr_override.R')

#Existing_variation in VEP

VEP_out %>% dplyr::select(1,Existing_variation,Extra) %>% head(100)

VEP_out %>% count(Existing_variation, sort = TRUE) #COSM: cosmic (catalogue of somatic mutations in cancer)

VEP_out %>% filter(Existing_variation!="-") %>% 
  filter(str_detect(Existing_variation,'rs')) %>% select(1,Existing_variation,everything()) %>% 
  count(Existing_variation,sort=TRUE)

########## CLINVAR ###########

#Rare SNPs in OMIM/ClinVar databases 

clinVar <- read_tsv("data/bedtools_output/clinVar_CLNDN.out", col_names = FALSE) %>% 
  unite("siteID" , c(X1,X2)) %>% select(-X3,-X7) %>% select(siteID,everything()) %>% 
  dplyr::rename(clinVarID=X4, clinVar_ref=X5, clinVar_alt=X6) %>% 
  separate(X8, into=c('heading','data'),sep="=") %>% spread(heading,data) %>% select(-`<NA>`)


clinVar %>% count(CLNDN)


####################################
####################################


#NB clinVar was created from an unfiltered list of putative editing sites (dt_siteStats.bed);
# Not all clinVar sites survive the filtering steps that produce dt_siteStats_TDjoin

dt_siteStats_TDjoin %>% left_join(clinVar,by='siteID') %>% filter(!is.na(CLNDN)) %>% 
  mutate(Darm_detect = ifelse(!is.na(pxID),1,NA)) %>% 
  select(-c(n_Cells_Darmanis, pxID,clinVarID,clinVar_ref,clinVar_alt )) %>% 
  filter(!str_detect(CLNDN,'not_')) %>% 
  distinct() %>% 
  left_join(dt_siteStats_filt %>% select(siteID,ENSGID,gene_name, description )) %>% 
  select(siteID, gene_name,description, CLNDN,everything()) %>% head()

dt_siteStats_TDjoin %>% left_join(clinVar,by='siteID') %>% filter(!is.na(CLNDN)) %>% 
  mutate(Darm_detect = ifelse(!is.na(pxID),1,NA)) %>% 
  filter(!str_detect(CLNDN,'not_')) %>% 
  distinct() %>% 
  left_join(dt_siteStats_filt %>% select(siteID,ENSGID,gene_name, description )) %>% count(CLNDN)

#6 rare clinical variants have genetic SNPs that overlap editing sites.



########################################
############## AUTISM ##################
########################################


#Differentially edited sites in ASD, reported by Tran et al Nat Neuro 2018

Tran_hg19 <- read_tsv('data/related_studies/Tran_etal_NatNeuro/Tran_hg19_ST2b_2c_hg19.txt')
Tran_hg38 <- read_tsv('data/related_studies/Tran_etal_NatNeuro/Tran_hg38.bed',col_names = FALSE)

#FC: BA9; TC: BA41-42-22. cf Tran et al STable 1 [41593_2018_287_MOESM3_ESM.xlsx]
Tran_FC_TC_sites <- readRDS("data/related_studies/Tran_etal_NatNeuro/Tran_FC_TC_sites.Rds")

dt_siteStats_TDjoin %>% left_join(Tran_FC_TC_sites,by='siteID') %>% 
  count( Tran_FC, Tran_TC, sort=TRUE) %>% 
  filter(!is.na(Tran_FC)|!is.na(Tran_TC)) %>% mutate(sum=sum(n))



########################################
################ SCZ ###################
########################################

Breen_FC_sites <- readRDS("data/related_studies/Breen_etal_NatNeuro/Breen_FC_sites.Rds")

dt_siteStats_TDjoin %>% select(siteID,site_type,site_status) %>% distinct() %>% 
  left_join(Breen_FC_sites, by='siteID') %>% 
  group_by_at(vars(contains('Breen_'))) %>% count(sort=TRUE) %>% 
  filter(!is.na(Breen_CM_ACC) | !is.na(Breen_CM_DLPFC)|!is.na(Breen_HBCC_DLPFC)) %>% 
  ungroup() %>% mutate(sum=sum(n))

dt_siteStats_filt %>% 
  select(siteID,site_type,site_status,gene_name) %>% distinct() %>% 
  left_join(Breen_FC_sites, by='siteID') %>% 
  filter_at(vars(contains('Breen_')), any_vars(!is.na(.))) %>% 
  filter(str_detect(gene_name,'SNHG|SNOR')) %>% count(gene_name)


##############################################################


VEP_out %>% filter(siteID=="9_20926380") %>% mutate(xtraNST = str_split(Extra,';')) %>% select(1, xtraNST) %>% 
  unnest(xtraNST) %>% separate(xtraNST, into = c('key','value'), sep = "=") %>% distinct() %>% head()

#FOCAD
VEP_out %>% filter(siteID=="9_20926380") %>% mutate(xtraNST = str_split(Extra,';')) %>% select(1,xtraNST) %>% 
  unnest(xtraNST) %>% separate(xtraNST, into = c('key','value'), sep = "=") %>% distinct() 

#GRIA2
VEP_out %>% filter(siteID=="4_157336723") %>% mutate(xtraNST = str_split(Extra,';')) %>% select(1,xtraNST) %>% 
  unnest(xtraNST) %>% separate(xtraNST, into = c('key','value'), sep = "=") %>% distinct() 


#### JOIN ALL sites #####


clin_join <- dt_siteStats_TDjoin %>% left_join(clinVar, by='siteID') %>% 
  mutate(Darm_detect = ifelse(!is.na(pxID),1,NA)) %>% 
  select(-c(n_Cells_Darmanis, pxID,clinVarID,clinVar_ref,clinVar_alt )) %>% 
  rename(ClinVar = CLNDN) %>% 
  left_join(Breen_FC_sites,by='siteID') %>% 
  left_join(Tran_FC_TC_sites,by='siteID') 

clin_join %>% select(siteID, contains('Breen'), contains('Tran_')) %>% distinct() %>% 
  group_by_at(vars(contains('Breen'), contains('Tran_'))) %>% count(sort=TRUE) 

clin_join %>% group_by_at(24:28) %>% count(sort=TRUE)

clin_join %>% filter(siteID=="1_54849678") %>% head() #Desmosterolosis


clin_join %>% select(siteID, Breen_CM_ACC) %>% distinct() %>% count(!is.na(Breen_CM_ACC)) 
clin_join %>% select(siteID, Breen_CM_DLPFC) %>% distinct() %>% count(!is.na(Breen_CM_DLPFC)) 
clin_join %>% select(siteID,Breen_HBCC_DLPFC) %>% distinct() %>% count(!is.na(Breen_HBCC_DLPFC)) 


saveRDS(clin_join,'data/join_clinical.Rds')


##### check ClinVar cell types #####

clin_join %>% filter(!is.na(ClinVar)) ## join cell type (BL_joinKey? dt_stats?) ## --> edDepth


dt_filt <- readRDS("data/phs000834/dt_filt.Rds")


dt_filt %>% filter(siteID %in% (clin_join %>% filter(!is.na(ClinVar)) %>% pull(siteID))) %>% 
  left_join(BL_joinKey %>% select(SRA,neuType,Group_ID,area), by=c('sample'='SRA')) %>% 
  group_by(siteID) %>% count(neuType) %>% 
  ggplot(aes(x=siteID,y=n,fill=neuType)) + geom_col(position='fill') +
  coord_flip()

dt_filt %>% filter(siteID %in% (clin_join %>% filter(!is.na(ClinVar)) %>% pull(siteID))) %>% 
  left_join(BL_joinKey %>% select(SRA,neuType,Group_ID,area), by=c('sample'='SRA')) %>% 
  group_by(siteID) %>% count(area) %>% 
  ggplot(aes(x=siteID,y=n,fill=area)) + geom_col(position='fill') +
  coord_flip()


clinVar_celltypes <- dt_filt %>% select(siteID, sample) %>% filter(siteID %in% (clin_join %>% filter(!is.na(ClinVar)) %>% pull(siteID))) %>% 
  left_join(clin_join %>% select(siteID,ClinVar), by='siteID') %>% 
  left_join(BL_joinKey %>% select(SRA,neuType,Group_ID,area), by=c('sample'='SRA')) %>% distinct() %>% 
  group_by(siteID,ClinVar) %>% count(neuType,area) %>% 
  ungroup() %>% 
  left_join(clinVar %>% select(siteID,clinVarID),by="siteID") %>% 
  select(siteID,clinVarID,ClinVar,neuType,area,n)


saveRDS(clinVar_celltypes, 'data/clinVar_cellTypes.Rds')

