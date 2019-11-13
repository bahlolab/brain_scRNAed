#Model Global Editing Index vs gene expression (TPM)

library(here)
here()

#devtools::install_github("tidyverse/multidplyr")

library(multidplyr)
library(tidyverse)

BL_joinKey <- readRDS("data/phs000834/BL_metadata_highQualCells.Rds")

cell_edProp <- readRDS("data/stat_tests/GEI.Rds")


tpm_reshape <- readRDS("data/expn_processing/BlueLake_tpmSCE_reshape.Rds") #omits cells where count == 0 (i.e. log(gene value) == -Inf ).


bighead <- function(x){x[1:5,1:5]}


nCell_geneObs <- tpm_reshape %>% group_by(ENSG) %>% summarize(nObs = n())
nCell_geneObs %>% ggplot(aes(x=nObs)) + geom_histogram(fill='dodger blue') #N cells in which each gene is detected (tpm>0)

saveRDS(nCell_geneObs,'data/stat_tests/nCell_geneObs.Rds')


tpm_reshape %>% head(n=100000) %>% ggplot(aes(x=log10(value))) + geom_density(aes(col=key))


#generate list of genes only expressed in one neuType:

nType_single_ENSG <- tpm_reshape %>% 
  left_join(BL_joinKey, by=c('key'='SRA')) %>% 
  filter(neuType!="NoN") %>%
  count(ENSG, neuType) %>% #arrange(ENSG)
  group_by(ENSG) %>% mutate(nType=n()) %>% ungroup() %>%
  filter(nType==1) %>% select(ENSG) %>% distinct() %>% pull(ENSG)


#generate list of genes expressed in ≥2 neuronal sub-types, with logTPM > -Inf (i.e. TPM > 0)

nType_GTE1_ENSG <- tpm_reshape %>% 
  left_join(BL_joinKey,by=c('key'='SRA')) %>% 
  filter(neuType!="NoN") %>%
  count(ENSG,neuType) %>% #arrange(ENSG)
  group_by(ENSG) %>% mutate(nType=n()) %>% ungroup() %>%
  filter(nType>1) %>% select(ENSG) %>% distinct() %>% pull(ENSG)


#generate list of genes expresed in ≥2 cortical areas with logTPM > -Inf

area_GTE1_ENSG <- tpm_reshape %>% 
  left_join(BL_joinKey, by=c('key'='SRA')) %>% 
  filter(neuType!="NoN") %>%
  count(ENSG,area) %>% 
  group_by(ENSG) %>% mutate(nArea=n()) %>% ungroup() %>%
  filter(nArea>1) %>% select(ENSG) %>% distinct() %>% pull(ENSG)


#NB multidplyr can't take gather(); need to supply gathered, grouped df.


#### CREATE CLUSTER ####

cluster <- new_cluster(24)


#### SITES in single neuronal type

#Partition 
tpm_single_nType <- tpm_reshape %>% filter(ENSG %in% nType_single_ENSG) %>% group_by(ENSG) %>% partition(cluster)

#Collect
lm_out_single_neuType <- tpm_single_nType %>% 
  dplyr::left_join(BL_joinKey %>% select(SRA,Group_ID,area,neuType), copy=TRUE, by=c('key'='SRA')) %>% 
  left_join(edProp, copy=TRUE, by=c('key'='sample')) %>% 
  filter(Group_ID!="NoN") %>% dplyr::rename(edPropn=edProp) %>% 
  do(broom::tidy(lm(edPropn ~ log10(value) , data = . ))) %>% #NB adding area fails for genes only expressed in a single area
  filter(term!="(Intercept)") %>% arrange(p.value) %>% 
  collect()

lm_out_single_neuType %>% head()


###### Contrl for NEURONAL TYPE (In vs Ex) ######

#Partition
tpm_multi_nType <- tpm_reshape %>% filter(ENSG %in% nType_GTE1_ENSG) %>% group_by(ENSG) %>% partition(cluster)

#Collect
# include interaction term

lm_out_INT_neuType <- tpm_multi_nType %>% 
  dplyr::left_join(BL_joinKey %>% select(SRA,Group_ID,area,neuType), copy=TRUE, by=c('key'='SRA')) %>% 
  left_join(edProp, copy=TRUE, by=c('key'='sample')) %>% 
  filter(Group_ID!="NoN") %>% dplyr::rename(edPropn=edProp) %>% 
  do(broom::tidy(lm(edPropn ~ log10(value) * neuType , data = . ))) %>% 
  filter(term!="(Intercept)") %>% arrange(p.value) %>% 
  collect()



lm_out_INT_neuType %>% filter(estimate<0) %>% head()
lm_out_INT_neuType %>% filter(estimate>0) %>% head()


######## ######## ######## ######## ######## ########
######## ######## ######## ######## ######## ########

# Write output

saveRDS(lm_out_single_neuType,'data/stat_tests/broom_lm/lm_out_single_neuType.Rds')
saveRDS(lm_out_INT_neuType,'data/stat_tests/broom_lm/lm_out_INT_neuType.Rds')

