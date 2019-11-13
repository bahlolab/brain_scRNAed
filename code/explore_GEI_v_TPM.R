

#Explore Global Editing Index v gene expression (tpm)

library(here)
here()


### READ DATA ###

anno <- readRDS("data/expn_processing/GRCh38_genes.Rds")

tpm_reshape <- readRDS("data/expn_processing/BlueLake_tpmSCE_reshape.Rds")

dt_siteStats_filt <- readRDS('data/phs000834/dt_siteStats_filt.Rds') 

nCell_geneObs <- readRDS("data/stat_tests/nCell_geneObs.Rds")

edProp <- readRDS("data/stat_tests/GEI.Rds")

lm_out_INT_neuType   <-  readRDS('data/stat_tests/broom_lm/lm_out_INT_neuType.Rds') %>% 
  ungroup() %>% mutate(FDR=p.adjust(p.value,method="BH")) %>% arrange(FDR)



##### ##### #####

#compare ADAR expn

tpm_reshape %>% filter(ENSG %in% c("ENSG00000160710", "ENSG00000197381","ENSG00000185736")) %>% 
  left_join(BL_joinKey,by=c('key'='SRA')) %>% filter(neuType!="NoN") %>% 
  group_by(ENSG) %>% 
  do(broom::tidy(lm(log10(value) ~ neuType , data = .))) %>% filter(term!="(Intercept)")

tpm_reshape %>% filter(ENSG %in% c("ENSG00000160710", "ENSG00000197381","ENSG00000185736")) %>% 
  left_join(edProp,by=c('key'='sample')) %>% 
  left_join(BL_joinKey,by=c('key'='SRA')) %>% filter(neuType!="NoN") %>% 
  group_by(ENSG) %>% 
  do(broom::tidy(lm(edProp ~ log10(value) * neuType, data = .))) %>% 
  filter(term!="(Intercept)") %>% 
  filter(p.value<0.05,term!="neuTypeIn")


tpm_reshape %>% filter(ENSG %in% c("ENSG00000160710", "ENSG00000197381","ENSG00000185736")) %>% 
  left_join(BL_joinKey,by=c('key'='SRA')) %>% filter(neuType!="NoN") %>% 
  mutate(gene_name=case_when(ENSG=="ENSG00000160710" ~ 'ADAR1',
                             ENSG=="ENSG00000197381"~ 'ADAR2',
                             ENSG=="ENSG00000185736"~ 'ADAR3')) %>% 
  mutate(gene_name=factor(gene_name, levels=c('ADAR1','ADAR2','ADAR3'))) %>% 
  ggplot(aes(x=neuType,y=log10(value),col=neuType)) + geom_jitter(width=0.25, cex=0.35) +
  geom_boxplot(alpha=0, col='black') + facet_wrap(~gene_name) + theme(legend.position = 'NULL')

tpm_reshape %>% filter(ENSG %in% c("ENSG00000160710", "ENSG00000197381","ENSG00000185736")) %>% 
  left_join(edProp,by=c('key'='sample')) %>% 
  left_join(BL_joinKey,by=c('key'='SRA')) %>% filter(neuType!="NoN") %>% 
  mutate(gene_name=case_when(ENSG=="ENSG00000160710" ~ 'ADAR1',
                             ENSG=="ENSG00000197381"~ 'ADAR2',
                             ENSG=="ENSG00000185736"~ 'ADAR3')) %>% 
  mutate(gene_name=factor(gene_name, levels=c('ADAR1','ADAR2','ADAR3'))) %>% 
  ggplot(aes(x=log10(value),y=edProp,col=neuType)) + geom_point(cex=0.25, alpha=0.5) +
  facet_wrap(~ gene_name, ncol=1,scales='free_x') + geom_smooth(method="lm", se = FALSE)

#################

lm_out_INT_neuType %>% left_join(anno, by=c('ENSG' = 'gene_id' )) %>% 
  filter(str_detect(gene_name,'ADAR')) %>% 
  filter(term == "log10(value)", FDR<0.05)


lm_out_INT_neuType %>% dim() #114942


#Significant correlations with all transcripts, including directly edited transcripts 

nType_sig_cor_allXscript <- lm_out_INT_neuType %>% filter(FDR<0.05) %>% filter(term=='log10(value)') %>% 
  left_join(nCell_geneObs, by='ENSG') %>% 
  left_join(anno, by=c('ENSG' = 'gene_id' )) %>% 
  select(1:FDR,nObs,region,gene_name,gene_biotype,description) %>% 
  arrange(desc(estimate))

#Significant interactions between transcript expression and neuronal type
nType_sig_int_allXscript <- lm_out_INT_neuType %>% filter(FDR<0.05) %>% filter(str_detect(term,':')) %>% 
  left_join(nCell_geneObs, by='ENSG') %>% 
  left_join(anno, by=c('ENSG' = 'gene_id' )) %>% 
  select(1:FDR,nObs,region,gene_name,gene_biotype,description) %>% 
  arrange(desc(estimate))

#Drop transcripts that are directly edited

dt_siteStats_filt %>% filter(ENSGID %in% c("ENSG00000160710", "ENSG00000197381","ENSG00000185736")) %>% count(ENSGID)
#detect 28 edited sites in ADARB2 (ADAR3)

nType_sig_cor <- lm_out_INT_neuType %>% filter(FDR<0.05) %>% filter(term=='log10(value)') %>% 
  filter(! ENSG %in% dt_siteStats_filt$ENSGID) %>% 
  left_join(nCell_geneObs, by='ENSG') %>% 
  left_join(anno, by=c('ENSG' = 'gene_id' )) %>% 
  select(1:FDR,nObs,region,gene_name,gene_biotype,description) %>% 
  arrange(desc(estimate))


nType_sig_int <- lm_out_INT_neuType %>% filter(FDR<0.05) %>% filter(term=='log10(value):neuTypeIn') %>% 
  filter(! ENSG %in% dt_siteStats_filt$ENSGID) %>% 
  left_join(nCell_geneObs,by='ENSG') %>% 
  left_join(anno, by=c('ENSG' = 'gene_id' )) %>% 
  select(1:FDR,nObs,region,gene_name,gene_biotype,description) %>% 
  arrange(desc(estimate))

########## ########## ##########

nType_sig_cor_allXscript %>% count(estimate<0) 

#NB ~6x more negative than positive correlations

########## ########## ##########


nType_sig_cor %>% count(estimate<0) #many more negative than positive correlations
nType_sig_cor %>% ggplot(aes(x=estimate)) + geom_histogram(bins=50, fill='dodger blue',col='grey')
nType_sig_cor %>% filter(abs(estimate) > 0.01) %>% count(estimate>0)


nType_sig_int %>% filter(abs(estimate) > 0.01) %>% count(estimate>0) %>% mutate(sum=sum(n))



##### GSEA for editing-correlated genes ######

library(limma)

entz_posCorfor_GSEA <- nType_sig_cor %>% 
  filter(estimate>0.01) %>% 
  left_join(anno, by=c('ENSG'='gene_id')) %>% dplyr::select(ENSG, entrezid) %>% 
  dplyr::select(ENSG,entrezid) %>% 
  tidyr::unnest(entrezid) %>%  #requires tidyr >= 0.8.9
  group_by(ENSG) %>% dplyr::slice(n=1) %>% ungroup()

entz_negCorfor_GSEA <- nType_sig_cor %>% 
  filter(estimate< (-0.01)) %>% 
  left_join(anno, by=c('ENSG'='gene_id')) %>% dplyr::select(ENSG, entrezid) %>% 
  dplyr::select(ENSG,entrezid) %>% 
  tidyr::unnest(entrezid) %>%  #requires tidyr >= 0.8.9
  group_by(ENSG) %>% dplyr::slice(n=1) %>% ungroup()


entz_negCorfor_GSEA %>% dim()
entz_negCorfor_GSEA %>% count(ENSG) #842 when abs(est.) > 0.01


lm_ctrlNeuType_GOANAposCor <- limma::goana(as.character(entz_posCorfor_GSEA$entrezid)) %>% rownames_to_column(var='term')
lm_ctrlNeuType_GOANAnegCor <- limma::goana(as.character(entz_negCorfor_GSEA$entrezid)) %>% rownames_to_column(var='term')


lm_ctrlNeuType_GOANAposCor %>% filter(Ont!="CC") %>% 
  arrange(P.DE) %>% filter(DE>10, P.DE < 0.05) %>% head()

lm_ctrlNeuType_GOANAnegCor %>% filter(Ont!="CC") %>% 
  arrange(P.DE) %>% filter(DE>10, P.DE < 0.05) %>% head()


############ ############ ############ ############ ############ ############
############ ############ ############ ############ ############ ############

#Annotate enriched GO terms with constituent correlated genes
#A T Lun solution

DB <- paste("org", "Hs", "eg", "db", sep = ".")
require(DB, character.only = TRUE)
GO2ALLEGS <- paste("org", "Hs", "egGO2ALLEGS", sep = ".")
EG.GO <- AnnotationDbi::toTable(get(GO2ALLEGS))
d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])
EG.GO <- EG.GO[!d, ]


source('code/dplyr_override.R')

############### ############### ############### ################
############### POSITIVE CORRELATION with EdProp ###############
############### ############### ############### ################

nType_sig_cor %>% 
  filter(estimate> (0.01)) %>% 
  left_join(anno, by=c('ENSG'='gene_id')) %>% head()

nType_sig_cor %>% 
  filter(estimate> (0.01)) %>% count(gene_biotype, sort=TRUE)

lm_ctrlNeuType_GOANAposCor %>% filter(DE>10, Ont!="CC") %>% filter(P.DE<0.05) %>%  arrange(P.DE) %>% head()



#GO RNA processing
entz_posCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:0006396') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% #count(region, sort=TRUE)
  count(gene_biotype)


#GO organic cyclic compound
entz_posCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:1901360') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% 
  head()


# GO:0010467 'gene expression'
entz_posCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:0010467') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% head()
  
nType_sig_cor %>% filter(estimate>0) %>% filter(str_detect(description,'RNA binding'))



#Output tables

GOANA_pos_annot <- lm_ctrlNeuType_GOANAposCor %>% filter(DE>10, Ont!="CC") %>% filter(P.DE<0.05) %>%  arrange(P.DE) %>% 
  left_join(EG.GO %>% filter(gene_id %in% as.character(entz_posCorfor_GSEA$entrezid)),
            by=c('term'='go_id')) %>% 
  left_join(entz_posCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)), 
            by=c('gene_id' = 'entrezid')) %>% 
  left_join(anno %>% dplyr::select(gene_id,gene_name,description), by=c('ENSG'='gene_id')) 


GOANA_neg_annot <- lm_ctrlNeuType_GOANAnegCor %>% filter(DE>10, Ont!="CC") %>% filter(P.DE<0.05) %>%  arrange(P.DE) %>% 
  left_join(EG.GO %>% filter(gene_id %in% as.character(entz_negCorfor_GSEA$entrezid)),
            by=c('term'='go_id')) %>% 
  left_join(entz_negCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)), 
            by=c('gene_id' = 'entrezid')) %>% 
  left_join(anno %>% dplyr::select(gene_id,gene_name,description), by=c('ENSG'='gene_id')) 

saveRDS(GOANA_pos_annot, 'data/stat_tests/broom_lm/GOANA_pos_annot.Rds')
saveRDS(GOANA_neg_annot, 'data/stat_tests/broom_lm/GOANA_neg_annot.Rds')




############### ############### ############### ################
############### NEGATIVE CORRELATION with EdProp ###############
############### ############### ############### ################

lm_ctrlNeuType_GOANAnegCor %>% filter(DE>10, Ont!="CC", P.DE<0.05) %>% arrange(P.DE) %>% head()

#GO:2001141 regulation of biosynthesis
entz_negCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:2001141') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% #count(region, sort=TRUE)
  head()


#DNA binding TF
entz_negCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:0043565') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% #count(region, sort=TRUE)
  head()

#signaling receptor activity GO:0038023 
entz_negCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:0038023') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% #count(region, sort=TRUE)
  head()

#neuron differentiation GO:0030182
entz_negCorfor_GSEA %>% mutate(entrezid=as.character(entrezid)) %>% left_join(EG.GO,by=c('entrezid'='gene_id')) %>% 
  filter(go_id=='GO:0030182') %>% 
  left_join(nType_sig_cor,by='ENSG') %>% #count(region, sort=TRUE)
  head()





###############################################################  
######## Related studies of RNA-editing correlations ##########
###############################################################  


#### Tan etal Nature 2017 correlates ####

Tan_posCor <- readxl::read_xlsx("data/related_studies/Tan_etal_Nature/Supplementary File 6.xlsx", 
                                sheet=1, col_names = TRUE) %>% 
  separate(EnsemblId,into=c('ENSGpref','ENSGsuff'))


Tan_negCor <- readxl::read_xlsx("data/related_studies/Tan_etal_Nature/Supplementary File 6.xlsx", 
                                sheet=2, col_names = TRUE) %>% 
  separate(EnsemblId,into=c('ENSGpref','ENSGsuff'))

### ### ### ### ### ###

Tan_posCor %>% mutate(tag='Tan_posCor') %>% select(1,2,Estimate,tag) %>% distinct() %>% 
  left_join(lm_out_INT_neuType %>% filter(term=='log10(value)'), by = c('ENSGpref'='ENSG')) %>% distinct() %>% 
  filter(FDR<0.05) %>% 
  count(Estimate>0,estimate>0)

edCor_pos_Tan <- Tan_posCor %>% mutate(tag='Tan_posCor') %>% select(1,2,Estimate,tag) %>% distinct() %>% 
  left_join(lm_out_INT_neuType %>% filter(term=='log10(value)'), by = c('ENSGpref'='ENSG')) %>% distinct() %>% 
  filter(FDR<0.05) %>% 
  filter(Estimate>0, estimate>0) %>% 
  left_join(anno %>% select(gene_id,description),by=c('ENSGpref'='gene_id'))

### ###

Tan_negCor %>% mutate(tag='Tan_negCor') %>% select(1,2,Estimate,tag) %>% distinct() %>% 
  left_join(lm_out_INT_neuType %>% filter(term=='log10(value)'), by = c('ENSGpref'='ENSG')) %>% distinct() %>% 
  filter(FDR<0.05) %>% 
  filter(Estimate<0,estimate<0) %>% 
  left_join(anno %>% select(gene_id,description),by=c('ENSGpref'='gene_id')) 

edCor_neg_Tan <- Tan_negCor %>% mutate(tag='Tan_negCor') %>% select(1,2,Estimate,tag) %>% distinct() %>% 
  left_join(lm_out_INT_neuType %>% filter(term=='log10(value)'), by = c('ENSGpref'='ENSG')) %>% distinct() %>% 
  filter(FDR<0.05) %>% 
  filter(Estimate<0, estimate<0) %>% 
  left_join(anno %>% select(gene_id,description),by=c('ENSGpref'='gene_id')) 


saveRDS(rbind(edCor_pos_Tan, edCor_neg_Tan), 'data/stat_tests/edCor_v_Tan.Rds')

#### Tan correlations v genetype

Tan_posCor %>% left_join(dt_siteStats_filt,by=c("ENSGpref"="ENSGID")) %>% count(ENSGbioType)
Tan_negCor %>% left_join(dt_siteStats_filt,by=c("ENSGpref"="ENSGID")) %>% count(ENSGbioType)


##### Ed-associated RNA-binding proteins #####

#[cf code/edAssoc_RBPs.R] Quinones-Valdez etal Nat Comms 2019

edAssoc_RBPs <- readRDS('data/related_studies/QuinonesValdez_etal_NatComms/edAssoc_RBPs.Rds')

RBP_QV_intersect <- lm_out_INT_neuType %>% filter(term=='log10(value)') %>% 
  left_join(anno %>% select(gene_id,gene_name,description),by=c('ENSG'='gene_id')) %>% 
  filter(FDR<0.05) %>% 
  left_join(edAssoc_RBPs,by=c('gene_name' = 'RBP')) %>% 
  filter(!is.na(tag)) %>% 
  mutate(concordance = case_when(estimate<0 & tag == 'inhibitor' ~ 'YES',
                                 estimate>0 & tag == 'activator' ~ 'YES',
                                 TRUE ~ 'NO' )) 

saveRDS(RBP_QV_intersect, 'data/stat_tests/edCor_v_QV.Rds')  

##################################

#Hudson & Ortlund DNA & RNA binding proteins

DRBPs <- read_tsv('data/related_studies/Hudson_Ortlund_NRMolCellBiol_2014/DNA_RNA_bindingProteins.txt',col_names = FALSE) %>% 
  magrittr::set_colnames('NA_BP') %>% 
  separate(NA_BP,into=c("gene_name",'description'),sep="\\(") %>% 
  mutate(description=str_remove_all(description,"\\)")) %>% 
  mutate(gene_name=str_remove_all(gene_name,' '))

RBP_HO_intersect <- lm_out_INT_neuType %>%  filter(term=='log10(value)') %>% 
  left_join(anno %>% select(gene_id,gene_name,description),by=c('ENSG'='gene_id')) %>% 
  #filter(FDR<0.05,abs(estimate>=0.01)) %>% 
  filter(FDR<0.05) %>% 
  filter(gene_name %in% DRBPs$gene_name) 

saveRDS(RBP_HO_intersect, 'data/stat_tests/edCor_v_HO.Rds')  

##################################

#Known correlates of RNA editing (interacting proteins)
lm_out_INT_neuType %>% left_join(anno %>% select(gene_id,gene_name, description), by=c('ENSG'='gene_id')) %>% 
  filter(term=='log10(value)', FDR<0.05) %>% 
  filter(str_detect(gene_name,'ADAR|FTO|METTL3|FTO|DROSHA|TARDBP|MATR3|XRCC6|ILF2|PABPC1|FXR1|AIMP2|FMRP|PIN1')) 




