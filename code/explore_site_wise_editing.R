#explore site-wise editing

library(tidyverse)

here::here()

source('code/dplyr_override.R')

BL_joinKey <- readRDS("BL_metadata_highQualCells.Rds")
dt_siteStats_filt <- readRDS("data/phs000834/dt_siteStats_filt.Rds")
dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")

anno <- readRDS('data/expn_processing/GRCh38_genes.Rds')


#created in test_site_wise_editing.R

areaXnType_edCount <- readRDS("data/stat_tests/binom_testing/area_neuType_edCount.Rds")
areaXnType_chisq <- readRDS("data/stat_tests/binom_testing/area_nType_chisq.Rds") #all sites with FDR < 0.1
areaXnType_FDR0.05_postHoc <- readRDS('data/stat_tests/binom_testing/area_neuType_FDR0.05_postHoc.Rds') #only sites with FDR < 0.05 were tested in post-hoc.

#Cluster marker genes
#STable S5

BL_markers <- readxl::read_xlsx("data/phs000834/aaf1204_Lake_SM_Tables-S1-S16_BA.xlsx",
                                col_names = TRUE, sheet=5, skip=5) %>% 
  filter(str_detect(cluster,'In|Ex')) %>% 
  left_join(anno,by=c('Gene'='gene_name')) %>% 
  filter(str_detect(gene_id,'ENSG'), !str_detect(region,'CHR')) 

BL_markers %>% count(cluster)

####################################

areaXnType_chisq %>% arrange(desc(FDR)) %>% filter(FDR<0.05) %>% count(siteID)

areaXnType_FDR0.05_postHoc %>% filter(p.adj.Fisher < 0.1) %>% count(siteID) # more lenient given larger n. groups (= smaller cell numbers)

areaXnType_FDR0.05_postHoc %>% filter(p.adj.Fisher<0.1) %>% count(siteID) %>% 
  left_join(dt_siteStats_filt %>% select(siteID,site_status,site_type,ENSGbioType) %>% distinct(),by='siteID') %>% 
  filter(str_detect(site_status,'single')) %>% 
  count(site_type,ENSGbioType,sort=TRUE)



#### Transform pair-wise output ####

areaXnType_ph_tbl <- areaXnType_FDR0.05_postHoc %>% mutate(comp = Comparison) %>% 
  separate(comp, into=c('group1','group2'),  sep=" : ") %>% 
  separate(group1, into=c('area1','nType1'), sep="_") %>% 
  separate(group2, into=c('area2','nType2'), sep="_") %>% 
  select(-Comparison) %>% 
  select(siteID,area1,nType1,area2,nType2,everything())

saveRDS(areaXnType_ph_tbl, 'data/stat_tests/binom_testing/areaXnType_ph_tbl.Rds')

areaXnType_ph_tbl %>% filter(p.adj.Fisher<0.1) %>% arrange(p.adj.Fisher) %>% count(siteID) 
areaXnType_ph_tbl %>% filter(p.adj.Fisher<0.1) %>% arrange(p.adj.Fisher) %>% 
  left_join(dt_siteStats_filt,by='siteID') %>% filter(str_detect(site_status,'single')) %>% count(siteID) 

#210 are in non-OL regions


#Check for presence of cluster marker genes in differentially edited (DEd) genes

#All sites incl. multi_cognate
areaXnType_dEd_sites <- areaXnType_ph_tbl %>% filter(p.adj.Fisher<0.1) %>% arrange(p.adj.Fisher) %>% 
  left_join(dt_siteStats_filt, by='siteID') %>% #filter(str_detect(site_status,'single')) %>% 
  left_join(BL_markers %>% select(cluster,gene_id), by=c('ENSGID'='gene_id')) 

areaXnType_dEd_sites %>% select(siteID,site_status) %>% distinct() %>% count(site_status)

#only non-OL sites (single_cognate)
areaXnType_dEd_nonOLsites <- areaXnType_ph_tbl %>% filter(p.adj.Fisher<0.1) %>% arrange(p.adj.Fisher) %>% 
  left_join(dt_siteStats_filt, by='siteID') %>% filter(str_detect(site_status,'single')) %>% 
  left_join(BL_markers %>% select(cluster, gene_id), by=c('ENSGID'='gene_id')) 




###################################
###################################


areaXnType_ph_tbl %>% count(p.adj.Fisher<0.05)

spread_p1 <- areaXnType_edCount %>% gather(key,value,-siteID) %>% 
  separate(key,into=c('area','nType','status'),sep="_") %>% 
  spread(status,value) 

areaXnType_ph_tbl_edCount <- areaXnType_ph_tbl %>% filter(p.adj.Fisher < 0.1) %>% 
  left_join(spread_p1,by=c('siteID','area1'='area','nType1'='nType')) %>% 
  rename(at1_edit=edit,at1_nonEdit=nonEdit) %>% 
  left_join(spread_p1,by=c('siteID','area2'='area','nType2'='nType')) %>% 
  rename(at2_edit=edit, at2_nonEdit=nonEdit) 


areaXnType_ph_tbl_edCount %>% count(area1,nType1,area2,nType2,sort=TRUE)
#Most commonly diff ed. groups are BA41 vs BA8 


#Tag cortical lobes and gross neuronal phenotypes for DEd results; 

areaXnType_diffType <- areaXnType_ph_tbl_edCount %>% arrange(p.adj.Fisher) %>% 
  filter(p.adj.Fisher < 0.05) %>% 
  left_join(dt_siteStats_filt %>% select(siteID,ENSGID), by='siteID') %>% 
  left_join(anno %>% select(gene_id,gene_name,description), by=c('ENSGID'='gene_id')) %>% 
  mutate(diffType=case_when(nType1=="Ex" & nType2=="Ex" ~ 'Ex_Ex',
                            nType1=="In" & nType2=="Ex" ~ 'Ex_In',
                            nType1=="Ex" & nType2=="In" ~ 'Ex_In',
                            nType1=="In" & nType2=="In" ~ 'In_In',
                            TRUE ~ 'other')) %>% 
  mutate(diffArea=case_when(str_detect(area1, 'BA10|BA8') & str_detect(area2,'BA21|BA22|BA41') ~ 'FC_TC',
                            str_detect(area2, 'BA10|BA8') & str_detect(area1,'BA21|BA22|BA41') ~ 'FC_TC',
                            str_detect(area1, 'BA10|BA8') & str_detect(area2,'BA17') ~ 'FC_VC',
                            str_detect(area2, 'BA10|BA8') & str_detect(area1,'BA17') ~ 'FC_VC',
                            str_detect(area1, 'BA21|BA22|BA41') & str_detect(area2,'BA17') ~ 'TC_VC',
                            str_detect(area2, 'BA21|BA22|BA41') & str_detect(area1,'BA17') ~ 'TC_VC',
                            str_detect(area1, 'BA21|BA22|BA41') & str_detect(area2,'BA21|BA22|BA41') ~ 'TC_TC',
                            str_detect(area1, 'BA10|BA8') & str_detect(area2,'BA10|BA8') ~ 'FC_FC',
                            area1==area2 ~ 'within_area',
                            TRUE ~ 'other')) 

saveRDS(areaXnType_diffType, 'data/stat_tests/areaXnType_diffType.Rds') #adj P  < 0.1

#--> [Figure5a_plot.R]

areaXnType_GOenrich.Rds

######### ######### ######### ######### ######### ######### #########

######### ######### ######### ######### ######### ######### #########


#Plot draft

Figure5b_table <- areaXnType_diffType %>% filter(description!="NULL") %>% 
  filter(p.adj.Fisher < 0.05) %>% 
  arrange(p.Fisher) %>% rownames_to_column(var='pair') %>% mutate(pair=as.numeric(pair)) %>% 
  select(pair, everything(), -contains('Fisher'), -ENSGID, -description, -contains('diff')) %>% 
  mutate(chr=str_extract(siteID,"[^_]+")) %>% 
  mutate(gp1=paste(area1,nType1,sep="_"),
         gp2=paste(area2,nType2,sep="_")) %>% 
  mutate(at1_prop=at1_edit/(at1_edit + at1_nonEdit),
         at2_prop=at2_edit/(at2_edit + at2_nonEdit)) %>% 
  select(pair,chr,gene_name,siteID,gp1,gp2,at1_prop,at2_prop) %>% 
  gather(key,value,-c(pair,chr,gene_name,siteID,gp1,gp2)) %>% arrange(pair) %>% 
  mutate(myKey=ifelse(key=='at1_prop', gp1,gp2)) %>% select(-contains('gp'),-key) %>% 
  mutate(chrSplit=ifelse(chr %in% 1:10,1,2)) %>% 
  mutate(site_name=paste(gene_name,siteID)) %>% 
  separate(myKey,into=c('area','nType'),sep="_") %>%
  mutate(area=factor(area,levels=paste0('BA',c(8,10,21,22,41,17)))) 

saveRDS(Figure5b_table, 'data/Figure5b_table.Rds')


# facet by highest proportionate editing

Figure5b_table %>% 
  group_by(siteID) %>% mutate(edRank=rank(desc(value),ties.method = 'random')) %>% 
  arrange(siteID) %>% 
  group_by(siteID) %>% 
  mutate(facetTag=ifelse(edRank==1,as.character(area),NA_character_)) %>% 
  arrange(siteID,facetTag) %>% 
  fill(facetTag, .direction='down') %>% ungroup() %>% 
  mutate(facetTag=factor(facetTag,levels=paste0('BA',c(8,10,17,21,22,41)))) %>%  
  arrange(site_name) %>% 
  mutate(site_name=factor(site_name,levels=rev(unique(site_name)))) %>% 
  ggplot(aes(x=site_name,y=value,group=pair)) + 
  geom_line(lwd=0.25, col='grey', position=position_dodge(width=0.75), 
            show.legend = FALSE) +
  geom_point(aes(color=area,shape=nType), alpha=0.75, position=position_dodge(width=0.75)) + 
  ylab('Proportion of detected sites edited') + xlab("") +
  coord_flip() + facet_wrap(~facetTag, scales='free_y',ncol=2,dir = 'v') +
  theme_classic()

#--> add ASD/SCZ DEd site overlay [Figure5b_plot.R]


#GSEA for differentially edited sites

dEd_entrez <- Figure5b_table %>% filter(value<0.05) %>% select(site_name) %>% distinct() %>% 
  separate(site_name, into=c('gene_name','siteID'),sep=' ') %>% 
  left_join(anno %>% select(gene_name,entrezid)) %>% unnest(cols=c(entrezid)) %>% 
  filter(!is.na(entrezid)) %>% 
  group_by(gene_name) %>% slice(1) 

dEd_GOenrich <- limma::goana(dEd_entrez$entrezid)

dEd_GOenrich_table <- dEd_GOenrich %>% 
  rownames_to_column(var='goid') %>% filter(P.DE < 0.05,DE>=5) %>% 
  arrange(Ont,P.DE)

saveRDS(dEd_GOenrich_table, 'data/stat_tests/binom_testing/areaXnType_GOenrich.Rds')



