#Plot Figure 5b 

#### READ DATA ####

Figure5b_table <- readRDS('data/Figure5b_table.Rds')

areaXnType_ph_tbl   <- readRDS('data/stat_tests/binom_testing/areaXnType_ph_tbl.Rds') 

clin_join <- readRDS('data/join_clinical.Rds')


source('code/dplyr_override.R')



###################### ###################### ######################
###################   PLOT w CLINICAL OVERLAY    ###################
###################### ###################### ######################


areaXnType_ph_tbl %>% filter(p.adj.Fisher < 0.1) %>% 
  left_join(clin_join,by='siteID') %>%   
  rowwise() %>% 
  filter(sum(is.na(Breen_CM_ACC),is.na(Breen_CM_DLPFC),is.na(Breen_HBCC_DLPFC),is.na(Tran_FC),is.na(Tran_TC))<5) %>% 
  count(siteID)
#29 sites are 'DEd' between area-neuTypes and in neuropsych. conditions


scz_asd_cross <- areaXnType_ph_tbl %>% filter(p.adj.Fisher< 0.1) %>% 
  left_join(clin_join,by='siteID') %>%   
  rowwise() %>% 
  filter(sum(is.na(Breen_CM_ACC),is.na(Breen_CM_DLPFC),is.na(Breen_HBCC_DLPFC),is.na(Tran_FC),is.na(Tran_TC))<5)


####### add clinical overlap to plot ########

disease_tag <- scz_asd_cross %>% select(siteID,contains('Breen'),contains('Tran', ignore.case =FALSE )) %>% 
  distinct() %>% 
  mutate(scz_tag = ifelse((!is.na(Breen_CM_ACC) | !is.na(Breen_CM_DLPFC) | !is.na(Breen_HBCC_DLPFC)),"*",""),
         asd_tag= ifelse((!is.na(Tran_FC) | !is.na(Tran_TC)),'°','')) %>% 
  mutate(disease_tag=paste0(asd_tag,scz_tag))



Figure5b_table %>% group_by(siteID) %>% mutate(edRank=rank(desc(value),ties.method = 'random')) %>% 
  arrange(siteID) %>% 
  mutate(nType = ifelse(nType =='In','Inhibitory','Excitatory')) %>% 
  group_by(siteID) %>% 
  mutate(facetTag=ifelse(edRank==1,as.character(area),NA_character_)) %>% 
  arrange(siteID,facetTag) %>% 
  fill(facetTag, .direction='down') %>% ungroup() %>% 
  mutate(facetTag=factor(facetTag,levels=paste0('BA',c(8,10,17,21,22,41)))) %>%  
  left_join(disease_tag,by='siteID') %>% 
  mutate(line_tag=ifelse(!is.na(disease_tag),disease_tag,'notag')) %>% 
  mutate(line_tag=forcats::fct_relevel(line_tag,'notag')) %>% 
  arrange(site_name) %>% 
  mutate(site_name=factor(site_name,levels=rev(unique(site_name)))) %>% 
  mutate(site_name_tag=ifelse(!is.na(disease_tag),paste0(disease_tag,site_name),as.character(site_name))) %>% 
  mutate(forRelevel = as.numeric(site_name)) %>% 
  mutate(site_name_tag = forcats::fct_reorder(site_name_tag, forRelevel)) %>% #count(siteID) #117
  rename(Area = area, `Neuronal type` = nType) %>% 
  ggplot(aes(x=site_name_tag,y=value,group=pair)) + 
  geom_line(lwd=0.25, aes(lty=line_tag), position=position_dodge(width=0.75), 
            show.legend = FALSE) +
  geom_point(aes(color=Area,shape=`Neuronal type`), alpha=0.75, position=position_dodge(width=0.75)) + 
  ylab('Proportion of detected sites edited') + xlab("") +
  coord_flip() + facet_wrap(~facetTag, scales='free_y',ncol=2,dir = 'v') +
  theme_classic()

ggsave('charts/areaXnType_facetArea_ASD_SCZ.pdf',width=10,height=10)


###### ###### ###### ###### ###### ######

#Write supp. table

areaXnType_diffType <- readRDS('data/stat_tests/binom_testing/areaXnType_diffType.Rds')

forSTable <- areaXnType_diffType %>% filter(description!="NULL") %>% 
  left_join(disease_tag,by='siteID') 


saveRDS(forSTable, 'data/stat_tests/binom_testing/fig5B_STable.Rds')



