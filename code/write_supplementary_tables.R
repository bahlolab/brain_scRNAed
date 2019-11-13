#create_supp_tables


library(tidyverse)
here::here()

########## ########## ########## ##########

#Ttable 1 (Lake)

BL_joinKey <- readRDS('data/phs000834/BL_metadata_highQualCells.Rds')

STable1 <- BL_joinKey %>% select(1:area) %>% rename(Neuronal_type = neuType, Brodmann_Area = area)

write_tsv(STable1,'data/Supplementary_Tables/STable1.tsv',col_names = TRUE)

########## ########## ########## ##########

#STable2 (Darmanis)

STable2 <- readRDS('data/GSE67835/Darmanis_Neurons_meta.Rds') %>% select(1:4,8) %>% 
  magrittr::set_colnames(c('SRA','Cell_type','Age_years','Donor_ID','Sex'))

write_tsv(STable2, 'data/Supplementary_Tables/STable2.tsv',col_names = TRUE)

########## ########## ########## ##########

#STable3 

dt_siteStats_filt <- readRDS('data/phs000834/dt_siteStats_filt.Rds')

for_STable3 <- dt_siteStats_filt %>%
  select(-strandMatch,-tag,-gene,-gene_version,-gene_source,-RM_strand) %>% 
  magrittr::set_colnames(c("Site_ID", "RADAR_type", 'base_change',
                           'chromosome','position','gene_source',"gene_strand",'ENSMBL_GeneID','gene_name','gene_biotype',
                           'N_overlapping_genes','Site_status','RepeatMasker_ID','RepeatMasker_type','Genic_feature','Site_type','N_cells_detected',
                           'Mean_AltAllele_Proportion','Median_AltAllele_Proportion','Gene_description')) %>% 
  mutate(Site_status = str_replace(Site_status,'Watson','Plus'),
         Site_status = str_replace(Site_status,'Crick','Minus')) %>% 
  select(1,chromosome,position,base_change,Site_type,Site_status,everything(),-gene_source,-RADAR_type) 


##

dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")

Darmanis_stats <- dt_siteStats_TDjoin %>% select(siteID,pxID,contains('Darmanis')) %>% distinct() %>% spread(pxID,n_Cells_Darmanis) %>% 
  select(-contains('Darmanis')) %>% select(-`<NA>`) %>% 
  gather(key,value,contains('AB_')) %>% filter(!is.na(value)) %>% 
  mutate(key=paste('Darmanis_2015',key,sep="_")) %>% 
  spread(key,value,fill=0) %>% 
  mutate(Total_nCells_Darmanis = Darmanis_2015_AB_S11 + Darmanis_2015_AB_S4 + Darmanis_2015_AB_S5 + Darmanis_2015_AB_S7)


Tan_stats <- dt_siteStats_TDjoin %>% select(siteID,Tan_etal) %>% distinct() %>% rename(Tan_etal_FI = Tan_etal )


STable3 <- for_STable3 %>% left_join(Darmanis_stats,by=c('Site_ID' = 'siteID')) %>% left_join(Tan_stats,by=c('Site_ID' = 'siteID'))

write_tsv(STable3, 'data/Supplementary_Tables/STable3.tsv', col_names = TRUE)


########## ########## ########## ##########

#STABLE 4 

# lm for GEI between neuronal types/cortical regions

lm_area_GEI <- readRDS('data/stat_tests/lm_area_GEI.Rds')

lm_Group_ID_GEI <- readRDS('data/stat_tests/lm_Group_ID_GEI.Rds')

STable4 <- rbind(lm_area_GEI %>% mutate(FDR= p.adjust(p.value,method = 'BH')) %>% mutate(test = 'CorticalRegion_v_GEI'),
                 lm_Group_ID_GEI %>% mutate(test = 'Neuronal_Subtype_v_GEI')) %>% 
  select(test,everything())

write_tsv(STable4, 'data/Supplementary_Tables/STable4.tsv')

########## ########## ########## ##########

#STABLE 5

anno <- readRDS('data/expn_processing/GRCh38_genes.Rds')

lm_out_INT_neuType   <-  readRDS('data/stat_tests/broom_lm/lm_out_INT_neuType.Rds') %>% 
  ungroup() %>% mutate(FDR=p.adjust(p.value,method="BH")) %>% arrange(FDR)


#reported DNA/RNA binding proteins & modulators of editing

edCor_v_Tan <- readRDS('data/stat_tests/edCor_v_Tan.Rds') %>% dplyr::rename(Tan_Estimate = Estimate)
edCor_v_QV <- readRDS('data/stat_tests/edCor_v_QV.Rds') %>% dplyr::rename(Quinones_Valdez_Mean_Editing_Change_on_Knockdown = Mean_Editing_Change,
                                                                          Quinones_Valdez_Tag = tag)
edCor_v_HO <- readRDS('data/stat_tests/edCor_v_HO.Rds')


STable5 <- lm_out_INT_neuType %>% left_join(anno %>% select(gene_id,gene_name,description), by = c('ENSG' = 'gene_id')) %>% 
  filter(term=='log10(value)') %>% mutate(term = 'log10(TPM)') %>% 
  left_join(edCor_v_Tan %>% select(ENSGpref,Tan_Estimate,tag),by=c('ENSG' = 'ENSGpref')) %>% rename(Tan_Corr = tag) %>% 
  left_join(edCor_v_QV %>% select(ENSG,contains('Quinones'),concordance), by='ENSG') %>% rename(Quinones_Valdez_concordance = concordance) %>% 
  mutate('Hudson_Ortlund_reported_DNA_RNA_binding_protein' = ifelse(ENSG %in% edCor_v_HO$ENSG, 'yes','no'))

write_tsv(STable5, 'data/Supplementary_Tables/STable5.tsv')


########## ########## ########## ##########

#STABLE 6

#GOANA for lm_out_INT_neuType [explore_edProp_v_TPM.R]

GOANA_pos_annot <-  readRDS('data/stat_tests/broom_lm/GOANA_pos_annot.Rds')

GOANA_neg_annot <-  readRDS('data/stat_tests/broom_lm/GOANA_neg_annot.Rds')

STable6 <- rbind(GOANA_pos_annot %>%  mutate(correlation = 'positive_on_GEI'),
                 GOANA_neg_annot %>%  mutate(correlation = 'negative_on_GEI')) %>% 
  dplyr::rename(goid = term, go_term = Term, N_background_genes = N, N_foreground_genes = DE, Probability_of_Enrichment = P.DE ,
                entrezid= gene_id) %>% select(-Ont) 


write_tsv(STable6, 'data/Supplementary_Tables/STable6.tsv')


########## ########## ########## ##########

#STABLE 7

#VEP + Nishikura known editing sites

STable7 <- readRDS('data/VEP_Nishikura_missense_tbl.Rds')

write_tsv(STable7 %>% select(-Refs), 'data/Supplementary_Tables/STable7.tsv')



########## ########## ########## ##########

#STABLE 8

#clinVar results

clinVar_celltypes <- readRDS('data/clinVar_cellTypes.Rds')

STable8 <- clinVar_celltypes %>% 
  rename(Site_ID= siteID, ClinVar_ID = clinVarID, ClinVar_desc = ClinVar, 
         Neuronal_type = neuType, Cortical_region = area, N_cells = n)

write_tsv(STable8, 'data/Supplementary_Tables/STable8.tsv')

########## ########## ########## ##########

#STABLE 9

#GO term enrichment in differentially edited genes

STable9 <- readRDS('data/stat_tests/binom_testing/areaXnType_GOenrich.Rds') %>% 
  rename(go_term = Term, Ontology = Ont, N_background_genes = N, N_foreground_genes = DE, Probability_of_Enrichment = P.DE )

write_tsv(STable9, path='data/Supplementary_Tables/STable9.tsv')

########## ########## ########## ##########

#STABLE 10

#site-wise differential editing (post-hoc)

STable10 <- readRDS('data/stat_tests/binom_testing/fig5B_STable.Rds') %>% select(-contains('tag'))

STable10 %>% filter(p.adj.Fisher<0.05) %>% count(siteID) 

STable10 %>% filter(!is.na(disease_tag)) %>% count(siteID) #28

write_tsv(STable10, 'data/Supplementary_Tables/STable10.tsv')

########## ########## ########## ##########



