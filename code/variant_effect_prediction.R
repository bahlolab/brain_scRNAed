#Predict bioligcal effects of edited sites (VEP)

library(here)
library(tidyverse)

###########################
########## VEP ############
###########################

#Create query file for VEP

dt_siteStats_bed <- read_tsv("data/bedtools_input/dt_siteStats.bed",
                             col_names = TRUE, col_types = cols('region'=col_character()))

#high quality filtered sites:
dt_siteStats_filt <- readRDS("data/phs000834/dt_siteStats_filt.Rds")
dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")

#sub-set on dt_siteStats_TDjoin
filtered_sites_toVEP <- dt_siteStats_bed %>% 
  mutate(siteID=paste(region,position,sep="_")) %>% 
  filter(siteID %in% dt_siteStats_TDjoin$siteID) %>% select(-siteID) %>% 
  mutate(basechange=str_replace(basechange,'_','/')) %>% 
  mutate(p2=position) %>% 
  select(region,position,p2,basechange,strand) 


write_delim(filtered_sites_toVEP,'data/VEP_input/filtered_sites_toVEP.txt', delim = " ", col_names = FALSE)

#run ensembl-vep on unix cluster

####################

VEP_out <- read_tsv('data/VEP_output/gatk_filtered_sites_vep_full.output.vcf', skip=87, col_names = TRUE,
                    col_types = cols(Location = col_character())) %>% 
  dplyr::rename(Uploaded_variation = `#Uploaded_variation` ) %>% 
  separate(Uploaded_variation, into=c('region','position','basechange'),sep="_") %>% 
  unite('siteID',region,position) %>% 
  filter(siteID %in% dt_siteStats_TDjoin$siteID)

saveRDS(VEP_out, file="data/VEP_output/VEP_out.Rds")

VEP_out %>% count(siteID) 40861


VEP_out %>% filter(!str_detect(Consequence,'intron|intergenic|downstream|upstream')) %>% 
  mutate(conSplit=str_split(Consequence,',')) %>% 
  unnest(conSplit) %>% 
  count(Consequence,sort=TRUE)


VEP_out %>% filter(Amino_acids!="-") %>% select(siteID, Consequence, Amino_acids) %>% distinct() %>% 
  filter(str_detect(Consequence,'missense|stop')) %>% 
  count(siteID,sort=TRUE) #327 coding missense, or stop mutations


VEP_out %>% filter(Amino_acids!="-") %>% select(siteID, Amino_acids) %>% distinct() %>% 
  left_join(dt_siteStats_filt,by='siteID') %>% 
  count(site_status,site_type,sort=TRUE)

VEP_out %>% filter(Amino_acids!="-") %>% select(siteID, Amino_acids) %>% distinct() %>% 
  left_join(dt_siteStats_filt,by='siteID') %>% 
  filter(str_detect(site_status,'single')) %>% count(gene_name,description,sort=TRUE)



######################################################
##################### EXPLORE VEP ###################
######################################################

#reverse complement

testCodon <- 'tTg'
chartr("ATGCatgc","TACGtacg", testCodon)


# codon context #

VEP_out %>% filter(Consequence=='missense_variant') %>% count(siteID,Amino_acids)

VEP_out %>% filter(Codons!="-") %>% count(Consequence,sort=TRUE) #only in coding regions

VEP_out %>% filter(Codons!="-") %>% select(siteID, basechange, Codons) %>% 
  separate(Codons, into=c('ref_codon',"edit_codon"),sep="/") %>% 
  count(ref_codon,edit_codon,sort=TRUE)

VEP_out %>% filter(Codons!="-") %>% select(siteID, basechange, Codons) %>% 
  separate(Codons, into=c('ref_codon',"alt_codon"),sep="/") %>% 
  mutate(refpos=case_when(str_detect(ref_codon,"^A|^T") ~ 1,
                          str_detect(ref_codon,"^.A|^.T") ~ 2,
                          str_detect(ref_codon,"^..A|^..T") ~ 3,
                          TRUE ~ 0)) %>% 
  filter((basechange == 'A/G' & str_detect(ref_codon,'A')) | basechange=='T/C' & str_detect(ref_codon,'T')) %>% 
  #count(siteID) #467 sites
  count(refpos, ref_codon,sort=TRUE) %>% 
  mutate(ref_c1=ifelse(str_detect(ref_codon,fixed('T',ignore_case = FALSE )), chartr("ATGCatgc","TACGtacg",ref_codon), ref_codon)) %>% 
  mutate(ref_c2=ifelse(str_detect(ref_c1,fixed('t', ignore_case=FALSE)), chartr('t','u',ref_c1),ref_c1)) %>% 
  group_by(refpos, ref_c2) %>% summarize(sum=sum(n)) %>% 
  ggplot(aes(x=ref_c2,y=sum)) + geom_col(col='white') +
  facet_wrap(~refpos, scales="free_x",ncol=1) +
  theme(axis.text.x=element_text(angle=90)) +
  xlab('Edited reference codon') + ylab('Sum of edited sites')

ggsave("charts/VEP_AG_edit_motifs.pdf",width=5,height=5)


####### Cross with Nishikura Nat Rev Mol Cell Biol 2015 #######

Nishikura_edited_codingGenes <- readxl::read_xlsx('data/related_studies/Nishikura_NatRevMolCellBiol/Table1.xlsx') %>% 
  fill(Gene,Protein) %>% separate(`Recoding ADAR responsible`, into=c('Amino_acids','ADARmember'), sep=" ") %>% 
  mutate(Amino_acids=str_replace(Amino_acids,'>','/')) %>% 
  mutate(Amino_acids=str_split(Amino_acids,",")) %>% unnest(Amino_acids)

Nishikura_edited_codingGenes %>% count(Gene) #19 edited coding genes reviewed by Nishikura

dt_siteStats_filt %>% filter(gene_name %in% Nishikura_edited_codingGenes$Gene) %>% count(gene_name, sort=TRUE) #10 of 19 are edited in phs000834

VEP_out %>% filter(Consequence=="missense_variant") %>% select(siteID, basechange,Gene,Amino_acids) %>% 
  left_join(anno %>% select(gene_id,gene_name),  by=c('Gene'='gene_id')) %>% 
  left_join(Nishikura_edited_codingGenes,  by=c('gene_name'='Gene'), suffix=c("_VEP",'_Nishikura')) %>% 
  filter(!is.na(Protein)) %>% distinct() %>% 
  filter(Amino_acids_VEP == Amino_acids_Nishikura) %>% 
  count(siteID, Amino_acids_VEP, gene_name,Function)

#6 of 19 exhibit identical edited site


#### Cross missense sites with Nishikura ###


VEP_Nishikura_missense_tbl <- VEP_out %>% filter(Amino_acids!="-") %>% select(siteID, Consequence, Amino_acids) %>% distinct() %>% 
  filter(str_detect(Consequence,'missense|stop')) %>% 
  left_join(dt_siteStats_filt,by='siteID') %>%  select(1:3,basechange,gene_name,description)  %>% 
  distinct() %>% arrange(siteID) %>% ungroup() %>% 
  left_join(Nishikura_edited_codingGenes, by=c('gene_name'='Gene'), suffix=c("_VEP",'_Nishikura')) 

saveRDS(VEP_Nishikura_missense_tbl, 'data/VEP_Nishikura_missense_tbl.Rds')



####### Intersect Tan etal Nature 2017 #######

#Which edited coding sites in phs000834 are also reported in Tan etal?

VEP_out %>% filter(Consequence=="missense_variant") %>% select(siteID, basechange,Gene,Amino_acids) %>% 
  left_join(dt_siteStats_TDjoin %>% select(-description),by='siteID') %>% 
  filter(!is.na(Tan_etal)) %>% 
  left_join(anno,by=c('Gene'='gene_id')) %>% 
  dplyr::select(siteID,Amino_acids,gene_name,description) %>% distinct()


######### ########## ######### ##########


#Distribution of n. cells in which missense sites are detected

dt_siteStats_TDjoin %>% select(siteID,n_Cells_Lake ) %>% 
  filter(siteID %in% (VEP_out %>% filter(Consequence=='missense_variant') %>% pull(siteID))) %>% 
  summary(n_Cells_Lake)
#median 15; mean 24; max 239.



