#p1 intersect genomic features, genic features, RepeatMasker tracks and stranding information


library(tidyverse)


setwd('~/Dropbox/Bahlo_Lab/scRNA_editing/phs000834_BlueLake_etal/gatk_output/for_GitHub/')

here::set_here()
here::here()

dt_siteStats <- readRDS("data/phs000834/dt_siteStats.Rds") 


########################
## Tag COGNATE STRAND ##
########################

#create bed file for bedtools intersect
dt_siteStats %>%  dim() #50K sites
dt_siteStats %>% count(siteID)

#check non-common variants (including putative RNAed sites intersection with cognate (same strand) bounding ensembl features.
#  write bed file:
write_tsv(dt_siteStats %>% 
            separate(siteID, into=c('region','position'), convert = TRUE) %>% 
            mutate(posPlus = position + 1,
                   val = 1,
                   strand = ifelse(basechange=="A_G","+","-")) %>% 
            arrange(region,position) %>% 
            select(1, 2, posPlus,basechange,val,strand) ,
          path="data/bedtools_input/dt_siteStats.bed")

# bedtools intersect -b <(grep -v "region" dt_siteStats.bed | sort | uniq)  -a ../../../hg38_merge_genes/hg38_reorder_genes.bed \
# | sed 's/"//g' | tr -d ';' > gatk_AGall_hg38_intersect.tsv

#genic features
# bedtools intersect -b <(grep -v region  dt_siteStats.bed) \
# -a ../../../hg38_merge_exons/Homo_sapiens_ProteinCodingGenes.GRCh38.91.bed |  \
# cut -f1-10 |tr -s ' ' '\t' | sed 's/"//g' | sed 's/;//g'  > dt_siteStats_genicFeatures.out

#bedtools intersect with RepeatMasker track (annotated repeat regions incl. Alu regions)

# bedtools intersect -b  <(grep -v "region" dt_siteStats.bed | sort | uniq) \
# -a hg38/UCSC_Repeats_reformat.bed.gz > gatk_AGall_UCSC_Repeats_intersect.tsv


########## ########## ########## ########## ########## ##########
########## ########## ########## ########## ########## ##########

#Tag sites by gene overlap status and cognate strand:


gatk_intersect <- read_tsv('data/bedtools_output/gatk_AGall_hg38_intersect.tsv', col_names = FALSE, col_types = cols('X1'=col_character()) ) %>% 
  mutate(siteID=paste(X1,X2,sep="_")) %>% select(-c(X3,X6,X8)) %>% 
  spread(X9,X10) %>% spread(X11,X12) %>% spread(X13,X14) %>% spread(X15,X16) %>% spread(X17,X18) %>% select(siteID,everything()) %>% 
  dplyr::rename(region=X1,position=X2,source=X4,gene=X5,strand=X7)
  

################## ################## ##################
################## ################## ##################

tag_cognate_p1 <- dt_siteStats %>% 
  select(1:basechange) %>% 
  left_join(gatk_intersect,by='siteID') %>% 
  group_by(siteID) %>% mutate(n_oLap=n()) %>% ungroup() %>% #filter(is.na(strand)) %>% select(siteID) %>% distinct() %>% dim()
  mutate(strandMatch = case_when(is.na(strand) ~ NA_character_, 
                               basechange=="A_G" & n_oLap==1 & strand=="+" ~ 'single_cogn_Watson',
                               basechange=="T_C" & n_oLap==1 & strand=="-" ~ 'single_cogn_Crick',
                               basechange=="A_G" & n_oLap==1 & strand=="-" ~ 'non_cognate',
                               basechange=="T_C" & n_oLap==1 & strand=="+" ~ 'non_cognate',
                               TRUE ~ 'other')) 
tag_cognate_p1 %>% count(strandMatch,sort=TRUE)

tag_cognate_p2 <- tag_cognate_p1 %>% filter(strandMatch=="other"|is.na(strandMatch)) %>% 
  group_by(siteID) %>% mutate(site_status = case_when(rdrType=='uncatalog' & is.na(strandMatch) ~ 'drop',
                                                      rdrType!='uncatalog' & strandMatch=='other' & basechange=="A_G" & all(strand=="+") ~ 'multi_cogn_Watson',
                                                      rdrType!='uncatalog' & strandMatch=='other' & basechange=='T_C' & all(strand=="-") ~ 'multi_cogn_Crick',
                                                      rdrType!='uncatalog' & strandMatch=='other' & basechange=="A_G" & any(strand=="+") ~ 'multi_cogn_mix',
                                                      rdrType!='uncatalog' & strandMatch=='other' & basechange=="T_C" & any(strand=="-") ~ 'multi_cogn_mix',
                                                      rdrType!='uncatalog' & strandMatch=='other' & basechange=="A_G" & !any(strand=="+") ~ 'multi_noncogn',
                                                      rdrType!='uncatalog' & strandMatch=='other' & basechange=="T_C" & !any(strand=="-") ~ 'multi_noncogn',
                                                      is.na(strandMatch) ~ 'no_ENSG',
                                                      TRUE ~ 'oLap_mix')) %>% ungroup() 
#NB oLap_mix indicates uncatalogued sites overlapped by genes on both plus ('Watson') and minus ('Crick') strands. 



tag_cognate_p2 %>% count(rdrType, site_status, sort=TRUE) 
#if all overlapping genes are cognate with the site, retain site;
#if there is a mix of cognate and non-cognate (i.e. overlapping gene regions), drop site

sites_tagged <- rbind(
  tag_cognate_p1 %>% filter(!is.na(strandMatch) , strandMatch!= 'other') %>% 
    mutate(site_status = strandMatch),
  tag_cognate_p2)




dt_siteStats_anno <- sites_tagged %>% filter(site_status!="drop", site_status!="non_cognate",site_status!='multi_noncogn',site_status!='oLap_mix') %>% 
  mutate(tag=ifelse(rdrType=='uncatalog' & site_status=='no_ENSG', 'drop','retain')) %>% 
  filter(tag!="drop") %>% 
  rename(ENSGID=gene_id,ENSGbioType=gene_biotype)
  
dt_siteStats_anno %>% select(siteID) %>% distinct() %>% dim()   #40.9K



# - retained reported sites from RADAR;
# - removed non_cognate sites; retained overlapping/no_ENSG RADAR sites only


################ ################ ################ ################
################ ################ ################ ################


dt_siteStats_anno %>% count(siteID) %>% dim()
dt_siteStats_anno %>% filter(str_detect(site_status,'single')) %>% count(siteID,site_status,sort=TRUE) 
dt_siteStats_anno %>% filter(!str_detect(site_status,'single')) %>% count(siteID,site_status,sort=TRUE) 

dt_siteStats_anno %>% count(countSingle = str_detect(site_status,'single')) %>% mutate(sum=sum(n))
#36K of 44K (83%) are bounded by a single non-overlapping ENSG feature region



####################  ######################
######### REPEAT MASKER Alu DATABASE #######
####################  ######################

#Alu features in ~/~/hg38/UCSC_Alu_reformat.bed.gz #subset of:
# Alu features in ~/~/hg38/UCSC_Repeats_reformat.bed.gz

RM_Repeats_intersect <- read_tsv("data/bedtools_output/gatk_AGall_UCSC_Repeats_intersect.tsv",col_names = FALSE,
                                 col_types=cols('X1'=col_character())) %>% 
  mutate(siteID=paste(X1,X2,sep="_")) %>% select(-X3,-X5,-X6) %>% 
  dplyr::rename(region=X1,position=X2,Repeat_name=X4,RM_strand=X7,Repeat_id=X8) %>% 
  mutate(Repeat_id=paste(Repeat_name,Repeat_id,sep="_")) %>%  
  select(siteID,Repeat_id,everything()) %>% 
  mutate(RM_repType=ifelse(str_detect(Repeat_name,'Alu'),'Alu','rep_nonAlu'))

RM_Repeats_intersect %>% count(siteID,sort=TRUE) #non-redundant
RM_Repeats_intersect %>% count(RM_repType)


dt_siteStats_anno %>% left_join(RM_Repeats_intersect %>% select(1,2,RM_strand,RM_repType),by='siteID') %>% 
  filter(rdrType=='uncatalog') %>% 
  select(siteID,basechange,strand,RM_repType,RM_strand) %>% 
  count(( RM_repType=="Alu" & basechange=="A_G" & RM_strand=="-" ) | (RM_repType=="Alu" & basechange=='T_C' & RM_strand=="+") ) 

#Uncatalogued sites are bounded by a cognate ensemble feature; the majority are bounded by a cognate Alu as well.
# Cannot distinguish Alu from non-Alu sites in cases where uncatalogued site is antisense to an Alu site. 
# These sites may be either Alu or (non-)rep non-Alu.


dt_siteStats_anno_RM <- dt_siteStats_anno %>% 
  left_join(RM_Repeats_intersect %>% 
                  filter(siteID %in% (dt_siteStats_anno %>% filter(rdrType=="uncatalog") %>% pull(siteID))) %>% 
                  select(1, Repeat_id, RM_strand, RM_repType),
            by='siteID')

dt_siteStats_anno_RM %>% filter(rdrType=="uncatalog") %>% 
  count(site_status,RM_repType, basechange, strand, RM_strand, sort=TRUE) %>% mutate(prop=n/sum(n))
  


################ ################ ################ ################
################ ################ ################ ################


#Report GENIC FEATURES (within coding genes) for simple sites (i.e., non-OL, cognate)


genFeatures <- read_tsv("data/bedtools_output/dt_siteStats_genicFeatures.out",
                        col_names = FALSE, col_types = cols("X1" = col_character())) %>% 
  magrittr::set_colnames(c('region','siteStart','siteEnd','feature','sense','strand','ENSGID')) %>% distinct() %>% 
  mutate(siteID=paste(region,siteStart,sep="_"))

genFeatures %>% filter(siteID %in% (dt_siteStats_anno_RM %>% filter(str_detect(site_status,'single')) %>% pull(siteID))) 
#NB sites in certain ~3000 genes can be in mutliple features depending on the transcript isoform.

genFeatures %>% count(siteID,sort=TRUE) 
genFeatures %>% count(feature)

# NB 3' UTR, 5' UTR start and stop codons are a sub-set of the exon (mature RNA transcript)

#cf:
#1       ensembl_havana  exon             100022386      100035637       .       +       .       gene_id "ENSG00000117620"; gene_version "14"; transcript_id "ENST00000533028"
#1       ensembl_havana  three_prime_utr 100022477       100035637       .       +       .       gene_id "ENSG00000117620"; gene_version "14"; transcript_id "ENST00000533028"

genFeatures %>% 
  filter(siteID %in% (dt_siteStats_anno_RM %>% filter(str_detect(site_status,'single')) %>% pull(siteID))) %>% 
   count(siteID,sort=TRUE) %>% count(n>1) 


genFeatures %>% filter(siteID %in% (dt_siteStats_anno_RM %>% filter(str_detect(site_status,'single')) %>% pull(siteID))) %>% 
  count(ENSGID,sort=TRUE) #4.3K genes represented


genFeatures %>% 
  select(siteID,ENSGID, feature) %>% 
  mutate(dummy=1) %>% spread(feature,dummy) %>% 
  filter(siteID=="1_100003296") #multi-cognate = overlapped by two ensembl genes

dt_siteStats_anno_RM %>% count(str_detect(site_status,'single')) %>% mutate(prop=n/sum(n))
#De-duplicate genic features (for the 83% of sites which are non-overlapping)



#SPREAD genic features

genFeatures_spread <- genFeatures %>% select(siteID, ENSGID, feature) %>% #filter(siteID=="1_70231974") #single
  group_by(siteID,feature) %>% slice(n=1) %>% ungroup() %>% #NB this is to ensure sites overlaping multiple cognate genes are represented only once.
  mutate(dummy=1) %>% spread(feature,dummy) %>% #count(gene) 34270
  mutate(intronic = ifelse( gene==1  & is.na(exon) & is.na(five_prime_utr) & is.na(stop_codon) & is.na(three_prime_utr),1,NA)) %>% 
  mutate(exonic   = ifelse(is.na(intronic),1,NA)) %>% 
  select(-c(gene,exon)) 

#all sites:
genFeatures_spread %>% count(intronic,exonic,five_prime_utr,start_codon, stop_codon , three_prime_utr,sort=TRUE)

#nonOL sites:
genFeatures_spread %>% filter(siteID %in% (dt_siteStats_anno_RM %>% filter(str_detect(site_status,'single')) %>% pull(siteID))) %>% 
  count(intronic,exonic,five_prime_utr,start_codon, stop_codon , three_prime_utr,sort=TRUE)


#deduplicate:
genFeatures_deDup <- genFeatures_spread %>% 
  filter(siteID %in% (dt_siteStats_anno_RM %>% filter(str_detect(site_status,'single')) %>% pull(siteID))) %>% 
  gather(key,value,-c(siteID,ENSGID)) %>% na.omit() %>% 
  group_by(siteID) %>% mutate(n_Obs = n()) %>% ungroup() %>% 
  arrange(desc(n_Obs),siteID) %>% #count(key,sort=TRUE) %>% #filter(siteID=="1_70231974") %>% 
  mutate(tag=ifelse(key=='exonic' & n_Obs > 1,'drop','retain')) %>% # Drops 'exonic' counts that are redundant (i.e. with 3'utr / 5'utr)
  #filter(key=='exonic' & n_Obs==1) %>% arrange(siteID) %>% dim() #915
  filter(tag=="retain") %>% #filter(siteID=="1_116389630") %>% 
  # Count genic features per site & identify sites that can occupy different transcriptomic features depending on isoform (i.e., n_Iso > 1)
  group_by(siteID) %>% mutate(n_Iso=n()) %>% ungroup() %>% 
  mutate(genicFeature = ifelse(n_Iso==1, key, paste0('multipl:_',n_Iso))) %>% 
  select(-tag,-n_Obs) 

genFeatures_deDup %>% count(genicFeature)
genFeatures_deDup %>% head()
genFeatures_deDup %>% select(1,2,genicFeature) %>% distinct() %>% dim() #34228; non-redundant



#join genic features with dt_siteStats_anno_RM & add description from GRCHg38

anno <- readRDS("data/expn_processing/GRCh38_genes.Rds")


dt_siteStats_filt <- dt_siteStats_anno_RM %>% #redundant on multi-cognate sites
  left_join(genFeatures_deDup %>% select(1,2,genicFeature) %>% distinct(), by=c('siteID','ENSGID')) %>% #NR
  mutate(site_type = case_when(rdrType!="uncatalog" ~ rdrType,
                               rdrType=='uncatalog' & RM_repType == 'Alu' ~ 'Alu_Novel',
                               rdrType=='uncatalog' & RM_repType == 'rep_nonAlu' ~ 'rep_nonAlu_Novel',
                               TRUE ~ 'nonRep_Novel')) %>%  
  left_join(dt_siteStats, by = c("siteID", "rdrType", "basechange")) %>% 
  left_join(anno %>% select(gene_id, description), by=c('ENSGID'='gene_id')) #NR


#40861 sites; NB this is redundant on overlapping genes
dt_siteStats_filt %>% count(siteID,site_status,sort=TRUE)

dt_siteStats_filt %>% 
  filter(ENSGbioType=='protein_coding') %>% 
  count(rdrType,RM_repType,genicFeature,sort=TRUE) %>% mutate(sum=sum(n))

dt_siteStats_filt %>% count(site_type,sort=TRUE)
#dt_siteStats_filt contains redundant overlapping cognate gene info

saveRDS(dt_siteStats_filt, file="data/phs000834/dt_siteStats_filt.Rds")



#Count sites in non-overlapping protein-coding regions

dt_siteStats_filt %>% count(siteID) %>% dim() #40.861K
dt_siteStats_filt  %>% filter(str_detect(site_status,'single')) %>% 
  filter(ENSGbioType=='protein_coding') %>% count(siteID) %>% dim() #34.218K
#34228 / 40861 = 83.7%.


#Calculate proportion of sites in protein coding genes that are documented in the RADAR database
dt_siteStats_filt  %>% filter(str_detect(site_status,'single')) %>% 
  filter(ENSGbioType=='protein_coding') %>% count(rdrType) %>% 
  mutate(prop=n/sum(n)) %>% mutate(cs=cumsum(prop))
