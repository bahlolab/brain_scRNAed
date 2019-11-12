library(tidyverse)


#Explore SRA and gene expression

BL_joinKey <- readRDS('BL_metadata.Rds')

source('dplyr_override.R')


##################################################
########### cross with expression data ###########
##################################################


cluster_data <- readRDS('expn_processing/cluster_data.Rds') %>% rename(apop_score=apoptotic) 
cluster_data %>% summary()


flag_records <- readRDS("expn_processing/flag_records.Rds") %>% dplyr::rename(sample=value,expn_flag=key)


BL_joinKey %>% left_join(cluster_data, by=c('SRA'='sample')) %>% 
               left_join(flag_records,by=c('SRA'='sample')) %>% 
  count(expn_flag)

BL_joinKey %>% left_join(cluster_data, by=c('SRA'='sample')) %>% 
  left_join(flag_records,by=c('SRA'='sample')) %>% 
  count(expn_flag, is.na(apop_score)) %>% mutate(sum = sum(n))

BL_joinKey %>% left_join(cluster_data, by=c('SRA'='sample')) %>% 
  left_join(flag_records,by=c('SRA'='sample')) %>% 
  filter(!is.na(expn_flag)) %>% 
  count(SRA, expn_flag) %>% count(SRA,sort=TRUE)


flag_spread <- flag_records %>% mutate(dummy=1) %>% spread(expn_flag,dummy) 
flag_spread %>% count(high_mitoCells, low_featureCells, low_libsizeCells) %>% mutate(sum=sum(n))


BL_joinKey %>% 
  left_join(flag_spread, by=c('SRA'='sample')) %>% 
  filter(is.na(low_featureCells), is.na( high_mitoCells))  %>% dim()
#3055 cells

BL_HQ_joinKey <- BL_joinKey %>% left_join(cluster_data, by=c('SRA'='sample')) %>% 
  left_join(flag_spread, by=c('SRA'='sample')) %>% 
  filter(is.na(low_featureCells), is.na( high_mitoCells)) %>% 
  mutate(area=factor(area,levels=c('BA8','BA10','BA17','BA21','BA22','BA41'))) 

saveRDS(BL_HQ_joinKey, "BL_metadata_highQualCells.Rds")



############################################################
######### cross with STAR 2-pass mapping statistics ########


mapping_stats <- read_tsv("data/mapping_output/starMapping_stats.out",col_names = FALSE) %>% 
  dplyr::rename(sample=X1,category=X2,value=X3)


mapping_stats %>% dplyr::count(sample,sort=TRUE)
mapping_stats %>% distinct() %>% dplyr::count(sample,sort=TRUE)
mapping_stats %>% dplyr::filter(sample=="SRR5900438") 

mapping_stats %>% dplyr::filter(sample %in% BL_HQ_joinKey$SRA) %>% dplyr::count(sample,sort=TRUE)
#some samples are duplicated (mapped twice)

mapping_stats_BL <- mapping_stats %>% dplyr::filter(sample %in% BL_HQ_joinKey$SRA) %>% group_by(sample) %>% 
  mutate(rowID = dplyr::row_number()) %>% filter(rowID %in% c(1:4)) %>% ungroup() 

mapping_stats_BL <- mapping_stats_BL %>% select(-rowID) %>% spread(category,value) %>% 
  rename(percent_unmapped_tooShort=`pct_reads_unmapped:_too_short`) %>% 
  mutate(pct_uniqMap = Uniquely_mapped_reads_number/Number_of_input_reads,
         pct_multiMap = Number_of_reads_mapped_to_multiple_loci/Number_of_input_reads,
         pct_unmapped = percent_unmapped_tooShort/100,
         map_sum = Uniquely_mapped_reads_number + Number_of_reads_mapped_to_multiple_loci)

mapping_stats_BL %>% 
  select(sample, contains('pct')) %>% 
  gather(key,value,-sample) %>% 
  ggplot(aes(x=value)) + geom_density(aes(col=key))
#NB the 'unmapped' includes a large proportion of ERCC spike-in reads.


mapping_stats_BL %>% summary()


mapping_stats_BL %>% 
  select(1,3,5) %>% mutate(sum=Number_of_reads_mapped_to_multiple_loci + Uniquely_mapped_reads_number) %>% 
  select(1,3,4) %>% gather(key,value,-sample) %>% 
  ggplot() + geom_histogram(aes(x=value,fill=key),position='dodge') 


saveRDS(mapping_stats_BL, 'data/mapping_output/mapping_stats_BL.Rds')


#chart mapped reads by cortical area
mapping_stats_BL %>% left_join(BL_HQ_joinKey,by=c('sample'='SRA')) %>% 
  filter(Group_ID!="NoN") %>% 
  select(area,Group_ID,Uniquely_mapped_reads_number) %>% 
  ggplot(aes(x=area,y=Uniquely_mapped_reads_number)) + geom_boxplot()

mapping_stats_BL %>% left_join(BL_HQ_joinKey,by=c('sample'='SRA')) %>% 
  filter(Group_ID!="NoN") %>% 
  dplyr::select(area,Group_ID,Uniquely_mapped_reads_number) %>% 
  do(broom::tidy(lm(Uniquely_mapped_reads_number ~ factor(area), data = .))) %>% 
  arrange(p.value)
#BA8 showed significantly more uniquely mapped reads than other cortical regions

#chart mapped reads by neuronal sub-group
mapping_stats_BL %>% left_join(BL_HQ_joinKey,by=c('sample'='SRA')) %>% 
  filter(Group_ID!="NoN") %>% 
  select(area,Group_ID,Uniquely_mapped_reads_number) %>% 
  ggplot(aes(x=Group_ID,y=Uniquely_mapped_reads_number)) + geom_boxplot()

mapping_stats_BL %>% left_join(BL_HQ_joinKey,by=c('sample'='SRA')) %>% 
  filter(Group_ID!="NoN") %>% 
  dplyr::select(area,Group_ID,Uniquely_mapped_reads_number) %>% 
  do(broom::tidy(lm(Uniquely_mapped_reads_number ~ factor(Group_ID), data = .))) %>% 
  arrange(p.value)


###### test independent assortment of neuronal phenotype and cortical area ######

neuType_X_area <- BL_HQ_joinKey %>% filter(neuType!="NoN") %>% 
  count(area,neuType,sort=TRUE) %>% 
  spread(neuType,n) %>% mutate(prop=In/(In + Ex))

BL_HQ_joinKey %>% filter(neuType!="NoN") %>% 
  count(area,Group_ID,sort=TRUE) %>% 
  spread(area,n) 


BL_HQ_joinKey %>% filter(neuType!="NoN") %>% 
  count(area,Group_ID) %>% 
  spread(Group_ID,n) %>% gather(key,value,-area) %>% 
  mutate(key=factor(key,levels=rev(unique(key)))) %>% 
  ggplot(aes(x=area,y=key)) + geom_tile(col='white',aes(fill=value)) +
  geom_text(aes(label=value),col='white') +
  scale_fill_viridis_c()  + ylab('Neuronal group') + xlab('Brodmann Area') 
ggsave('charts/GroupID_x_area_count.pdf')


BL_HQ_joinKey %>% filter(neuType!="NoN") %>% 
  count(area,neuType) %>% 
  ggplot(aes(x=area,y=n,fill=neuType)) + geom_col() 
ggsave("charts/neuType_x_area_count.pdf",width=4,height=3)


BL_HQ_joinKey %>% select(SRA,Group_ID) %>% distinct() %>% count(Group_ID) %>% summary() #mean n:179
BL_HQ_joinKey %>% select(SRA,area) %>% distinct() %>% count(area) %>% summary() #mean n:509

neu_area_mat <- neuType_X_area %>% select(Ex,In) %>% as.matrix()
rownames(neu_area_mat) <- neuType_X_area$area
chisq.test(neu_area_mat) #p 1.3e-15


