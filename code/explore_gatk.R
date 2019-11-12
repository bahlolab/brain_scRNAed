#explore all GATK

library(tidyverse)

dt_label <- readRDS("phs000834_aggregateSNPs_inclUncatalogd.Rds")

BL_HQ_joinKey <- readRDS("BL_metadata_highQualCells.Rds")

head(dt_label)
snp_dp5_alt2 <- dt_label %>% filter(totalDP>=5, altDP>=2) 


#retain variants covered with a depth of at least 5 reads, with an alternate allele count >=2.

snp_gt5_count <- snp_dp5_alt2 %>% select(chromosome,position,ref_al,alt_al) %>% 
   distinct() %>% mutate(complex=nchar(alt_al)) %>% filter(complex<2) %>% 
  count(ref_al,alt_al,sort=TRUE)


snp_gt5_count  #NB only A and T ref alleles are retained in the aggregation script [filter_vcf.R]
#AG and TC SNPs outweight other alternative alleles by at least 5:1

dt_filt <- snp_dp5_alt2 %>% 
  mutate(siteID=paste(chromosome,position,sep="_")) %>% 
  mutate(basechange=paste(ref_al,alt_al,sep="_")) %>% 
    filter(basechange %in% c('A_G','T_C')) 


dt_filt_BL <- dt_filt %>% filter(sample %in% BL_HQ_joinKey$SRA) %>% 
  filter(dbSNP_status!="commonSNP") %>% filter(chromosome!="Y") %>% 
  group_by(siteID) %>% mutate(n_Cells=n()) %>% 
  ungroup()


#retain SNPs detected in at least 10 individual cells:
dt_filt_gte10 <- dt_filt_BL %>% filter(n_Cells>=10) %>% 
  mutate(rdrType=ifelse(rdrType=='drnd_Only','uncatalog',rdrType))

dt_filt_gte10 %>% count(siteID) #50,472 sites; 54.2% previously reported
dt_filt_gte10 %>% count(rdrType,sort=TRUE) %>% mutate(prop=n/sum(n))

#distribution of vcf site statistics
dt_filt_gte10 %>% select(dbSNP_status,rdrType,totalDP:SOR) %>% head(n=10000) %>% 
  select(-c(AN,ClippingRankSum,DS,ExcessHet,FS,MQ,MQRankSum)) %>% 
  mutate(status=case_when(dbSNP_status=="commonSNP" ~ dbSNP_status,
                          ( QD >= 2 & ReadPosRankSum > -4) ~ "filt",
                          TRUE ~ "drop")) %>% 
  select(-dbSNP_status) %>% 
  gather(key,value,-c(rdrType,status)) %>%  
  ggplot(aes(x=value)) + 
  geom_density(aes(col=rdrType), alpha=0.5) + 
  facet_wrap(~ key, scales="free")



dt_filt_gte10 <- dt_filt_gte10 %>% filter(!str_detect(siteID,'Y'))

dt_filt_gte10 %>% dim() #1.14m records of ~50K SNPs


dt_siteStats <- dt_filt_gte10 %>% group_by(siteID,rdrType,basechange,n_Cells) %>% 
  summarize(meanAP=mean(altProp),medAP=median(altProp)) %>% ungroup()



##### ##### ##### #####
##### ##### ##### #####



dt_filt_gte10 %>% select(siteID,rdrType,dbSNP_status) %>% distinct() %>% 
  count(rdrType,dbSNP_status) %>% mutate(sum=sum(n), pct=n/sum)
#More than 50% of variants are previously reported in the RADAR database



dt_filt_gte10 %>% 
  select(siteID,rdrType,dbSNP_status) %>% distinct() %>% 
  count(rdrType,dbSNP_status) %>% 
  ggplot(aes(x=rdrType,y=n,fill=dbSNP_status))+ geom_col() +
  coord_flip()

dt_filt_gte10 %>% 
  select(siteID,rdrType,dbSNP_status,basechange) %>% distinct() %>% 
  count(rdrType,basechange) %>% 
  ggplot(aes(x=basechange, y=n, fill=rdrType)) + geom_col(col="white",position="dodge") 
ggsave('charts/basechange_frequency.pdf',width=5,height=4)


#### GOLD STANDARD GLUR2 (GRIA2) editing site; Nishikura 2015 Nat Rev Mol Cell Biol ####

dt_filt_gte10 %>% filter(dbSNP_status!="commonSNP") %>% 
  filter(siteID=="4_157336723") #59 cells  

dt_filt_gte10 %>% filter(siteID=="4_157336723") %>% 
  ggplot(aes(x=altDP,y=totalDP)) + 
  geom_abline() +
  geom_point() + coord_equal()
#100% penetrance in all 59 cells.

#What proportion of sites have 100% penetrance but are not known RADAR / dbSNP sites?

dt_siteStats %>% ggplot(aes(x=meanAP)) + geom_density(aes(col=rdrType))

dt_siteStats %>% ggplot(aes(x=medAP)) + geom_density(aes(col=rdrType))


#Mean alternate allele frequency (FI) is close to 1.


####### ####### ####### #######
  ###     WRITE OUTPUT    ###
####### ####### ####### #######

#Common genomic SNPs dropped:

saveRDS(dt_filt_BL, file="data/phs000834/dt_filt.Rds") 
saveRDS(dt_siteStats %>% ungroup(), file="data/phs000834/dt_siteStats.Rds")


