##locus 15q11 (Prader-Willi snoRNA host genes)


library(here)
library(tidyverse)

#Explore snoRNA cluster in 15q11


####### ####### ####### ####### ####### #######
####### #######    READ DATA    ####### ####### 
####### ####### ####### ####### ####### #######

#snoRNA targets (in 5HT receptor transcripts)

# from Jorjani et al NAR 2016 doi: 10.1093/nar/gkw386; NB original data is Hg19

snoRNA_db <- readRDS("data/related_studies/Jorjani_etal_NAR/snoRNA_db.Rds")

BL_joinKey <- readRDS('data/phs000834/BL_metadata_highQualCells.Rds')

anno <- readRDS('data/expn_processing/GRCh38_genes.Rds')

dt_siteStats_filt <- readRDS('data/phs000834/dt_siteStats_filt.Rds')

#lm test results
lm_out_INT_neuType   <-  readRDS('data/stat_tests/broom_lm/lm_out_INT_neuType.Rds') %>% 
  ungroup() %>% mutate(FDR=p.adjust(p.value,method="BH")) %>% arrange(FDR)


####### ####### ####### ####### ####### #######
####### ####### ####### ####### ####### #######


snoRNA_db %>% left_join(anno,by=c('Name'='gene_name')) %>% 
  count(is.na(gene_biotype)) #only 50% can be joined to Hg38.

snoRNA_db %>% filter(Chromosome=='chr15') 

snoRNA_db %>% filter(Chromosome=='chr15') %>% 
  filter(str_detect(Name,'SNORD115|SNORD116')) %>% 
  mutate(base_name=str_remove(Name,'SNORD'),
         base_name=str_extract(base_name,'[^-]+')) %>% 
  count(`Host gene`,base_name) 


#########################
### test 5HT editing ###
#########################

#5HT2C : ENSG00000147246

dt_siteStats_filt %>% filter(ENSGID=="ENSG00000147246")  #no edited sites detected in this transcript

anno %>% filter(gene_id=="ENSG00000147246") %>% select(region,start,end) #sense strand (+)



########## ########## ########## ########## ########## ########## ##########

########## ########## ########## ########## ########## ########## ##########


#consider snoRNAs

snoRNA_db %>% filter(Chromosome=='chr15') %>% 
  filter(str_detect(Name,'SNORD115|SNORD116')) %>% 
  mutate(base_name=str_remove(Name,'SNORD'),
         base_name=str_extract(base_name,'[^-]+')) %>% 
  count(`Host gene`,base_name) 


snoRNA_db %>% filter(Chromosome=='chr15') %>% 
  filter(str_detect(Name,'SNORD115|SNORD116')) %>% select(1:6) %>% mutate(Chromosome=str_remove(Chromosome,'chr')) %>% 
  left_join(anno, by=c('Chromosome'='region',"Name" = 'gene_name')) %>% 
  mutate(diff = `Gene start` - start) %>% ggplot(aes(diff)) + geom_histogram(bins=50)
#must be hg19 vs hg38.

#which are edited vs non-edited in current data; and related studies (RADAR, Tan et al, Breen et al, Tran et al)

ed_snord_snhg <- dt_siteStats_filt %>% 
  filter(str_detect(gene_name,'SNRPN|SNORD|SNHG14|UBE3A')) 

ed_snord_snhg %>% select(siteID,region,gene_name,site_status) %>% distinct() %>% count(gene_name,region,site_status)

ed_snord_snhg %>% filter(gene_name=='SNHG14') %>% count(site_type,sort=TRUE)


anno %>% filter(region=="15", gene_name %in% c('SNRPN','UBE3A','SNHG14','SNHG15')) %>%
  rowwise() %>% mutate(forY=rnorm(1)) %>% 
  select(gene_name,region,start,end,forY) %>% 
  gather(key,value,start,end) %>% 
  ggplot(aes(x=value,y=forY, group=gene_name)) + geom_line(aes(col=gene_name))



#extract beta for y val
lm_out_INT_neuType %>% 
  filter(term == 'log10(value)') %>% 
  left_join(anno %>% select(gene_id,gene_name) %>% distinct(), 
            by=c('ENSG'='gene_id')) %>% 
  mutate(gene_family = case_when(str_detect(gene_name,'SNORD115') ~ 'SNORD115 cluster',
                                 str_detect(gene_name,'SNORD116') ~ 'SNORD116 cluster',
                                 TRUE ~ 'Background')) %>% 
  ggplot(aes(x=estimate,fill=factor(gene_family), col=factor(gene_family))) +
  geom_density(alpha = 0.25) + # facet_wrap(~goi,scales='free_y',ncol=1)
  xlim(-0.03,0.03) + scale_y_reverse() + coord_flip()

ggsave('charts/SNORD115_116_beta_density.pdf',width=5,height=5)

SNORD115_beta <- lm_out_INT_neuType %>% 
  filter(term == 'log10(value)') %>% 
  left_join(anno %>% select(gene_id,gene_name,region,start,end) %>% distinct(), 
            by=c('ENSG'='gene_id')) %>% 
  filter(str_detect(gene_name,'SNORD115'))

SNORD115_116_beta <- lm_out_INT_neuType %>% 
  filter(term == 'log10(value)') %>% 
  left_join(anno %>% select(gene_id,gene_name,region,start,end) %>% distinct(), 
            by=c('ENSG'='gene_id')) %>% 
  filter(str_detect(gene_name,'SNORD115|SNORD116'))


#what is the beta for SNHG14 and SNRPN?
SNRPN_SNHG14_beta <- lm_out_INT_neuType %>% 
  filter(term == 'log10(value)') %>% 
  left_join(anno %>% select(gene_id,gene_name,region,start,end) %>% distinct(), 
            by=c('ENSG'='gene_id')) %>% 
  filter(gene_name %in% c('SNHG14','SNRPN','UBE3A'))



#plot SNORDs relative to genes. All coords are Hg38
#y for SNORDs = ed effect (beta); add editing cluster in snhg14

anno %>% filter(region=="15", str_detect(gene_name, 'SNRPN|UBE3A|SNHG14|SNHG15')) %>% arrange(end) %>% 
  rowwise() %>% mutate(forY=rnorm(1)) %>% 
  select(gene_name,region,start,end,forY) %>% 
  gather(key,value,start,end) %>% 
  ggplot() + 
  geom_histogram(data = ed_snord_snhg %>% filter(position<2.9e7 & position>2.4e7), bins=50, aes(x=position, fill=site_type)) +
  geom_line(aes(x=value,y=0, lwd=forY, group=gene_name, col=gene_name)) +
  geom_point(data = SNORD115_beta %>% filter(start<2.9e7) ,
             aes(x=start,y=estimate), cex=2)


ed_snord_snhg %>% 
  filter(position > min(SNORD115_beta$start), position< max(SNORD115_beta$start)) %>% 
  ggplot() + geom_histogram(aes(x=position,fill=site_type) , bins=50) +
  geom_point(data = SNORD115_beta %>% filter(start<2.9e7) ,
             aes(x=start, y=estimate), cex=2)


gene_boundaries_beta <- anno %>% filter(region=="15", str_detect(gene_name, 'SNRPN|UBE3A|SNHG14|SNHG15')) %>% arrange(end) %>% 
  rowwise() %>% mutate(forY=rnorm(1)) %>% 
  select(gene_name,region,strand,start,end,forY) %>% 
  gather(key,value,start,end) %>% 
  left_join(SNRPN_SNHG14_beta,by='gene_name')


####### ####### ####### #######
####### ####### ####### #######

#Expand width

ed_snord_snhg %>% 
  filter(position>2.4e7, position<25.4e6) %>% 
  ggplot() + 
  geom_hline(yintercept=0, lty=2, col='grey') + 
  geom_segment(aes(x=position, y=0.0025, xend=position, yend = -0.0025, col = site_type)) +
  geom_line(data=gene_boundaries_beta  %>% filter(gene_name!="UBE3A"), 
            col='grey',
            aes(x=value, y=estimate, group=gene_name), lwd=2, alpha=0.75,
            arrow = arrow(angle = 45, length=unit(0.35,'cm'), ends = "last", type = "closed")) + 
  geom_point(data = SNORD115_116_beta %>% filter(start<2.9e7) %>% mutate(stat = ifelse(FDR<0.05, 1, 0),
                                                                         shape = ifelse(str_detect(gene_name,'115'),1,0)),
             aes(x=start,y=estimate,shape=factor(shape)), cex=1, show.legend = FALSE) + 
  geom_text(data = tribble(~gene_name, ~ xVal, ~ yVal, 
                           "SNRPN", 24.9e6,-0.007, 
                           #"SNHG14", 25.15e6, -0.005,
                           "SNHG14", 25.0e6, 0.0075,
                           'SNORD115 cluster',25.21e6, 0.0325,
                           'SNORD116 cluster',25.075e6, 0.0325),
            aes(x=xVal,y=yVal,label=gene_name)) + 
  ylab('Effect size of correlation \n with global editing (beta)') +
  xlab('Locus 15q11.2 (BP) ') + ylim(-0.01,0.04) + 
  scale_color_manual(values=c('blue','light blue', 'orange','pink','green','dark green'))

ggsave('charts/15q11.2_chart.pdf', width=10, height=4)


######### ######### ######### ######### ######### ######### #########
######### ######### ######### ######### ######### ######### #########





