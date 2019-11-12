#cell-wise editing proportion testing

#devtools::install_github("tidymodels/broom")

library(tidyverse)

here::here()

BL_joinKey <- readRDS('data/phs000834/BL_metadata_highQualCells.Rds')

mapping_stats_BL <- readRDS("data/mapping_output/mapping_stats_BL.Rds")

dt_filt <- readRDS("data/phs000834/dt_filt.Rds")

dt_siteStats_TDjoin <- readRDS("data/phs000834/dt_siteStats_TDjoin.Rds")


#samtools bqsr depth; contains candidate editing sites in each neuron covered by at least 5 reads.
gte5_DPsites <- read_tsv("data/samtools_depth_output/samDepth_BQSR_sitesGTE5.out", col_names = FALSE,
                         col_types = cols('X2'=col_character())) %>%
  magrittr::set_colnames(c('sample','siteID','depth')) %>%
  filter(sample %in% BL_joinKey$SRA) %>% 
  filter(siteID %in% dt_siteStats_TDjoin$siteID)


#plot covered sites per cell
gte5_DPsites %>% count(sample) %>% ggplot(aes(x=n)) + geom_histogram(fill='dodger blue')
gte5_DPsites %>% count(sample) %>% summary()

#plot n. cells per covered site
gte5_DPsites %>% count(siteID) %>% ggplot(aes(x=n)) + geom_histogram(fill='dodger blue')
gte5_DPsites %>% count(siteID) %>% summary()


#### ### ### ### ### ### ### ### ### ### ###
### ### ### GLOBAL EDITING INDEX ### ### ###
#### ### ### ### ### ### ### ### ### ### ###

### Calculate cell-wise editing as a proportion of covered sites ###

cell_edProp <- left_join(gte5_DPsites %>% count(sample),  #n sites covered ≥ 5 high-quality (de-duplicated BQSR) reads per cell.
                         dt_filt %>% filter(siteID %in% dt_siteStats_TDjoin$siteID) %>% count(sample), by='sample',
                         suffix=c('_totalCov',"_ed")) %>% 
  mutate(edProp = n_ed/n_totalCov) 


saveRDS(cell_edProp, "data/stat_tests/GEI.Rds")


cell_edProp %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  filter(neuType!="NoN") %>% 
  ggplot(aes(x=area,y=edProp, fill=neuType,col=neuType)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8), cex=0.25) +
  geom_boxplot(alpha=0) +
  ylab('Proportion of covered sites edited') + xlab('Brodmann area')


#### Neuronal GroupID ####

#uniquely mapped reads
mapping_stats_BL %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  filter(neuType!="NoN") %>% 
  ggplot(aes(x=Group_ID,y=Uniquely_mapped_reads_number, fill=neuType,col=neuType)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8), cex=0.25) +
  geom_boxplot(alpha=0) 

#GEI
# add combined Ex and In groups


## FIG 4c
cell_edProp %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  filter(neuType!="NoN") %>% 
  mutate("Collapsed" = neuType) %>% 
  gather(key,value,Group_ID,`Collapsed`) %>% 
  mutate(key=factor(key,levels=c('Group_ID','Collapsed'))) %>% 
  ggplot(aes(x=value,y=edProp, fill=neuType,col=neuType)) + 
  geom_point(position=position_jitterdodge(jitter.width = 0.2,dodge.width = 0.8), cex=0.25) +
  geom_boxplot(alpha=0)  +
  #facet_wrap(~key)
  ggforce::facet_row(vars(key), scales='free_x',space='free') +
  xlab('Neuronal phenotype') + ylab('Proportion of sites edited') +
  theme_grey() + theme(legend.position = "NONE")
ggsave('charts/Figure_4c.pdf', width=7,height=3.5)




## FIG 4b
cell_edProp %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  filter(neuType!="NoN") %>% mutate(`Neuronal type`=neuType) %>% 
  mutate(area=factor(area,levels=paste0('BA', c(8,10,21,22,41,17)))) %>% 
  ggplot(aes(x=area,y=edProp, fill=`Neuronal type`,col=`Neuronal type`)) + 
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width = 0.15,dodge.width = 0.8),alpha=0.5, cex=0.25) +
  geom_boxplot(alpha=0) +
  xlab('Brodmann area') + ylab('Proportion of sites edited') +
  theme(legend.position='Null')

ggsave('charts/4b.pdf',width=7,height=3.5)


#Calculate mean cell type proportions across cortical regions

BL_joinKey %>% group_by(area) %>% filter(neuType !="NoN") %>% count(neuType) %>% 
  mutate(prop=n/sum(n)) %>% group_by(neuType) %>% summarize(meanProp = mean(prop))


#Count neuronal type per cortical region


p1 <- cell_edProp %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  filter(neuType!="NoN") %>% mutate(`Neuronal type`=neuType) %>% 
  mutate(area=factor(area,levels=paste0('BA', rev(c(8,10,21,22,41,17))))) %>% 
  ggplot(aes(x=area,y=edProp, fill=`Neuronal type`,col=`Neuronal type`)) + 
  geom_point(position=position_jitterdodge(jitter.height=0, jitter.width = 0.15,dodge.width = 0.8),alpha=0.5, cex=0.25) +
  geom_boxplot(alpha=0) +
  xlab('Brodmann area') + ylab('Proportion of sites edited') +
  theme(legend.position='bottom') +
  coord_flip()


source('code/colour_palettes.R')
col_mat_full <- matrix(rep(c(drsimonj_pal('hot')(8),drsimonj_pal('cool')(8)),16),ncol=16,nrow=16)

p2 <- cell_edProp %>% 
  left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  mutate(area=factor(area,levels = paste0('BA', c(8,10,21,22,41,17)))) %>% 
  filter(neuType!="NoN") %>% mutate(`Neuronal type`=neuType) %>% 
  group_by(area,neuType,Group_ID) %>% summarize("N. nuclei"=n()) %>% 
  ggplot(aes(x=neuType,y=`N. nuclei`)) + geom_col(aes(fill=Group_ID),position='stack',col='white') +
  xlab('') +
  coord_flip() + 
  theme(legend.position = 'bottom') +
  scale_fill_manual(values=col_mat_full) +
  ggforce::facet_col(vars(area), scales='free',space='free') 


cowplot::plot_grid(p1 + theme_grey(), p2 + theme_grey()  + theme(legend.position='bottom'),
                   labels='AUTO',align='h',axis='b')

Fig4a <- p1 + theme_grey() + theme(legend.position='bottom')

ggsave(plot=Fig4a, 'charts/Fig4a.pdf',width=4,height=7)

Fig4b <- p2 + theme_grey() + theme(legend.position='bottom')

ggsave(plot=Fig4b, 'charts/Fig4b.pdf',width=4,height=7)


###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ######



broom::tidy(summary(lm(map_sum ~ factor(neuType), 
                       data = mapping_stats_BL %>% 
                         left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
                         filter(neuType!="NoN") )))

#No significant difference in read mapping numbers for excitatory vs inhibitory neurons

broom::tidy(summary(
  lm(map_sum ~ factor(area), 
     data = mapping_stats_BL %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN") %>% 
       mutate(map_sum = Uniquely_mapped_reads_number + Number_of_reads_mapped_to_multiple_loci)))) %>% 
  arrange(p.value)
#Significant differences in read total mapping per cortical region.


broom::tidy(summary(
  lm(map_sum ~ factor(Group_ID), 
     data = mapping_stats_BL %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN") %>% 
       mutate(map_sum = Uniquely_mapped_reads_number + Number_of_reads_mapped_to_multiple_loci)))) %>% 
  arrange(p.value)
#Significant differences in total read mapping per neuronal sub-group: Ex2, In6, In7 and In8 are different to Ex1.


#### TEST GEI by neuronal sub-type (Group_ID)  ####

broom::tidy(summary(
  lm(edProp ~ factor(Group_ID), 
     data = cell_edProp %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN")))) %>% 
  arrange(p.value) %>% mutate(FDR=p.adjust(p.value,method="BH")) %>% 
  filter(p.value > 0.05)

#all sub-groups EXCEPT for Ex8 and In3 have GEI significantly different to Ex1.


#some relationship between mapping stats / total coverage and editing proportion
mapping_stats_BL %>% left_join(cell_edProp,by='sample') %>%
  ggplot(aes(x=map_sum,y=edProp)) + geom_hex(bins=100) +
  geom_smooth(method="lm")

mapping_stats_BL %>% left_join(cell_edProp,by='sample') %>%
  ggplot(aes(x=n_totalCov, y=edProp)) + geom_hex(bins=100) +
  geom_smooth(method="lm")


#test collapsed neuType controlling for read coverage
broom::tidy(summary(
  lm(edProp ~ factor(neuType) + n_totalCov, #  + map_sum + Number_of_reads_mapped_to_multiple_loci
     data = cell_edProp %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN") %>% 
       left_join(mapping_stats_BL,by='sample')))) %>% 
  arrange(p.value) %>% mutate(FDR=p.adjust(p.value,method="BH"))


#control for n sites covered [and/or sum of mapped reads]
#test Group_ID
lm_Group_ID_GEI <- broom::tidy(summary(
  lm(edProp ~ factor(Group_ID) + n_totalCov , # + map_sum
     data = cell_edProp %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN") %>% 
       left_join(mapping_stats_BL,by='sample')))) %>% 
  arrange(estimate) %>% mutate(FDR=p.adjust(p.value,method="BH")) 

lm_Group_ID_GEI %>% filter(FDR>0.05)

#GEI differences remain after controlling for site coverage.


#test cortical area ± total coverage
lm_area_GEI <- broom::tidy(summary(
  lm(edProp ~ factor(area), 
     data = cell_edProp %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN") %>% 
       left_join(mapping_stats_BL,by='sample')))) %>% 
  arrange(estimate)

lm_area_GEI_ctrlCov <- broom::tidy(summary(
  lm(edProp ~ factor(area)  + n_totalCov, 
     data = cell_edProp %>% 
       left_join(BL_joinKey, by = c('sample'='SRA')) %>% 
       filter(neuType!="NoN") %>% 
       left_join(mapping_stats_BL,by='sample')))) %>% 
  arrange(estimate)


rbind(lm_area_GEI %>% mutate(covar='none'),
      lm_area_GEI_ctrlCov %>% mutate(covar = 'n_totalCov')) %>% 
  select(term,estimate,covar) %>% 
  spread(covar,estimate) %>% filter(term!='(Intercept)', term!='n_totalCov') %>% 
  ggplot(aes(x=none,y=n_totalCov,col=term)) + geom_point() + geom_smooth(method='lm') +
  geom_hline(yintercept = 0)+ geom_vline(xintercept = 0)

#only BA22 sign changes.


#for STables

saveRDS(lm_area_GEI, 'data/stat_tests/lm_area_GEI.Rds')

saveRDS(lm_Group_ID_GEI, 'data/stat_tests/lm_Group_ID_GEI.Rds')



###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ######
###### ###### ###### ###### ###### ###### ###### ######


#Do inhibitory neurons have greater editing in shared genes; or editing in e.g. inhib.-specific genes?



dt_depth <- left_join(gte5_DPsites, dt_filt %>% filter(siteID %in% gte5_DPsites$siteID), by=c('siteID','sample')) %>%
  left_join(BL_joinKey %>% select(SRA,neuType,area), by=c('sample'='SRA'))

# saveRDS(dt_depth, 'data/phs000834/dt_depth.Rds')
# dt_depth <- readRDS('data/phs000834/dt_depth.Rds')

dt_depth_summ <- dt_depth %>% filter(neuType!="NoN") %>% group_by(sample, neuType) %>% summarize(meanDP = mean(depth),medDP = median(depth)) 


dt_depth_summ %>% ggplot(aes(x=medDP)) + geom_histogram(aes(fill=neuType)) + 
  facet_wrap(~ neuType, ncol = 1, scales='free_y') +
  xlab('Median editing site coverage')

dt_depth_summ %>% group_by(neuType) %>% summarize(meanMED = mean(medDP))
broom::tidy(lm(medDP ~ factor(neuType), data = dt_depth_summ))

# --> Inhibitory neurons have an averaeg median edSite coverage 0.5 units greater than excitatory neurons.

#sample 500 cellsof each neuronal type
set.seed(1234); sample_Ex_In <- dt_depth %>% 
  select(neuType, sample) %>% distinct() %>% 
  filter(neuType!="NoN") %>% 
  group_by(neuType) %>% sample_n(500)  


dt_depth %>% filter(str_detect(siteID,"1_")) %>%  #sites on chromosome 1  
  filter(!is.na(alt_al)) %>% 
  filter(sample %in% sample_Ex_In$sample) %>% 
  group_by(siteID) %>% mutate(nCells = n()) %>% ungroup() %>% 
  select(nCells, everything()) %>% arrange(desc(nCells)) %>% #--> exclude marker genes
  filter(nCells > 200) %>% count(neuType)

#For inhibitory neurons there are 82 sites in chr1 in > 200 cells; excitatory have 160 sites in chr1 in >200 cells.



#Do In. neurons have more sites than Ex., or greater editing per site?

dt_depth %>% count(sample, neuType)  %>% 
  filter(neuType!="NoN") %>% 
  ggplot(aes(x=n)) + 
  geom_histogram(aes(fill=neuType),bins=50)

dt_depth %>% 
  filter(!is.na(alt_al)) %>% 
  count(sample, neuType)  %>% 
  filter(neuType!="NoN") %>% 
  ggplot(aes(x=n)) + 
  geom_histogram(aes(fill=neuType),bins=50)

#Ex have more editing sites covered, and more edited sites, than Inhibitory...


dt_depth %>% select(sample,siteID,depth,alt_al,neuType) %>% 
  filter(neuType!="NoN") %>% 
  mutate(status = ifelse(is.na(alt_al),0,1)) %>% 
  group_by(neuType) %>% summarize(meanEd = mean(status))

#But In neurons have more editing as a _proportion_ of sites covered.

#N. and % of sites are restricted to In / Ex neurons

dt_depth %>% filter(neuType!="NoN") %>% count(neuType) 
dt_depth %>% filter(neuType!="NoN") %>% group_by(neuType) %>% count(siteID) %>% 
  arrange(siteID) %>%
  ungroup() %>% count(siteID,sort=TRUE) %>% count(n)

#Only 253 sites were unique one neuronal type. 99.4% sites detected in both cell types.
1-(253/sum(253,40608))

### ### ### ### ### ### ### ### ###

# If ensemble averaging is suppressing novel sites, we expect those sites to be expressed in small numbers of cells.


dt_siteStats_filt <- readRDS('data/phs000834/dt_siteStats_filt.Rds')

dt_siteStats_filt %>% select(site_type,n_Cells) %>% 
  mutate(site_status = ifelse(str_detect(site_type,'Novel'),'Novel','Reported')) %>% 
  ggplot(aes(x=site_type,y=log10(n_Cells))) + 
  geom_boxplot(aes(col=site_type),outlier.alpha = 0) +
  geom_jitter(width=0.05,cex=0.1,alpha=0.25, aes(col=site_type)) + coord_flip() +
  scale_color_brewer(palette='Paired') + xlab('Site type') +
  theme(legend.position = 'None')
ggsave('charts/Ncells_per_site_bySiteType.pdf',width=6,height=5)  


#### #### #### #### #### #### #### #### 

dt_siteStats_filt %>% select(site_type,n_Cells) %>% 
  mutate(site_status = ifelse(str_detect(site_type,'Novel'),'Novel','Reported')) %>% 
  mutate(log_nCells = log10(n_Cells),
         site_status = factor(site_status)) %>% 
  do(broom::tidy(lm(log_nCells ~ site_status, data = .)))

dt_siteStats_filt %>% select(site_type,n_Cells) %>% 
  mutate(site_status = ifelse(str_detect(site_type,'Novel'),'Novel','Reported')) %>% 
  mutate(log_nCells = log10(n_Cells),
         site_status = factor(site_status)) %>% 
  do(broom::tidy(lm(log_nCells ~ site_type, data = .)))

#Indeed, previously reported nonRep and rep-nonAlu sites are edited in significantly more cells than Alu sites; \
# whereas novel sites in each class are detected in significantly fewer cells on average than reported sites.



