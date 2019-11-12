
#assisted by Jacob Munro

library(tidyverse)
library(here)

here()

#tables created in [Fig1_code.R]

F1a_table <- readRDS('data/figure_data/F1a_table.Rds')
F1b_table <- readRDS('data/figure_data/F1b_table.Rds')
F1c_alt_table <- readRDS('data/figure_data/F1c_alt_table.Rds')
F1d_table <- readRDS('data/figure_data/F1d_table.Rds')


#Figure 1.


F1a <- F1a_table %>% 
  ggplot(aes(x=n)) + geom_histogram(fill='dodger blue',bins=50) +
  xlab('N. sites') + ylab('N. cells')


F1b <- F1b_table %>% 
  ggplot(aes(x=n_Cells_Lake)) + geom_histogram(bins=50, fill="dodgerblue") +
  xlim(0, 250) +
  xlab('N. cells') + ylab('N. sites')


cowplot::plot_grid(F1a,F1b,nrow=1, labels='AUTO')

ggsave(plot=F1a,'charts/F1a.pdf', width=3,height=2.5)
ggsave(plot=F1b,'charts/F1b.pdf', width=3,height=2.5)


###### ###### ###### ###### ###### ######

#Figure 1c & d.

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


#These charts are be limited to 83% of sites covered by a single non-OL gene region;

F1c_alt <- F1c_alt_table %>% ggplot(aes(x=`ENSG Bio-type`)) + 
  geom_bar(col='black', lwd=0.1, aes(fill=site_type)) + 
  geom_bar(alpha=0, aes(col=lineBound),lwd=0.5, show.legend = FALSE ) + 
  ylab('N. sites') + 
  ggforce::facet_col(vars(gene_type), scales='free',space='free')+
  coord_flip() +
  guides(fill = guide_legend(reverse = T)) +
  scale_fill_brewer(palette='Paired') +
  scale_color_manual(values=rep('black', 4)) +
  theme(legend.position = 'bottom', legend.key.size = unit(0.3, "cm"), 
        legend.title = element_text(size=9),legend.text=element_text(size=8))



F1d <- F1d_table %>%   
  ggplot(aes(x=site_type)) + 
  geom_bar(col='black', lwd=0.2, aes(fill=`Genic feature`)) + 
  xlab('Site location') + ylab('N. sites')  + 
  guides(fill = guide_legend(reverse=T)) + 
  scale_fill_brewer() +
  theme(legend.position = 'bottom', legend.key.size = unit(0.3, "cm"), 
        legend.title = element_text(size=9),legend.text=element_text(size=8)) + 
  coord_flip() +
  facet_wrap(~ tag, scales='free_y', ncol=1)

## ############################## ##

## JM Suggestion ##


F1c_fill <- F1c_alt_table %>% 
  group_by(`ENSG Bio-type`) %>% mutate(nType = n()) %>% 
  filter(nType>=5) %>% 
  ggplot(aes(x=reorder(`ENSG Bio-type`,nType))) + 
  geom_bar(col = 'black', lwd = 0.1, aes(fill = site_type), position='fill') + 
  geom_bar(alpha = 0, aes(col = lineBound),lwd=0.5, show.legend = FALSE, position='fill' ) + 
  ylab('Proportion of sites') + xlab('Biotype') + 
  coord_flip() +
  guides(fill = guide_legend(reverse = T)) +
  scale_fill_brewer(palette='Paired') +
  scale_color_manual(values=rep('black', 4)) +
  theme(legend.position = 'bottom', legend.key.size = unit(0.3, "cm"), 
        legend.title = element_text(size=9),legend.text=element_text(size=8))


F1c_count <- F1c_alt_table %>% 
  group_by(`ENSG Bio-type`) %>% summarize(nType = n()) %>% 
  filter(nType>=5) %>% 
  mutate(lnType = log10(nType)) %>% 
  ggplot(aes(x = reorder(`ENSG Bio-type`,lnType), y=lnType)) + geom_col() + coord_flip() +
  ylab('log10(N. sites)') + xlab('') + theme(axis.text.y = element_blank())

cowplot::plot_grid(F1c_fill, F1c_count, align = 'hv', axis='b')

ggsave('charts/F1c_forSplit.pdf',width=10,height=3)


####### ####### ####### #######



F1d_fill <- F1d_table %>% 
  group_by(`Genic feature`) %>% mutate(nFeature = n()) %>% 
  mutate(lineBound=case_when(site_type=='Alu' | site_type=='Alu_Novel' ~ 'A',
                             site_type=='nonRep' | site_type=='nonRep_Novel'~ 'B',
                             site_type=='rep_nonAlu' | site_type=='rep_nonAlu_Novel'~ 'C',
                             TRUE ~ 'NULL')) %>% 
  ungroup() %>% 
  mutate(`Genic feature` = factor(`Genic feature`, levels = c('three_prime_utr','stop_codon','intronic','exonic','five_prime_utr'))) %>% 
  ggplot(aes(x=`Genic feature`)) + 
  geom_bar(col = 'black', lwd = 0.1, aes(fill = site_type), position='fill') + 
  geom_bar(alpha = 0, aes(col = lineBound),lwd=0.5, show.legend = FALSE, position='fill' ) + 
  ylab('Proportion of sites') + xlab('Genic feature') +  coord_flip() +
  guides(fill = guide_legend(reverse = T)) +
  scale_fill_brewer(palette='Paired') +
  scale_color_manual(values=rep('black', 4)) +
  theme(legend.position = 'bottom', legend.key.size = unit(0.3, "cm"), 
        legend.title = element_text(size=9),legend.text=element_text(size=8))


F1d_count <- F1d_table %>% 
  group_by(`Genic feature`) %>% mutate(nFeature = n()) %>% 
  ungroup() %>% 
  mutate(`Genic feature` = factor(`Genic feature`, levels = c('three_prime_utr','stop_codon','intronic','exonic','five_prime_utr'))) %>% 
  mutate(lnFeature = log10(nFeature)) %>% 
  ggplot(aes(x=`Genic feature`,y=lnFeature)) + geom_col(position = 'identity') + 
  ylab('log10(N. sites)') +  xlab('') +
  coord_flip() +
  theme(axis.text.y = element_blank())

cowplot::plot_grid(F1d_fill, F1d_count, align = 'hv', axis='b') #
ggsave('charts/F1d_forSplit.pdf',width=10,height=3)


########### ########### ########### ########### ###########
########### ########### ########### ########### ###########

cowplot::plot_grid(F1c_fill, F1c_count, F1d_fill, F1d_count, align='hv', axis='b')



