#Figure 5a plot

library(tidyverse)
library(circlize)

here::here()

areaXnType_diffType <- readRDS('data/stat_tests/binom_testing/areaXnType_diffType.Rds')

dim(areaXnType_diffType) #853

areaXnType_diffType %>% count(siteID) #239
areaXnType_diffType %>% filter(p.adj.Fisher < 0.05) %>% count(siteID) #124

areaXnType_diffType %>% filter(!is.na(gene_name), p.adj.Fisher< 0.05) %>% count(siteID)

areaXnType_diffType %>% arrange(desc(p.adj.Fisher))

areaXnType_diffType %>% count(siteID)

areaXnType_diffType %>% filter(!is.na(gene_name), p.adj.Fisher< 0.05) %>% 
  group_by(siteID) %>% count(gene_name,sort=TRUE) %>% count(siteID,sort=TRUE)


areaXnType_diffType %>% filter(description!="NULL") %>% 
  filter(p.adj.Fisher < 0.05) %>% select(1:p.adj.Fisher) %>% distinct() #353


areaXnType_diffType %>% filter(description!="NULL", p.adj.Fisher< 0.05) %>% 
  select(-ENSGID,-gene_name) %>% distinct() %>% arrange(p.adj.Fisher) %>% 
  count(diffType,diffArea, sort=TRUE) %>% 
  mutate(sumN= sum(n)) %>% 
  mutate(propN = n/sum(n))

#### chorddiag edges indicate 'flow strength between each pair of nodes'

myCols <- c('blue','red','orange','yellow','pink','light blue')

areaXnType_diffType %>% filter(p.adj.Fisher < 0.05) %>% count(diffType)

areaXnType_diffType %>% filter(p.adj.Fisher < 0.05) %>% count(area1,area2,diffType,sort=TRUE)

#Ex-Ex dEd
areaXnType_diffType %>% filter(description!="NULL", p.adj.Fisher< 0.05) %>% 
  select(-ENSGID,-gene_name) %>% distinct() %>% 
  filter(diffType=='Ex_Ex') %>% select(-diffType) %>% 
  count(area1,area2) %>% 
  chordDiagramFromDataFrame(grid.col=myCols)

#Ex-In dEd
areaXnType_diffType %>% filter(description!="NULL", p.adj.Fisher< 0.05) %>% 
  select(-ENSGID,-gene_name) %>% distinct() %>% 
  filter(diffType=='Ex_In') %>% select(-diffType) %>% # dim() #194
  count(area1,area2) %>% 
  chordDiagramFromDataFrame(grid.col=myCols,self.link = TRUE) 

#In-In dEd
areaXnType_diffType %>% filter(description!="NULL", p.adj.Fisher < 0.05) %>% 
  select(-ENSGID,-gene_name) %>% distinct() %>% 
  filter(diffType=='In_In') %>% select(-diffType) %>% # dim() 33
  count(area1,area2) %>% 
  chordDiagramFromDataFrame(grid.col=myCols)



