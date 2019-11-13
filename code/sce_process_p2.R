
#############################
######## SCE object #########
#############################

library(tidyverse)

library(caret)
library(scater)
library(scran) 


library(iSEE)
theme_set(theme_classic())


library(annotate)
library(AnnotationHub)
library(org.Hs.eg.db)

##### READ IN DATA ######

BL_joinKey <- readxl::read_xlsx("fromBL_Suppl_Tables_v5.xlsx",
                                sheet=2, skip=5) %>% 
  rename_all(funs(str_replace_all(., ' ','_'))) %>% 
  rename_all(funs(str_replace_all(., '-','_'))) %>% 
  dplyr::mutate(neuType=ifelse(str_detect(Sub_Group_ID,'^Ex'),'Ex','In'),
         neuType=ifelse(str_detect(Group_ID, "NoN"),'NoN', neuType)) %>% 
  dplyr::mutate(area=str_extract(Sub_Group_ID,"BA.*")) %>% 
  dplyr::filter(!is.na(SRA))
names(BL_joinKey) <- str_remove(names(BL_joinKey),'†')

saveRDS(BL_joinKey, 'BL_metadata.Rds')


#created in sce_p1.R

sce <- readRDS(file='expn_processing/Lake_sce.Rds')


## CLUSTERING ##
sce <- runPCA(sce)
reducedDim(sce,'PCA')[1:5, 2] 


pca1 <- reducedDim(sce,'PCA')[,1]
pca2 <- reducedDim(sce,'PCA')[,2]

#plot PCA

tibble('pc1'=pca1,'pc2'=pca2) %>% 
  ggplot(aes(x=pc1,y=pc2)) + geom_point()

library(uwot)

sce <- scater::runUMAP(sce, use_dimred = NULL, n_neighbors = 20)

plot_dat <- tibble('sample'=rownames(reducedDim(sce,"UMAP")),
                   'dim1'=reducedDim(sce, "UMAP")[, 1],
                   'dim2'=reducedDim(sce, "UMAP")[, 2])

ggplot(plot_dat, aes(x = dim1, y = dim2)) +
  geom_point( alpha = 0.5,cex=0.5) +
  scale_colour_brewer(palette = "Set1")
#2989 cells plotted 



#### COLOR by Blue Lake cell phenotype designations [BL_joinKey] ####


#join assigned phenotype

plot_dat %>% 
  left_join(BL_joinKey,by=c('sample'='SRA')) %>% #count(is.na(`C1 Protocol`))
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point( aes(col=C1_Protocol), alpha = 0.5,cex=0.5) 


plot_dat %>% 
  left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point( aes(col=Group_ID), alpha = 0.5,cex=0.5) 

plot_dat %>% 
  left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  dplyr::filter(Group_ID!="NoN") %>% 
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point( aes(col=Group_ID), alpha = 0.5,cex=0.5, show.legend=TRUE) + 
  theme_grey() + theme(legend.position = 'bottom')


##### ##### ##### ##### #####

plot_dat %>% 
  left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point( aes(col=neuType), alpha = 0.5,cex=0.5) +
  theme_grey()


plot_dat %>% 
  left_join(BL_joinKey,by=c('sample'='SRA')) %>% dplyr::rename(BA=area) %>% 
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point( aes(col=BA),cex=0.5) + 
  scale_color_manual(values=c('blue','orange','grey','yellow','navy','dark green')) +
  theme_grey() 
ggsave("charts/sceCluster_byArea.pdf",width=6.5,height=4)

#Clustering by inhibitory vs exitatory neuron Ids.



# Assess apoptotic nuclei
## calculates mean log expression of genes in GO apoptosis pathway.


go_entrez_id <- get("GO:0006915", revmap(org.Hs.egGO))
apoptotic_genes <- unlist(mget(go_entrez_id , org.Hs.egENSEMBL))

index_apoptotic <- match(apoptotic_genes, rownames(sce))
index_apoptotic <- index_apoptotic[!is.na(index_apoptotic)]
plot_dat$apoptotic <- colMeans(logcounts(sce)[index_apoptotic, ])

ggplot(plot_dat, aes(x = dim1, y = dim2)) +
  geom_point(aes(col = apoptotic), cex=0.5, alpha = 0.5) +
  scale_color_viridis_c()

plot_dat %>% dplyr::count(apoptotic>3)


# Cluster annotation
#BiocManager::install("scran", dependencies=TRUE)

library(scran)


snn_gr <- buildSNNGraph(sce, use.dimred = "PCA", k = 50) 

clusters <- igraph::cluster_louvain(snn_gr)
table(clusters$membership) 

sce$Cluster <- factor(clusters$membership)
plot_dat$Cluster <- sce$Cluster

plot_dat %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point(aes(col = Cluster), cex=0.5, alpha = 0.5)

plot_dat %>% left_join(BL_joinKey,by=c('sample'='SRA')) %>% 
  ggplot(aes(x = dim1, y = dim2)) +
  geom_point(aes(col = Group_ID), cex=0.5, alpha = 0.5)


saveRDS(plot_dat, "expn_processing/cluster_data.Rds")








