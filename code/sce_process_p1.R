


#############################
######## SCE object #########
#############################

library(caret)
library(scater)
library(scran) 

library(DropletUtils) 
library(batchelor) #
library(umap) #

library(iSEE)

theme_set(theme_classic())

library(annotate)
library(AnnotationHub)
library(org.Hs.eg.db)

# if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
# 
#BiocManager::install("scran", dependencies = TRUE)
#BiocManager::install("ensembldb")
#BiocManager::install("uwot")
#BiocManager::install("edgeR")
#BiocManager::install("AnnotationHub")


library(SingleCellExperiment)


####### READ DATA #########


#created in sce_join_tables.R
BL_expn_table <- readRDS('expn_processing/BL_expn_table.Rds')


BL_expn_table <- as.data.frame(BL_expn_table)
rownames(BL_expn_table) <- BL_expn_table[,1]

samples <- colnames(BL_expn_table[,-1])


sce <- SingleCellExperiment(assays=list(counts = as.matrix(BL_expn_table[,-1])))


colData(sce)
rowData(sce)

sce$Sample_ID <- as.factor(rownames(colData(sce)))
rowData(sce)$geneID <- BL_expn_table[,1]

######### ######### #########
######### ANNOTATE ##########
######### ######### #########

#anno <- readRDS("data/ENSG_annotation.Rds")
#anno %>% dplyr::filter(gene_id=="ENSG00000284676") #AnnotationHub lacks information for several ENSG

library(tidyverse)

GRCh38_genes <- readRDS('data/expn_processing/GRCh38_genes.Rds')
GRCh38_genes$gene_name %>% head()


SCEgeneTable <- tibble('id' = rowData(sce)$geneID) %>% 
  left_join(GRCh38_genes, by = c('id' = 'gene_id')) 


SCEgeneTable %>% dplyr::filter(!complete.cases(.)) %>% dplyr::count(id)



length(SCEgeneTable$gene_name)  

rowData(sce)$Symbol <- SCEgeneTable$gene_name
#rowData(sce)$width <- SCEgeneTable$width #for feature_length in calculateTPM()

#counts matrix: assay(sce)


# IDENTIFY mitochondrial genes

is.mito <- grepl("^MT-", rowData(sce)$Symbol)


# Specify ERCC spikes

table(grepl("^ERCC", rownames(sce)))

is.spike <- grepl("^ERCC", rownames(sce))

isSpike(sce,"ERCC") <- which(is.spike) #58303:58394


# Library sizes
sizeFactors(sce) <- colSums(assay(sce))


sizeFactors(sce, "ERCC") <- colSums(assay(sce)[isSpike(sce, "ERCC"), ])

head(sizeFactors(sce))
head(sizeFactors(sce, "ERCC"))

## FILTERING

sce <- calculateQCMetrics(sce, feature_controls = list(Mt = is.mito, ERCC=is.spike))

as_tibble(colData(sce)$pct_counts_Mt) %>% 
  ggplot(aes(x=value)) + geom_histogram(bins=100,col="grey",fill="dodger blue")

high_mito <- sapply(levels(sce$Sample_ID), function(x)
  sum((colData(sce)$pct_counts_Mt > 15)[sce$Sample_ID == x]))

table(high_mito)
high_mitoCells <- tibble('high_mitoCells' = names(high_mito[high_mito==1]) ) 


rownames(colData(sce))[colData(sce)$pct_counts_Mt > 15]

#flag these cells:

toDrop <- colData(sce)$pct_counts_Mt > 15
table(is.na(toDrop))
toDrop[is.na(toDrop)] <- TRUE #replace 2 NA values with 'TRUE', which will result in them dropping out from data
table(toDrop)

#drop samples with mito RNA >15 or mitoRNA = 'NA' 
sce <- sce[, !toDrop]



#Filter by total RNA content and diversity
as_tibble(sce$total_counts) %>% 
  ggplot(aes(x=value)) + 
  geom_histogram(bins=100,col="grey",fill="dodger blue")


libsize_flag <- isOutlier(sce$total_counts,
                          nmads = 3,type = "lower", log = TRUE
)
feature_flag <- isOutlier(sce$total_features_by_counts,
                          nmads = 3, type = "lower", log = TRUE
)
table(feature_flag) #drops 30 cells due to low variance in RNA type ( < 3x median absolute SD)



low_libsize <- sapply(levels(sce$Sample_ID), function(x)
  sum(libsize_flag[sce$Sample_ID == x]))
low_feature <- sapply(levels(sce$Sample_ID), function(x)
  sum(feature_flag[sce$Sample_ID == x]))

table(low_feature) #30 cells
table(low_libsize) #72 cells

low_libsizeCells <- tibble('low_libsizeCells' = names(low_libsize[low_libsize==1]) ) 
low_featureCells <- tibble('low_featureCells' = names(low_feature[low_feature==1]) ) 

##### record dropped cells #####

flag_records <- rbind(gather(high_mitoCells), 
                      gather(low_libsizeCells), 
                      gather(low_featureCells)) 


flag_records %>% dplyr::count(value,sort=TRUE)

flag_records %>% dplyr::count(key) %>% dplyr::mutate(sum=sum(n))

saveRDS(flag_records, "expn_processing/flag_records.Rds")


#remove offending cells
dim(sce) #3085
sce <- sce[, !(libsize_flag | feature_flag)]
dim(sce) #2989


#Filter genes by expression
keep_genes <- apply(counts(sce), 1, function(x)
  sum(x > 0) > (dim(sce)[2] * 0.01))

table(keep_genes)

sce <- sce[keep_genes, ]

paste0("Low-abundance genes removed: ", sum(!keep_genes))
paste0("Genes retained: ", sum(keep_genes))

### NORMALIZE ###

sce <- normalize(sce)


### CPM ###
# https://bioconductor.org/packages/devel/bioc/vignettes/scater/inst/doc/vignette-intro.html

cpm(sce) <- calculateCPM(sce); cpm(sce)[1:4,1:4]


colSums(cpm(sce)[ ,-1 ])[1:10]

cpmSCE <- tidyGE::tableTOdf(cpm(sce)) %>% dplyr::rename(ENSG=key)
cpmSCE[ 1:5,1:5]


saveRDS(sce, file='expn_processing/Lake_sce.Rds')

saveRDS(cpmSCE, file="expn_processing/BlueLake_cpmSCE.Rds")



#### TPM #####
#requires effective length - using gene length as many introns are present in nuclear transcripts

#sce <- readRDS("expn_processing/Lake_sce.Rds") #normalized above

GRCh38_genes <- readRDS("expn_processing/GRCh38_genes.Rds")

gene_len <- tibble('id' = rowData(sce)$geneID) %>% left_join(GRCh38_genes,by=c('id'='gene_id')) %>% 
  select(id,width)

length(gene_len$width)
table(is.na(gene_len$width)) 
table(sizeFactors(sce)==0); table(sizeFactors(sce)>0)

sce_noERCC <- sce[1:38324]
dim(sce_noERCC); length(gene_len$width)

tpm(sce_noERCC) <- calculateTPM(sce_noERCC, effective_length = gene_len$width[1:38324]) 

colSums(tpm(sce_noERCC)[ ,-1 ]) %>% head()

tpmSCE <- tidyGE::tableTOdf(tpm(sce_noERCC)) %>% dplyr::rename(ENSG=key)

#compare:
tpmSCE[ 1:5,1:5] ; cpmSCE[1:5,1:5]

saveRDS(tpmSCE, file="expn_processing/BlueLake_tpmSCE.Rds")

#### --> sce_p2.R




