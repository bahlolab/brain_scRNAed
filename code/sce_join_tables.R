#Lake et al 2016 Science
#integrate counts data for human genes and ERCC spike-in 


library(tidyverse)


# metaData from Dr Blue Lake ; communicated via email 15 Jun 2019
BL_joinKey <- readxl::read_xlsx("fromBL_Suppl_Tables_v5.xlsx",
                                sheet=2, skip=5) %>% 
  rename_all(funs(str_replace_all(., ' ','_'))) %>% 
  rename_all(funs(str_replace_all(., '-','_'))) %>% 
  mutate(neuType=ifelse(str_detect(Sub_Group_ID,'^Ex'),'Ex','In'),
         neuType=ifelse(str_detect(Group_ID, "NoN"),'NoN', neuType)) %>% 
  mutate(area=str_extract(Sub_Group_ID,"BA.*")) %>% 
  filter(!is.na(SRA))

names(BL_joinKey) <- str_remove(names(BL_joinKey),'†')

BL_joinKey %>% select(SRA) %>% distinct() %>% dim() #3127 unique SRA IDs only
BL_joinKey %>% filter(str_detect(SRA,"SRR1764085|SRR3230612"))
#SRR1764085 excluded: not in BL_joinKey
#SRR3230612 excluded: did not finish for unknown reason

############ ############ ############
# Integrate counts data
############ ############ ############


bigTable1 <- read_tsv("~/hg_ercc_matrix.tsv", col_names = TRUE)
hg2 <- read_tsv("~/remaining1432_fc_matrix.tsv", col_names = TRUE)
ercc2 <- read_tsv("~/remaining1432_ercc_fc_matrix.tsv", col_names = TRUE)
bigTable2 <- rbind(hg2, ercc2) 

#add 5 remaining cells
r6 <- read_tsv("~/remaining6_fCount.txt",col_names = FALSE) %>% 
  filter(X1!="SRR1764085")
r6_ercc <- read_tsv("~/remaining6_ercc_fCount.txt",col_names = FALSE) %>% 
  filter(X1!="SRR1764085")


###### ###### ###### ###### ###### ###### 

r6_bound <- rbind(r6,r6_ercc) %>% spread(X1,X3) %>% dplyr::rename(Geneid = X2)
dim(r6_bound)

bigTable1 %>% select(-c(SRR1764085,SRR3230612)) %>% filter(!complete.cases(.))
bigTable2 %>% filter(!complete.cases(.))

p1 <- bigTable1 %>% select(Geneid, which(colnames(bigTable1) %in% BL_joinKey$SRA))
p2 <- bigTable2 %>% select(Geneid, which(colnames(bigTable2) %in% BL_joinKey$SRA))

p3 <- left_join(p1,p2 , by="Geneid")  
BL_expn_table <- left_join(p3 %>% select(-SRR3230612), r6_bound, by='Geneid') 

which(! colnames(BL_expn_table) %in% BL_joinKey$SRA) 
table(colnames(BL_expn_table) %in% BL_joinKey$SRA )

BL_expn_table %>% filter(!complete.cases(.))

saveRDS(BL_expn_table,'expn_processing/BL_expn_table.Rds')

### --> continued in sce_p1.R
