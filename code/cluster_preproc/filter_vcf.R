
#filter vcf

#vcf_filter_loop

#This script is run for each sample var_call.vcf file, to extract variants of interest.


#options(repos = c(CRAN = "http://cran.rstudio.com"))

#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#install.packages("SeqArray")


library(tidyverse)
library(gdsfmt)
library(SeqArray)


vcfDir="~/scRNA_editing/study_SRA/gatk/"
dataDir="~/scRNA_editing/study_SRA/vcf_processed/"
writeDir=dataDir 

dbSNP_loci <- readRDS("~/SNP_RNAed_reference/dbSNP_lociOnly.Rds")
DRcomb <- readRDS("~/SNP_RNAed_reference/DARNED_RADARall.Rds")

sampleID <- "MYSAMPLE"


seqVCF2GDS(paste0(dataDir,sampleID,"_var_call.vcf"),
            paste0(writeDir,sampleID,"_var_call.gds"),
            storage.option = "ZIP_RA")


## read vcf header

VCF_head <- seqVCF_Header(paste0(vcfDir, sampleID, "_var_call.vcf"))


# load gds file

gds <- paste0(dataDir,sampleID,"_var_call.gds")

genofile <- seqOpen(gds, allow.duplicate = TRUE)

vcf_filt_join <- c(1:22,"X","Y") %>% 

    map(function(chrmsm, QD=10, DP=5){
    
  mychr <- chrmsm  
  myQD <- 10;    myDP <- 5
    
  seqSetFilterChrom(genofile, mychr) 
    
  mytbl <- str_c(pull(VCF_head$info[1]))[c(3:11,14:18)] %>% 
      map(function(x){
        seqGetData(genofile,paste0("annotation/info/", x)) %>% unlist()
      }) %>% bind_cols() 
    
    names(mytbl) <- unlist(VCF_head$info[1])[c(3:11,14:18)]
    
    
    mytbl_2 <- mytbl %>% mutate("chromosome"= seqGetData(genofile,"chromosome"), 
                                "position"  =  seqGetData(genofile,"position"),
                                "allele"    =  seqGetData(genofile,"allele")) %>% 
      select(chromosome,position,allele,everything()) 
    
    tbl_tagged <- mytbl_2 %>% 
      #dplyr::filter(QD > myQD) %>% dplyr::filter(DP >= myDP) %>% 
      mutate(dbSNP_status=ifelse(position %in% dbSNP_loci[dbSNP_loci$chr==mychr , ]$pos , 
                                 'commonSNP','other')) %>% 
      mutate(chrJoin=paste0('chr',chromosome)) %>% 
      left_join(DRcomb %>% select(chr_hg38,pos_hg38,rdrType,annot1,annot2,inchr,inrna,source),
                by=c('chrJoin'="chr_hg38", "position"='pos_hg38')) %>% 
      mutate(rdrType = ifelse(is.na(rdrType),'uncatalog',rdrType)) 
    
    
    tbl_tagged %>% 
      select(position) %>% pull() %>% map(function(mypos){
        myindex <- which(seqGetData(genofile, "position") == mypos)
        
        myAD_mat <- seqGetData(genofile,"annotation/format/AD")
        NalleleVars <- myAD_mat$length[myindex]
        refCol <- sum(myAD_mat$length[1:myindex]) + 1 - NalleleVars
        altCol <- sum(myAD_mat$length[1:myindex])
        
        if(refCol <= length(myAD_mat$data)){
          refDP <- myAD_mat$data[ , refCol]} 
          else{refDP <- NA}
          altDP <- myAD_mat$data[ , altCol] 

        alleles <- seqGetData(genofile,"allele")[myindex]
        
        tibble("posn"=mypos,
               "NalleleVars" = NalleleVars,
               'refDP'=refDP,
               'altDP'=altDP,
               'alleles'=alleles) %>% 
          mutate(totalDP = sum(altDP,refDP),
                 altProp=altDP/totalDP) %>% 
          separate(alleles,into=c('ref_al','alt_al'),sep=",") %>% 
          filter(ref_al=="A") %>% 
          select(1,2,ref_al,alt_al,refDP,altDP,totalDP,altProp) #,error = function(e) print(NA))
    
      }) %>% bind_rows() -> depth_ratio_table
    
    
	if(nrow(depth_ratio_table>0)){
   	 outTable <- left_join(depth_ratio_table,tbl_tagged,by=c('posn'='position')) %>% 
     		      mutate(position=posn) %>% select(posn:chromosome,position,everything()) %>% filter(ref_al %in% c('A','T'))
    } else{outTable <- as_tibble(matrix(nrow=0,ncol=24), .name_repair = 'minimal')
          names(outTable) <- c("chromosome","position","posn","rdrType","dbSNP_status","NalleleVars","ref_al","alt_al","refDP","altDP","totalDP","altProp","AN","BaseQRankSum","ClippingRankSum","DP","DS","ExcessHet","FS","MQ","MQRankSum","QD","ReadPosRankSum","SOR")}
        
    return(outTable) 
    
  }) %>% bind_rows() %>% 
  select(chromosome,position,posn,rdrType,dbSNP_status,2:8,AN:SOR) %>% 
  select(-c(HaplotypeScore, InbreedingCoeff)) 

seqClose(genofile)
rm(genofile)

write_tsv(vcf_filt_join,
         path=paste0(writeDir,sampleID,"_processed.txt"))
