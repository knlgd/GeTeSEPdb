library(dplyr)
library(stringr)
library(stringi)
library(tibble)
library(readxl)
library(readr)
library(data.table)

###数据预处理，对于转录组数据计算TPM值，芯片数据已经取了对数，不用再做处理。如果有生物学重复，制作重复信息表
PreProcess <- function(sampleinfo,rslt){
  species <- unique(sampleinfo$Organism) %>% tolower() %>% gsub(" ","_",.) %>% str_to_title(.)
  sym_enz_ens <- fread(paste0('~/mqf/Temporal/workflow/src/',species,'_sym_enz_ens.txt')) %>% dplyr::select(.,c('GeneID','Symbol'))
  sym_enz_ens <- sym_enz_ens[!duplicated(sym_enz_ens),]
  sample <- sampleinfo$Sample
  datasource <- unique(sampleinfo$RNAseq_other)  ###2023.06.17 新增一列判断数据是RNA-seq流程跑出 还是自己手动整理
  # strategy <- unique(sampleinfo$Strategy)
  GSE <- unique(sampleinfo$Accession)
  #读取原始数据
  if(datasource == "RNA-Seq"){
    #读取GSM-SRR对应关系表
    combineresult <- read.delim('~/mqf/Temporal/workflow/src/combineresult.txt')
    srr <- combineresult[combineresult$GSM %in% sample,]$SRR
    data_dir <- '~/RNA-seq/temporal/combine'
    all_profile <- fread(paste0(data_dir,'/',species,'_allprofile_count.csv'),sep = ',') %>% as.data.frame()
    raw_data <- all_profile[,which(colnames(all_profile) %in% c('Geneid',srr,"Length"))]
    raw_data <- raw_data[,c('Geneid',srr,"Length")]
    ##替换列名
    names(raw_data) <- combineresult$GSM[match(names(raw_data),combineresult$SRR)]
    colnames(raw_data)[1] <- 'Geneid'
    colnames(raw_data)[ncol(raw_data)] <- 'Length'
    raw_data <- raw_data[,c('Geneid',sample,'Length')]  ###确保列名的顺序与时序一致
    #2023.06.17，先全部转成Entrezid，后续再转回来
    raw_data <- left_join(sym_enz_ens,raw_data,by=c('Symbol' = 'Geneid')) %>% drop_na()
    raw_data <- raw_data[,-'Symbol'] %>% as.data.frame()
    
    mycounts <- raw_data[,c(-ncol(raw_data))] %>% column_to_rownames(.,'GeneID')
    #计算TPM
    len <- raw_data$Length
    kb <- len / 1000
    RPKM <- mycounts / kb
    exp <- t(t(RPKM)/colSums(RPKM) * 1000000) %>% as.data.frame()
    exp$GeneID <- rownames(exp)
    exp <- exp[,c(ncol(exp),1:(ncol(exp)-1))]
    write.table(exp,file = paste0(rslt,'/tpm_profile.txt'),row.names = F,sep = '\t',quote = F)
    write.table(raw_data,file = paste0(rslt,'/raw_profile.txt'),row.names = F,sep = '\t',quote = F)
  }else{
    #读取数据
    all_profile <- read.delim(paste0('~/mqf/Temporal/workflow/Microarray/',GSE,'_exp.txt'))
    raw_data <- all_profile[,which(colnames(all_profile) %in% c('genename',sample))]
    raw_data <- raw_data[,c('genename',sample)]  ###确保列名的顺序与时序一致
    colnames(raw_data)[1] <- 'ID'
    #2023.06.17，先全部转成Entrezid，后续再转回来
    raw_data <- left_join(sym_enz_ens,raw_data,by=c('Symbol' = 'ID')) %>% drop_na()
    raw_data <- raw_data[,-'Symbol'] %>% as.data.frame()
    
    write.table(raw_data,file = paste0(rslt,'/tpm_profile.txt'),row.names = F,sep = '\t',quote = F)##这里实际不是tpm值，只是为了保持跟RNA-seq的名称一致便于流程，所以叫tpm
  }
  
  ###制作样本信息表
  if(length(unique(sampleinfo$Time)) != length(sampleinfo$Time)){ ##Time去重不等于Time的行数，说明有生物学重复
    ReplicatesFile <- data.frame()
    time <- unique(sampleinfo$Time)
    for(i in 1:length(time)){
      time_info <- paste(time[i],unique(sampleinfo$Unit),sep = '_')
      time_sample <- paste(sampleinfo[sampleinfo$Time == time[i],]$Sample,collapse = ',')
      ReplicatesFile <- rbind(ReplicatesFile,data.frame('tpm_profile.txt',time_info,time_sample))
    }
    write.table(ReplicatesFile,file = paste0(rslt,'/ReplicatesFile.txt'),row.names = F,col.names = F,sep = '\t',quote = F)
  }
}