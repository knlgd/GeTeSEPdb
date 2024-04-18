# 2023.04.07 主流程
# 2023.05.21 更新，加入富集分析模块
args <- commandArgs(T)
library(tidyverse)
library(Biobase)
library(dplyr)
library(dbplyr)
library(data.table)
library(stringr)
library(readxl)
library(readr)
library(tidyr)
library(ggplot2)
library(enrichplot)
library(clusterProfiler)

#SampleInfo <- read_excel("~/mqf/Temporal/workflow/src/SampleInfo.xlsx")
#SampleInfo <- as.data.frame(read.table("~/mqf/Temporal/workflow/src/SampleInfo1.txt",header=T,sep="\t",fill=TRUE,encoding='UTF-8'))
##########2023.6.19修改sampleinfo数据读取/2023.6.21加入手动处理RNAseq数据

SampleInfo <- read.delim("~/mqf/Temporal/workflow/src/SampleInfo1.txt")

######2023.6.20删除手动下载的RNAseq数据，重新读取新的数据
#SampleInfo <- read.delim("~/mqf/Temporal/workflow/src/SampleInfo2.txt")


Studyid <- args[1]#'TCD1681' #'TCD1718' 'TCD1719' 'TCD1796'
sampleinfo <- SampleInfo[SampleInfo$Studyid == Studyid,]
rslt <- '~/mqf/Temporal/workflow/result_RNAseq'
if(!dir.exists(paste0(rslt,'/',Studyid))){dir.create(paste0(rslt,'/',Studyid))}
rslt <- paste0(rslt,'/',Studyid)

##1.先进行数据预处理，计算TPM，如有生物学重复，制作重复列表文件。
source('~/mqf/Temporal/workflow/lib/PreProcess.R')
PreProcess(sampleinfo = sampleinfo,rslt = rslt)

##2.运行Clust
source('~/mqf/Temporal/workflow/lib/Clust.R')
Clust(sampleinfo = sampleinfo,rslt = rslt)

##3.回归分析
source('~/mqf/Temporal/workflow/lib/regression.R')
Regression(sampleinfo = sampleinfo,rslt = rslt)

##4.富集分析
source('~/mqf/Temporal/workflow/lib/enrich_fun.R')
enrich_fun(sampleinfo = sampleinfo,rslt = rslt)



