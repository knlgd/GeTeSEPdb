# 2023.04.07 主流程
# 2023.05.21 更新，加入富集分析模块
args <- commandArgs(T)
library(tidyverse)
library(Biobase)
library(enrichplot)
library(dplyr)
library(dbplyr)
library(data.table)
library(stringr)
library(readxl)
library(readr)
library(tidyr)
library(ggplot2)
library(clusterProfiler)


SampleInfo <- read.table("F:\\project\\数据\\SampleInfo_mysql20230624.txt",head=T,fill=TRUE)

Studyid <- 'TCD0342'
sampleinfo <- SampleInfo[SampleInfo$Studyid == Studyid,]
rslt <- 'F:/project/timecourse/timecourse_app/static/study_detail'
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



