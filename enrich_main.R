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
setwd("F:\\microarray_analysis")
resultdir = "F:\\microarray_analysis\\result"
source('F:\\microarray_analysis\\Microarray\\lib\\enrich_fun.R')
# reference group
group_list =  read.table("F:\\microarray_analysis\\Microarray\\src\\array.txt",header=T,sep = '\t',fill=T)

#enrich
all_study = list.files('./result/')
for(i in 1:length(all_study)){
  study = group_list %>% filter(Studyid == all_study[i])
  #sample = unlist(str_split(study$Sample,";"))
  enrich_fun(all_study[i])
  print(i)
}
