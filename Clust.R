###Clust流程
# library(reticulate)
# use_condaenv("/public/home/panjianbo/.conda/envs/scanpy38_mqf") #激活虚拟环境

Clust <- function(sampleinfo,rslt){
  #clust结果
  if(!dir.exists(paste0(rslt,'/clust_result')))(dir.create(paste0(rslt,'/clust_result')))
  #先判断平台
  strategy <- unique(sampleinfo$Strategy)
  if(strategy == "RNA-Seq"){
    if(file.exists(paste0(rslt,'/ReplicatesFile.txt'))){
      clust_cmd <- paste0("source /public/software/apps/anaconda3/2021.05/bin/activate scanpy38_mqf && clust ",rslt,'/tpm_profile.txt -o ',rslt,'/clust_result/ -r ',rslt,'/ReplicatesFile.txt -n 101 3 4 -t 1 -q3s 1 -K 10 12 14 16 18 20')
    }else{
      clust_cmd <- paste0("source /public/software/apps/anaconda3/2021.05/bin/activate scanpy38_mqf && clust ",rslt,'/tpm_profile.txt -o ',rslt,'/clust_result/ -n 3 4 -t 1 -q3s 1 -K 10 12 14 16 18 20')
    }
  }else{
    if(file.exists(paste0(rslt,'/ReplicatesFile.txt'))){
      clust_cmd <- paste0("source /public/software/apps/anaconda3/2021.05/bin/activate scanpy38_mqf && clust ",rslt,'/tpm_profile.txt -o ',rslt,'/clust_result/ -r ',rslt,'/ReplicatesFile.txt -n 101 4 -t 1 -q3s 1 -K 10 12 14 16 18 20')
    }else{
      clust_cmd <- paste0("source /public/software/apps/anaconda3/2021.05/bin/activate scanpy38_mqf && clust ",rslt,'/tpm_profile.txt -o ',rslt,'/clust_result/ -n 101 4 -t 1 -q3s 1 -K 10 12 14 16 18 20')
    }
  }
  system(clust_cmd)
}