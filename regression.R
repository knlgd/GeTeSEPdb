###提取每个cluster的Eigengenes做回归分析
library(minpack.lm)
library(ggplot2)
library(forecast)
library(segmented)
library(splines)
library(Hmisc)
library(rms)
library(mgcv)
library(caret)
library(tidyverse)
library(MASS)
library(drc)
library(nlme)
library(aomisc)
library(dplyr)
library(tibble)
library(data.table)

Regression <- function(sampleinfo,rslt){
  studyid <- unique(sampleinfo$Studyid)
  species <- unique(sampleinfo$Organism) %>% tolower() %>% gsub(" ","_",.) %>% str_to_title(.)
  #读取Eigengens进行拟合
  Eigengenes <- read.table(paste0(rslt,'/clust_result/Eigengenes.tsv'))
  colnames(Eigengenes) <- unique(as.numeric(sampleinfo$phototime)+1)######2023.6.21修改

  #读取每个cluster的基因信息
  Clusters <- read.delim(paste0(rslt,'/clust_result/Clusters_Objects.tsv')) %>% dplyr::slice(-1)
  #读取处理后的表达谱
  processed_tpm <- read.delim(paste0(rslt,'/clust_result/Processed_Data/tpm_profile.txt_processed.tsv'),check.names = F)
  #各个cluster的评价指标
  EvaIndicators <- data.frame()
  #总的各个cluster的基因表达
  ResGeneDf <- data.frame()
  #各个cluster模型，参数以及表达式
  ModelDf <- data.frame()

  for(clu in 1:nrow(Eigengenes)){
    tryCatch({
      gene_expr <- as.numeric(t(Eigengenes[clu,]))
      #colnames(gene_expr) <- 'gene_expr'
      time <- as.numeric(names(Eigengenes))
      data <- data.frame(time = time, gene_expr = gene_expr)
      
      Linear <- fit_linear(data)
      exponential_growth <- fit_exp_growth(data)
      exponential_decay <- fit_exp_decay(data)
      Logarithmic <- fit_log(data)
      Power <- fit_power(data)
      Logis <- fit_logis(data)
      Sin <- fit_sin(data)
      Quadratic <- fit_Gaussian(data)
    
      Modles <- list(Linear,exponential_growth,exponential_decay,Logarithmic,Power,Logis,Sin,Quadratic)
      names(Modles) <- c('Linear','exponential_growth','exponential_decay','Logarithmic','Power','Logis','Sin','Quadratic')
      #比较各个参数，通过rank并给定权重，选择最优模型
      Goodness <- data.frame()
      for(modle in names(Modles)){
        df <- Modles[[modle]]$result
        rownames(df) <- modle
        colnames(df) <- c('rmse','rsq','adjrsq','aicc')
        Goodness <- rbind(Goodness,df)
      }
      Goodness$cluster <- paste0('cluster_',clu)
      if(TRUE %in% (Goodness$adjrsq > 0.7)){  
        best_model <- rownames(Goodness[Goodness$adjrsq > 0.7,])
        cluster_gene <- Clusters[,clu][nchar(Clusters[,clu])>0]  ##去除了有些空的字符串
        cluster_gene_tpm <- processed_tpm[processed_tpm$Genes %in% cluster_gene,] %>% remove_rownames %>%
          column_to_rownames(.,'Genes')
        for(m in best_model){
          bm <- Modles[[m]]$modle
          if(m=='Quadratic'){
            p<- 3
          }else if(m=='Sin' | m=='Logis'){
            p<- 4
          }else{
            p<- 2
          }
        
          
          #再对各个cluster中的基因使用最佳模型进行拟合，控制R2阈值，保留R2>0.9的基因做后续分析
          res_gene_rsq <- apply(cluster_gene_tpm, 1, function(y) {
            pred <- predict(bm)
            tss <- sum((y - mean(y))^2)
            rss <- sum((y - pred)^2)
            ssr <- tss - rss
            rsq <- ssr / tss
            n <- length(y)
            ajr2 <- 1 - ((1 - rsq) * (n - 1)) / (n - p - 1)
            return(ajr2)
          })
          res_gene <- names(res_gene_rsq[res_gene_rsq > 0.7])
          res_gene_df <- processed_tpm[processed_tpm$Genes %in% res_gene,]
          res_gene_df$cluster <- paste0('cluster_',clu)
          res_gene_df$model <- m
          ResGeneDf <- rbind(ResGeneDf,res_gene_df)
          #输出ModelDf
          res_modeldf <- switch(m,
                                "Linear" = data.frame(cluster = paste0('cluster_',clu),
                                                      modle = 'Linear',
                                                      Parameters = paste0('a = ',coefficients(bm)[1],'; b = ',coefficients(bm)[2]), 
                                                      exp = 'y = a + b*x'),
                                "exponential_growth" = data.frame(cluster = paste0('cluster_',clu),
                                                                  modle = 'Exponential',
                                                                  Parameters = paste0('a = ',coefficients(bm)[1],'; k = ',coefficients(bm)[2]), 
                                                                  exp = 'y = a*exp(k*x)'),
                                "exponential_decay" = data.frame(cluster = paste0('cluster_',clu),
                                                                 modle = 'Exponential',
                                                                 Parameters = paste0('a = ',coefficients(bm)[1],'; k = ',coefficients(bm)[2]), 
                                                                 exp = 'y = a*exp(k*x)'),
                                "Logarithmic" = data.frame(cluster = paste0('cluster_',clu),
                                                           modle = 'Logarithmic',
                                                           Parameters = paste0('a = ',coefficients(bm)[1],'; b = ',coefficients(bm)[2]), 
                                                           exp = 'y = a + b * log(x)'),
                                "Power" = data.frame(cluster = paste0('cluster_',clu),
                                                     modle = 'Power',
                                                     Parameters = paste0('a = ',coefficients(bm)[1],'; b = ',coefficients(bm)[2]), 
                                                     exp = 'y = a * x ^ b'),
                                "Logis" = data.frame(cluster = paste0('cluster_',clu),
                                                     modle = 'Logis',
                                                     Parameters = paste0('b = ',coefficients(bm)[1],'; c = ',coefficients(bm)[2],'; d = ',coefficients(bm)[3],'; e = ',coefficients(bm)[4]), 
                                                     exp = 'y = c + (d - c)/(1 + exp(b * (x - e)))'),
                                "Sin" = data.frame(cluster = paste0('cluster_',clu),
                                                   modle = 'Sin',
                                                   Parameters = paste0('A = ',coefficients(bm)[1],'; freq = ',coefficients(bm)[2], '; phi = ',coefficients(bm)[3], '; C = ',coefficients(bm)[4]), 
                                                   exp = 'y = A * sin(2*pi*freq * x + phi) + C'),
                                "Quadratic" = data.frame(cluster = paste0('cluster_',clu),
                                                        modle = 'Quadratic',
                                                        Parameters = paste0('a = ',coefficients(bm)[1],'; b = ',coefficients(bm)[2],'; c = ',coefficients(bm)[3]), 
                                                        exp = 'y = a + b*x + c*x^2')
          )
          ModelDf <- rbind(ModelDf,res_modeldf)
        }
      }else{
        sink(paste0('~/mqf/Temporal/workflow/result/regression.log'), append = T)
        print(paste0(studyid,' cluster_',clu,' fit failed'))
        sink()
      }
      Goodness <- rownames_to_column(Goodness,'model')
      EvaIndicators <- rbind(EvaIndicators,Goodness)
    },
    error = function(e){
      sink(paste0('~/mqf/Temporal/workflow/result/regression.log'), append = T)
      print(paste0(studyid,' cluster_',clu,' fit failed'))
      sink()
      
    })
  }
  
  #对housekeeping基因拟合
  tpm_profile <- read.delim(paste0(rslt,'/tpm_profile.txt'),check.names = F)
  all_cluster_gene <- Clusters[,1]
  if(ncol(Clusters)>1){
    for(k in 2:ncol(Clusters)){
      all_cluster_gene <- union(all_cluster_gene,Clusters[,k])
    }
  }
  rest_tpm <- tpm_profile[!tpm_profile$GeneID %in% all_cluster_gene,] %>% as.data.frame()
  #select HKG
  #1, 剔除所有样本中表达为0的基因
  #2, 评估每个基因的变异系数（CV），标准差和平均值的比值
  #3, 设定阈值，变异系数越低，说明越不表达
  rest_tpm <- aggregate(.~GeneID, mean, data = rest_tpm)
  rownames(rest_tpm) <- rest_tpm$GeneID
  rest_tpm <- rest_tpm[,-1]
  #删除包含0值的行
  rest_tpm[rest_tpm == 0] <- NA
  rest_tpm <- na.omit(rest_tpm)
  #计算每个基因的CV
  cal_cv = function(x){
    y=na.omit(x)
    return(sd(y)/mean(y))
  }
  rest_tpm$CV <- apply(rest_tpm,1,cal_cv)
  rest_tpm <- rest_tpm[order(rest_tpm$CV),]
  if(nrow(rest_tpm[rest_tpm$CV < 0.1,])>0 && nrow(rest_tpm[rest_tpm$CV < 0.1,])-nrow(Clusters) < 100){  ####2023.06.20 这里再加一个判断 有的研究所有基因的cv都大于0.1
    constantDF <- rest_tpm[rest_tpm$CV < 0.1,]
  }else if((nrow(rest_tpm) * 0.05) - nrow(Clusters) < 100){
    constantDF <- rest_tpm[1:(nrow(rest_tpm) * 0.05),]
  }else{
    constantDF <- rest_tpm[1:nrow(Clusters),]
  }
  constantDF$mean <- rowMeans(constantDF[,1:(ncol(constantDF)-1)])
  constantDF <- rownames_to_column(constantDF,var = 'Genes')
  if(file.exists(paste0(rslt,'/ReplicatesFile.txt'))){
    group_info <- read.delim(paste0(rslt,'/ReplicatesFile.txt'),check.names = F,header=FALSE)
    constantDF_mean <-constantDF[,c("Genes","CV","mean")]
    for(g in 1:nrow(group_info)){
      name_group <- group_info[g,]$V2
      sample_group <- unlist(strsplit(group_info[g,]$V3,','))
      if(length(sample_group) == 1){
        constantDF_group <- constantDF[,sample_group] %>% as.data.frame()
        colnames(constantDF_group) <- sample_group
      }else{
        constantDF_group <- constantDF[,sample_group]
      }
      
      constantDF_group$mean <- apply(constantDF_group,1,mean)
      colnames(constantDF_group) <- c(sample_group,name_group)
      constantDF_group <- as.data.frame(constantDF_group[,name_group])
      colnames(constantDF_group) <- name_group
      constantDF_mean <- cbind(constantDF_mean,constantDF_group)
      # constantDF_mean <- constantDF_mean[,c(1,4:ncol(constantDF_mean),2:3)]
    }
  }else{
    constantDF_mean <- constantDF
  }
  
  constantDF_mean$Studyid <- studyid
  constantDF_mean$cluster <- 'cluster_0'
  constantDF_mean$model <- 'constant'
  constantDF_mean$Genes <- as.numeric(constantDF_mean$Genes)
  
  if(nrow(ResGeneDf) > 0){   ###2023.06.21 添加判断，解决没有拟合出来结果，ResGeneDf为空的报错
    ResGeneDf$Studyid <- studyid
    EvaIndicators$Studyid <- studyid
    ModelDf$Studyid <- studyid
    
    #2023.06.17 将替换的Entrezid替换回Symbol
    sym_enz_ens <- fread(paste0('~/mqf/Temporal/workflow/src/',species,'_sym_enz_ens.txt')) %>% dplyr::select(.,c('GeneID','Symbol'))
    ResGeneDf <- left_join(sym_enz_ens,ResGeneDf,by = c('GeneID'='Genes')) %>% drop_na() %>% dplyr::select(.,-'GeneID')
    colnames(ResGeneDf)[1] <- 'Genes'
    constantDF_mean <- left_join(sym_enz_ens,constantDF_mean,by = c('GeneID'='Genes')) %>% drop_na() %>% dplyr::select(.,-'GeneID')
    colnames(constantDF_mean)[1] <- 'Genes'
    
    
    #和表，基因表格
    browse_gene_1 <- ResGeneDf[,c('Studyid','cluster','model','Genes')]
    browse_gene_1$CV <- NA
    browse_gene_1$mean <- NA
    browse_gene_2 <- constantDF_mean[,c('Studyid','cluster','model','Genes','CV','mean')]
    browse_gene <- rbind(browse_gene_2,browse_gene_1)
    
    write.csv(ResGeneDf, file = paste0(rslt,'/ResGeneDF0409.csv'),row.names = F)
    write.csv(ModelDf,file = paste0(rslt,'/ModelDF0409.csv'),row.names = F)
    write.csv(EvaIndicators,paste0(rslt,'/EvaIndicators0409.csv'),row.names = F)
    write.csv(constantDF_mean,paste0(rslt,'/ConstantDF0409.csv'),row.names = F)
    write.csv(browse_gene,paste0(rslt,'/BrowseGene0409.csv'),row.names = F)
  }else{
    sym_enz_ens <- fread(paste0('~/mqf/Temporal/workflow/src/',species,'_sym_enz_ens.txt')) %>% dplyr::select(.,c('GeneID','Symbol'))
    constantDF_mean <- left_join(sym_enz_ens,constantDF_mean,by = c('GeneID'='Genes')) %>% drop_na() %>% dplyr::select(.,-'GeneID')
    colnames(constantDF_mean)[1] <- 'Genes'
    browse_gene <- constantDF_mean[,c('Studyid','cluster','model','Genes','CV','mean')]
    
    write.csv(constantDF_mean,paste0(rslt,'/ConstantDF0409.csv'),row.names = F)
    write.csv(browse_gene,paste0(rslt,'/BrowseGene0409.csv'),row.names = F)
  }
  
}

#封装拟合模型
#线性模型
fit_linear <- function(data){
  x <- data$time
  y <- data$gene_expr
  lm_fit <- lm(y ~ x,data = data)
  lm_pred <- predict(lm_fit)
  lm_rmse <- sqrt(mean((lm_pred - y)^2))
  lm_r2 <- summary(lm_fit)$r.squared
  lm_ajr2 <- summary(lm_fit)$adj.r.squared
  lm_aicc <- AIC(lm_fit)
  lm_result=data.frame(lm_rmse=lm_rmse,lm_r2=lm_r2,lm_ajr2=lm_ajr2,lm_aicc=lm_aicc)
  return(list(modle = lm_fit, result = lm_result))
}
#指数增长模型
fit_exp_growth <- function(data){
  x <- data$time
  y <- data$gene_expr
  exp_fit <- drm(y ~ x, fct = DRC.expoGrowth(), data = data)######2023.6.21增加了参数logDose=2
  # exp_fit <- nls(y ~ NLS.expoGrowth(x,a,k),data = data)
  exp_pred <- predict(exp_fit)
  exp_tss <- sum((y - mean(y))^2)# Calculate TSS
  exp_rss <- sum((y - exp_pred)^2)# Calculate RSS
  exp_ssr <- exp_tss - exp_rss # Calculate SSR
  exp_rsq <- exp_ssr / exp_tss# Calculate R^2
  exp_rmse <- sqrt(mean((exp_pred - y)^2))
  exp_aicc <- AIC(exp_fit)
  n <- length(y)
  p <- 2 # Assuming 2 parameters in the model
  exp_ajr2 <- 1 - ((1 - exp_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  exp_result <- data.frame(exp_rmse=exp_rmse, exp_rsq=exp_rsq, exp_ajr2=exp_ajr2, exp_aicc=exp_aicc)
  #exp_result=data.frame(exp_rmse=exp_rmse,exp_rsq=exp_rsq,exp_aicc=exp_aicc)
  return(list(modle = exp_fit,result = exp_result))
}
#指数递减模型
fit_exp_decay <- function(data){
  x <- data$time
  y <- data$gene_expr
  exp_fit <- drm(y ~ x, fct = DRC.expoDecay(), data = data)
  # exp_fit <- nls(y ~ NLS.expoGrowth(x,a,k),data = data)
  exp_pred <- predict(exp_fit)
  exp_tss <- sum((y - mean(y))^2)# Calculate TSS
  exp_rss <- sum((y - exp_pred)^2)# Calculate RSS
  exp_ssr <- exp_tss - exp_rss # Calculate SSR
  exp_rsq <- exp_ssr / exp_tss# Calculate R^2
  exp_rmse <- sqrt(mean((exp_pred - y)^2))
  exp_aicc <- AIC(exp_fit)
  n <- length(y)
  p <- 2 # Assuming 2 parameters in the model
  exp_ajr2 <- 1 - ((1 - exp_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  exp_result <- data.frame(exp_rmse=exp_rmse, exp_rsq=exp_rsq, exp_ajr2=exp_ajr2, exp_aicc=exp_aicc)
  #exp_result=data.frame(exp_rmse=exp_rmse,exp_rsq=exp_rsq,exp_aicc=exp_aicc)
  return(list(modle = exp_fit,result = exp_result))
}
#对数函数
fit_log <- function(data){
  x <- data$time
  y <- data$gene_expr
  #loglm_fit <- lm(y ~ log(x), data = data)
  loglm_fit <- drm(y ~ x, fct = DRC.logCurve(),data = data) ##################修改
  loglm_pred <- predict(loglm_fit)
  loglm_tss <- sum((y - mean(y))^2)# Calculate TSS
  loglm_rss <- sum((y - loglm_pred)^2)# Calculate RSS
  loglm_ssr <- loglm_tss - loglm_rss # Calculate SSR
  loglm_rsq <- loglm_ssr / loglm_tss# Calculate R^2
  loglm_rmse <- sqrt(mean((loglm_pred - y)^2))
  loglm_aicc <- AIC(loglm_fit)
  n <- length(y)
  p <- 2 # Assuming 2 parameters in the model
  loglm_ajr2 <- 1 - ((1 - loglm_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  loglm_result <- data.frame(loglm_rmse=loglm_rmse,loglm_rsq=loglm_rsq, loglm_ajr2=loglm_ajr2,loglm_aicc=loglm_aicc)
  #loglm_result=data.frame(loglm_rmse=loglm_rmse,loglm_rsq=loglm_rsq,loglm_aicc=loglm_aicc)
  return(list(modle = loglm_fit, result = loglm_result))
}
#幂函数
fit_power <- function(data){
  x <- data$time
  y <- data$gene_expr
  power_fit <- drm(y ~ x, fct = DRC.powerCurve(),data = data)
  power_pred <- predict(power_fit)
  power_tss <- sum((y - mean(y))^2)# Calculate TSS
  power_rss <- sum((y - power_pred)^2)# Calculate RSS
  power_ssr <- power_tss - power_rss # Calculate SSR
  power_rsq <- power_ssr / power_tss# Calculate R^2
  power_rmse <- sqrt(mean((power_pred - y)^2))
  power_aicc <- AIC(power_fit)
  n <- length(y)
  p <- 2 # Assuming 2 parameters in the model
  power_ajr2 <- 1 - ((1 - power_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  power_result <- data.frame(power_rmse=power_rmse,power_rsq=power_rsq, power_ajr2=power_ajr2,power_aicc=power_aicc)
  #power_result=data.frame(power_rmse=power_rmse,power_rsq=power_rsq,power_aicc=power_aicc)
  return(list(modle = power_fit,result = power_result))
}
#逻辑回归
fit_logis <- function(data){
  x <- data$time
  y <- data$gene_expr
  Logis_fit <- drm(y ~  x,fct = L.4(), data = data)
  Logis_pred <- predict(Logis_fit)
  Logis_tss <- sum((y - mean(y))^2)# Calculate TSS
  Logis_rss <- sum((y - Logis_pred)^2)# Calculate RSS
  Logis_ssr <- Logis_tss - Logis_rss # Calculate SSR
  Logis_rsq <- Logis_ssr / Logis_tss# Calculate R^2
  Logis_rmse <- sqrt(mean((Logis_pred - y)^2))
  Logis_aicc <- AIC(Logis_fit)
  n <- length(y)
  p <- 4 # Assuming 2 parameters in the model
  Logis_ajr2 <- 1 - ((1 - Logis_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  Logis_result <- data.frame(Logis_rmse=Logis_rmse,Logis_rsq=Logis_rsq, Logis_ajr2=Logis_ajr2,Logis_aicc=Logis_aicc)
  #Logis_result=data.frame(Logis_rmse=Logis_rmse,Logis_rsq=Logis_rsq,Logis_aicc=Logis_aicc)
  return(list(modle = Logis_fit, result = Logis_result))
}
#三角函数
SStrig <- function(x, A, freq, phi, C) {
  A * sin(2*pi*freq*x + phi) + C
}
fit_sin <- function(data){
  x <- data$time
  # x <- data$time * 2 * pi / max(data$time)
  y <- data$gene_expr
  # 循环结构拟合每个数据
  nlc <- nls.control(maxiter = 1000)
  ssp <- spectrum(y)
  per <- 1/ssp$freq[ssp$spec==max(ssp$spec)]
  start_values <- c(A = (max(y) - min(y))/2,
                    freq = 1/per,
                    phi = 0,
                    C = (max(y) + min(y))/2)
  sin_fit <- nlsLM(gene_expr ~ SStrig(x, A, freq, phi, C), 
                   start = start_values, control = nlc, data = data)
  
  #sin_fit <- nlsLM(y ~ a * sin(b * x + c), start = list(a = 1, b = 2, c = 0))
  # sin_fit <- lm(y ~ sin(x),data = data)
  sin_pred <- predict(sin_fit)
  sin_rmse <- sqrt(mean((sin_pred - y)^2))
  sin_tss <- sum((y - mean(y))^2)# Calculate TSS
  sin_rss <- sum((y - sin_pred)^2)# Calculate RSS
  sin_ssr <- sin_tss - sin_rss # Calculate SSR
  sin_rsq <- sin_ssr / sin_tss# Calculate R^2
  sin_aicc <- AIC(sin_fit)
  n <- length(y)
  p <- 4 # Assuming 2 parameters in the model
  sin_ajr2 <- 1 - ((1 - sin_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  sin_result <- data.frame(sin_rmse=sin_rmse,sin_rsq=sin_rsq, sin_ajr2=sin_ajr2,sin_aicc=sin_aicc)
  #sin_result=data.frame(sin_rmse=sin_rmse,sin_rsq=sin_rsq,sin_aicc=sin_aicc)
  return(list(modle = sin_fit,result = sin_result))
}
#高斯函数
fit_Gaussian <- function(data){
  x <- data$time
  y <- data$gene_expr
  # Gaussian_fit <- drm(y ~ x, fct = DRC.bragg.3())
  # Gaussian_fit <- nls(y ~ a * exp(-((x-b)/c)^2),data = data, start = list(a = max(y),b=mean(x),c = sd(x)))
  #Gaussian_fit <- lm(y~poly(x,2))
  Gaussian_fit <- drm(gene_expr ~ time, fct = DRC.poly2(),data = data) #######################################################修改
  
  Gaussian_pred <- predict(Gaussian_fit)
  Gaussian_tss <- sum((y - mean(y))^2)# Calculate TSS
  Gaussian_rss <- sum((y - Gaussian_pred)^2)# Calculate RSS
  Gaussian_ssr <- Gaussian_tss - Gaussian_rss # Calculate SSR
  Gaussian_rsq <- Gaussian_ssr / Gaussian_tss# Calculate R^2
  Gaussian_rmse <- sqrt(mean((Gaussian_pred - y)^2))
  Gaussian_aicc <- AIC(Gaussian_fit)
  n <- length(y)
  p <- 3 # Assuming 2 parameters in the model
  Gaussian_ajr2 <- 1 - ((1 - Gaussian_rsq) * (n - 1)) / (n - p - 1) # Calculate adjusted R^2
  Gaussian_result <- data.frame(Gaussian_rmse=Gaussian_rmse,Gaussian_rsq=Gaussian_rsq, Gaussian_ajr2=Gaussian_ajr2,Gaussian_aicc=Gaussian_aicc)
  #Gaussian_result=data.frame(Gaussian_rmse=Gaussian_rmse,Gaussian_rsq=Gaussian_rsq,Gaussian_aicc=Gaussian_aicc)
  return(list(modle = Gaussian_fit,result = Gaussian_result))
}
