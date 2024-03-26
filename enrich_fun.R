# library(org.Mm.eg.db)
# library(org.Hs.eg.db)
# library(org.At.tair.db)
# library(org.Rn.eg.db)
# library(org.Sc.sgd.db)
# library(org.Dm.eg.db)
# library(org.Dr.eg.db)
# library(org.EcK12.eg.db)
# library(org.Ce.eg.db)
# library(org.Bt.eg.db)
# library(org.Ss.eg.db)
# library(org.Mmu.eg.db)
library(tidyr)
library(ggplot2)

enrich_GO_plot = function(cluster_gene,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,sym_enz_ens,cluster_name){
  x <- clusterProfiler::enricher(cluster_gene,TERM2GENE = GO_TERM2GENE,pvalueCutoff = 1,qvalueCutoff = 1)
  #x = clusterProfiler::setReadable(x,orgDb,keyType="ENTREZID") %>% dplyr::select(-Description)
  x <- x %>% dplyr::select(-Description)
  GOenrich_result = x@result %>% left_join(GO_TERM2NAME,by = c('ID'='PATHID')) %>% 
    left_join(GO_TERM2ONT,by = c('ID'='PATHID')) %>% dplyr::select(ID,NAME,ONT,everything()) %>% dplyr::rename(Description = 'NAME')
  for(k in 1:nrow(GOenrich_result)){
    idlist <- GOenrich_result[k,]$geneID
    idlist <- unlist(str_split(idlist,'/'))
    tmp <- sym_enz_ens[sym_enz_ens$GeneID %in% idlist,]
    genelist = unique(tmp$Symbol)
    GOenrich_result[k,]$geneID <- paste(genelist,collapse = '/')
  }
  GOenrich_result$cluster <- cluster_name
  GOenrich_result$Studyid <- studyid
  return(GOenrich_result)
  # write_tsv(GOenrich_result,paste0(rslt,'/enrich/',cluster_name,"_GO.tsv"))
  # GObubble(GOenrich_result,regulate,DEG_length)
}
enrich_KEGG_plot = function(cluster_gene,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,sym_enz_ens,cluster_name){
  x <- clusterProfiler::enricher(cluster_gene,TERM2GENE = KEGG_TERM2GENE,pvalueCutoff =1,qvalueCutoff = 1)
  #x = clusterProfiler::setReadable(x,orgDb,keyType="ENTREZID") %>% dplyr::select(-Description)
  x <- x %>% dplyr::select(-Description)
  KEGGenrich_result = x@result %>% left_join(KEGG_TERM2NAME,by = c('ID'='KEGGPATHID')) %>% dplyr::select(ID,NAME,everything()) %>% dplyr::rename(Description = 'NAME')
  for(k in 1:nrow(KEGGenrich_result)){
    idlist <- KEGGenrich_result[k,]$geneID
    idlist <- unlist(str_split(idlist,'/'))
    tmp <- sym_enz_ens[sym_enz_ens$GeneID %in% idlist,]
    genelist = unique(tmp$Symbol)
    KEGGenrich_result[k,]$geneID <- paste(genelist,collapse = '/')
  }
  KEGGenrich_result$cluster <- cluster_name
  KEGGenrich_result$Studyid <- studyid
  return(KEGGenrich_result)
  # write_tsv(KEGGenrich_result,paste0(rslt,'/enrich/',cluster_name,"_KEGG.tsv"))
  # KEGGbubble(KEGGenrich_result,regulate)
}

BUBBLE = function(enrich,type){
  if(toupper(type) == "GO" | toupper(type) == "KEGG"){
    enrich = enrich %>% drop_na(ID)
    enrich = enrich %>% separate(GeneRatio,c("enrichcount","genecount"),"/",remove = TRUE)
    cluster_gene_length = enrich$genecount[1] %>% as.numeric()
    enrich$GeneRatio = apply(enrich %>% dplyr::select(enrichcount,genecount),1,function(x)as.numeric(x[1])/as.numeric(x[2]))
    enrich = enrich %>% arrange(desc(GeneRatio))
    x=enrich$GeneRatio
    y<-factor(enrich$Description,levels = rev(enrich$Description))
    ggplot(enrich,aes(x,y))+
      geom_point(aes(size=Count,color=p.adjust))+ 
      guides(shape = guide_legend(order=1,override.aes=list(size=3.5)),color = guide_colourbar(order=2),size = guide_legend(order=3))+
      scale_color_gradient(low = "red", high = "blue")+ 
      labs(color=expression(p.adjust),size="Count",x="GeneRatio",y="",title="")+
      theme_bw()+
      theme(plot.title = element_text(hjust = 0.5),
            #去网格去背景色
            axis.line = element_line(colour = "black"), axis.text = element_text(color = "black",size = 14),
            #刻度字体大小
            legend.text = element_text(size = 14),legend.title=element_text(size=14),
            axis.title.x = element_text(size = 14))+scale_size_continuous(range=c(4,8))
  }else{
    print("Invalid data type input!")
  }
}


######
# organism_orgDb = data.frame(organism = c("Mus_musculus","Homo_sapiens","Arabidopsis_thaliana",
#                                          "Rattus_norvegicus","Saccharomyces_cerevisiae","Drosophila_melanogaster",
#                                          "Danio_rerio","Escherichia_coli","Caenorhabditis_elegans","Bos_taurus",
#                                          "Sus_scrofa","Macaca_mulatta","Oryza_sativa"),
#                             OrgDbs = c('org.Mm.eg.db','org.Hs.eg.db','org.At.tair.db','org.Rn.eg.db',
#                                        'org.Sc.sgd.db','org.Dm.eg.db','org.Dr.eg.db','org.EcK12.eg.db',
#                                        'org.Ce.eg.db','org.Bt.eg.db','org.Ss.eg.db','org.Mmu.eg.db',
#                                        'org.Os.eg.db'))

### GO ### KEGG ###
localdbpath = '~/mqf/Temporal/workflow/lib/lib1/'
enrich_fun = function(sampleinfo,rslt){
  if(!dir.exists(paste0(rslt,'/enrich')))(dir.create(paste0(rslt,'/enrich')))
  #if(!dir.exists(paste0(rslt,'/enrich/pics')))(dir.create(paste0(rslt,'/enrich/pics')))
  studyid <- unique(sampleinfo$Studyid)
  strategy <- unique(sampleinfo$Strategy)
  if(file.exists(paste0('~/mqf/Temporal/workflow/result/',studyid,'/','clust_result','/','Clusters_Objects.tsv'))){
    organism = gsub(' ','_',unique(sampleinfo$Organism))
    # if(organism == 'Oryza_sativa'){
    #   ####Oryza_sativa_orgDb
    #   library(AnnotationHub)
    #   hub <- AnnotationHub()
    #   #query(hub, "oryza sativa")
    #   # org.Oryza_sativa_Japonica_Group.eg.sqlite <- hub[['AH96212']]
    #   org.Os.eg.db <- hub[['AH96212']]
    #   # org.Os.eg.db <- readRDS('~/mqf/Temporal/workflow/lib/lib1/org.Os.eg.db.Rds')
    # }
    sym_enz_ens <- fread(paste0('~/mqf/Temporal/workflow/src/',organism,'_sym_enz_ens.txt'))
    # orgDb = organism_orgDb[match(organism,organism_orgDb$organism),2]
    if(organism == 'Homo_sapiens'){
      load(paste0(localdbpath,'GO_KEGG_Hs_ENTREZID_2023May20.RData'))
    }else if(organism == 'Mus_musculus'){
      load(paste0(localdbpath,'GO_KEGG_Mm_ENTREZID_2023May20.RData'))
    }else if(organism == 'Arabidopsis_thaliana'){
      load(paste0(localdbpath,'GO_KEGG_At_ENTREZID_2023May20.RData'))
    }else if(organism == 'Rattus_norvegicus'){
      load(paste0(localdbpath,'GO_KEGG_Rn_ENTREZID_2023May20.RData'))
    }else if(organism == 'Saccharomyces_cerevisiae'){
      load(paste0(localdbpath,'GO_KEGG_Sc_ENTREZID_2023May20.RData'))
    }else if(organism == 'Drosophila_melanogaster'){
      load(paste0(localdbpath,'GO_KEGG_Dm_ENTREZID_2023May20.RData'))
    }else if(organism == 'Danio_rerio'){
      load(paste0(localdbpath,'GO_KEGG_Dr_ENTREZID_2023May20.RData'))
    }else if(organism == 'Escherichia_coli'){
      load(paste0(localdbpath,'GO_KEGG_EcK12_ENTREZID_2023May20.RData'))
    }else if(organism == 'Caenorhabditis_elegans'){
      load(paste0(localdbpath,'GO_KEGG_Ce_ENTREZID_2023May20.RData'))
    }else if(organism == 'Bos_taurus'){
      load(paste0(localdbpath,'GO_KEGG_Bt_ENTREZID_2023May20.RData'))
    }else if(organism == 'Sus_scrofa'){
      load(paste0(localdbpath,'GO_KEGG_Ss_ENTREZID_2023May20.RData'))
    }else if(organism == 'Macaca_mulatta'){
      load(paste0(localdbpath,'GO_KEGG_Mmu_ENTREZID_2023May20.RData'))
    }else{
      load(paste0(localdbpath,'GO_KEGG_osa_ENTREZID_2023May20.RData'))
    }
    # all_cluster_gene = read.table(paste0(rslt,'/','clust_result','/','Clusters_Objects.tsv'),sep = '\t',header = T) %>% 
      # as.data.frame() %>% dplyr::slice(-1)
    GOenrich <- data.frame()
    KEGGenrich <- data.frame()
    if(file.exists(paste0(rslt,'/','ResGeneDF.csv'))){   ###2023.06.21mqf 加一步判断，如果拟合报错的只对管家基因富集
      all_cluster_gene = as.data.frame(read.table(paste0(rslt,'/','ResGeneDF.csv'),sep = ',',header = T,check.names=F))
      cluster_group=unique(all_cluster_gene$cluster)
      # df_exp=read.table(paste0(rslt,'/','clust_result/Processed_Data','/','tpm_profile.txt_processed.tsv'),sep = '\t',header = F)
      for(clu in cluster_group){
        genename=all_cluster_gene[which(all_cluster_gene$cluster==clu),1]
        if(length(genename) != 0 ){
          cluster_name <- gsub('_','',clu)
          # heatmapdata=df_exp[match(genename,df_exp[,1]),]
          heatmapdata=all_cluster_gene[which(all_cluster_gene$cluster==clu),-c((ncol(all_cluster_gene)-2):ncol(all_cluster_gene))]
          write_tsv(heatmapdata,paste0(rslt,'/',studyid,'_',cluster_name,"_heatmapdata.tsv"))
          
          eg <- sym_enz_ens[sym_enz_ens$Symbol %in% genename,]
          cluster_gene=eg$GeneID
          # if(strategy == 'RNA-Seq'){
          #   eg <- sym_enz_ens[sym_enz_ens$Symbol %in% genename,]
          #   cluster_gene=eg$GeneID
          # }else{
          #   cluster_gene=genename
          # }
          GOenrich_result <- data.frame()   ####2023.06.22mqf 先预定义两个数据框，后面判断是否存在 cluster0也添加
          KEGGenrich_result <- data.frame()
          tryCatch({
            # if(organism == 'Oryza_sativa'){
            #   enrich_GO_plot(cluster_gene,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,org.Os.eg.db,cluster_name)
            # }else{
            #   enrich_GO_plot(cluster_gene,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,orgDb,cluster_name)
            # }
            GOenrich_result <- enrich_GO_plot(cluster_gene,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,sym_enz_ens,cluster_name)
            # GO = read_tsv(paste0(rslt,'/enrich/',cluster_name,"_GO.tsv"))
            # BP = GO %>% filter(ONT == "BP");CC = GO %>% filter(ONT == "CC");MF = GO %>% filter(ONT == "MF");
            # p = BUBBLE(BP[1:10,],"GO")
            # ggsave(paste0(rslt,'/enrich/',"pics/",cluster_name,"_BP.png"),plot = p,width = 12,height = 8)
            # p = BUBBLE(CC[1:10,],"GO")
            # ggsave(paste0(rslt,'/enrich/',"pics/",cluster_name,"_CC.png"),plot = p,width = 12,height = 8)
            # p = BUBBLE(MF[1:10,],"GO")
            # ggsave(paste0(rslt,'/enrich/',"pics/",cluster_name,"_MF.png"),plot = p,width = 12,height = 8)
          },error=function(e){
            sink('~/mqf/Temporal/workflow/result/enrich.log',append = T)
            print(paste0(studyid,'_',cluster_name,'_GO.fail'))
            sink()
          })
          tryCatch({
            # if(organism == 'Oryza_sativa'){
            #   enrich_KEGG_plot(cluster_gene,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,org.Os.eg.db,cluster_name)
            # }else{
            #   enrich_KEGG_plot(cluster_gene,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,orgDb,cluster_name)
            # }
            KEGGenrich_result <- enrich_KEGG_plot(cluster_gene,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,sym_enz_ens,cluster_name)
            # KEGG = read_tsv(paste0(rslt,'/enrich/',cluster_name,'_KEGG.tsv'))
            # p = BUBBLE(KEGG[1:10,],"KEGG")
            # ggsave(paste0(rslt,'/enrich/',"pics/",cluster_name,"_KEGG.png"),plot =p,width = 12,height = 8)
          },error=function(e){
            sink('~/mqf/Temporal/workflow/result/enrich.log',append = T)
            print(paste0(studyid,'_',cluster_name,'_KEGG.fail'))
            sink()
          })
          if(nrow(GOenrich_result)>0){  ### 2023.06.22mqf 添加判断，如果某个cluster富集失败就跳过
            GOenrich <- rbind(GOenrich, GOenrich_result)
          }
          if(nrow(KEGGenrich_result)>0){
            KEGGenrich <- rbind(KEGGenrich, KEGGenrich_result)
          }
          
        }else{
          print('enrich fail: the most diff fail!')
        }
      }
    }
    
    #对housekeeping基因富集
    ConstantDF = as.data.frame(read.table(paste0(rslt,'/','ConstantDF.csv'),sep = ',',header = T,check.names=F))
    ConstantDF <- dplyr::select(ConstantDF,-c('Studyid','CV','mean','cluster','model'))
    #标准化，取对数，Z-score
    ConstantDF <- aggregate(.~Genes, mean, data = ConstantDF) #2023.06.17 重复基因取均值
    ConstantDF <- tibble::column_to_rownames(ConstantDF,'Genes')
    ConstantDF <- log2(ConstantDF) %>% scale(.) %>% as.data.frame()
    ConstantDF <- tibble::rownames_to_column(ConstantDF,'Genes')
    write_tsv(ConstantDF,paste0(rslt,'/',studyid,'_','cluster0',"_heatmapdata.tsv"))
    genename <- ConstantDF$Genes
    eg <- sym_enz_ens[sym_enz_ens$Symbol %in% genename,]
    cluster_gene=eg$GeneID
    
    cluster_name <- 'cluster0'
    GOenrich_result <- data.frame()   ####2023.06.22mqf 先预定义两个数据框，后面判断是否存在
    KEGGenrich_result <- data.frame()
    tryCatch({
      GOenrich_result <- enrich_GO_plot(cluster_gene,studyid,GO_TERM2GENE,GO_TERM2NAME,GO_TERM2ONT,sym_enz_ens,cluster_name)
    },error=function(e){
      sink('~/mqf/Temporal/workflow/result/enrich.log',append = T)
      print(paste0(studyid,'_','cluster0','_GO.fail'))
      sink()
    })
    tryCatch({
      KEGGenrich_result <- enrich_KEGG_plot(cluster_gene,studyid,KEGG_TERM2GENE,KEGG_TERM2NAME,sym_enz_ens,cluster_name)
    },error=function(e){
      sink('~/mqf/Temporal/workflow/result/enrich.log',append = T)
      print(paste0(studyid,'_','cluster0','_KEGG.fail'))
      sink()
    })
    
    if(nrow(GOenrich_result)>0){  ### 2023.06.22mqf 添加判断，如果某个cluster富集失败就跳过
      GOenrich <- rbind(GOenrich, GOenrich_result)
    }
    if(nrow(KEGGenrich_result)>0){
      KEGGenrich <- rbind(KEGGenrich, KEGGenrich_result)
    }

    write_tsv(GOenrich,paste0(rslt,'/enrich/',studyid,"_GO.tsv"))
    write_tsv(KEGGenrich,paste0(rslt,'/enrich/',studyid,"_KEGG.tsv"))
  }
}
    
    
