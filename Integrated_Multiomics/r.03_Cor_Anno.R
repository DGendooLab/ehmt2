#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/Integrated_Multiomics/')
dir_name<-'03_Cor_Anno';if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

library(corrplot)
library(ggplot2)
library(ggpubr)
library(readr)
library(dplyr)

#OpenUP
{
  #读取候选基因
  candg<-readRDS('../01_CommonGenes/3.Union_DAAOpen-DEGsUP.RDS')
  #读取表达矩阵
  td<-read_tsv('../../RNASeq_data_Saad/results/star_salmon/salmon.merged.gene_tpm.tsv') %>%
    dplyr::select(c('gene_id','gene_name',
                    paste0(rep(c("A_","B_","C_"),each=3),rep(seq(1:3),3))))
  #筛选矩阵行
  td<-rbind(
    td %>% 
      dplyr::filter(gene_name %in% c('EHMT2')) %>% 
      dplyr::select(-'gene_id') %>%
      column_to_rownames('gene_name'),
    td %>% 
      dplyr::filter(gene_id %in% candg) %>% 
      dplyr::select(-'gene_name') %>%
      column_to_rownames('gene_id')
  ) %>% t() %>% as.data.frame()
  #相关性系数矩阵
  tdc <- cor (td, method="pearson")
  #显著性
  testRes <- cor.mtest(td, method="pearson",conf.level = 0.95)
  
  #热图
  corrplot(tdc, method = "ellipse", type = "upper",p.mat = testRes$p,
           tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
           insig = 'label_sig')
  corrplot(tdc, method = "number", type = "lower",
           tl.col = "n", tl.cex = 0.8, tl.pos = "n",
           add = T)
  corPlot1<-recordPlot()
  export::graph2pdf(corPlot1,'1.CorHeatmap.pdf_OpenUP',height=9,width=12)
  
  #与EHMT2显著相关的基因
  SigGeneCorEHMT2<-rownames(testRes$p %>%
                              as.data.frame() %>%
                              dplyr::filter(EHMT2<0.05))
  SigGeneCorEHMT2<-SigGeneCorEHMT2[which(SigGeneCorEHMT2!='EHMT2')]
  #转换ID
  library(biomaRt)
  # 建立数据库链接，ensembl
  mart<-useMart("ensembl")
  mydataset <- useDataset("hsapiens_gene_ensembl",mart = mart) # 人数据库
  ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name",
                              "hgnc_symbol","uniprot_gn_symbol"),
               filters = "ensembl_gene_id_version",
               values = SigGeneCorEHMT2,
               mart = mydataset)  # attributes为输出ID类型
  write.csv(ids,'2.SigCorGenesWithEHMT2_OpenUP.csv',row.names = F)
  
  # 建立数据库链接，ensembl
  mart<-useMart("ensembl")
  mydataset <- useDataset("hsapiens_gene_ensembl",mart = mart) # 人数据库
  ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name",
                              "hgnc_symbol","uniprot_gn_symbol"),
               filters = "ensembl_gene_id_version",
               values = candg,
               mart = mydataset)  # attributes为输出ID类型
  write.csv(ids,'3.CandgSymbol_OpenUP.csv',row.names = F)
  
}

#CloseDOWN
{
  #读取候选基因
  candg<-readRDS('../01_CommonGenes/6.Union_DAAClose-DEGsDOWN.RDS')
  #读取表达矩阵
  td<-read_tsv('../../RNASeq_data_Saad/results/star_salmon/salmon.merged.gene_tpm.tsv') %>%
    dplyr::select(c('gene_id','gene_name',
                    paste0(rep(c("A_","B_","C_"),each=3),rep(seq(1:3),3))))
  #筛选矩阵行
  td<-rbind(
    td %>% 
      dplyr::filter(gene_name %in% c('EHMT2')) %>% 
      dplyr::select(-'gene_id') %>%
      column_to_rownames('gene_name'),
    td %>% 
      dplyr::filter(gene_id %in% candg) %>% 
      dplyr::select(-'gene_name') %>%
      column_to_rownames('gene_id')
  ) %>% t() %>% as.data.frame()
  #相关性系数矩阵
  tdc <- cor (td, method="pearson")
  #显著性
  testRes <- cor.mtest(td, method="pearson",conf.level = 0.95)
  
  #热图
  corrplot(tdc, method = "ellipse", type = "upper",p.mat = testRes$p,
           tl.col = "black", tl.cex = 0.8, tl.srt = 45,tl.pos = "lt",
           sig.level = c(0.001, 0.01, 0.05), pch.cex = 1.2,
           insig = 'label_sig')
  corrplot(tdc, method = "number", type = "lower",
           tl.col = "n", tl.cex = 0.8, tl.pos = "n",
           add = T)
  corPlot1<-recordPlot()
  export::graph2pdf(corPlot1,'1.CorHeatmap.pdf_CloseDOWN',height=9,width=12)
  
  #与EHMT2显著相关的基因
  SigGeneCorEHMT2<-rownames(testRes$p %>%
                              as.data.frame() %>%
                              dplyr::filter(EHMT2<0.05))
  SigGeneCorEHMT2<-SigGeneCorEHMT2[which(SigGeneCorEHMT2!='EHMT2')]
  #转换ID
  library(biomaRt)
  # 建立数据库链接，ensembl
  mart<-useMart("ensembl")
  mydataset <- useDataset("hsapiens_gene_ensembl",mart = mart) # 人数据库
  ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name",
                              "hgnc_symbol","uniprot_gn_symbol"),
               filters = "ensembl_gene_id_version",
               values = SigGeneCorEHMT2,
               mart = mydataset)  # attributes为输出ID类型
  write.csv(ids,'2.SigCorGenesWithEHMT2_CloseDOWN.csv',row.names = F)
  
  # 建立数据库链接，ensembl
  mart<-useMart("ensembl")
  mydataset <- useDataset("hsapiens_gene_ensembl",mart = mart) # 人数据库
  ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name",
                              "hgnc_symbol","uniprot_gn_symbol"),
               filters = "ensembl_gene_id_version",
               values = candg,
               mart = mydataset)  # attributes为输出ID类型
  write.csv(ids,'3.CandgSymbol_CloseDOWN.csv',row.names = F)
  
}