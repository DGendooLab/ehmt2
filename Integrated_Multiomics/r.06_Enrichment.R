#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/Integrated_Multiomics/')
folder_name<-"06_Enrichment";if (!dir.exists(folder_name)) {dir.create(folder_name)}
setwd(folder_name)
# 加载R包
suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(org.Hs.eg.db))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GOplot))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
library(export)
library(circlize)
library(grid)
library(graphics)
library(ComplexHeatmap)
select=dplyr::select
#候选基因 Candidate Genes
candg_up<-unique(read.csv('../03_Cor_Anno/3.CandgSymbol_OpenUP.csv')$ensembl_gene_id_version)
candg_up<-candg_up[which(candg_up!='')]
candg_down<-unique(read.csv('../03_Cor_Anno/3.CandgSymbol_CloseDOWN.csv')$ensembl_gene_id_version)
candg_down<-candg_down[which(candg_down!='')]
candg<-c(candg_up,candg_down) %>% str_sub(1,15)
## id转换 Convert Gene ID
symbol2entrezid <- bitr(geneID = candg,
                        fromType = 'ENSEMBL',
                        toType = 'ENTREZID',
                        OrgDb = 'org.Hs.eg.db')
# GO富集分析 GO Enrichment####---------------------
pvalueCutoff <- 1
qvalueCutoff <- 1

if(!file.exists('ego.rds')){
  ego <- enrichGO(gene = as.numeric(symbol2entrezid$ENTREZID),
                  keyType = "ENTREZID",
                  OrgDb = org.Hs.eg.db, 
                  pvalueCutoff = pvalueCutoff, 
                  qvalueCutoff = qvalueCutoff,
                  ont="ALL",
                  readable =T)
  saveRDS(ego,'ego.rds')
}else{
  ego<-readRDS('ego.rds')
}
as.data.frame(ego) %>% filter(pvalue <= 0.05) %>% 
  group_by(ONTOLOGY) %>% dplyr::count()
# # A tibble: 3 × 2
# # Groups:   ONTOLOGY [3]
# ONTOLOGY     n
# <chr>    <int>
#   1 BP         116
# 2 CC          18
# 3 MF          21
write_csv(as.data.frame(ego) %>% filter(pvalue <= 0.05),"1.Rich_GO_enrich_sig.csv")
go.df<-as.data.frame(ego) %>% filter(pvalue <= 0.05)
go.df <- go.df %>% group_by(ONTOLOGY) %>% slice_head(n=5)
# 使画出的GO term的顺序与输入一致
# Make sure the order in plot is same with frame
go.df$Description <- factor(go.df$Description,levels = rev(go.df$Description))
# 绘图 plot
GO_Plot<-
  ggplot(data = go.df)+ # 绘图使用的数据
  geom_point(aes(x = Description, y=(-log10(pvalue)), 
                 size = Count,color = ONTOLOGY))+
  scale_color_manual(values = c("#0000CD","orange","#43CD80"))+
  coord_flip()+# 横纵坐标反转
  theme_bw()+ #去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "GO terms",y = paste0('-log10(P)'))+ # 设置坐标轴标题及标题
  ggtitle('GO Enrichment')+
  ggplot2::theme(axis.title = element_text(size = 13), # 坐标轴标题大小
                 axis.text = element_text(size = 11), # 坐标轴标签大小
                 plot.title = element_text(size = 10,hjust = 0.5,face = "bold"), # 标题设置
                 legend.title = element_text(size = 10), # 图例标题大小
                 legend.text = element_text(size = 10) # 图例标签大小
  )
ggsave('2.GO_bubble.pdf',GO_Plot,width = 10,height=6)
ggsave('2.GO_bubble.svg',GO_Plot,width = 10,height=6)

# KEGG富集分析 KEGG enrichment####---------------------
if(!file.exists('ego.rds')){
  ekegg <- enrichKEGG(gene = symbol2entrezid$ENTREZID ,
                      keyType = "kegg",
                      organism = "hsa",
                      pAdjustMethod = "BH",
                      pvalueCutoff = 1,
                      qvalueCutoff = 1)
  saveRDS(ekegg,'ekegg.rds')
}else{
  ekegg <- readRDS('ekegg.rds')
}

ekegg2 <- setReadable(ekegg, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
write_csv(as.data.frame(ekegg2),"3.KEGG_enrich.csv")

kegg.df <- ekegg2[order(ekegg2$pvalue),]
kegg.df$Y_Axis_Value<-kegg.df[,which(colnames(kegg.df)=='pvalue')]
#计算数值型GeneRatio Colculate GeneRatio
kegg.df$GeneRatio_Number<-NA
for (i in 1:nrow(kegg.df)) {
  kegg.df$GeneRatio_Number[i]<-
    as.numeric(
      str_sub(kegg.df$GeneRatio[i],1,str_locate_all(kegg.df$GeneRatio[i],'/')[[1]][1,1]-1)
    )/as.numeric(
      str_sub(kegg.df$GeneRatio[i],
              str_locate_all(kegg.df$GeneRatio[i],'/')[[1]][1,1]+1,
              nchar(kegg.df$GeneRatio[i]))
    )
}

kegg_df_top<-kegg.df %>% 
  slice_head(n=10)
# 使画出的kegg term的顺序与输入一致
kegg_df_top$Description <- factor(kegg_df_top$Description,levels = rev(kegg_df_top$Description))
# 绘图
KEGG_Plot<-
  ggplot(data = kegg_df_top)+ # 绘图使用的数据
  geom_point(aes(x = Description, y=GeneRatio_Number, size = Count,color = (-log10(pvalue))))+
  scale_color_gradient(low="blue",high="red")+
  coord_flip()+# 横纵坐标反转
  theme_bw()+ #去除背景色
  scale_x_discrete(labels = function(x) str_wrap(x,width = 50))+ # 设置term名称过长时换行
  labs(x = "KEGG terms",y = "GeneRatio", 
       color = paste0('-log10(P)'))+ # 设置坐标轴标题及标题
  ggtitle('KEGG Enrichment')+
  ggplot2::theme(axis.title = element_text(size = 13), # 坐标轴标题大小
                 axis.text = element_text(size = 11), # 坐标轴标签大小
                 plot.title = element_text(size = 10,hjust = 0.5,face = "bold"), # 标题设置
                 legend.title = element_text(size = 10), # 图例标题大小
                 legend.text = element_text(size = 10) # 图例标签大小
  )
ggsave('4.KEGG_bubble.pdf',KEGG_Plot,width = 10,height=6)
ggsave('4.KEGG_bubble.svg',KEGG_Plot,width = 10,height=6)