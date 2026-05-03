#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/Integrated_Multiomics/')
dir_name<-'01_CommonGenes';if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

library(ggVennDiagram)
library(ggplot2)

#bluebear的X11选项加载
options(bitmapType='cairo')

#Common
{
  DAA_Genes_Open<-readRDS('../../ATAC_Seq/DAAs/MethodCross/2.MethodCross_Open.RDS')
  DEGs_UP<-readRDS('../../RNA_Seq/DOWNSTREAM_ANALYSIS/01_DEGs/4.CommonDEGs.RDS')$UP
  CommonVenn<-
  ggVennDiagram(list(`ATAC-Seq_Open`=DAA_Genes_Open,`RNA-Seq_Up`=DEGs_UP), color = 1, lwd = 0.7,
                label_alpha = 0) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    labs(title = 'ATAC-RNA (Common)')+
    theme(legend.position = "none",
          plot.title = element_text(face='bold',hjust=0.5))+
    coord_flip()
  ggsave('1.CommonVenn_OpenUP.pdf',CommonVenn,width = 5,height = 3.75,units = 'in')
  ggsave('1.CommonVenn_OpenUP.svg',CommonVenn,width = 5,height = 3.75,units = 'in')
  
  DAA_Genes_Close<-readRDS('../../ATAC_Seq/DAAs/MethodCross/5.MethodCross_Close.RDS')
  DEGs_DOWN<-readRDS('../../RNA_Seq/DOWNSTREAM_ANALYSIS/01_DEGs/4.CommonDEGs.RDS')$DOWN
  CommonVenn<-
    ggVennDiagram(list(`ATAC-Seq_Close`=DAA_Genes_Close,`RNA-Seq_Down`=DEGs_DOWN), color = 1, lwd = 0.7,
                  label_alpha = 0) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    labs(title = 'ATAC-RNA (Common)')+
    theme(legend.position = "none",
          plot.title = element_text(face='bold',hjust=0.5))+
    coord_flip()
  ggsave('2.CommonVenn_CloseDOWN.pdf',CommonVenn,width = 5,height = 3.75,units = 'in')
  ggsave('2.CommonVenn_CloseDOWN.svg',CommonVenn,width = 5,height = 3.75,units = 'in')
}

#Union
{
  DAA_Genes_Open<-readRDS('../../ATAC_Seq/DAAs/MethodCross/3.UnionMethodsDAA_Open.RDS')
  DEGs_UP<-readRDS('../../RNA_Seq/DOWNSTREAM_ANALYSIS/01_DEGs/5.UnionDEGs.RDS')$UP
  UnionVenn<-
    ggVennDiagram(list(`ATAC-Seq_Open`=DAA_Genes_Open,`RNA-Seq`=DEGs_UP), color = 1, lwd = 0.7,
                  label_alpha = 0) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    labs(title = 'ATAC-RNA (Union)')+
    theme(legend.position = "none",
          plot.title = element_text(face='bold',hjust=0.5))+
    coord_flip()
  ggsave('3.UnionVenn_OpenUp.pdf',UnionVenn,width = 5,height = 3.75,units = 'in')
  ggsave('3.UnionVenn_OpenUp.svg',UnionVenn,width = 5,height = 3.75,units = 'in')
  
  saveRDS(intersect(DAA_Genes_Open,DEGs_UP),'3.Union_DAAOpen-DEGsUP.RDS')
  write.csv(data.frame(gene=intersect(DAA_Genes_Open,DEGs_UP)),
            '4.Union_DAAOpen-DEGsUP.csv',row.names = F)
  library(biomaRt)
  # 建立数据库链接，ensembl
  mart<-useMart("ensembl")
  mydataset <- useDataset("hsapiens_gene_ensembl",mart = mart) # 人数据库
  ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name",
                              "hgnc_symbol","uniprot_gn_id","uniprot_gn_symbol"),
               filters = "ensembl_gene_id_version",
               values = intersect(DAA_Genes_Open,DEGs_UP),
               mart = mydataset)  # attributes为输出ID类型
  write.csv(ids,'5.Union_DAAopen-DEGsUP_MultiIDs.csv',row.names = F)
  
  
  DAA_Genes_Close<-readRDS('../../ATAC_Seq/DAAs/MethodCross/6.UnionMethodsDAA_Close.RDS')
  DEGs_DOWN<-readRDS('../../RNA_Seq/DOWNSTREAM_ANALYSIS/01_DEGs/5.UnionDEGs.RDS')$DOWN
  UnionVenn<-
    ggVennDiagram(list(`ATAC-Seq_Close`=DAA_Genes_Close,`RNA-Seq`=DEGs_DOWN), color = 1, lwd = 0.7,
                  label_alpha = 0) + 
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    labs(title = 'ATAC-RNA (Union)')+
    theme(legend.position = "none",
          plot.title = element_text(face='bold',hjust=0.5))+
    coord_flip()
  ggsave('6.UnionVenn_CloseDown.pdf',UnionVenn,width = 5,height = 3.75,units = 'in')
  ggsave('6.UnionVenn_CloseDown.svg',UnionVenn,width = 5,height = 3.75,units = 'in')
  saveRDS(intersect(DAA_Genes_Close,DEGs_DOWN),'6.Union_DAAClose-DEGsDOWN.RDS')
  write.csv(data.frame(gene=intersect(DAA_Genes_Close,DEGs_DOWN)),
            '7.Union_DAAClose-DEGsDOWN.csv',row.names = F)
  library(biomaRt)
  # 建立数据库链接，ensembl
  mart<-useMart("ensembl")
  mydataset <- useDataset("hsapiens_gene_ensembl",mart = mart) # 人数据库
  ids <- getBM(attributes = c("ensembl_gene_id_version","external_gene_name",
                              "hgnc_symbol","uniprot_gn_id","uniprot_gn_symbol"),
               filters = "ensembl_gene_id_version",
               values = intersect(DAA_Genes_Close,DEGs_DOWN),
               mart = mydataset)  # attributes为输出ID类型
  write.csv(ids,'8.Union_DAAClose-DEGsDOWN_MultiIDs.csv',row.names = F)
}


