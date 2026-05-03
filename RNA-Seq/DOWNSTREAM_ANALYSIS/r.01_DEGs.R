#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/DOWNSTREAM_ANALYSIS/')
dir_name<-'01_DEGs';if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

#加载R包
library(stringr)#用于处理字符串
library(readr)
library(DESeq2)
library(tidyverse)
library(biomaRt)
library(ggrepel)
library(psych)
library(corrplot)
library(ComplexHeatmap)
library(ggsci)
library(export)
library(ggvenn)
library(GenomicFeatures)
library(tximport)
select=dplyr::select

#读取数据
countData<-read_tsv('../../results/star_salmon/salmon.merged.gene_counts.tsv')
# 检查gene_name列是否有重复的基因
length(countData$gene_name)==length(unique(countData$gene_name))
# [1] FALSE
# 检查gene_name列是否有重复的基因
length(countData$gene_id)==length(unique(countData$gene_id))
# [1] TRUE
# 调整表达矩阵
countData <- countData %>%
  dplyr::select(-'gene_name') %>%
  column_to_rownames('gene_id') #使用EnsemblID作为行名
countData <- countData[which(rowSums(countData)>1),] # 过滤，剔除低表达基因

# PCA
{
  countData_t<-t(countData)
  pca <- prcomp(countData_t, center = T, scale. = T)
  pca.data <- data.frame(pca$x)
  pca.variance <- pca$sdev^2 / sum(pca$sdev^2)
  pca.data.plot<-pca.data %>% 
    dplyr::select(c('PC1','PC2'))
  pca.data.plot$Class1<-rep(c('A','B','C'),each=6)
  pca.data.plot$Class2<-rep(c('Control','EHMT2_KD1','EHMT2_KD2',
                              'NoT_Control','NoT_KD1','NoT_KD2'),3)
  PCA_ABC<-
    ggplot(pca.data.plot, aes(x = PC1, y = PC2, color = Class1)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95) + 
    theme_bw()+
    labs(title = 'PCA - Group')+
    theme(axis.title=element_text(face = 'bold'),
          plot.title = element_text(face = 'bold',hjust=0.5))
  ggsave('1.PCA_ABCGroup.pdf',PCA_ABC,width=5,height = 5,units = 'in')
  
  PCA_Traetment<-
    ggplot(pca.data.plot, aes(x = PC1, y = PC2, color = Class2)) +
    geom_point(size = 3) +
    geom_hline(yintercept = 0) + 
    geom_vline(xintercept = 0) +
    stat_ellipse(aes(x = PC1, y = PC2), linetype = 2, size = 0.5, level = 0.95) + 
    theme_bw()+
    labs(title = 'PCA - Group')+
    theme(axis.title=element_text(face = 'bold'),
          plot.title = element_text(face = 'bold',hjust=0.5))
  ggsave('2.PCA_Traetment.pdf',PCA_Traetment,width=5,height = 5,units = 'in')
}

#分组信息
sample_group<-data.frame(
  sample=c(paste0('A_',seq(1:6)),paste0('B_',seq(1:6)),paste0('C_',seq(1:6))),
  group=rep(c('A','B','C'),each=6),
  treatment=rep(c('Control','EHMT2_KD1','EHMT2_KD2',
                  'NoT_Control','NoT_KD1','NoT_KD2'),3)
)
saveRDS(countData,'countData.RDS')
saveRDS(sample_group,'sample_group.RDS')

#------------------------------导入数据-----------------------------------------
txd <- makeTxDbFromGFF(file = "/rds/projects/g/gendood-preclinomics/EHMT2/References_Index/references/gencode.v45.annotation.gtf.gz", 
                       format = "gtf", dataSource = "gencode.v45.annotation.gtf.gz", 
                       organism = "Homo sapiens")
# 保存 TxDb 对象到本地
saveDb(txd, file="./Gencode_V45_TxDb.sqlite")
txTogene <- AnnotationDbi::select(
  txd, keys=keys(txd, "TXNAME"), keytype="TXNAME", 
  columns=c("TXNAME", "GENEID"))
sampleList <- c(paste0(rep(c('A_','B_','C_'),each=6),rep(seq(1:6),3)))
fileList <- file.path("../../results/star_salmon", 
                      sampleList, "quant.sf")
names(fileList) <- sampleList
txi <- tximport(fileList, type = "salmon", tx2gene = txTogene)
saveRDS(txi,'txi.RDS')

#-------------------------------KD1_vs_C-----------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/DOWNSTREAM_ANALYSIS/01_DEGs/')
#------------
Control_group<-'Control'
Case_Group<-'EHMT2_KD1'
#------------
dir_name<-paste0(Case_Group,'_vs_',Control_group)
if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

#提取对应数据
sample_group<-readRDS('../sample_group.RDS') %>%
  dplyr::filter(treatment %in% c(Control_group,Case_Group)) %>%
  dplyr::select(c('sample','treatment')) %>%
  clusterProfiler::rename('group'='treatment')
txi<-readRDS('../txi.RDS')

#提取对应样本的数据
{
  txi$abundance<-txi$abundance %>% 
    as.data.frame() %>% 
    dplyr::select(sample_group$sample) %>% 
    as.matrix()
  txi$counts<-txi$counts %>% 
    as.data.frame() %>% 
    dplyr::select(sample_group$sample) %>% 
    as.matrix()
  txi$length<-txi$length %>% 
    as.data.frame() %>% 
    dplyr::select(sample_group$sample) %>% 
    as.matrix()
}

#因子化分组信息
condition <- factor(c(sample_group$group),#与列进行对应
                    levels = c(Control_group,Case_Group))#对照组在前，实验组在后
#表达矩阵中列名对应的样本信息
colData <- data.frame(row.names=sample_group$sample, condition)#样本信息与分组信息匹配
#差异分析
dds <- DESeqDataSetFromTximport(txi = txi, #表达矩阵
                              colData = colData, #样本信息
                              design = ~ condition)#分组信息
#标准化
dds1 <- DESeq(dds, fitType = 'mean', #使用均值作为拟合模型
              minReplicatesForReplace = 3, #在进行差异分析时，至少需要3个样本
              parallel = T) #使用多线程运行
#将结果用result()函数来获取
res <- results(dds1)
#res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
#依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
if(length(which(is.na(res1$pvalue)))>0){#判断是否有行P值不存在
  res1<-res1[-which(is.na(res1$pvalue)),]#如果有，去除该行
}
#判断变化
logFC_cutoff<-(log2(2))
res1$Change<-NA#添加变化的列
for (i in 1:nrow(res1)) {#循环和判断语句逐个加变化标签
  if(res1$pvalue[i]<0.05){
    if(res1$log2FoldChange[i]>logFC_cutoff){
      res1$Change[i]<-'Up'
    } else if (res1$log2FoldChange[i]<(-logFC_cutoff)) { 
      res1$Change[i]<-'Down'
    } else {
      res1$Change[i]<-
        'Stable'
    }
  } else {
    res1$Change[i]<-
      'Stable'
  }
}
table(res1$Change)
# Down Stable     Up 
# 133  29152    229 
sig_diff <- subset(res1,
                   res1$pvalue < 0.05 & 
                     abs(res1$log2FoldChange) > logFC_cutoff)
write.csv(res1,'1.DESeq2_res1.csv',row.names = T)
saveRDS(res1,'1.DESeq2_res1.RDS')

res1_sig<-res1[which(res1$Change!='Stable'),]
gene_sig<-rownames(res1_sig)
write.csv(res1_sig,'2.DESeq2_res1_sig.csv',row.names = T)
saveRDS(res1_sig,'2.DESeq2_res1_sig.RDS')

#火山图
#挑选变化倍数前10
dat_rep <- rbind(head(res1_sig[order(res1_sig$log2FoldChange,decreasing = T),],10),
                 head(res1_sig[order(res1_sig$log2FoldChange,decreasing = F),],10)) %>%
  rownames_to_column('gene_name')
dat_rep <- res1_sig %>% rownames_to_column('gene_name') %>% 
  filter(gene_name %in% dat_rep$gene_name)
row.names(dat_rep) <- dat_rep$gene_name
#绘图
volcano_plot <- 
  ggplot(data = res1,aes(x = log2FoldChange,y = -log10(pvalue),color=Change)) +
  scale_color_manual(values = c("blue", "darkgray","darkorange")) +
  scale_x_continuous(breaks = c((-10),(-5),(-1),0,1,5,10),limits = c(-10,10)) +
  scale_y_continuous(trans = "log1p",limits = c(0,30),breaks=c(0,1,3,10,30,70)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1), lty = 4, col = "darkgray", lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "darkgray", lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",color="black",family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",color = "black",size = 15),
        axis.text.y = element_text(face = "bold",color = "black",size = 15),
        axis.title.x = element_text(face = "bold",color = "black",size = 15),
        axis.title.y = element_text(face = "bold",color = "black",size = 15)) +
  geom_label_repel(data = dat_rep, aes(label = gene_name),
                   max.overlaps = 20, size = 4,
                   box.padding = unit(0.5, "lines"),
                   min.segment.length = 0,
                   point.padding = unit(0.8, "lines"), segment.color = "black",
                   show.legend = FALSE )+
  labs(x = "log2(Fold Change)",y = "-log10 (P.Value)")
#输出文件
export::graph2png(volcano_plot,
                  file='./3.DEGs_volcano.png',
                  width = 8.78, height = 7.43)
export::graph2pdf(volcano_plot,
                  file='./3.DEGs_volcano.pdf',
                  width = 8.78, height = 7.43)

#密度热图
#前十FC的上下调差异基因
top10_gene_down <- rownames(sig_diff %>% 
                              arrange(log2FoldChange) %>% 
                              head(10))
top10_gene_up <- rownames(sig_diff %>% 
                            arrange(desc(log2FoldChange)) %>% 
                            head(10))
top_gene <- c(top10_gene_down,top10_gene_up )
#分组信息
countData<-txi$counts
rt <- countData[top_gene,sample_group$sample]
group<-sample_group[order(sample_group$group),]
x <- rt
mat <- t(scale(t(x)))#归一化
df1 <- as.data.frame(mat)
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2
identical(group$sample, colnames(mat))
mat <- mat[, group$sample]
identical(group$sample, colnames(mat))
#绘图
ppp <- densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",
                      height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("EHMT2_KD1" = "#B72230",
                                                              "Control" = "#104680"))) %v%
  Heatmap(mat,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression",
          ###cluster_cols = F,
          cluster_rows = F,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
#输出文件
export::graph2png(ppp, 
                  file="./4.DEGs_deg_top_heatmap.png",
                  width=6.51,height=5.4)
export::graph2pdf(ppp, 
                  file="./4.DEGs_deg_top_heatmap.pdf",
                  width=6.51,height=5.4)

DEGs_list<-list(
  UP=rownames(res1 %>%
                dplyr::filter(Change=='Up')),
  DOWN=rownames(res1 %>%
                  dplyr::filter(Change=='Down'))
)
saveRDS(DEGs_list,'5.DEGs_list.RDS')


#-------------------------------KD2vs_C-----------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/DOWNSTREAM_ANALYSIS/01_DEGs/')
#------------
Control_group<-'Control'
Case_Group<-'EHMT2_KD2'
#------------
dir_name<-paste0(Case_Group,'_vs_',Control_group)
if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

#提取对应数据
sample_group<-readRDS('../sample_group.RDS') %>%
  dplyr::filter(treatment %in% c(Control_group,Case_Group)) %>%
  dplyr::select(c('sample','treatment')) %>%
  clusterProfiler::rename('group'='treatment')
txi<-readRDS('../txi.RDS')

#提取对应样本的数据
{
  txi$abundance<-txi$abundance %>% 
    as.data.frame() %>% 
    dplyr::select(sample_group$sample) %>% 
    as.matrix()
  txi$counts<-txi$counts %>% 
    as.data.frame() %>% 
    dplyr::select(sample_group$sample) %>% 
    as.matrix()
  txi$length<-txi$length %>% 
    as.data.frame() %>% 
    dplyr::select(sample_group$sample) %>% 
    as.matrix()
}

#因子化分组信息
condition <- factor(c(sample_group$group),#与列进行对应
                    levels = c(Control_group,Case_Group))#对照组在前，实验组在后
#表达矩阵中列名对应的样本信息
colData <- data.frame(row.names=sample_group$sample, condition)#样本信息与分组信息匹配
#差异分析
dds <- DESeqDataSetFromTximport(txi = txi, #表达矩阵
                                colData = colData, #样本信息
                                design = ~ condition)#分组信息
#标准化
dds1 <- DESeq(dds, fitType = 'mean', #使用均值作为拟合模型
              minReplicatesForReplace = 3, #在进行差异分析时，至少需要3个样本
              parallel = T) #使用多线程运行
#将结果用result()函数来获取
res <- results(dds1)
#res格式转化：用data.frame转化为表格形式
res1 <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)
#依次按照pvalue值log2FoldChange值进行排序
res1 <- res1[order(res1$pvalue, res1$log2FoldChange, decreasing = c(FALSE, TRUE)), ]
if(length(which(is.na(res1$pvalue)))>0){#判断是否有行P值不存在
  res1<-res1[-which(is.na(res1$pvalue)),]#如果有，去除该行
}
#判断变化
logFC_cutoff<-(log2(2))
res1$Change<-NA#添加变化的列
for (i in 1:nrow(res1)) {#循环和判断语句逐个加变化标签
  if(res1$pvalue[i]<0.05){
    if(res1$log2FoldChange[i]>logFC_cutoff){
      res1$Change[i]<-'Up'
    } else if (res1$log2FoldChange[i]<(-logFC_cutoff)) { 
      res1$Change[i]<-'Down'
    } else {
      res1$Change[i]<-
        'Stable'
    }
  } else {
    res1$Change[i]<-
      'Stable'
  }
}
table(res1$Change)
# Down Stable     Up 
# 141  28946    156 
sig_diff <- subset(res1,
                   res1$pvalue < 0.05 & 
                     abs(res1$log2FoldChange) > logFC_cutoff)
write.csv(res1,'1.DESeq2_res1.csv',row.names = T)
saveRDS(res1,'1.DESeq2_res1.RDS')

res1_sig<-res1[which(res1$Change!='Stable'),]
gene_sig<-rownames(res1_sig)
write.csv(res1_sig,'2.DESeq2_res1_sig.csv',row.names = T)
saveRDS(res1_sig,'2.DESeq2_res1_sig.RDS')

#火山图
#挑选变化倍数前10
dat_rep <- rbind(head(res1_sig[order(res1_sig$log2FoldChange,decreasing = T),],10),
                 head(res1_sig[order(res1_sig$log2FoldChange,decreasing = F),],10)) %>%
  rownames_to_column('gene_name')
dat_rep <- res1_sig %>% rownames_to_column('gene_name') %>% 
  filter(gene_name %in% dat_rep$gene_name)
row.names(dat_rep) <- dat_rep$gene_name
#绘图
volcano_plot <- 
  ggplot(data = res1,aes(x = log2FoldChange,y = -log10(pvalue),color=Change)) +
  scale_color_manual(values = c("blue", "darkgray","darkorange")) +
  scale_x_continuous(breaks = c((-10),(-5),(-1),0,1,5,10),limits = c(-10,10)) +
  scale_y_continuous(trans = "log1p",limits = c(0,30),breaks=c(0,1,3,10,30,70)) +
  geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept = c(-1,1), lty = 4, col = "darkgray", lwd = 0.6)+
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "darkgray", lwd = 0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",color="black",family = "Times",
                                   size=13),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(face = "bold",color = "black",size = 15),
        axis.text.y = element_text(face = "bold",color = "black",size = 15),
        axis.title.x = element_text(face = "bold",color = "black",size = 15),
        axis.title.y = element_text(face = "bold",color = "black",size = 15)) +
  geom_label_repel(data = dat_rep, aes(label = gene_name),
                   max.overlaps = 20, size = 4,
                   box.padding = unit(0.5, "lines"),
                   min.segment.length = 0,
                   point.padding = unit(0.8, "lines"), segment.color = "black",
                   show.legend = FALSE )+
  labs(x = "log2(Fold Change)",y = "-log10 (P.Value)")
#输出文件
export::graph2png(volcano_plot,
                  file='./3.DEGs_volcano.png',
                  width = 8.78, height = 7.43)
export::graph2pdf(volcano_plot,
                  file='./3.DEGs_volcano.pdf',
                  width = 8.78, height = 7.43)

#密度热图
#前十FC的上下调差异基因
top10_gene_down <- rownames(sig_diff %>% 
                              arrange(log2FoldChange) %>% 
                              head(10))
top10_gene_up <- rownames(sig_diff %>% 
                            arrange(desc(log2FoldChange)) %>% 
                            head(10))
top_gene <- c(top10_gene_down,top10_gene_up )
#分组信息
countData<-txi$counts
rt <- countData[top_gene,sample_group$sample]
group<-sample_group[order(sample_group$group),]
x <- rt
mat <- t(scale(t(x)))#归一化
df1 <- as.data.frame(mat)
mat[mat < (-2)] <- (-2)
mat[mat > 2] <- 2
identical(group$sample, colnames(mat))
mat <- mat[, group$sample]
identical(group$sample, colnames(mat))
#绘图
ppp <- densityHeatmap(mat ,title = "Distribution as heatmap", ylab = " ",
                      height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = group$group, col = list(Group = c("EHMT2_KD2" = "#B72230",
                                                              "Control" = "#104680"))) %v%
  Heatmap(mat,
          row_names_gp = gpar(fontsize = 7),
          show_column_names = F,
          show_row_names = T,
          ###show_colnames = FALSE,
          name = "expression",
          ###cluster_cols = F,
          cluster_rows = F,
          height = unit(6, "cm"),
          #cluster_columns = FALSE,
          ###cluster_rows = FALSE,
          col = colorRampPalette(c("#0A878D", "white","#D80305"))(100))
#输出文件
export::graph2png(ppp, 
                  file="./4.DEGs_deg_top_heatmap.png",
                  width=6.51,height=5.4)
export::graph2pdf(ppp, 
                  file="./4.DEGs_deg_top_heatmap.pdf",
                  width=6.51,height=5.4)

DEGs_list<-list(
  UP=rownames(res1 %>%
                dplyr::filter(Change=='Up')),
  DOWN=rownames(res1 %>%
                  dplyr::filter(Change=='Down'))
)
saveRDS(DEGs_list,'5.DEGs_list.RDS')


#-------------------------------Combine-----------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/RNASeq_data_Saad/DOWNSTREAM_ANALYSIS/01_DEGs/')

library(ggVennDiagram)

CommonDEGs<-list()
CommonDEGs_plot<-list()
UnionDEGs<-list()

for (i in 1:length(c('UP','DOWN'))) {
  KD1_sigGenes<-readRDS('./EHMT2_KD1_vs_Control/5.DEGs_list.RDS')
  KD1_sigGenes<-KD1_sigGenes[[which(names(KD1_sigGenes)==c('UP','DOWN')[i])]]
  KD2_sigGenes<-readRDS('./EHMT2_KD2_vs_Control/5.DEGs_list.RDS')
  KD2_sigGenes<-KD2_sigGenes[[which(names(KD2_sigGenes)==c('UP','DOWN')[i])]]
  CommonDEGs_plot[[i]]<-
  ggVennDiagram(x=list(KD1=KD1_sigGenes,KD2=KD2_sigGenes),
                color = 1, lwd = 0.7,
                label_font = "Times") + 
    labs(title = paste0(c('UP','DOWN')[i],'-Regulated Genes'))+
    scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
    theme(legend.position = "none",
          plot.title = element_text(face='bold',hjust=0.5))+
    coord_flip()
  CommonDEGs[[i]]<-intersect(KD1_sigGenes,KD2_sigGenes)
  UnionDEGs[[i]]<-c(KD1_sigGenes,KD2_sigGenes)
  names(CommonDEGs)[[i]]<-c('UP','DOWN')[i]
  names(UnionDEGs)[[i]]<-c('UP','DOWN')[i]
}

ggsave('3.CommonDEGs_KDvsC.pdf',
       patchwork::wrap_plots(CommonDEGs_plot,nrow=1),
       width = 7.5,height = 3,units = 'in')
ggsave('3.CommonDEGs_KDvsC.svg',
       patchwork::wrap_plots(CommonDEGs_plot,nrow=1),
       width = 7.5,height = 3,units = 'in')
saveRDS(CommonDEGs,'4.CommonDEGs.RDS')
saveRDS(UnionDEGs,'5.UnionDEGs.RDS')



