#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls())
setwd('/home/liuzhaoshuo/UoB_Work/HCC_scRNA/')

library(tidyverse);library(data.table);library(vroom);library('Seurat');library(uwot);
library(Matrix);library(cowplot);library(data.table);library(magrittr);
library(RCurl);library(purrr);library(pheatmap);library(RColorBrewer);library(grid);
library(gtable);library(scater);library(rsvd);library(Rtsne);
library(cowplot);library(SingleCellExperiment);library(harmony);library(clustree);
library(Seurat);library(SeuratData);library(cowplot);library(ggsci);library(export);library(patchwork)
library(SingleR);library(devEMF)
library(celldex)
library(reshape2)
library(IOBR)
library(colorspace)
library(openxlsx)
library(CellChat)
library(SeuratData)
library(RcppML)
library(showtext)
library(progress)
select=dplyr::select
set.seed(123)  #设置随机数种子，使结果可重复
options(timeout=999999)

#绘图风格
theme.set = theme(
  axis.title = element_text(size = 20, face = "bold", family = "Arial"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Arial"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Arial"),
  legend.text = element_text(size = 16, face = "bold", family = "Arial"),
  legend.title = element_text(size = 18,face='bold',family = "Arial"),
  text = element_text(family = "Arial"))


sample_cols <- c(
  "#E31A1C", "#1F78B4", "#B2DF8A", "#FF7F00", "#6A3D9A", "#FFFF99", 
  "#33A02C", "#FB9A99", "#A6CEE3", "#FDBF6F", "#CAB2D6", "#000000", 
  "#8DD3C7", "#E6550D", "#9E9AC8", "#31A354", "#FD8D3C", "#BDBDBD", 
  "#756BB1", "#A65628", "#F781BF", "#3182BD", "#FFFF33", "#D9D9D9", 
  "#636363", "#00FFFF", "#843C39", "#9C9EDE", "#8CA252", "#E7969C", 
  "#393B79", "#E7BA52", "#CE6DBD", "#8C564B", "#17BECF", "#CEDB9C", 
  "#A55194", "#BCBD22", "#C6DBEF", "#D6616B", "#637939", "#FDD0A2", 
  "#5254A3", "#AD494A", "#C7E9C0", "#BD9E39", "#8C6D31", "#DADAEB", 
  "#7B4173", "#B5CF6B", "#969696", "#FF00FF", "#008080", "#FFD700", 
  "#800000", "#00FF00", "#000080", "#FF1493", "#0000FF", "#808000"  
)

save.image('./env.RData')

#---------------------------------读取数据GSE149614--------------------------------------
#整理数据
##建立文件夹
datadir_name<-'Data';if (!dir.exists(datadir_name)) {dir.create(datadir_name)}
setwd(datadir_name)

#创建单细胞对象
##下载数据
if(!file.exists('GSE149614_HCC.metadata.updated.txt.gz')){
  download.file(url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2Emetadata%2Eupdated%2Etxt%2Egz',
                destfile = 'GSE149614_HCC.metadata.updated.txt.gz')}
if(!file.exists('GSE149614_HCC.scRNAseq.S71915.count.txt.gz')){
  download.file(url = 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE149614&format=file&file=GSE149614%5FHCC%2EscRNAseq%2ES71915%2Ecount%2Etxt%2Egz',
                destfile = 'GSE149614_HCC.scRNAseq.S71915.count.txt.gz')}
##读取meta
meta <- fread('./GSE149614_HCC.metadata.updated.txt.gz')
##读取表达矩阵
scdata<-fread("./GSE149614_HCC.scRNAseq.S71915.count.txt.gz", data.table = FALSE)
rownames(scdata)<-scdata[,1]
scdata<-scdata[,-1]
scdata<- as.matrix(scdata)
scdata <- as(scdata, "dgCMatrix")

##创建seurat对象
sce1 <- Seurat::CreateSeuratObject(
  counts = scdata, 
  assay = "RNA",
  min.cells = 3,
  # min.features = 400,
  meta.data = meta,
  project = 'GSE149614')
sce1[["Batch"]] <-  'GSE149614'
##保存
saveRDS(sce1,'scRNA.RDS')

#---------------------------------------质控QC----------------------------------------
setwd('/home/liuzhaoshuo/UoB_Work/HCC_scRNA/')
if (!dir.exists("1_scRNA_QC")) {dir.create("1_scRNA_QC")}; setwd("1_scRNA_QC")
All_Data <- readRDS('../Data/scRNA.RDS')
#评价
# Mitochondrial 线粒体基因
All_Data <- PercentageFeatureSet(All_Data, "^MT-", col.name = "percent_mito")
# Ribosomal 核糖体基因
All_Data <- PercentageFeatureSet(All_Data, "^RP[SL]", col.name = "percent_ribo")
# Percentage hemoglobin genes - includes all genes starting with HB except HBP.
##血红蛋白基因百分比-包括除HBP外的所有以HB开头的基因。
All_Data <- PercentageFeatureSet(All_Data, "^HB[^(P|E|S)]", col.name = "percent_hb")
# Percentage for some platelet markers
#某些血小板标志物的百分比
All_Data <- PercentageFeatureSet(All_Data, "PECAM1|PF4", col.name = "percent_plat")
#评价项目
feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")

#小提琴图
p1<-VlnPlot(All_Data, group.by = "sample", features = feats, 
            pt.size = 0, ncol = 3)
p1
pdf('1-1.QC_violin_before.pdf',
    width = 15,height = 6)
p1
dev.off()
#散点图
p2<-FeatureScatter(All_Data, "nCount_RNA", "nFeature_RNA", 
                   group.by = "sample", pt.size = .5)
p2
pdf('1-2.QC_point_before.pdf',
    width = 10,height = 8)
p2
dev.off()

#过滤
#基于检测的过滤
#至少在3个细胞表达的基因和至少有200个基因表达的细胞
alldata<-All_Data
selected_c <- WhichCells(alldata, expression = nFeature_RNA >= 200)
selected_f <- rownames(alldata)[Matrix::rowSums(alldata[["RNA"]]$counts) >= 3]
data.filt <- subset(alldata, features = selected_f, cells = selected_c)
dim(data.filt)
# [1] 25479 71915

#细胞内基因数量 ≤ 4000、基因表达count数 ≤ 60000
data.filt <- subset(data.filt, cells=WhichCells(data.filt, expression = nFeature_RNA <= 4000))
data.filt <- subset(data.filt, cells=WhichCells(data.filt, expression = nCount_RNA <= 60000))

#基于核糖体和线粒体基因的过滤
#线粒体含量百分比 < 25%
data.filt <- subset(data.filt, percent_mito < 25)
dim(data.filt)
# [1] 25479 63778

feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito")
#小提琴图
p3<-VlnPlot(data.filt, group.by = "orig.ident", features = feats, pt.size = 0,
            ncol = 3) + NoLegend()
p3
pdf('1-3.QC_violin_post.pdf',
    width = 15,height = 6)
p3
dev.off()

saveRDS(data.filt,'data.filt.RDS')


#-----------------------------------降维、聚类和注释----------------------------
rm(list = ls());gc()
setwd('/home/liuzhaoshuo/UoB_Work/HCC_scRNA/')
load('./env.RData')
if (!dir.exists("2_scRNA_Dimension-Cluster_Anno")) {
  dir.create("2_scRNA_Dimension-Cluster_Anno")}
setwd("2_scRNA_Dimension-Cluster_Anno")

#读取过滤后数据
alldata<-readRDS('../Data/scRNA.RDS')

# 标准化，HVG
alldata <- alldata %>% 
  NormalizeData(normalization.method = "LogNormalize",scale.factor = 10000) %>% 
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000, 
                       verbose = FALSE, assay = "RNA") %>% 
  ScaleData()

#top20
top20 <- head(VariableFeatures(alldata), 20)
p4<-LabelPoints(plot = VariableFeaturePlot(alldata), points = top20, 
                repel = TRUE,max.overlaps=100)
p4
ggsave('2-1.HVG.pdf',width = 10,height = 8,p4)

#Z分数转换缩放基因
alldata <- ScaleData(alldata, 
                     vars.to.regress = c("percent_mito", "nFeature_RNA"), 
                     assay = "RNA")

#PCA
alldata <- RunPCA(alldata, npcs = 100, verbose = F)
print(alldata[["pca"]], dims = 1:5, nfeatures = 5)
#绘制第一个主分量
PCA_Plot<-DimPlot(alldata, reduction = "pca", 
                  group.by = "orig.ident", dims = 1:2)+
  theme.set
PCA_Plot
export::graph2pdf(PCA_Plot,'2-2.PCA.pdf',w=7.35,h=5.9)

#碎石图
ElbowPlot(alldata,reduction="pca",ndims = 100) +theme.set #根据ElbowPlot确定PC数量
ElbowPlot <- recordPlot()
export::graph2pdf(ElbowPlot,'2-3.ElbowPlot.pdf',w=8,h=6)

# if(!file.exists('alldata.JS.RData')){
#   alldata.JS <- JackStraw(alldata, num.replicate = 100, dims = 50) #运行速度慢，影响JackStrawPlot
#   save(alldata.JS,file='alldata.JS.RData')
# }else{
#   load('alldata.JS.RData')
# }
# alldata <- ScoreJackStraw(alldata.JS, dims = 1:50)

# #JackStrawPlot
# JackStrawPlot(object = alldata, dims = 1:50,reduction = "pca",xmax = 0.005) + 
#   theme(
#     axis.title = element_text(size = 20, face = "bold", family = "Arial"),
#     axis.text.x = element_text(size = 14,  face = "bold", family = "Arial"),
#     axis.text.y = element_text(size = 14,  face = "bold", family = "Arial"),
#     legend.text = element_text(size = 9, face = "bold", family = "Arial"),
#     legend.title = element_text(size = 12,face='bold',family = "Arial"),
#     text = element_text(family = "Arial"))+
#   guides(color = guide_legend(ncol = 3))+  # 将图例分为两列
#   showtext_auto()
# ggsave("2-4.JackStrawPlot.pdf", width = 12, height = 7.24)

# 聚类
pcSelect <- 25
alldata <- alldata %>% FindNeighbors(dims = 1:pcSelect) 
# 选择合适的分辨率
# 从0.1-4的resolution结果均运行一遍
seq = seq(0.1,4,by=0.1)
for (res in seq) {
  alldata = FindClusters(alldata, resolution = res)
}
# 画图
clustree <- clustree(alldata,prefix = 'RNA_snn_res.') + 
  coord_flip()+
  guides(color = guide_legend(ncol = 2))  # 将图例分为两列
clustree
ggsave('2-5.clustree.pdf',clustree,w=30,h=14)


# 根据clustree图选择res，在此res上增加并不会明显增加分群数目
alldata <- FindClusters(object = alldata, resolution = 0.5)   

# UMAP TSNE ----
alldata <- RunUMAP(alldata, dims = 1:pcSelect)
alldata <- RunTSNE(alldata, dims = 1:pcSelect) 
num_dir <- '2-6.umap_tsne'; if (!dir.exists(num_dir)) {dir.create(num_dir)}

# UMAP
## 按簇
UMAP_p <- DimPlot(alldata, reduction = "umap")+
  scale_color_manual(values = sample_cols)
UMAP_p
export::graph2pdf(UMAP_p,"2-6.umap_tsne/1.DimPlot_umap_by_cluster.pdf",width = 7.72,height = 6.08)
export::graph2svg(UMAP_p,"2-6.umap_tsne/1.DimPlot_umap_by_cluster.svg",width = 7.72,height = 6.08)
emf("2-6.umap_tsne/1.DimPlot_umap_by_cluster.emf", width = 7.72,height = 6.08)
print(UMAP_p)
dev.off()

## 按样本
UMAP_p2 <- DimPlot(alldata, group.by = "sample", 
                   reduction = "umap")+
  scale_color_manual(values = sample_cols)+
  labs(title=NULL)
UMAP_p2
export::graph2pdf(UMAP_p2,"2-6.umap_tsne/2.DimPlot_umap_by_sample.pdf",width = 7.72,height = 6.08)
export::graph2svg(UMAP_p2,"2-6.umap_tsne/2.DimPlot_umap_by_sample.svg",width = 7.72,height = 6.08)
emf("2-6.umap_tsne/2.DimPlot_umap_by_sample.emf", width = 7.72,height = 6.08)
print(UMAP_p2)
dev.off()


## 按组织类型
UMAP_p3 <- DimPlot(alldata, group.by = "site", reduction = "umap",label = TRUE, 
                   repel = TRUE)+
  scale_color_manual(values = sample_cols)+
  labs(title=NULL)
UMAP_p3
export::graph2pdf(UMAP_p3,"2-6.umap_tsne/3.DimPlot_umap_by_site.pdf",width = 7.72,height = 6.08)
export::graph2svg(UMAP_p3,"2-6.umap_tsne/3.DimPlot_umap_by_site.svg",width = 7.72,height = 6.08)
emf("2-6.umap_tsne/3.DimPlot_umap_by_site.emf", width = 7.72,height = 6.08)
print(UMAP_p3)
dev.off()

## 按细胞类型
UMAP_p4 <- DimPlot(alldata, group.by = "celltype", reduction = "umap")+
  scale_color_manual(values = sample_cols)+
  labs(title=NULL)
UMAP_p4
export::graph2pdf(UMAP_p4,"2-6.umap_tsne/4.DimPlot_umap_by_cellType.pdf",width = 7.72,height = 6.08)
export::graph2svg(UMAP_p4,"2-6.umap_tsne/4.DimPlot_umap_by_cellType.svg",width = 7.72,height = 6.08)
emf("2-6.umap_tsne/4.DimPlot_umap_by_cellType.emf", width = 7.72,height = 6.08)
print(UMAP_p4)
dev.off()

Idents(alldata) <- 'celltype'
library(ggsci)
library("scales")
colors <- c(pal_ucscgb()(26))
show_col(colors)

library(ggpubr)
x <- table(alldata@meta.data$celltype,alldata@meta.data$sample);x
x <- t(t(x)/rowSums(t(x))) %>% as.table() %>% as.data.frame();x
colnames(x) <- c('celltype','sample','Freq')
x$celltype <- factor(x$celltype,
                     levels = unique(alldata@meta.data$celltype)[
                       order(unique(alldata@meta.data$celltype))])
x$group<-NA
for (i in 1:nrow(x)) {
  if(stringr::str_detect(x$sample[i],'T')){
    x$group[i] <- 'Tumor'
  }else if(stringr::str_detect(x$sample[i],'N')){
    x$group[i] <- 'Normal'
  }else if(stringr::str_detect(x$sample[i],'L')){
    x$group[i] <- 'Lymph'
  }else if(stringr::str_detect(x$sample[i],'P')){
    x$group[i] <- 'PVTT'
  }
}
x$group <- factor(x$group,levels = c('Normal','Tumor','Lymph','PVTT'))
x<-arrange(x,celltype,group)#重新排序
#绘图风格
theme.set = theme(
  axis.title = element_text(size = 20, face = "bold", family = "Arial"),
  axis.text.x = element_text(size = 14, angle = 45, hjust = 1,face = "bold", family = "Arial"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Arial"),
  #legend.position = 'top',
  legend.text = element_text(size = 16, face = "bold", family = "Arial"),
  legend.title = element_text(size = 18,face='bold',family = "Arial"),
  text = element_text(family = "Arial"))

compare_cell_anno <- 
  ggplot(data = x, aes(x = group, y = Freq, fill = group)) +
  geom_boxplot(varwidth = 0.3) + 
  scale_fill_aaas()+
  ggpubr::stat_compare_means(comparisons = 
                        list(c('Normal','Tumor'),
                             c('Normal','Lymph'),
                             c('Normal','PVTT'),
                             c('Tumor','Lymph'),
                             c('Tumor','PVTT'),
                             c('Lymph','PVTT')),
                      test = "wilcox.test",
                      hide.ns = T,
                      label = 'p.signif')+
  labs(x='',y='Proportion')+
  facet_wrap(~celltype,nrow=1)+
  theme_bw() + theme.set
compare_cell_anno
graph2pdf(compare_cell_anno,'2-9.compare_cell_anno.pdf', w=13,h=8)

# 堆叠柱状图
library('RColorBrewer')
cells_barplot <- 
  ggplot(data = x, aes(x = sample, y = Freq, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill",width = 0.7)+
  scale_fill_manual(values = brewer.pal(10,"BrBG"))+
  scale_y_continuous(expand = c(0,0))+
  theme_minimal() + 
  labs(x='Sample',y='Freq',fill='Cell Type')+
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"),  # 上、右、下、左
        axis.title = element_text(face = 'bold'),
        axis.text.x = element_text(angle=90,vjust=0.5),
        legend.title = element_text(face = 'bold',hjust=0.5))
graphics.off()
cells_barplot
graph2pdf(cells_barplot,'2-10.cells_barplot.pdf', w=12,h=5.86)

saveRDS(alldata,'./scRNA_anno.RDS')

#Endo
{
  if(!dir.exists('Endo')){dir.create('Endo')};setwd('Endo')
  # 读取带标签的原始数据
  alldata <- readRDS('../../Data/scRNA.RDS')
  # 提取正常和肿瘤组织中的内皮细胞
  Idents(alldata) = "celltype" # 设置 细胞类型 为分群信息
  Endothelial_cells <- alldata[, Idents(alldata) == 'Endothelial']#提取内皮细胞亚群
  Idents(Endothelial_cells) = "site" # 设置 组织类型 为分群信息
  Endothelial_cells <- Endothelial_cells[,Idents(Endothelial_cells) %in% 
                                           c('Normal','Tumor')]#提取内皮细胞亚群
  
  # 标准化
  Endothelial_cells <- Endothelial_cells %>% 
    NormalizeData(normalization.method = "LogNormalize",scale.factor = 10000) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, 
                         verbose = FALSE, assay = "RNA") %>% 
    ScaleData()
  # Z分数转换缩放基因
  Endothelial_cells <- ScaleData(Endothelial_cells, 
                                 vars.to.regress = c("percent_mito", "nFeature_RNA"), 
                                 assay = "RNA")
  #PCA
  Endothelial_cells <- RunPCA(Endothelial_cells, npcs = 100, verbose = F)
  
  #碎石图
  ElbowPlot(Endothelial_cells,reduction="pca",ndims = 100) +theme.set #根据ElbowPlot确定PC数量
  ElbowPlot <- recordPlot()
  export::graph2pdf(ElbowPlot,'1.ElbowPlot.pdf',w=8,h=6)
  
  # 聚类
  pcSelect <- 32
  Endothelial_cells <- Endothelial_cells %>% FindNeighbors(dims = 1:pcSelect) 
  # 选择合适的分辨率
  # 从0.1-4的resolution结果均运行一遍
  seq = seq(0.1,4,by=0.1)
  for (res in seq) {
    Endothelial_cells = FindClusters(Endothelial_cells, resolution = res)
  }
  # 画图
  clustree <- clustree(Endothelial_cells,prefix = 'RNA_snn_res.') + 
    coord_flip()+
    guides(color = guide_legend(ncol = 2))  # 将图例分为两列
  clustree
  ggsave('2.clustree.pdf',clustree,w=30,h=14)
  
  # 根据clustree图选择res，在此res上增加并不会明显增加分群数目
  Endothelial_cells <- FindClusters(object = Endothelial_cells, resolution = 0.5)   
  
  # 正常和肿瘤组织中的内皮细胞亚群UMAP TSNE ----
  Endothelial_cells <- RunUMAP(Endothelial_cells, dims = 1:pcSelect)
  Endothelial_cells <- RunTSNE(Endothelial_cells, dims = 1:pcSelect) 
  
  #手动注释
  #来源：CellMarker2.0
  genes_markers<-openxlsx::read.xlsx('./EndoMarkers.xlsx')
  # 对 cell 名自动换行，width 可根据你的最长名字适当调整
  genes_markers$cell_wrapped <- str_wrap(genes_markers$cell, width = 10)
  
  Idents(Endothelial_cells) <- 'seurat_clusters'
  DotPlot(Endothelial_cells, 
          features = split(genes_markers$marker_gene,
                           genes_markers$cell_wrapped), 
          assay = "RNA") + 
    scale_color_gradientn(colours = c('#330066','#336699','white','#FFCC33','red')) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
          # 可选：让 facet strip 的文字稍小一点，防止换行后仍然拥挤
          strip.text.x = element_text(size = 12, margin = ggplot2::margin(2,0,2,0)))
  ggsave('3.cells_anno_clusters_dotplot.pdf', w = 9, h = 6)
  ggsave('3.cells_anno_clusters_dotplot.svg', w = 9, h = 6)
  
  #anno结果
  Endothelial_cells$celltype <-
    case_when(Endothelial_cells$seurat_clusters %in% c(3,5,10,14,17)~'Precursor',#
              Endothelial_cells$seurat_clusters %in% c(2,4,7,8,9,12,13,16)~'LSECs',#
              Endothelial_cells$seurat_clusters %in% c(1,6,11)~'Vascular',#
              TRUE~'Other')#
  
  table(Endothelial_cells$celltype)
  # LSECs     Other Precursor  Vascular 
  # 1419       520       839       759 
  table(Endothelial_cells$seurat_clusters) %>% sum() == table(Endothelial_cells$celltype) %>% sum()
  # [1] TRUE
  
  num_dir <- '4.umap_tsne_Endothelial'; if (!dir.exists(num_dir)) {dir.create(num_dir)}
  # 内皮细胞亚群UMAP
  ## 按簇
  UMAP_p <- DimPlot(Endothelial_cells, reduction = "umap")+
    scale_color_manual(values = sample_cols)
  UMAP_p
  export::graph2pdf(UMAP_p,"4.umap_tsne_Endothelial/1.DimPlot_umap_by_cluster.pdf",width = 7.72,height = 6.08)
  emf("4.umap_tsne_Endothelial/1.DimPlot_umap_by_cluster.emf", width = 7.72,height = 6.08)
  print(UMAP_p)
  dev.off()
  
  ## 按样本
  UMAP_p2 <- DimPlot(Endothelial_cells, group.by = "sample", 
                     reduction = "umap")+
    scale_color_manual(values = sample_cols)+
    labs(title=NULL)
  UMAP_p2
  export::graph2pdf(UMAP_p2,"4.umap_tsne_Endothelial/2.DimPlot_umap_by_sample.pdf",width = 7.72,height = 6.08)
  emf("4.umap_tsne_Endothelial/2.DimPlot_umap_by_sample.emf", width = 7.72,height = 6.08)
  print(UMAP_p2)
  dev.off()
  
  ## 按组织类型
  UMAP_p3 <- DimPlot(Endothelial_cells, group.by = "site", reduction = "umap")+
    scale_color_manual(values = sample_cols)+
    labs(title=NULL)
  UMAP_p3
  export::graph2pdf(UMAP_p3,"4.umap_tsne_Endothelial/3.DimPlot_umap_by_site.pdf",width = 7.72,height = 6.08)
  emf("4.umap_tsne_Endothelial/3.DimPlot_umap_by_site.emf", width = 7.72,height = 6.08)
  print(UMAP_p3)
  dev.off()
  
  ## 按亚型类型
  UMAP_p4 <- DimPlot(Endothelial_cells, group.by = "celltype", reduction = "umap")+
    scale_color_manual(values = sample_cols)+
    labs(title=NULL)
  UMAP_p4
  export::graph2pdf(UMAP_p4,"4.umap_tsne_Endothelial/4.DimPlot_umap_by_celltype.pdf",width = 7.72,height = 6.08)
  emf("4.umap_tsne_Endothelial/4.DimPlot_umap_by_celltype.emf", width = 7.72,height = 6.08)
  print(UMAP_p4)
  dev.off()
  
  
  #标志基因热图
  {
    library(Seurat)
    library(pheatmap)
    
    # 1. 提取感兴趣的标记基因列表
    genes_to_plot <- unique(genes_markers$marker_gene)
    
    # 2. 确保这些基因存在于 Seurat 对象中
    genes_to_plot <- intersect(genes_to_plot, rownames(Endothelial_cells))
    if(length(genes_to_plot) == 0) stop("None of the marker genes found in the Seurat object.")
    
    # 3. 计算每个细胞类型的平均表达量（默认使用 data 插槽的归一化值）
    avg_expr <- AverageExpression(
      Endothelial_cells,
      features = genes_to_plot,
      group.by = "celltype",
      slot = "data"    # 使用 log-normalized 数据；若想用 counts 可改为 "counts"
    )
    avg_matrix <- avg_expr$RNA   # 如果 assay 是 RNA，否则换成对应 assay 名称
    # 若 Seurat 版本不同，可能需要用 avg_expr[[1]]
    
    # 4. 按行（基因）进行 Z-score 标准化，使热图更清晰展示相对表达模式
    avg_matrix_scaled <- scale(t(avg_matrix))
    genes_ordered<- c("FCGR2B","GPR182","LYVE1","ITGA4","FLT1","ITGA2","ITGA3")
    cells_ordered<-c('LSECs','Precursor','Vascular','Other')
    avg_matrix_scaled <- avg_matrix_scaled[cells_ordered, genes_ordered, drop = FALSE]
    
    # 5. 绘制热图
    pheatmap(
      avg_matrix_scaled,
      cluster_rows = F,      # 对基因聚类
      cluster_cols =F,      # 对细胞类型聚类
      scale = "none",           # 已手动标准化
      row_names_side = "left",   # 行名放在右侧
      color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
      main = "Marker gene expression per cell type",
      angle_col = '45',
      fontsize_row = 10,
      fontsize_col = 10
    )
    marker_plot<-recordPlot()
    export::graph2pdf(marker_plot,'5.marker_plot.pdf',height=3,width=8)
    export::graph2svg(marker_plot,'5.marker_plot.svg',height=3,width=8)
    export::graph2png(marker_plot,'5.marker_plot.png',height=3,width=8,dpi=600)
  }
  
  saveRDS(Endothelial_cells,'scRNA_Endothelial.RDS')
  
  # 提取正常和肿瘤组织中的内皮细胞
  Idents(Endothelial_cells) = "celltype" # 设置 细胞类型 为分群信息
  LSECs_cells <- Endothelial_cells[, Idents(Endothelial_cells) == 'LSECs']#提取内皮细胞亚群
  # 标准化
  LSECs_cells <- LSECs_cells %>% 
    NormalizeData(normalization.method = "LogNormalize",scale.factor = 10000) %>% 
    FindVariableFeatures(selection.method = "vst", nfeatures = 2000, 
                         verbose = FALSE, assay = "RNA") %>% 
    ScaleData()
  # Z分数转换缩放基因
  LSECs_cells <- ScaleData(LSECs_cells, 
                           vars.to.regress = c("percent_mito", "nFeature_RNA"), 
                           assay = "RNA")
  saveRDS(LSECs_cells,'scRNA_LSECs_nromalized.RDS')
  setwd('../')
}

#--------------------------关键细胞鉴定-----------------------------------------
setwd('/home/liuzhaoshuo/UoB_Work/HCC_scRNA/')
if (!dir.exists("3_scRNA_KeyCell")) {dir.create("3_scRNA_KeyCell")}
setwd("3_scRNA_KeyCell")

rm(list = ls());gc()
scedata<-readRDS('../2_scRNA_Dimension-Cluster_Anno/Endo/scRNA_LSECs_nromalized.RDS')
load('../env.RData')

Idents(scedata) <- 'site'
dim(scedata)



# 预后基因在关键细胞中的表达
#UMAP Dot
{
  EHMT2_umap<-FeaturePlot(scedata, features = 'EHMT2',by.col = F, 
                          split.by = "site", 
                          cols = c("lightgrey", "blue"), 
                          pt.size = 1)
  export::graph2pdf(EHMT2_umap,'1.EHMT2_umap.pdf',w=4,h=6)
  export::graph2svg(EHMT2_umap,'1.EHMT2_umap.svg',w=4,h=6)
  
  #读取候选基因
  hubgenes <- c('PAX6', 'ANXA10', 'DNAH7', 'OTOGL', 'ZP3', 'ATF7-NPFF', 
                'FAT2', 'SH3PXD2A', 'ADGRB2', 'CFAP161', 'ERICH6B', 'LINGO1', 
                'KSR2', 'GARIN1A', 'MUC19', 'EML6')
  
  hubgene_umap_grid<-list()
  for(i in 1:length(hubgenes)){
    
    # 检查基因是否存在于默认 assay 的 scale.data/data 中
    if (! hubgenes[i] %in% rownames(scedata)) {
      message("Skipping missing gene: ", hubgenes[i])
      next  # 跳过本次循环
    }
    
    hubgene_umap_grid[[i]] <- FeaturePlot(scedata, features = hubgenes[i],by.col = F, 
                                          split.by = "site", 
                                          cols = c("lightgrey", "blue"), 
                                          pt.size = 1)
  }
  # Skipping missing gene: ATF7-NPFF
  # Skipping missing gene: FAT2
  # Skipping missing gene: ERICH6B
  # Skipping missing gene: GARIN1A
  # Skipping missing gene: MUC19
  
  # 去除列表中的 NULL 元素
  hubgene_umap_grid <- hubgene_umap_grid[!sapply(hubgene_umap_grid, is.null)]
  hubgene_umap_grid <- patchwork::wrap_plots(hubgene_umap_grid)
  hubgene_umap_grid
  export::graph2pdf(hubgene_umap_grid,'2.hubgene_umap.pdf',w=15,h=18)
  export::graph2svg(hubgene_umap_grid,'2.hubgene_umap.svg',w=15,h=18)
  
  EHMT2_dotplot <- DotPlot(scedata, features = 'EHMT2',
                           assay = "RNA",
                           group.by = "stage", 
                           split.by = "site",
                           cols = c("lightgrey", "blue")) + 
    RotatedAxis()
  EHMT2_dotplot
  export::graph2pdf(EHMT2_dotplot,'3.EHMT2_dotplot_site-stage.pdf',w=4,h=4.5)
  export::graph2svg(EHMT2_dotplot,'3.EHMT2_dotplot_site-stage.svg',w=4,h=4.5)
  
  hubgene_dotplot <- DotPlot(scedata, features = hubgenes,
                             assay = "RNA",
                             group.by = "stage", 
                             split.by = "site",
                             cols = c("lightgrey", "blue")) + 
    RotatedAxis()
  hubgene_dotplot
  export::graph2pdf(hubgene_dotplot,'4.hubgene_dotplot_site-stage.pdf',w=6,h=4.5)
  export::graph2svg(hubgene_dotplot,'4.hubgene_dotplot_site-stage.svg',w=6,h=4.5)
}


# 绘制小提琴图观察不同celltype中不同时期的EHMT2表达差异，Wilcoxon检验
{
  # 一次性提取基因表达与分组信息
  plot_data <- FetchData(scedata, vars = c("EHMT2", "site", "stage"))
  
  # 此时 plot_data 中必然有 "EHMT2" 列，且行名是细胞条码
  head(plot_data)
  
  library(ggpubr)
  library(ggsci)
  
  # 小提琴图：不同 stage
  function1 <- function(gene){
    p <- ggplot(plot_data, aes(x = stage, y = .data[[gene]], fill = stage)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
      scale_fill_npg() +
      stat_compare_means(
        mapping = aes(group = stage),
        comparisons = list(c("I","II"), c("I","IIIA"), c("I","IIIB"),
                           c("I","IV"),  c("IIIA","II"), c("IIIB","II"),
                           c("II","IV"), c("IIIB","IIIA"), c("IIIA","IV"),
                           c("IIIB","IV")),
        label = "p.signif", method = "wilcox.test",
        paired = FALSE, label.x = 1.4, size = 3
      ) +
      labs(title = gene, y = "Expression Level") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank()
      )
    return(p)
  }
  
  p1 <- function1("EHMT2")
  p1
  export::graph2pdf(p1, '5.EHMT2_expr_sc_stage.pdf', w = 4, h = 5)
  export::graph2svg(p1, '5.EHMT2_expr_sc_stage.svg', w = 4, h = 5)
  
  # 小提琴图：Tumor vs Normal
  function2 <- function(gene){
    p <- ggplot(plot_data, aes(x = site, y = .data[[gene]], fill = site)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
      scale_fill_npg() +
      stat_compare_means(
        mapping = aes(group = site),
        comparisons = list(c("Tumor", "Normal")),
        label = "p.signif", method = "wilcox.test",
        paired = FALSE, label.x = 1.4, size = 3
      ) +
      labs(title = gene, y = "Expression Level") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank()
      )
    return(p)
  }
  
  p2 <- function2("EHMT2")
  p2
  export::graph2pdf(p2, '6.EHMT2_expr_sc_tissue.pdf', w = 4, h = 5)
  export::graph2svg(p2, '6.EHMT2_expr_sc_tissue.svg', w = 4, h = 5)
}

# 绘制小提琴图观察不同celltype中不同时期的hubgene表达差异，Wilcoxon检验
{
  # 一次性提取基因表达与分组信息
  plot_data <- FetchData(scedata, vars = c(hubgenes, "site", "stage"))
  
  # 此时 plot_data 中必然有候选基因列，且行名是细胞条码
  head(plot_data)
  
  library(ggpubr)
  library(ggsci)
  
  # 小提琴图：不同 stage
  function1 <- function(gene){
    p <- ggplot(plot_data, aes(x = stage, y = .data[[gene]], fill = stage)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
      scale_fill_npg() +
      stat_compare_means(
        mapping = aes(group = stage),
        comparisons = list(c("I","II"), c("I","IIIA"), c("I","IIIB"),
                           c("I","IV"),  c("IIIA","II"), c("IIIB","II"),
                           c("II","IV"), c("IIIB","IIIA"), c("IIIA","IV"),
                           c("IIIB","IV")),
        label = "p.signif", method = "wilcox.test",
        paired = FALSE, label.x = 1.4, size = 3
      ) +
      labs(title = gene, y = "Expression Level") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank()
      )
    return(p)
  }
  
  p1_list<-list()
  for(i in 1:length(hubgenes)){
    # 检查基因是否存在于默认 assay 的 scale.data/data 中
    if (! hubgenes[i] %in% rownames(scedata)) {
      message("Skipping missing gene: ", hubgenes[i])
      next  # 跳过本次循环
    }
    
    p1_list[[i]] <- function1(hubgenes[i])
  }
  # Skipping missing gene: ATF7-NPFF
  # Skipping missing gene: FAT2
  # Skipping missing gene: ERICH6B
  # Skipping missing gene: GARIN1A
  # Skipping missing gene: MUC19
  
  # 去除列表中的 NULL 元素
  p1_list <- p1_list[!sapply(p1_list, is.null)]
  p1 <- patchwork::wrap_plots(p1_list,nrow=3)
  p1
  export::graph2pdf(p1, '7.hubgenes_expr_sc_stage.pdf', w = 12, h = 10)
  export::graph2svg(p1, '7.hubgenes_expr_sc_stage.svg', w = 12, h = 10)
  
  # 小提琴图：Tumor vs Normal
  function2 <- function(gene){
    p <- ggplot(plot_data, aes(x = site, y = .data[[gene]], fill = site)) +
      geom_violin(scale = "width", trim = TRUE) +
      geom_boxplot(width = 0.1, position = position_dodge(0.9), outlier.shape = NA) +
      scale_fill_npg() +
      stat_compare_means(
        mapping = aes(group = site),
        comparisons = list(c("Tumor", "Normal")),
        label = "p.signif", method = "wilcox.test",
        paired = FALSE, label.x = 1.4, size = 3
      ) +
      labs(title = gene, y = "Expression Level") +
      theme_bw() +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"),
        axis.title.x = element_blank(),
        legend.position = "top", legend.title = element_blank(),
        panel.grid = element_blank()
      )
    return(p)
  }
  
  p2_list<-list()
  for(i in 1:length(hubgenes)){
    # 检查基因是否存在于默认 assay 的 scale.data/data 中
    if (! hubgenes[i] %in% rownames(scedata)) {
      message("Skipping missing gene: ", hubgenes[i])
      next  # 跳过本次循环
    }
    
    p2_list[[i]] <- function2(hubgenes[i])
  }
  # Skipping missing gene: ATF7-NPFF
  # Skipping missing gene: FAT2
  # Skipping missing gene: ERICH6B
  # Skipping missing gene: GARIN1A
  # Skipping missing gene: MUC19
  
  # 去除列表中的 NULL 元素
  p2_list <- p2_list[!sapply(p2_list, is.null)]
  p2 <- patchwork::wrap_plots(p2_list,nrow=3)
  p2
  export::graph2pdf(p2, '8.hubgenes_expr_sc_tissue.pdf', w = 12, h = 10)
  export::graph2svg(p2, '8.hubgenes_expr_sc_tissue.svg', w = 12, h = 10)
}

# save.image('Image.RData')

#--------------------------差异基因分析-----------------------------------------
setwd('/home/liuzhaoshuo/UoB_Work/HCC_scRNA/')
if (!dir.exists("4_scRNA_DEGs")) {dir.create("4_scRNA_DEGs")}
setwd("4_scRNA_DEGs")

rm(list = ls());gc()
scedata<-readRDS('../2_scRNA_Dimension-Cluster_Anno/Endo/scRNA_LSECs_nromalized.RDS')
load('../env.RData')

# 设置分组标识为 site 列
Idents(scedata) <- "site"

# Tumor vs Normal 差异表达基因（Tumor 相对 Normal 高/低表达）
degs <- FindMarkers(
  object = scedata,
  ident.1 = "Tumor",      # 实验组
  ident.2 = "Normal",     # 对照组
  only.pos = FALSE,       # FALSE 表示同时保留上调和下调基因
  logfc.threshold = 0.25, # log2 fold change 最小阈值
  min.pct = 0.1,          # 至少在 ident.1 或 ident.2 中表达细胞的百分比
  test.use = "wilcox"     # 默认秩和检验，适合单细胞数据
)

# 导出结果
write.csv(degs, "1.Tumor_vs_Normal_DEGs.csv")

# 筛选显著差异基因（使用原始p值）
library(dplyr)

sig_degs <- degs %>%
  filter(p_val < 0.05, abs(avg_log2FC) > 0.5) %>%
  arrange(desc(avg_log2FC))

# 火山图（使用原始p值）
library(ggplot2)
library(ggrepel)

# 添加区分列（基于原始p值）
degs$significance <- "Not significant"
degs$significance[degs$p_val < 0.05 & degs$avg_log2FC > 0.5] <- "Up"
degs$significance[degs$p_val < 0.05 & degs$avg_log2FC < -0.5] <- "Down"
table(degs$significance)
# Down Not significant              Up 
# 243            2550            1059 

# 标记前几个差异最显著的基因（按原始p值排序）
top_genes <- rbind(
  degs %>% 
    filter(significance == "Up") %>%
    slice_min(order_by = p_val, n = 10),
  degs %>% 
    filter(significance == "Down") %>%
    slice_min(order_by = p_val, n = 10)
)

ggplot(degs, aes(x = avg_log2FC, y = -log10(p_val), color = significance)) +
  geom_point(alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Down" = "blue", "Not significant" = "grey", "Up" = "red")) +
  geom_text_repel(data = top_genes, aes(label = rownames(top_genes)), size = 3) +
  geom_hline(yintercept = (-log10(0.05)), linetype = 'dashed')+
  geom_vline(xintercept = c((-0.25), 0.25), linetype = 'dashed')+
  theme_bw() +
  labs(title = "Tumor vs Normal", x = "log2 Fold Change", y = "-log10 P-value")+
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5),
    axis.title = element_text(face = 'bold', hjust = 0.5),
    legend.position = 'none'
  )
ggsave('2.Volcano.pdf', width = 6, height = 6)
ggsave('2.Volcano.svg', width = 6, height = 6)


# 定义需要高亮的基因集（hubgenes）
hubgenes <- c('PAX6', 'ANXA10', 'DNAH7', 'OTOGL', 'ZP3', 'ATF7-NPFF', 
              'FAT2', 'SH3PXD2A', 'ADGRB2', 'CFAP161', 'ERICH6B', 'LINGO1', 
              'KSR2', 'GARIN1A', 'MUC19', 'EML6')

# 将基因名作为一列
degs$gene <- rownames(degs)
highlight_genes <- intersect(hubgenes, degs$gene)

# 提取这些基因的子集数据框
degs_highlight <- degs %>% dplyr::filter(gene %in% highlight_genes)

# 高亮火山图（使用原始p值）
ggplot(degs_highlight, aes(x = avg_log2FC, y = -log10(p_val))) +
  geom_point(aes(color = significance), alpha = 0.5, size = 1) +
  scale_color_manual(values = c("Down" = "blue", 
                                "Not significant" = "grey", 
                                "Up" = "red")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = degs_highlight, 
                  aes(label = gene), 
                  size = 3.5, 
                  max.overlaps = 20,
                  fontface = "italic") +
  theme_bw() +
  labs(title = "Tumor vs Normal (hub genes)",
       x = "log2 Fold Change",
       y = "-log10 P-value") +
  guides(color = guide_legend(override.aes = list(size = 3)))+
  theme(
    plot.title = element_text(face = 'bold', hjust = 0.5),
    axis.title = element_text(face = 'bold', hjust = 0.5),
    legend.position = 'none'
  )+
  scale_x_continuous(limits = c((-1),1))+
  scale_y_continuous(limits = c((-0.1),(-log10(0.05))*2))
ggsave('3.Volcano_hubgenes.pdf', width = 6, height = 6)
ggsave('3.Volcano_hubgenes.svg', width = 6, height = 6)


