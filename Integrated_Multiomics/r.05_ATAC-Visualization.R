#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/Integrated_Multiomics/')
dir_name<-'05_ATAC-Visualization';if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

source("https://raw.githubusercontent.com/XLions/MyOwnTools/c4d651dcc5fd92c3106a4bd5888b7985c7da25e1/R/ATACCombined.R")

#主体函数是我自己写的自定义函数
targetGenes<-c("PAX6", "ANXA10", "DNAH7", "OTOGL", "ZP3",
               "ATF7-NPFF", "FAT2", "SH3PXD2A", "ADGRB2", 
               "CFAP161", "ERICH6B", "LINGO1", "KSR2", "GARIN1A", 
               "MUC19", "EML6")

Genome<-rtracklayer::import('D:/OneDrive/PROJECTS/UoB/M6_GroupProject_EHMT2/Homo_sapiens.GRCh38.115.chr.gtf.gz')
Genome<-as.data.frame(Genome)

if(!dir.exists('Double')){dir.create('Double')};setwd('Double')
for(i1 in 1:length(targetGenes)){
  
  GeneID<-targetGenes[i1]
  #B2vsB1 B1vsB2
  {
    macs2_dir_case<-'D:/OneDrive/PROJECTS/UoB/M6_GroupProject_EHMT2/ATAC_Seq/CallPeak_MACS2/B2_vs_B1'
    macs2_dir_ctl<-'D:/OneDrive/PROJECTS/UoB/M6_GroupProject_EHMT2/ATAC_Seq/CallPeak_MACS2/B1_vs_B2'
    regionData1<-getLocationByID(GeneID=GeneID,
                                 ID_Type = 'gene_name',
                                 Genome,
                                 macs2_dir_case,
                                 sig_level=0.05,
                                 FC_cutoff=2)
    regionData2<-getLocationByID(GeneID=GeneID,
                                 ID_Type = 'gene_name',
                                 Genome,
                                 macs2_dir_ctl,
                                 sig_level=0.05,
                                 FC_cutoff=2)
    regionData<-list(
      gene_region=rbind(regionData1$gene_region,regionData2$gene_region) %>% dplyr::distinct(),
      macs2_peaks=rbind(regionData1$macs2_peaks,regionData2$macs2_peaks) %>% dplyr::distinct()
    )
    p1<-doubleBedGraphCtlCasePlot_Break(
      bedGraphData_case=read.table(paste0(macs2_dir_case,'/B2_vs_B1_summits.bed')), #实验组数据,read.table读取（X轴上方）
      bedGraphData_ctl=read.table(paste0(macs2_dir_ctl,'/B1_vs_B2_summits.bed')), #对照组数据,read.table读取（X轴下方）
      chrom_A=unique(regionData$gene_region$seqnames), #染色体序号
      plot_location=paste0(regionData$gene_region$start,'-',regionData$gene_region$end), #绘图X轴区域
      highlight_location_case=paste0(regionData1$macs2_peaks$start,'-',regionData1$macs2_peaks$end), #高亮区域
      highlight_location_ctl=paste0(regionData2$macs2_peaks$start,'-',regionData2$macs2_peaks$end), #高亮区域
      sample_name_ctl='B1 vs B2',
      sample_name_case='B2 vs B1', #样本名称
      DataRange=c(0,5), #绘图Y轴范围
      x_axis_label=T, #默认要不要X轴的标度标签
      expand_times=1.5
    )
  }
  
  #C2vsC1 C1vsC2
  {
    macs2_dir_case<-'D:/OneDrive/PROJECTS/UoB/M6_GroupProject_EHMT2/ATAC_Seq/CallPeak_MACS2/C2_vs_C1'
    macs2_dir_ctl<-'D:/OneDrive/PROJECTS/UoB/M6_GroupProject_EHMT2/ATAC_Seq/CallPeak_MACS2/C1_vs_C2'
    regionData1<-getLocationByID(GeneID=GeneID,
                                 ID_Type = 'gene_name',
                                 Genome,
                                 macs2_dir_case,
                                 sig_level=0.05,
                                 FC_cutoff=2)
    regionData2<-getLocationByID(GeneID=GeneID,
                                 ID_Type = 'gene_name',
                                 Genome,
                                 macs2_dir_ctl,
                                 sig_level=0.05,
                                 FC_cutoff=2)
    regionData<-list(
      gene_region=rbind(regionData1$gene_region,regionData2$gene_region) %>% dplyr::distinct(),
      macs2_peaks=rbind(regionData1$macs2_peaks,regionData2$macs2_peaks) %>% dplyr::distinct()
    )
    p2<-doubleBedGraphCtlCasePlot(
      bedGraphData_case=read.table(paste0(macs2_dir_case,'/C2_vs_C1_summits.bed')), #实验组数据,read.table读取（X轴上方）
      bedGraphData_ctl=read.table(paste0(macs2_dir_ctl,'/C1_vs_C2_summits.bed')), #对照组数据,read.table读取（X轴下方）
      chrom_A=unique(regionData$gene_region$seqnames), #染色体序号
      plot_location=paste0(regionData$gene_region$start,'-',regionData$gene_region$end), #绘图X轴区域
      highlight_location=paste0(regionData$macs2_peaks$start,'-',regionData$macs2_peaks$end), #高亮区域
      sample_name_ctl='C1 vs C2',
      sample_name_case='C2 vs C1', #样本名称
      DataRange=c(0,5), #绘图Y轴范围
      x_axis_label=T #默认要不要X轴的标度标签
    )
  }
  p_double<-
    patchwork::wrap_plots(
      list(p1,p2),ncol=1
    )+patchwork::plot_annotation(
      title=GeneID,
      theme=theme(
        plot.title = element_text(face = 'bold', hjust =0.5)
      )
    )
  ggsave(paste0(GeneID,'_p_double.svg'),height = 6,width = 10)
  ggsave(paste0(GeneID,'_p_double.png'),height = 6,width = 10,dpi=600)
  # setwd('../')
}





