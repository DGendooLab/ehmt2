
#设置工作目录，清除内存缓存
rm(list=ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DAAs/')

#bluebear的X11选项加载
options(bitmapType='cairo')

#加载R包
library(openxlsx)
library(tidyverse)
library(readr)
library(progress)
library(ggVennDiagram)
library(ggplot2)
library(ggforce)
library(patchwork)
library(ggbio)
library(GenomicRanges)
library(Gviz)

#1.准备工作
if(!file.exists('PreLoadEnv.RData')){
  #1.1.样本信息
  callpeak_info<-data.frame(
    sample=c('B2_vs_B1','B3_vs_B1','B4_vs_B2','B4_vs_B3',
             'C2_vs_C1','C3_vs_C1','C4_vs_C2','C4_vs_C3',
             'B1_vs_B2','B1_vs_B3','B2_vs_B4','B3_vs_B4',
             'C1_vs_C2','C1_vs_C3','C2_vs_C4','C3_vs_C4'),
    group=c(rep(c('KD-C','CT-C','KDT-KD','KDT-CT'),2),
            rep(c('C-KD','C-CT','KD-KDT','CT-KDT'),2)),
    direction=rep(c('Open','Close'),each=8)
  )
  #1.2.注释所需数据
  {
    #读取来源于GENCODE的基因组注释数据
    #Read genome annotation data from GENCODE
    ann_data_gff3<-read.table('../../References_Index/references/gencode.v49.chr_patch_hapl_scaff.annotation.gff3.gz')
    #按照GFF3标准格式重命名列
    #Rename the column according to GFF3 standard format
    colnames(ann_data_gff3)<-c(
      'seqid','source','type','start','end','score','strand','phase','attributes'
    )
  }
  ##1.3.定义函数：拆分attributes
  ##1.3.Define function: Split attributes
  ##输出：包含向量Item和Value及对应数量的数据框
  #Output: A data box containing vectors Item, Value, and their corresponding quantities
  separateAttributesCharToListIncludingItemValue<-function(attribute_char){
    singleChars<-stringr::str_split_1(attribute_char,';')
    assign('Items',c())
    assign('Values',c())
    for(i in 1:length(singleChars)){
      Items<-c(Items,
               stringr::str_split_1(singleChars[i],'=')[1])
      Values<-c(Values,
                stringr::str_split_1(singleChars[i],'=')[2])
    }
    return(data.frame(Item=Items,Value=Values,Length=length(singleChars)))
  }
  #1.4.保存预加载数据
  save.image('PreLoadEnv.RData')
}else{
  #读取预加载数据
  load('PreLoadEnv.RData')
}

#2.通过结果excel表格筛选
#定义筛选函数
getSigPeak<-function(dataadd,pcutoff,foldcutoff){
  data<-read.table(dataadd,header = T)
  library(dplyr)
  filteredData<-
    data %>%
    filter(`X.log10.pvalue.`>(-log10(pcutoff))) %>%
    filter(`fold_enrichment`>(foldcutoff))
  return(filteredData)
}
#建立新变量
SigPeaks<-data.frame(matrix(ncol=12,nrow=0))
colnames(SigPeaks)<-  
  c("chr","start","end","length","abs_summit","pileup","X.log10.pvalue.",
    "fold_enrichment","X.log10.qvalue.","name","sample",'group')
#设置参数
pcutoff<-0.05
foldcutoff<-2
for (i in 1:nrow(callpeak_info)) {
  add<-list.files(
    path = paste0('../CallPeak_MACS2/',callpeak_info$sample[i]),
    pattern = '.xls',
    full.names = T
  )
  filtered_data<-getSigPeak(add,pcutoff,foldcutoff)
  filtered_data$sample<-callpeak_info$sample[i]
  filtered_data$group<-callpeak_info$group[which(callpeak_info$sample==callpeak_info$sample[i])]
  SigPeaks<-rbind(
    SigPeaks,
    filtered_data
  )
}
saveRDS(SigPeaks,paste0('./RawSigPeaks_p',pcutoff,'_fold',foldcutoff,'.RDS'))

#Method1：先注释在对基因取交集
{
  #设置工作目录，清除内存缓存
  rm(list=ls());gc()
  setwd('/rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DAAs/')
  if(!dir.exists('Method1')){dir.create('Method1')};setwd('Method1')
  
  load('../PreLoadEnv.RData')
  SigPeaks<-readRDS(list.files(path = '..', pattern = '.RDS', full.names = T))
  
  #注释【【【【十分费时间】】】】
  #Annotate (computationally intensive process)
  SigPeaks$GENE_ID<-NA
  pb <- progress_bar$new(total = nrow(SigPeaks))#创建一个进度条对象 #Create a progress bar object
  for(i in 1:nrow(SigPeaks)){
    selectedAnnoData<-
      ann_data_gff3 %>%
      dplyr::filter(seqid==SigPeaks$chr[i]) %>%
      dplyr::filter(start<=SigPeaks$start[i]) %>%
      dplyr::filter(end>=SigPeaks$end[i]) %>%
      dplyr::filter(type %in% c('gene'))
    if(nrow(selectedAnnoData)==0){
      selectedAttributes<-NA
    }else{
      selectedAttributes<-
        (separateAttributesCharToListIncludingItemValue(selectedAnnoData$attributes[1]) %>%
           dplyr::filter(Item=='gene_id'))$Value
    }
    SigPeaks$GENE_ID[i]<-selectedAttributes
    if((i %% 1000)==0){
      print(i)
    }
    pb$tick()  #更新进度条#Update progress bar
  }
  saveRDS(SigPeaks,'./SigPeaks.RDS')
  
  #取交集-KD后开放区域
  DAA_gene_list<-list()
  DAA_gene_list_plot<-list()
  for(i in 1:length(c('KD-C','CT-C','KDT-KD','KDT-CT'))){
    group_i<-c('KD-C','CT-C','KDT-KD','KDT-CT')[i]
    samples_i<-callpeak_info$sample[which(callpeak_info$group==group_i)]
    DAA_gene_i1<-unique((SigPeaks %>%
                           dplyr::filter(group==group_i,
                                         sample==samples_i[1]))$GENE_ID)
    DAA_gene_i1<-DAA_gene_i1[!is.na(DAA_gene_i1)] %>%
      str_sub(1,15)
    DAA_gene_i2<-unique((SigPeaks %>%
                           dplyr::filter(group==group_i,
                                         sample==samples_i[2]))$GENE_ID)
    DAA_gene_i2<-DAA_gene_i2[!is.na(DAA_gene_i2)] %>%
      str_sub(1,15)
    DAA_gene_list[[i]]<-intersect(DAA_gene_i1,DAA_gene_i2)
    
    x<-list(DAA_gene_i1,DAA_gene_i2)
    names(x)<-samples_i
    DAA_gene_list_plot[[i]]<-
      ggVennDiagram(x, , color = 1, lwd = 0.7,label_alpha = 0,
                    label_font = "Times") + 
      labs(title = group_i)+
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      theme(legend.position = "none",
            plot.title = element_text(face='bold',hjust=0.5))+
      coord_flip()
  }
  names(DAA_gene_list)<-c('KD-C','CT-C','KDT-KD','KDT-CT')
  names(DAA_gene_list_plot)<-c('KD-C','CT-C','KDT-KD','KDT-CT')
  
  saveRDS(DAA_gene_list,'1.DAA_Open_gene_list.RDS')
  saveRDS(DAA_gene_list_plot,'2.Venns_Open.RDS')
  ggsave('2.Venns_Open.pdf',wrap_plots(DAA_gene_list_plot),width = 7.5,height = 6,units = 'in')
  ggsave('2.Venns_Open.svg',wrap_plots(DAA_gene_list_plot),width = 7.5,height = 6,units = 'in')
  
  KD_C_info<-SigPeaks %>%
    dplyr::filter(group=='KD-C',
                  str_sub(GENE_ID,1,15) %in% DAA_gene_list$`KD-C`)
  head(
    KD_C_info %>% dplyr::select(c('chr','start','end','name','sample'))
  )
  
  write.csv(KD_C_info,'3.DAAgenes_KD-C.csv',row.names = F)
  
  #取交集-KD后关闭区域
  DAA_gene_list<-list()
  DAA_gene_list_plot<-list()
  for(i in 1:length(c('C-KD','C-CT','KD-KDT','CT-KDT'))){
    group_i<-c('C-KD','C-CT','KD-KDT','CT-KDT')[i]
    samples_i<-callpeak_info$sample[which(callpeak_info$group==group_i)]
    DAA_gene_i1<-unique((SigPeaks %>%
                           dplyr::filter(group==group_i,
                                         sample==samples_i[1]))$GENE_ID)
    DAA_gene_i1<-DAA_gene_i1[!is.na(DAA_gene_i1)] %>%
      str_sub(1,15)
    DAA_gene_i2<-unique((SigPeaks %>%
                           dplyr::filter(group==group_i,
                                         sample==samples_i[2]))$GENE_ID)
    DAA_gene_i2<-DAA_gene_i2[!is.na(DAA_gene_i2)] %>%
      str_sub(1,15)
    DAA_gene_list[[i]]<-intersect(DAA_gene_i1,DAA_gene_i2)
    
    x<-list(DAA_gene_i1,DAA_gene_i2)
    names(x)<-samples_i
    DAA_gene_list_plot[[i]]<-
      ggVennDiagram(x, , color = 1, lwd = 0.7,label_alpha = 0,
                    label_font = "Times") + 
      labs(title = group_i)+
      scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
      theme(legend.position = "none",
            plot.title = element_text(face='bold',hjust=0.5))+
      coord_flip()
  }
  names(DAA_gene_list)<-c('C-KD','C-CT','KD-KDT','CT-KDT')
  names(DAA_gene_list_plot)<-c('C-KD','C-CT','KD-KDT','CT-KDT')
  
  saveRDS(DAA_gene_list,'4.DAA_Close_gene_list.RDS')
  saveRDS(DAA_gene_list_plot,'5.Venns_Close.RDS')
  ggsave('5.Venns_Close.pdf',wrap_plots(DAA_gene_list_plot),width = 7.5,height = 6,units = 'in')
  ggsave('5.Venns_Close.svg',wrap_plots(DAA_gene_list_plot),width = 7.5,height = 6,units = 'in')
  
  C_KD_info<-SigPeaks %>%
    dplyr::filter(group=='C-KD',
                  str_sub(GENE_ID,1,15) %in% DAA_gene_list$`C-KD`)
  head(
    C_KD_info %>% dplyr::select(c('chr','start','end','name','sample'))
  )
  
  write.csv(C_KD_info,'6.DAAgenes_C-KD.csv',row.names = F)
  
  setwd('../')
}

#Method2：先取峰值区域交集再注释
{
  #设置工作目录，清除内存缓存
  rm(list=ls());gc()
  setwd('/rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DAAs/')
  if(!dir.exists('Method2')){dir.create('Method2')};setwd('Method2')
  
  load('../PreLoadEnv.RData')
  SigPeaks<-readRDS(list.files(path = '..', pattern = '.RDS', full.names = T))
  
  #将一堆数值转换为连续的区间开始结束位点的函数
  #输入：所有数值
  #输出：列表：开始位点列和结束位点列
  Number_PointsToSegment<-function(numbers){
    if(length(numbers)>1){
      numbers<-sort(numbers)
      starts<-c(numbers[1])
      ends<-c()
      for(i in 2:length(numbers)){
        if((numbers[i]-numbers[i-1])>1){
          ends<-c(ends,numbers[i-1])
          starts<-c(starts,numbers[i])
        }
      }
      ends<-c(ends,numbers[length(numbers)])
      if(length(ends)==length(starts)){
        result<-list(
          start=starts,
          end=ends,
          length=length(ends)
        )
      }else{
        result<-list(
          start=starts,
          end=ends,
          length='Ends and Starts have different lengths')
      }
      return(result) 
    }else if(length(numbers)==1){
      return(list(
        start=numbers,
        end=numbers,
        length=1))
    }
  }
  
  #KD后关闭区域
  #定义比较的两组
  Comparison<-'KD-C'
  #提取数据
  Peak_Comparison<-SigPeaks %>%
    dplyr::filter(group==Comparison)
  #获取BC两组比较的样本
  samples<-unique(Peak_Comparison$sample)
  sample_B<-samples[which(str_detect(samples,'vs_B'))]
  sample_C<-samples[which(str_detect(samples,'vs_C'))]
  
  #有两组峰的染色体
  Chrom_Peak<-unique(intersect(
    unique((SigPeaks %>%
              dplyr::filter(group==Comparison) %>%
              dplyr::filter(sample==sample_B))$chr),
    unique((SigPeaks %>%
              dplyr::filter(group==Comparison) %>%
              dplyr::filter(sample==sample_C))$chr)
  ))
  
  #建立结果变量
  CommonPeakRegion<-data.frame(matrix(ncol=3,nrow=0))
  colnames(CommonPeakRegion)<-c('chr','start','end')
  
  #根据染色体建立循环
  for(i in 1:length(Chrom_Peak)){
    #建立B组所有碱基位点的位置数的变量
    Peak_Comparison_B<-SigPeaks %>%
      dplyr::filter(group==Comparison) %>%
      dplyr::filter(sample==sample_B) %>%
      dplyr::filter(chr==Chrom_Peak[i])
    #建立B组所有位点位置
    assign('pos_B',c())
    for(i1 in 1:nrow(Peak_Comparison_B)){
      pos_B<-c(
        pos_B,
        seq(Peak_Comparison_B$start[i1],Peak_Comparison_B$end[i1])
      )
    }
    
    #建立C组所有碱基位点的位置数的变量
    Peak_Comparison_C<-SigPeaks %>%
      dplyr::filter(group==Comparison) %>%
      dplyr::filter(sample==sample_C) %>%
      dplyr::filter(chr==Chrom_Peak[i])
    #建立B组所有位点位置
    assign('pos_C',c())
    for(i2 in 1:nrow(Peak_Comparison_C)){
      pos_C<-c(
        pos_C,
        seq(Peak_Comparison_C$start[i2],Peak_Comparison_C$end[i2])
      )
    }
    
    #取交集位点
    common_pos<-intersect(pos_B,pos_C)
    
    if(length(common_pos)!=0){
      getCommonRegion<-Number_PointsToSegment(common_pos)
      common_pos_frame<-data.frame(matrix(ncol=3,nrow=getCommonRegion$length))
      colnames(common_pos_frame)<-c('chr','start','end')
      common_pos_frame$chr<-Chrom_Peak[i]
      common_pos_frame$start<-getCommonRegion$start
      common_pos_frame$end<-getCommonRegion$end
      CommonPeakRegion<-rbind(CommonPeakRegion,common_pos_frame)
    }
    
    if((i%% 10) ==0){
      print(i)
    }
  }
  
  #注释
  #Annotate
  CommonPeakRegion$GENE_ID<-NA
  pb2 <- progress_bar$new(total = nrow(CommonPeakRegion))#创建一个进度条对象 #Create a progress bar object
  for(i in 1:nrow(CommonPeakRegion)){
    selectedAnnoData<-
      ann_data_gff3 %>%
      dplyr::filter(seqid==CommonPeakRegion$chr[i]) %>%
      dplyr::filter(start<=CommonPeakRegion$start[i]) %>%
      dplyr::filter(end>=CommonPeakRegion$end[i]) %>%
      dplyr::filter(type %in% c('gene'))
    if(nrow(selectedAnnoData)==0){
      selectedAttributes<-NA
    }else{
      selectedAttributes<-
        (separateAttributesCharToListIncludingItemValue(selectedAnnoData$attributes[1]) %>%
           dplyr::filter(Item=='gene_id'))$Value
    }
    CommonPeakRegion$GENE_ID[i]<-selectedAttributes
    if((i %% 100)==0){
      print(i)
    }
    pb2$tick()  #更新进度条#Update progress bar
  }
  saveRDS(CommonPeakRegion,'./CommonPeakRegion_Open.RDS')
  
  CommonPeakRegion_withID<-CommonPeakRegion %>%
    dplyr::filter(!is.na(GENE_ID))
  
  library(ggplot2)
  library(ggforce)
  
  pieplot<-
    ggplot()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.title=element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
    xlab("")+ylab('')+#添加颜色
    scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
                                 '#D6E7A3', '#57C3F3', '#476D87',
                                 '#E59CC4', '#AB3282', '#23452F', '#BD956A'))+
    geom_arc_bar(data=data.frame(V1=c('With Ensembl ID',
                                      'No Ensembl ID'),
                                 ratio=c(
                                   (length(which(!is.na(CommonPeakRegion$GENE_ID)))/nrow(CommonPeakRegion)),
                                   (length(which(is.na(CommonPeakRegion$GENE_ID)))/nrow(CommonPeakRegion))
                                 )),
                 stat = "pie",
                 aes(x0=0,y0=0,r0=1,r=2,,
                     amount=ratio,fill=V1)
    )+#饼图
    annotate("text",x=-1.6,y=1.5,label="37.44%",angle=45)+
    annotate("text",x=1.6,y=-1.5,label="62.56%",angle=45)#手动注释，还是很麻烦
  ggsave('1.piePlot_ifID_Open.pdf',pieplot,width=6.5,height = 5,units = 'in')
  ggsave('1.piePlot_ifID_Open.svg',pieplot,width=6.5,height = 5,units = 'in')
  
  #染色体频数图
  chr_freq_data<-table(CommonPeakRegion$chr) %>%
    data.frame() %>%
    dplyr::filter(str_detect(Var1,'chr'))
  chr_freq_data$Var1<-factor(chr_freq_data$Var1,
                             levels = paste0('chr',c(seq(1:22),'X','Y')))
  chr_freq_plot<-
  ggplot(chr_freq_data)+
    geom_col(aes(x=Var1,y=Freq),fill='#FFFACD',color='black')+
    theme_bw()+
    labs(x='Chromosome',y='Frequence')+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5),
          axis.title = element_text(face='bold'))
  ggsave(filename = '2.Chrom_Freq_Open.pdf',chr_freq_plot,
         height = 6,width = 8,units = 'in')
  
  #非染色体频数图
  Notchr_freq_data<-table(CommonPeakRegion$chr) %>%
    data.frame() %>%
    dplyr::filter(!str_detect(Var1,'chr'))
  Notchr_freq_plot<-
    ggplot(Notchr_freq_data)+
    geom_col(aes(x=Var1,y=Freq),fill='#87CEEB',color='black')+
    theme_bw()+
    labs(x='SeqID',y='Frequence')+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5),
          axis.title = element_text(face='bold'))
  ggsave(filename = '3.NotChrom_Freq_Open.pdf',Notchr_freq_plot,
         height = 6,width = 8,units = 'in')
  
  write.csv(CommonPeakRegion,'4.CommonPeakRegion_Open.csv')
  
  
  #KD后关闭区域
  #定义比较的两组
  Comparison<-'C-KD'
  #提取数据
  Peak_Comparison<-SigPeaks %>%
    dplyr::filter(group==Comparison)
  #获取BC两组比较的样本
  samples<-unique(Peak_Comparison$sample)
  sample_B<-samples[which(str_detect(samples,'vs_B'))]
  sample_C<-samples[which(str_detect(samples,'vs_C'))]
  
  #有两组峰的染色体
  Chrom_Peak<-unique(intersect(
    unique((SigPeaks %>%
              dplyr::filter(group==Comparison) %>%
              dplyr::filter(sample==sample_B))$chr),
    unique((SigPeaks %>%
              dplyr::filter(group==Comparison) %>%
              dplyr::filter(sample==sample_C))$chr)
  ))
  
  #建立结果变量
  CommonPeakRegion<-data.frame(matrix(ncol=3,nrow=0))
  colnames(CommonPeakRegion)<-c('chr','start','end')
  
  #根据染色体建立循环
  for(i in 1:length(Chrom_Peak)){
    #建立B组所有碱基位点的位置数的变量
    Peak_Comparison_B<-SigPeaks %>%
      dplyr::filter(group==Comparison) %>%
      dplyr::filter(sample==sample_B) %>%
      dplyr::filter(chr==Chrom_Peak[i])
    #建立B组所有位点位置
    assign('pos_B',c())
    for(i1 in 1:nrow(Peak_Comparison_B)){
      pos_B<-c(
        pos_B,
        seq(Peak_Comparison_B$start[i1],Peak_Comparison_B$end[i1])
      )
    }
    
    #建立C组所有碱基位点的位置数的变量
    Peak_Comparison_C<-SigPeaks %>%
      dplyr::filter(group==Comparison) %>%
      dplyr::filter(sample==sample_C) %>%
      dplyr::filter(chr==Chrom_Peak[i])
    #建立B组所有位点位置
    assign('pos_C',c())
    for(i2 in 1:nrow(Peak_Comparison_C)){
      pos_C<-c(
        pos_C,
        seq(Peak_Comparison_C$start[i2],Peak_Comparison_C$end[i2])
      )
    }
    
    #取交集位点
    common_pos<-intersect(pos_B,pos_C)
    
    if(length(common_pos)!=0){
      getCommonRegion<-Number_PointsToSegment(common_pos)
      common_pos_frame<-data.frame(matrix(ncol=3,nrow=getCommonRegion$length))
      colnames(common_pos_frame)<-c('chr','start','end')
      common_pos_frame$chr<-Chrom_Peak[i]
      common_pos_frame$start<-getCommonRegion$start
      common_pos_frame$end<-getCommonRegion$end
      CommonPeakRegion<-rbind(CommonPeakRegion,common_pos_frame)
    }
    
    if((i%% 10) ==0){
      print(i)
    }
  }
  
  #注释
  #Annotate
  CommonPeakRegion$GENE_ID<-NA
  pb2 <- progress_bar$new(total = nrow(CommonPeakRegion))#创建一个进度条对象 #Create a progress bar object
  for(i in 1:nrow(CommonPeakRegion)){
    selectedAnnoData<-
      ann_data_gff3 %>%
      dplyr::filter(seqid==CommonPeakRegion$chr[i]) %>%
      dplyr::filter(start<=CommonPeakRegion$start[i]) %>%
      dplyr::filter(end>=CommonPeakRegion$end[i]) %>%
      dplyr::filter(type %in% c('gene'))
    if(nrow(selectedAnnoData)==0){
      selectedAttributes<-NA
    }else{
      selectedAttributes<-
        (separateAttributesCharToListIncludingItemValue(selectedAnnoData$attributes[1]) %>%
           dplyr::filter(Item=='gene_id'))$Value
    }
    CommonPeakRegion$GENE_ID[i]<-selectedAttributes
    if((i %% 100)==0){
      print(i)
    }
    pb2$tick()  #更新进度条#Update progress bar
  }
  saveRDS(CommonPeakRegion,'./CommonPeakRegion_Close.RDS')
  
  CommonPeakRegion_withID<-CommonPeakRegion %>%
    dplyr::filter(!is.na(GENE_ID))
  
  library(ggplot2)
  library(ggforce)
  
  pieplot<-
    ggplot()+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(), 
          axis.text.y = element_blank(),
          axis.text.x = element_blank(),
          legend.title=element_blank(), 
          panel.border = element_blank(),
          panel.background = element_blank())+#去除没用的ggplot背景，坐标轴
    xlab("")+ylab('')+#添加颜色
    scale_fill_manual(values = c('#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', 
                                 '#D6E7A3', '#57C3F3', '#476D87',
                                 '#E59CC4', '#AB3282', '#23452F', '#BD956A'))+
    geom_arc_bar(data=data.frame(V1=c('With Ensembl ID',
                                      'No Ensembl ID'),
                                 ratio=c(
                                   (length(which(!is.na(CommonPeakRegion$GENE_ID)))/nrow(CommonPeakRegion)),
                                   (length(which(is.na(CommonPeakRegion$GENE_ID)))/nrow(CommonPeakRegion))
                                 )),
                 stat = "pie",
                 aes(x0=0,y0=0,r0=1,r=2,,
                     amount=ratio,fill=V1)
    )+#饼图
    annotate("text",x=-1.6,y=1.5,label="37.44%",angle=45)+
    annotate("text",x=1.6,y=-1.5,label="62.56%",angle=45)#手动注释，还是很麻烦
  ggsave('5.piePlot_ifID_Close.pdf',pieplot,width=6.5,height = 5,units = 'in')
  ggsave('5.piePlot_ifID_Close.svg',pieplot,width=6.5,height = 5,units = 'in')
  
  #染色体频数图
  chr_freq_data<-table(CommonPeakRegion$chr) %>%
    data.frame() %>%
    dplyr::filter(str_detect(Var1,'chr'))
  chr_freq_data$Var1<-factor(chr_freq_data$Var1,
                             levels = paste0('chr',c(seq(1:22),'X','Y')))
  chr_freq_plot<-
    ggplot(chr_freq_data)+
    geom_col(aes(x=Var1,y=Freq),fill='#FFFACD',color='black')+
    theme_bw()+
    labs(x='Chromosome',y='Frequence')+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5),
          axis.title = element_text(face='bold'))
  ggsave(filename = '6.Chrom_Freq_Close.pdf',chr_freq_plot,
         height = 6,width = 8,units = 'in')
  
  #非染色体频数图
  Notchr_freq_data<-table(CommonPeakRegion$chr) %>%
    data.frame() %>%
    dplyr::filter(!str_detect(Var1,'chr'))
  Notchr_freq_plot<-
    ggplot(Notchr_freq_data)+
    geom_col(aes(x=Var1,y=Freq),fill='#87CEEB',color='black')+
    theme_bw()+
    labs(x='SeqID',y='Frequence')+
    theme(axis.text.x = element_text(angle = 90,vjust=0.5),
          axis.title = element_text(face='bold'))
  ggsave(filename = '7.NotChrom_Freq_Close.pdf',Notchr_freq_plot,
         height = 6,width = 8,units = 'in')
  
  write.csv(CommonPeakRegion,'8.CommonPeakRegion_Close.csv')
  
  setwd('../')
}

#取交集看看情况
#KD后开放
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DAAs/')
if(!dir.exists('MethodCross')){dir.create('MethodCross')}
x<-list(
  Method1=unique((read.csv('./Method1/3.DAAgenes_KD-C.csv') %>% 
                    dplyr::filter(!is.na(GENE_ID)))$GENE_ID),
  Method2=unique((read.csv('./Method2/4.CommonPeakRegion_Open.csv') %>% 
                    dplyr::filter(!is.na(GENE_ID)))$GENE_ID)
)
Method_intersect_Venn<-
ggVennDiagram(x, , color = 1, lwd = 0.7,label_alpha = 0,
              label_font = "Times") + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = 'KD vs C')+
  theme(legend.position = "none",
        plot.title = element_text(face='bold',hjust=0.5))+
  coord_flip()
ggsave('./MethodCross/1.Venn_Open.pdf',Method_intersect_Venn,height = 3,width=4,
       units = 'in')
ggsave('./MethodCross/1.Venn_Open.svg',Method_intersect_Venn,height = 3,width=4,
       units = 'in')
write.csv(
  data.frame(GeneID=Reduce(intersect,x)),
  './MethodCross/2.MethodCross_Open.csv',row.names = F
)
saveRDS(Reduce(intersect,x),'./MethodCross/2.MethodCross_Open.RDS')

saveRDS(
  unique(c(
    unique((read.csv('./Method1/3.DAAgenes_KD-C.csv') %>% 
              dplyr::filter(!is.na(GENE_ID)))$GENE_ID),
    unique((read.csv('./Method2/4.CommonPeakRegion_Open.csv') %>% 
              dplyr::filter(!is.na(GENE_ID)))$GENE_ID)
  )),
  './MethodCross/3.UnionMethodsDAA_Open.RDS'
)

#KD后关闭
x<-list(
  Method1=unique((read.csv('./Method1/6.DAAgenes_C-KD.csv') %>% 
                    dplyr::filter(!is.na(GENE_ID)))$GENE_ID),
  Method2=unique((read.csv('./Method2/8.CommonPeakRegion_Close.csv') %>% 
                    dplyr::filter(!is.na(GENE_ID)))$GENE_ID)
)
Method_intersect_Venn<-
  ggVennDiagram(x, , color = 1, lwd = 0.7,label_alpha = 0,
                label_font = "Times") + 
  scale_fill_gradient(low = "#F4FAFE", high = "#4981BF") +
  labs(title = 'C vs KD')+
  theme(legend.position = "none",
        plot.title = element_text(face='bold',hjust=0.5))+
  coord_flip()
ggsave('./MethodCross/4.Venn_Close.pdf',Method_intersect_Venn,height = 3,width=4,
       units = 'in')
ggsave('./MethodCross/4.Venn_Close.svg',Method_intersect_Venn,height = 3,width=4,
       units = 'in')
write.csv(
  data.frame(GeneID=Reduce(intersect,x)),
  './MethodCross/5.MethodCross_Close.csv',row.names = F
)
saveRDS(Reduce(intersect,x),'./MethodCross/5.MethodCross_Close.RDS')

saveRDS(
  unique(c(
    unique((read.csv('./Method1/6.DAAgenes_C-KD.csv') %>% 
              dplyr::filter(!is.na(GENE_ID)))$GENE_ID),
    unique((read.csv('./Method2/8.CommonPeakRegion_Close.csv') %>% 
              dplyr::filter(!is.na(GENE_ID)))$GENE_ID)
  )),
  './MethodCross/6.UnionMethodsDAA_Close.RDS'
)

#通路富集：初步探索
{
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
  
  #KD后开放
  {
    setwd('/rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DAAs/MethodCross/')
    if(!dir.exists('Enrichment_Open')){dir.create('Enrichment_Open')};setwd('Enrichment_Open')
    x<-list(
      Method1=unique((read.csv('../../Method1/3.DAAgenes_KD-C.csv') %>% 
                        dplyr::filter(!is.na(GENE_ID)))$GENE_ID),
      Method2=unique((read.csv('../../Method2/4.CommonPeakRegion_Open.csv') %>% 
                        dplyr::filter(!is.na(GENE_ID)))$GENE_ID)
    )
    #候选基因 Candidate Genes
    candg <- Reduce(intersect,x) %>%
      str_sub(1,15) %>%
      unique()
    ## id转换 Convert Gene ID
    symbol2entrezid <- bitr(geneID = candg,
                            fromType = 'ENSEMBL',
                            toType = 'ENTREZID',
                            OrgDb = 'org.Hs.eg.db')
    # GO富集分析 GO Enrichment####---------------------
    pvalueCutoff <- 1
    qvalueCutoff <- 1
    
    ego <- enrichGO(gene = as.numeric(symbol2entrezid$ENTREZID),
                    keyType = "ENTREZID",
                    OrgDb = org.Hs.eg.db,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    ont="ALL",
                    readable =T)
    saveRDS(ego,'ego.rds')
    as.data.frame(ego) %>% filter(pvalue < 0.05) %>%
      group_by(ONTOLOGY) %>% dplyr::count()
    # # A tibble: 3 × 2
    # # Groups:   ONTOLOGY [3]
    # ONTOLOGY     n
    # <chr>    <int>
    #   1 BP         424
    # 2 CC          67
    # 3 MF          72
    write_csv(as.data.frame(ego) %>% filter(pvalue < 0.05),"1.Rich_GO_enrich_sig.csv")
    go.df<-as.data.frame(ego) %>% filter(pvalue < 0.05)
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
    ggsave('2.GO_bubble.pdf',GO_Plot,width = 10,height=6,units = 'in')
    
    # KEGG富集分析 KEGG enrichment####---------------------
    ekegg <- enrichKEGG(gene = symbol2entrezid$ENTREZID ,
                        keyType = "kegg",
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
    saveRDS(ekegg,'ekegg.rds')
    
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
    ggsave('4.KEGG_bubble.pdf',KEGG_Plot,width = 10,height=6,units = 'in')
  }
  
  
  #KD后关闭
  {
    setwd('/rds/projects/g/gendood-preclinomics/EHMT2/ATACSeq_data_Saad/DAAs/MethodCross/')
    if(!dir.exists('Enrichment_Close')){dir.create('Enrichment_Close')};setwd('Enrichment_Close')
    x<-list(
      Method1=unique((read.csv('../../Method1/6.DAAgenes_C-KD.csv') %>% 
                        dplyr::filter(!is.na(GENE_ID)))$GENE_ID),
      Method2=unique((read.csv('../../Method2/8.CommonPeakRegion_Close.csv') %>% 
                        dplyr::filter(!is.na(GENE_ID)))$GENE_ID)
    )
    #候选基因 Candidate Genes
    candg <- Reduce(intersect,x) %>%
      str_sub(1,15) %>%
      unique()
    ## id转换 Convert Gene ID
    symbol2entrezid <- bitr(geneID = candg,
                            fromType = 'ENSEMBL',
                            toType = 'ENTREZID',
                            OrgDb = 'org.Hs.eg.db')
    # GO富集分析 GO Enrichment####---------------------
    pvalueCutoff <- 1
    qvalueCutoff <- 1
    
    ego <- enrichGO(gene = as.numeric(symbol2entrezid$ENTREZID),
                    keyType = "ENTREZID",
                    OrgDb = org.Hs.eg.db,
                    pvalueCutoff = pvalueCutoff,
                    qvalueCutoff = qvalueCutoff,
                    ont="ALL",
                    readable =T)
    saveRDS(ego,'ego.rds')
    as.data.frame(ego) %>% filter(pvalue < 0.05) %>%
      group_by(ONTOLOGY) %>% dplyr::count()
    # # A tibble: 3 × 2
    # # Groups:   ONTOLOGY [3]
    # ONTOLOGY     n
    # <chr>    <int>
    #   1 BP         522
    # 2 CC         113
    # 3 MF         126
    write_csv(as.data.frame(ego) %>% filter(pvalue < 0.05),"1.Rich_GO_enrich_sig.csv")
    go.df<-as.data.frame(ego) %>% filter(pvalue < 0.05)
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
    ggsave('2.GO_bubble.pdf',GO_Plot,width = 10,height=6,units = 'in')
    
    # KEGG富集分析 KEGG enrichment####---------------------
    ekegg <- enrichKEGG(gene = symbol2entrezid$ENTREZID ,
                        keyType = "kegg",
                        organism = "hsa",
                        pAdjustMethod = "BH",
                        pvalueCutoff = 1,
                        qvalueCutoff = 1)
    saveRDS(ekegg,'ekegg.rds')
    
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
    ggsave('4.KEGG_bubble.pdf',KEGG_Plot,width = 10,height=6,units = 'in')
  }
}


