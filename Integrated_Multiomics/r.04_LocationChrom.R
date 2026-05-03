#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/Integrated_Multiomics/')
dir_name<-'04_LocationChrom';if (!dir.exists(dir_name)) {dir.create(dir_name)}
setwd(dir_name)

#加载R包
#Load R Packages
library(stringr)
library(progress)
library(tidyverse)
library(ggrepel)
library(RCircos)
library(rtracklayer)
library(dplyr)
library(ggplot2)
data("UCSC.HG38.Human.CytoBandIdeogram")

#读取数据
candg_up<-unique(read.csv('../03_Cor_Anno/3.CandgSymbol_OpenUP.csv')$external_gene_name)
candg_up<-candg_up[which(candg_up!='')]
candg_down<-unique(read.csv('../03_Cor_Anno/3.CandgSymbol_CloseDOWN.csv')$external_gene_name)
candg_down<-candg_down[which(candg_down!='')]
candg<-c(candg_up,candg_down,'EHMT2')

#-------------------------------染色体定位--------------------------------------
#参数设置
RCircos.Set.Core.Components(UCSC.HG38.Human.CytoBandIdeogram,#这是上面load的基因组文件
                            chr.exclude<- NULL, #无排除的染色体
                            tracks.inside=10, #染色体圆圈内部一共要画10个圆圈
                            tracks.outside=0)#在外部画0个圆圈
settings<-RCircos.Get.Plot.Parameters()
settings$text.size<-1
RCircos.Reset.Plot.Parameters(settings)

#骨架
RCircos.Set.Plot.Area()#建立一个画板
RCircos.Chromosome.Ideogram.Plot()#在当前的画板上画基因组的圆圈骨架

#候选基因位置信息
if(!file.exists('01.candg_location.csv')){
  if(!file.exists('Homo_sapiens.GRCh38.115.chr.gtf.gz')){
    download.file(
      url='https://ftp.ensembl.org/pub/release-115/gtf/homo_sapiens/Homo_sapiens.GRCh38.115.chr.gtf.gz',
      destfile = 'Homo_sapiens.GRCh38.115.chr.gtf.gz'
    )
  }
  #读取基因组注释文件，时间略长
  Hs_GRCh38.115<-import('./Homo_sapiens.GRCh38.115.chr.gtf.gz')
  #转换成数据框
  Hs_GRCh38.115_frame<-data.frame(Hs_GRCh38.115)
  #提取候选基因位置信息
  candg_location<-Hs_GRCh38.115_frame %>% 
    filter(type=='gene') %>% #保留的是完整基因的位置
    filter(gene_name %in% candg) %>% #筛选出候选基因
    dplyr::select('seqnames','start','end','gene_name') %>% #保留有用的列
    dplyr::rename('Chromosome'='seqnames',
                  'chromStart'='start',
                  'chromEnd'='end',
                  'Gene'='gene_name')#重命名列
  #染色体字符加上chr
  candg_location$Chromosome<-paste0('chr',candg_location$Chromosome)
  #保存
  write.csv(candg_location,'01.candg_location.csv',row.names = F)
}else{
  candg_location<-read.csv('01.candg_location.csv')
}

name.col <- 4 #数据是4列
side <- "in" #画在基因组骨架的内侧
track.num <- 1 #基因组骨架内侧的第一个track位置上画图
RCircos.Gene.Connector.Plot(candg_location,
                            + track.num, side)#画connector（连接基因名称和基因组位置）
track.num <- 2
RCircos.Gene.Name.Plot(candg_location,
                       + name.col,track.num, side)#加基因名称
geneChrLocation<-recordPlot()

png('02.geneChrLocation.png',width=8,height=8,unit='in',res=600)
geneChrLocation
dev.off()
pdf('02.geneChrLocation.pdf',width=8,height=8)
geneChrLocation
dev.off()


