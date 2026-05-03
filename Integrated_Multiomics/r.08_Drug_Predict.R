# 清除环境变量，设置工作目录
rm(list = ls());gc()
num_dir <- '08_Drug_Predict';if (!dir.exists(num_dir)){dir.create(num_dir)};setwd(num_dir)
system('ls -l')


library(enrichR)
library(tidyverse)

dbs <- listEnrichrDbs()

##DsigDB
dbs$libraryName
dbs <- c("DSigDB")
symbol <-  c('PAX6', 'ANXA10', 'DNAH7', 'OTOGL', 'ZP3', 'ATF7-NPFF', 'FAT2', 
             'SH3PXD2A', 'ADGRB2', 'CFAP161', 'ERICH6B', 'LINGO1', 'KSR2', 
             'GARIN1A', 'MUC19', 'EML6')
enrichr <- enrichr(symbol,dbs)
result <- data.frame(enrichr$DSigDB)
result <- result[result$P.value < 0.05, ]
# result <- result[result$Combined.Score > 8000,]

write_csv(result,'1.result.csv')

#建立Cytoscape绘图文件
Cytoscape_edge<-data.frame(
  Drug=NA,
  Gene=NA
)
for(i in 1:nrow(result)){
  TargetGenes<-str_split(result$Genes[i],';')[[1]]
  Cytoscape_edge<-rbind(
    Cytoscape_edge,
    data.frame(
      Drug=rep(result$Term[i],length(TargetGenes)),
      Gene=TargetGenes
    )
  )
}
Cytoscape_edge<-Cytoscape_edge[which(!is.na(Cytoscape_edge$Drug)),]

Cytoscape_node<-data.frame(
  node=c(unique(Cytoscape_edge$Drug),unique(Cytoscape_edge$Gene)),
  type=c(rep('Drug',length(unique(Cytoscape_edge$Drug))),
         rep('Gene',length(unique(Cytoscape_edge$Gene))))
)

write.table(Cytoscape_edge,'Cytoscape_edge.txt',row.names = F,quote = F,sep = '\t')
write.table(Cytoscape_node,'Cytoscape_node.txt',row.names = F,quote = F,sep = '\t')
