#00.清除环境变量，设置工作路径--------------------------------------------------
rm(list = ls());gc()
setwd('/rds/projects/g/gendood-preclinomics/EHMT2/Integrated_Multiomics/')
folder_name<-"07_Factor";if (!dir.exists(folder_name)) {dir.create(folder_name)}
setwd(folder_name)

#加载安装包
library(decoupleR);library(dplyr);library(tibble);library(tidyr);
library(ggplot2);library(pheatmap);library(ggrepel)

#自定义函数
{
  #转换行名到Symbol
  convert_ensembl_to_symbol <- function(expr_mat,
                                        species = c("human", "mouse"),
                                        method = "mean",
                                        drop_unmapped = TRUE) {
    # 参数匹配
    species <- match.arg(species)
    method <- match.arg(method, choices = "mean")
    
    # 确保输入为 matrix
    if (!is.matrix(expr_mat)) {
      expr_mat <- as.matrix(expr_mat)
    }
    
    # 1. 提取并清理 Ensembl ID (去除版本号)
    ensembl_ids <- rownames(expr_mat)
    ensembl_ids_clean <- sub("\\..*$", "", ensembl_ids)  # 去掉 .数字 后缀
    
    # 2. 加载对应物种的注释包
    if (species == "human") {
      if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
        stop("请先安装 Bioconductor 包: org.Hs.eg.db")
      }
      db <- org.Hs.eg.db::org.Hs.eg.db
    } else if (species == "mouse") {
      if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) {
        stop("请先安装 Bioconductor 包: org.Mm.eg.db")
      }
      db <- org.Mm.eg.db::org.Mm.eg.db
    }
    
    # 3. 映射 Ensembl ID 到 Symbol
    # 使用 AnnotationDbi::select 获取对应关系
    # 注意: select 函数来自 AnnotationDbi，通过 db 对象调用
    suppressMessages({
      mapping <- AnnotationDbi::select(db,
                                       keys = ensembl_ids_clean,
                                       columns = c("ENSEMBL", "SYMBOL"),
                                       keytype = "ENSEMBL")
    })
    
    # 去重 (同一个 Ensembl ID 可能有多个 Symbol，我们保留第一个非 NA 的 Symbol)
    mapping <- mapping[!is.na(mapping$SYMBOL), ]
    mapping <- mapping[!duplicated(mapping$ENSEMBL), ]
    
    # 4. 将原始表达矩阵的行名替换为 Symbol
    # 建立 Ensembl clean -> Symbol 的查找表
    symbol_lookup <- setNames(mapping$SYMBOL, mapping$ENSEMBL)
    new_rownames <- symbol_lookup[ensembl_ids_clean]
    
    # 处理无法映射的行
    if (drop_unmapped) {
      keep_idx <- !is.na(new_rownames)
      if (sum(keep_idx) == 0) {
        stop("没有能够映射到 Symbol 的 Ensembl ID，请检查物种或表达矩阵。")
      }
      expr_mat <- expr_mat[keep_idx, , drop = FALSE]
      new_rownames <- new_rownames[keep_idx]
    } else {
      # 若保留未映射行，则用原 Ensembl ID (清理版本号) 作为行名
      new_rownames[is.na(new_rownames)] <- ensembl_ids_clean[is.na(new_rownames)]
    }
    
    rownames(expr_mat) <- new_rownames
    
    # 5. 处理重复 Symbol: 按 Symbol 分组，对每个样本求均值
    # 检查是否有重复
    if (anyDuplicated(rownames(expr_mat))) {
      # 按行名分组，计算每个样本的均值
      # 使用 rowsum 按组求和，再除以每组计数得到均值
      groups <- rownames(expr_mat)
      unique_groups <- unique(groups)
      
      # 初始化结果矩阵
      result <- matrix(NA_real_, 
                       nrow = length(unique_groups),
                       ncol = ncol(expr_mat),
                       dimnames = list(unique_groups, colnames(expr_mat)))
      
      for (i in seq_along(unique_groups)) {
        grp <- unique_groups[i]
        idx <- which(groups == grp)
        if (length(idx) == 1) {
          result[i, ] <- expr_mat[idx, ]
        } else {
          # 按列求均值 (处理可能存在的 NA 值)
          result[i, ] <- colMeans(expr_mat[idx, , drop = FALSE], na.rm = TRUE)
        }
      }
      return(result)
    } else {
      return(expr_mat)
    }
  }
}

#读取数据
##表达矩阵
counts<-read.table('../../RNA_Seq/results/salmon.merged.gene_counts.tsv',
                   header = T) %>%
  column_to_rownames('gene_id') %>%
  dplyr::select(-c('gene_name')) %>%
  dplyr::select(paste0(rep(c('A','B','C'),each=3),rep(c('_1','_2','_3'),3)))
counts<-convert_ensembl_to_symbol(counts,'human')
##样本信息
design <- data.frame(
  sample=paste0(rep(c('A','B','C'),each=3),rep(c('_1','_2','_3'),3)),
  condition=rep(c('Control','KD1','KD2'),3)
)
##差异基因结果
deg_KD1<-readRDS('../../RNA_Seq/DOWNSTREAM_ANALYSIS/01_DEGs/EHMT2_KD1_vs_Control/1.DESeq2_res1.RDS') %>%
  dplyr::select(c('log2FoldChange','stat','pvalue')) %>%
  dplyr::rename('logFC'='log2FoldChange',t='stat','P.Value'='pvalue') %>%
  dplyr::filter(!is.na(t)) %>%
  as.matrix()
deg_KD2<-readRDS('../../RNA_Seq/DOWNSTREAM_ANALYSIS/01_DEGs/EHMT2_KD2_vs_Control/1.DESeq2_res1.RDS') %>%
  dplyr::select(c('log2FoldChange','stat','pvalue')) %>%
  dplyr::rename('logFC'='log2FoldChange',t='stat','P.Value'='pvalue') %>%
  dplyr::filter(!is.na(t)) %>%
  as.matrix()

#加载知识网络
if(!file.exists('net.RDS')){
  net <- decoupleR::get_collectri(organism = 'human',
                                  split_complexes = FALSE)
  saveRDS(net,'net.RDS')
}else{
  net<-readRDS('./net.RDS')
}

#使用ULM（Univariate Linear Model）活性推断
sample_acts <- decoupleR::run_ulm(mat = counts,
                                  net = net,
                                  .source = 'source',
                                  .target = 'target',
                                  .mor = 'mor',
                                  minsize = 5)

#可视化
#前30名TF活性热图
{
  n_tfs <- 30
  
  # Transform to wide matrix
  sample_acts_mat <- sample_acts %>%
    tidyr::pivot_wider(
      id_cols = 'condition',
      names_from = 'source',
      values_from = 'score'
    ) %>%
    tibble::column_to_rownames('condition') %>%
    as.matrix()
  
  # Get top tfs with more variable means across clusters
  tfs <- sample_acts %>%
    dplyr::group_by(source) %>%
    dplyr::summarise(std = stats::sd(score)) %>%
    dplyr::arrange(-abs(std)) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  sample_acts_mat <- sample_acts_mat[, tfs]
  
  # Scale per sample
  sample_acts_mat <- scale(sample_acts_mat)
  
  # Choose color palette
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
  colors.use <- grDevices::colorRampPalette(colors = colors)(100)
  
  my_breaks <- c(
    seq(-2, 0, length.out = ceiling(100 / 2) + 1),
    seq(0.05, 2, length.out = floor(100 / 2))
  )
  
  # Plot
  pheatmap::pheatmap(
    mat = sample_acts_mat,
    color = colors.use,
    border_color = "white",
    breaks = my_breaks,
    cellwidth = 15,
    cellheight = 15,
    treeheight_row = 20,
    treeheight_col = 20
  )
  p1<-recordPlot()
  export::graph2svg(p1,'1.Top30TF_heatmap.svg',width=8,height=4)
}

#基于差异分析的 TF 活性对比
##Limma分析差异基因
{
  if(!dir.exists('limma')){dir.create('limma')};setwd('limma')
  #阈值
  logFC_cutoff<-1
  PValue_cutoff<-0.05
  #加载R包
  library(limma)
  library(DESeq2)
  library(tidyverse)
  library(ggrepel)
  library(psych)
  library(corrplot)
  library(ComplexHeatmap)
  select=dplyr::select
  #KD1
  {
    if(!dir.exists('KD1')){dir.create('KD1')};setwd('KD1')
    
    # 输入表达矩阵和分组文件
    dat<-counts %>% as.data.frame() %>% 
      dplyr::select(paste0(rep(c('A','B','C'),each=2),rep(c('_1','_2'),3)))
    group<-data.frame(
      sample=paste0(rep(c('A','B','C'),each=2),rep(c('_1','_2'),3)),
      group=rep(c('Control','KD1'),3)
    )
    
    #分组因子化
    group$group <- factor(group$group, levels = c("Control", "KD1"))
    #建立分组矩阵
    design.mat <- cbind(Control = ifelse(group$group == "Control", 1, 0), 
                        KD1 = ifelse(group$group == "Control", 0, 1))
    
    #差异分析
    contrast.mat <- makeContrasts(contrasts="KD1-Control", levels=design.mat)
    fit <- lmFit(dat, design.mat)
    fit <- contrasts.fit(fit, contrast.mat)
    fit <- eBayes(fit)
    fit <- topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
    DEG <- na.omit(fit)
    res1<-DEG
    
    #依次按照P.Value值logFC值进行排序
    res1 <- res1[order(res1$P.Value, res1$logFC, decreasing = c(FALSE, TRUE)), ]
    if(length(which(is.na(res1$P.Value)))>0){#判断是否有行P值不存在
      res1<-res1[-which(is.na(res1$P.Value)),]#如果有，去除该行
    }
    
    #判断变化
    res1$Change<-NA#添加变化的列
    for (i in 1:nrow(res1)) {#循环和判断语句逐个加变化标签
      if(res1$P.Value[i]<PValue_cutoff){
        if(res1$logFC[i]>logFC_cutoff){
          res1$Change[i]<-'Up'
        } else if (res1$logFC[i]<(-logFC_cutoff)) { 
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
    # 188  36458    276
    res1$Sig<-NA
    for (i in 1:nrow(res1)) {#循环和判断语句逐个加显著性标签
      if(res1$P.Value[i]<PValue_cutoff){
        res1$Sig[i]<-'Sig'
      } else {
        res1$Sig[i]<-
          'Stable'
      }
    }
    
    write.csv(res1,'01.KD1_limma_res1.csv',row.names = T)
    
    #03.火山图----------------------------------------------------------------------
    sig_diff<-res1[which(!res1$Change=='Stable'),]
    #挑选变化倍数前10
    dat_rep <- rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                     head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))
    
    #绘图
    volcano_plot <- ggplot(data = res1,aes(x = logFC,y = -log10(P.Value),color =Change)) +
      scale_color_manual(values = c("blue", "darkgray","darkorange")) +
      scale_x_continuous(breaks = seq(-5,5,1),limits = c(-2,2)) +
      scale_y_continuous(trans = "log1p",limits = c(0,15),breaks=c(0,1,3,10,30)) +
      geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
      theme_bw(base_size = 12, base_family = "Times") +
      geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff), lty = 4, col = "darkgray", lwd = 0.6)+
      geom_hline(yintercept = -log10(PValue_cutoff), lty = 4, col = "darkgray", lwd = 0.6)+
      theme(legend.position = "right",
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(face="bold",color="black",family = "Times",size=13),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(face = "bold",color = "black",size = 15),
            axis.text.y = element_text(face = "bold",color = "black",size = 15),
            axis.title.x = element_text(face = "bold",color = "black",size = 15),
            axis.title.y = element_text(face = "bold",color = "black",size = 15)) +
      geom_label_repel(data = dat_rep, aes(label = rownames(dat_rep)),
                       max.overlaps = 20, size = 4,
                       box.padding = unit(0.5, "lines"),
                       min.segment.length = 0,
                       point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
      labs(x = "log2(Fold Change)",y = "-log10 (P.Value)")
    #输出文件
    export::graph2png(volcano_plot,
                      file='./02.KD1_DEGs_volcano.png',
                      width = 8.78, height = 7.43)
    export::graph2pdf(volcano_plot,
                      file='./02.KD1_DEGs_volcano.pdf',
                      width = 8.78, height = 7.43)
    
    
    #04.密度热图--------------------------------------------------------------------
    #前十FC的上下调差异基因
    top10_gene_down <- rownames(sig_diff %>% 
                                  arrange(logFC) %>% 
                                  head(10))
    top10_gene_up <- rownames(sig_diff %>% 
                                arrange(desc(logFC)) %>% 
                                head(10))
    top_gene <- c(top10_gene_down,top10_gene_up)
    #分组信息
    rt <- dat[top_gene,group$sample]
    group<-group[order(group$group),]
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
      HeatmapAnnotation(Group = group$group, col = list(Group = c("KD1" = "#B72230",
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
    ppp
    #输出文件
    export::graph2png(ppp, 
                      file="./03.KD1_DEGs_deg_top_heatmap.png",
                      width=6.51,height=5.4)
    export::graph2pdf(ppp, 
                      file="./03.KD1_DEGs_deg_top_heatmap.pdf",
                      width=6.51,height=5.4)
    
    setwd('../')
  }
  #KD2
  {
    if(!dir.exists('KD2')){dir.create('KD2')};setwd('KD2')
    
    # 输入表达矩阵和分组文件
    dat<-counts %>% as.data.frame() %>% 
      dplyr::select(paste0(rep(c('A','B','C'),each=2),rep(c('_1','_3'),3)))
    group<-data.frame(
      sample=paste0(rep(c('A','B','C'),each=2),rep(c('_1','_3'),3)),
      group=rep(c('Control','KD2'),3)
    )
    
    #分组因子化
    group$group <- factor(group$group, levels = c("Control", "KD2"))
    #建立分组矩阵
    design.mat <- cbind(Control = ifelse(group$group == "Control", 1, 0), 
                        KD2 = ifelse(group$group == "Control", 0, 1))
    
    #差异分析
    contrast.mat <- makeContrasts(contrasts="KD2-Control", levels=design.mat)
    fit <- lmFit(dat, design.mat)
    fit <- contrasts.fit(fit, contrast.mat)
    fit <- eBayes(fit)
    fit <- topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
    DEG <- na.omit(fit)
    res1<-DEG
    
    #依次按照P.Value值logFC值进行排序
    res1 <- res1[order(res1$P.Value, res1$logFC, decreasing = c(FALSE, TRUE)), ]
    if(length(which(is.na(res1$P.Value)))>0){#判断是否有行P值不存在
      res1<-res1[-which(is.na(res1$P.Value)),]#如果有，去除该行
    }
    
    #判断变化
    res1$Change<-NA#添加变化的列
    for (i in 1:nrow(res1)) {#循环和判断语句逐个加变化标签
      if(res1$P.Value[i]<PValue_cutoff){
        if(res1$logFC[i]>logFC_cutoff){
          res1$Change[i]<-'Up'
        } else if (res1$logFC[i]<(-logFC_cutoff)) { 
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
    # 145  36516    261 
    res1$Sig<-NA
    for (i in 1:nrow(res1)) {#循环和判断语句逐个加显著性标签
      if(res1$P.Value[i]<PValue_cutoff){
        res1$Sig[i]<-'Sig'
      } else {
        res1$Sig[i]<-
          'Stable'
      }
    }
    
    write.csv(res1,'01.KD2_limma_res1.csv',row.names = T)
    
    #03.火山图----------------------------------------------------------------------
    sig_diff<-res1[which(!res1$Change=='Stable'),]
    #挑选变化倍数前10
    dat_rep <- rbind(head(sig_diff[order(sig_diff$logFC,decreasing = T),],10),
                     head(sig_diff[order(sig_diff$logFC,decreasing = F),],10))
    
    #绘图
    volcano_plot <- ggplot(data = res1,aes(x = logFC,y = -log10(P.Value),color =Change)) +
      scale_color_manual(values = c("blue", "darkgray","darkorange")) +
      scale_x_continuous(breaks = seq(-5,5,1),limits = c(-2,2)) +
      scale_y_continuous(trans = "log1p",limits = c(0,15),breaks=c(0,1,3,10,30)) +
      geom_point(size = 1.2, alpha = 0.4, na.rm=T) +
      theme_bw(base_size = 12, base_family = "Times") +
      geom_vline(xintercept = c(-logFC_cutoff,logFC_cutoff), lty = 4, col = "darkgray", lwd = 0.6)+
      geom_hline(yintercept = -log10(PValue_cutoff), lty = 4, col = "darkgray", lwd = 0.6)+
      theme(legend.position = "right",
            panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(face="bold",color="black",family = "Times",size=13),
            plot.title = element_text(hjust = 0.5),
            axis.text.x = element_text(face = "bold",color = "black",size = 15),
            axis.text.y = element_text(face = "bold",color = "black",size = 15),
            axis.title.x = element_text(face = "bold",color = "black",size = 15),
            axis.title.y = element_text(face = "bold",color = "black",size = 15)) +
      geom_label_repel(data = dat_rep, aes(label = rownames(dat_rep)),
                       max.overlaps = 20, size = 4,
                       box.padding = unit(0.5, "lines"),
                       min.segment.length = 0,
                       point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
      labs(x = "log2(Fold Change)",y = "-log10 (P.Value)")
    #输出文件
    export::graph2png(volcano_plot,
                      file='./02.KD2_DEGs_volcano.png',
                      width = 8.78, height = 7.43)
    export::graph2pdf(volcano_plot,
                      file='./02.KD2_DEGs_volcano.pdf',
                      width = 8.78, height = 7.43)
    
    
    #04.密度热图--------------------------------------------------------------------
    #前十FC的上下调差异基因
    top10_gene_down <- rownames(sig_diff %>% 
                                  arrange(logFC) %>% 
                                  head(10))
    top10_gene_up <- rownames(sig_diff %>% 
                                arrange(desc(logFC)) %>% 
                                head(10))
    top_gene <- c(top10_gene_down,top10_gene_up)
    #分组信息
    rt <- dat[top_gene,group$sample]
    group<-group[order(group$group),]
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
      HeatmapAnnotation(Group = group$group, col = list(Group = c("KD2" = "#B72230",
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
    ppp
    #输出文件
    export::graph2png(ppp, 
                      file="./03.KD2_DEGs_deg_top_heatmap.png",
                      width=6.51,height=5.4)
    export::graph2pdf(ppp, 
                      file="./03.KD2_DEGs_deg_top_heatmap.pdf",
                      width=6.51,height=5.4)
    
    setwd('../')
  }
  
  setwd('../')
}

##KD1
{
  deg_KD1<-read.csv('./limma/KD1/01.KD1_limma_res1.csv',row.names = 1)
  contrast_acts_KD1 <- decoupleR::run_ulm(
    mat = deg_KD1[, 't', drop = FALSE],
    net = net,
    .source = 'source',
    .target = 'target',
    .mor = 'mor',
    minsize = 5
  )
  
  f_contrast_acts_KD1 <- contrast_acts_KD1 %>%
    dplyr::mutate(rnk = NA)
  
  msk_KD1 <- f_contrast_acts_KD1$score > 0
  
  f_contrast_acts_KD1[msk_KD1, 'rnk'] <- rank(-f_contrast_acts_KD1[msk_KD1, 'score'])
  f_contrast_acts_KD1[!msk_KD1, 'rnk'] <- rank(-abs(f_contrast_acts_KD1[!msk_KD1, 'score']))
  
  tfs_KD1 <- f_contrast_acts_KD1 %>%
    dplyr::arrange(rnk) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  f_contrast_acts_KD1 <- f_contrast_acts_KD1 %>%
    filter(source %in% tfs_KD1)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  p_KD1 <- ggplot2::ggplot(
    data = f_contrast_acts_KD1,
    mapping = ggplot2::aes(
      x = stats::reorder(source, score),
      y = score
    )
  ) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(fill = score),
      color = "black",
      stat = "identity"
    ) +
    ggplot2::scale_fill_gradient2(
      low = colors[1],
      mid = "whitesmoke",
      high = colors[2],
      midpoint = 0
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ggplot2::xlab("TFs")
  ggsave('2.KD1_TFs-Score.svg',p_KD1,height = 5,width = 7.5)
}

##KD2
{
  deg_KD2<-read.csv('./limma/KD2/01.KD2_limma_res1.csv',row.names = 1)
  contrast_acts_KD2 <- decoupleR::run_ulm(
    mat = deg_KD2[, 't', drop = FALSE],
    net = net,
    .source = 'source',
    .target = 'target',
    .mor = 'mor',
    minsize = 5
  )
  
  f_contrast_acts_KD2 <- contrast_acts_KD2 %>%
    dplyr::mutate(rnk = NA)
  
  msk_KD2 <- f_contrast_acts_KD2$score > 0
  
  f_contrast_acts_KD2[msk_KD2, 'rnk'] <- rank(-f_contrast_acts_KD2[msk_KD2, 'score'])
  f_contrast_acts_KD2[!msk_KD2, 'rnk'] <- rank(-abs(f_contrast_acts_KD2[!msk_KD2, 'score']))
  
  tfs_KD2 <- f_contrast_acts_KD2 %>%
    dplyr::arrange(rnk) %>%
    head(n_tfs) %>%
    dplyr::pull(source)
  
  f_contrast_acts_KD2 <- f_contrast_acts_KD2 %>%
    filter(source %in% tfs_KD2)
  
  colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
  p_KD2 <- ggplot2::ggplot(
    data = f_contrast_acts_KD2,
    mapping = ggplot2::aes(
      x = stats::reorder(source, score),
      y = score
    )
  ) +
    ggplot2::geom_bar(
      mapping = ggplot2::aes(fill = score),
      color = "black",
      stat = "identity"
    ) +
    ggplot2::scale_fill_gradient2(
      low = colors[1],
      mid = "whitesmoke",
      high = colors[2],
      midpoint = 0
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      axis.title = element_text(face = "bold", size = 12),
      axis.text.x = ggplot2::element_text(
        angle = 45,
        hjust = 1,
        size = 10,
        face = "bold"
      ),
      axis.text.y = ggplot2::element_text(
        size = 10,
        face = "bold"
      ),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    ) +
    ggplot2::xlab("TFs")
  ggsave('2.KD2_TFs-Score.svg',p_KD2,height = 5,width = 7.5)
}

#候选基因
candg<-c('PAX6', 'ANXA10', 'DNAH7', 'OTOGL', 'ZP3', 'NA', 'ATF7-NPFF', 'NA', 
         'FAT2', 'SH3PXD2A', 'ADGRB2', 'CFAP161', 'ERICH6B', 'LINGO1', 'KSR2', 
         'GARIN1A', 'MUC19', 'EML6')
{
  #仅有PAX6有数据
  
  #候选基因里TF活性热图
  {
    sample_acts_mat_candg<-sample_acts %>% 
      dplyr::filter(source %in% candg) %>%
      tidyr::pivot_wider(
        id_cols = 'condition',
        names_from = 'source',
        values_from = 'score'
      ) %>%
      tibble::column_to_rownames('condition') %>%
      as.matrix()
    
    tfs_candg <- sample_acts$source[which(sample_acts$source %in% candg)] %>% unique()
    
    sample_acts_mat_candg <- sample_acts_mat_candg[, tfs_candg]
    
    # Scale per sample
    sample_acts_mat_candg <- scale(sample_acts_mat_candg)
    
    # Choose color palette
    colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu"))
    colors.use <- grDevices::colorRampPalette(colors = colors)(100)
    
    my_breaks <- c(
      seq(-2, 0, length.out = ceiling(100 / 2) + 1),
      seq(0.05, 2, length.out = floor(100 / 2))
    )
    
    # Plot
    sample_acts_mat_candg_t <- t(sample_acts_mat_candg)
    rownames(sample_acts_mat_candg_t)<-tfs_candg
    pheatmap::pheatmap(
      mat = sample_acts_mat_candg_t,
      color = colors.use,
      border_color = "white",
      breaks = my_breaks,
      cellwidth = 30,
      cellheight = 30,
      treeheight_row = 20,
      treeheight_col = 20,
      cluster_rows = FALSE,   # 关闭行聚类
      cluster_cols = FALSE   # 关闭列聚类
    )
    p_candg<-recordPlot()
    export::graph2svg(p_candg,'3.TF_heatmap_heatmap.svg',width=8,height=4)
  }
  #PAX6下游调节基因
  {
    volPlotByKD<-function(KD){
      if(KD=='KD1'){
        deg<-deg_KD1
      }else if(KD=='KD2'){
        deg<-deg_KD2
      }
      tf <- 'PAX6'
      df <- net %>%
        dplyr::filter(source == tf) %>%
        dplyr::arrange(target) %>%
        dplyr::mutate(ID = target, color = "3") %>%
        tibble::column_to_rownames('target')
      inter <- sort(dplyr::intersect(rownames(deg), rownames(df)))
      df <- df[inter, ]
      df[,c('logfc', 't_value', 'p_value')] <- deg[inter, c('logFC', 't', 'P.Value')]
      df <- df %>%
        dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value > 0, '1', color)) %>% 
        dplyr::mutate(color = dplyr::if_else(mor > 0 & t_value < 0, '2', color)) %>%
        dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value > 0, '2', color)) %>% 
        dplyr::mutate(color = dplyr::if_else(mor < 0 & t_value < 0, '1', color))
      colors <- rev(RColorBrewer::brewer.pal(n = 11, name = "RdBu")[c(2, 10)])
      p <- ggplot2::ggplot(data = df,
                           mapping = ggplot2::aes(x = logfc,
                                                  y = -log10(p_value),
                                                  color = color, 
                                                  size = abs(mor))) +
        ggplot2::geom_point(size = 2.5, 
                            color = "black") +
        ggplot2::geom_point(size = 1.5) +
        ggplot2::scale_colour_manual(values = c(colors[2], colors[1], "grey")) + 
        ggrepel::geom_label_repel(mapping = ggplot2::aes(label = ID,
                                                         size = 1)) +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = "none",
                       plot.title = ggplot2::element_text(face = 'bold',hjust = 0.5),
                       axis.title = ggplot2::element_text(face = 'bold')) +
        ggplot2::geom_vline(xintercept = 0, linetype = 'dotted') +
        ggplot2::geom_hline(yintercept = 0, linetype = 'dotted') +
        ggplot2::labs(title = KD,x='logFC',y='-log10(P.Value)')
      return(
        list(
          data=df,
          plot=p
        )
      )
    }
    KD1_downGenes<-volPlotByKD('KD1')
  }

}
