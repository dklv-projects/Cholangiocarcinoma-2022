##### 细胞亚群markers生存分析
## 表达矩阵预处理
mtx <- readRDS("/home/rlhua/project/cholangiocarcinoma/result/bulkRNAseq/FPKM_gencode.v22.annotation.gtf_AllSampleCountMatrix.RDS") ##读入表达矩阵
geneID2symbol <- fread("/home/rlhua/project/cholangiocarcinoma/result/bulkRNAseq/geneID2symbol.tsv",sep = "\t",header = F,col.names = c("ID","symbol")) ## 读入ensenbl ID和gene symbol对照表
idx <- match(rownames(mtx), geneID2symbol$ID) 
mtx <- mtx[!is.na(idx),]
mtx[["symbol"]] <- geneID2symbol$symbol[idx[!is.na(idx)]] ## 添加gene symbol
mtx <- aggregate(mtx[,-ncol(mtx)],list(mtx$symbol),sum) ## 相同gene symbol的ID的表达值相加
rownames(mtx) <- mtx$Group.1
mtx <- mtx[,-1]
mtx <- as.matrix(mtx)
mtx <- mtx[rowMeans(mtx)>0,] ## 去除在所有样本中表达量为0的基因
mtx <- log1p(mtx) ## 表达值取对数

## 72个患者生存数据预处理
survival_info <- fread("/home/rlhua/project/cholangiocarcinoma/result/bulkRNAseq/survival_info.tsv",sep = "\t",header = T)
colnames(survival_info) <- c("sample","time","status","type","histology")
survival_info$status <- ifelse(survival_info$status=="Alive",0,1)

## M3_clustermarker2
# add RNA markers to integrated markers
lst <- list()
for (i in unique(df_integrated$celltype)) {
  df1 <- df_RNA[df_RNA$celltype==i]
  df2 <- df_integrated[df_integrated$celltype==i]
  df1 <- df1[!df1$gene %in% df2$gene]
  df <- rbind(df1[,c("gene","cluster","celltype")],df2[,c("gene","cluster","celltype")])
  lst[[i]] <- df
}
df <- rbindlist(lst)

## 生存曲线和箱线图
clusters <- unique(df$cluster)
lst <- list()
lst1 <- list()
lst2 <- list()
for (i in (1:length(clusters))[24:53]) {
  cat(clusters[i],"\n")
  genes <- df[cluster==clusters[i],gene] ## 取marker基因
  genes <- intersect(genes,rownames(mtx)) ## marker基因与表达矩阵基因取交集
  gsva_score <- GSVA::gsva(mtx,list(genes)) ## GSVA
  dt <- survival_info
  dt$gsva_score <- gsva_score[,dt$sample]
  dt$type <- factor(dt$type,levels = c("normal","paratumor","cancer"))
  dt_n <- dt[type=="normal"]
  dt_p <- dt[type=="paratumor"]
  dt_t <- dt[type=="cancer"]
  dt_n <- dt_n[!gsva_score %in% boxplot(dt_n$gsva_score)$out] 
  dt_p <- dt_p[!gsva_score %in% boxplot(dt_p$gsva_score)$out]
  dt_t <- dt_t[!gsva_score %in% boxplot(dt_t$gsva_score)$out]
  ymax <- max(c(dt_n$gsva_score,dt_p$gsva_score,dt_t$gsva_score))*1.2
  ## 箱线图
  lst[[clusters[i]]] <- ggboxplot(dt,x="type",y="gsva_score",fill = "type",xlab = F,ylab = "GSVA score",
                                  title = clusters[i],
                                  bxp.errorbar=T,outlier.shape=NA,ylim=c(min(dt$gsva_score)+0.0001,ymax))+
    theme_classic(base_size = 10)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none",axis.title.x = element_blank())+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.x.npc = 0.5,label.y = ymax*0.95,size=3)
  dt <- dt[type=="cancer"] 
  dt$group <- ifelse(dt$gsva_score > median(dt$gsva_score),"high","low")
  dt$group <- factor(dt$group,levels = c("high","low"))
  fit <- survfit(Surv(time,status)~group,data = dt)
  data.survdiff <- survdiff(Surv(time,status)~group,data = dt)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  HR = (data.survdiff$obs[1]/data.survdiff$exp[1])/(data.survdiff$obs[2]/data.survdiff$exp[2])
  up95 = exp(log(HR) + qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
  low95 = exp(log(HR) - qnorm(0.975)*sqrt(1/data.survdiff$exp[1]+1/data.survdiff$exp[2]))
  HR <- paste("Hazard Ratio = ", round(HR,2), sep = "")
  CI <- paste("95% CI: ", paste(round(low95,2), round(up95,2), sep = " - "), sep = "")
  ## 生存曲线图
  lst1[[clusters[i]]] <- ggsurvplot(fit,
                                    conf.int = F, #不画置信区间，想画置信区间就把F改成T
                                    #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                                    censor = F, #不显示观察值所在的位置
                                    palette = c("#D95F02","#1B9E77"), #线的颜色对应高、低
                                    xlab="Months",
                                    legend=c(0.85,0.95),
                                    title=clusters[i],
                                    legend.title = "",#基因名写在图例题目的位置
                                    font.legend = 11,#图例的字体大小
                                    #font.title = 12,font.x = 10,font.y = 10,#设置其他字体大小
                                    legend.labs=c("high","low"),
                                    #在左下角标出pvalue、HR、95% CI
                                    #太小的p value标为p < 0.001
                                    pval = paste(pval = ifelse(p.val < 0.001, "p < 0.001", 
                                                               paste("p = ",round(p.val,3), sep = "")),
                                                 HR, CI, sep = "\n"),
                                    pval.size=4,pval.coord=c(0,0.2),
                                    ggtheme = theme_survminer()+theme(plot.title = element_text(hjust = 0.5)))
  dt <- lapply(genes,function(x){
    cor(as.vector(mtx[x,]),as.vector(gsva_score),method="pearson")
  }) %>% unlist %>% data.table(cor=.)
  dt$gene <- genes
  dt$cluster <- clusters[i]
  dt <- dt[order(dt$cor,decreasing = T)]
  lst2[[i]] <- dt
  lst3 <- lapply(dt$gene, function(x){
    data_cor <- round(cor(as.vector(mtx[x,]),as.vector(gsva_score),method="pearson"),4)
    ggplot(mapping = aes(x=as.vector(mtx[x,]),y=as.vector(gsva_score)))+geom_point(size=.5)+
      ggtitle(paste0("cor=",data_cor))+xlab(x)+ylab("GSVA score")+theme(aspect.ratio = 1)
  })
  p <- cowplot::plot_grid(plotlist = lst3,ncol = 8,nrow = ceiling(length(lst3)/8),align = "v")
  fn <- paste0("/home/rlhua/project/cholangiocarcinoma/result/new/",clusters[i],"_markers_GSVAscore_correlation_dotplot.pdf")
  save_plot(fn,p,ncol = 8,nrow = ceiling(length(lst3)/8),base_width=2,base_height=2,limitsize = FALSE)
}
p1 <- cowplot::plot_grid(plotlist = lst,ncol = 8,nrow = 7,align = "v")
# save_plot("/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_integrated_assay_boxplot.pdf",p1,ncol = 8,nrow = 8,base_width=2,base_height=2)
# save_plot("/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_RNA_assay_boxplot.pdf",p1,ncol = 8,nrow = 8,base_width=2,base_height=2)
save_plot("/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_boxplot.pdf",p1,ncol = 8,nrow = 8,base_width=2,base_height=2)
p2 <- arrange_ggsurvplots(lst1, print = T,ncol = 8, nrow = 7)
# ggsave("/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_integrated_assay_survplot.pdf",p2,width = 24,height = 24)
# ggsave("/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_RNA_assay_survplot.pdf",p2,width = 24,height = 24)
ggsave("/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_survplot.pdf",p2,width = 24,height = 24)
dt <- rbindlist(lst2)
fwrite(dt,"/home/rlhua/project/cholangiocarcinoma/result/new/celltype_subcluster_markers2_correlation.tsv",sep="\t",quote=F)
