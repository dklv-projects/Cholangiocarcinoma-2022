library(scCancer)
library(ggplot2)
library(magrittr)
library(Seurat)
library(cowplot)
library(data.table)
library(survival)
library(survminer)

single.savePaths <- list.files("/home/rlhua/scRNApip/M1scCancer/anno","190065",full.names=T)
expr.list <- list()
for(i in 1:length(single.savePaths)){
  sampleName <- gsub(".*/","",single.savePaths[i])
  print(sampleName)
  expr.list[[sampleName]] <- readRDS(paste0(single.savePaths[i], "/expr.RDS"))
} # 读取单样本数据

npc=50
expr.anchors <- Seurat::FindIntegrationAnchors(object.list = expr.list,dims = 1:npc) # 找anchor
anchors <- expr.anchors@anchors
anchors$cellType1 <- "NULL"
anchors$cellType2 <- "NULL"
anchors$malignType1 <- "NULL"
anchors$malignType2 <- "NULL"
anchors$malignScore1 <- -1
anchors$malignScore2 <- -1
for(oi in expr.anchors@reference.objects){
  cur.ix <- which(anchors$dataset1 == oi)
  anchors$cellType1[cur.ix] <- expr.list[[oi]]@meta.data$Cell.Type[anchors$cell1[cur.ix]] # 给anchor添加细胞类型信息
  anchors$malignType1[cur.ix] <- expr.list[[oi]]@meta.data$Malign.type[anchors$cell1[cur.ix]] # 给anchor添加恶性程度信息
  anchors$malignScore1[cur.ix] <- expr.list[[oi]]@meta.data$Malign.score[anchors$cell1[cur.ix]]
  cur.ix <- which(anchors$dataset2 == oi)
  anchors$cellType2[cur.ix] <- expr.list[[oi]]@meta.data$Cell.Type[anchors$cell2[cur.ix]]
  anchors$malignType2[cur.ix] <- expr.list[[oi]]@meta.data$Malign.type[anchors$cell2[cur.ix]]
  anchors$malignScore2[cur.ix] <- expr.list[[oi]]@meta.data$Malign.score[anchors$cell2[cur.ix]]
}
anchors.new <- subset(anchors, cellType1 != "Epithelial" & cellType1 != "Unknown" & cellType2 != "Epithelial" & cellType2 != "Unknown") # 根据细胞类型筛选新的anchor
expr.anchors@anchors <- anchors.new # 使用新的anchor
genes <- c()
for (i in 1:length(expr.list)) {
  genes <- c(genes,rownames(expr.list[[i]]))
}
#genes <- table(genes)
#genes <- names(genes)[genes == length(expr.list)]
genes <- unique(genes)
expr <- Seurat::IntegrateData(anchorset = expr.anchors,dims = 1:npc,features.to.integrate =genes) # 整合数据
expr <- Seurat::ScaleData(expr,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")) # 中心化
expr <- Seurat::RunPCA(expr,npcs=100) # pca
expr <- Seurat::FindNeighbors(expr, reduction = "pca", dims = 1:npc)
expr <- Seurat::FindClusters(expr, resolution = 0.8) # 聚类
expr <- Seurat::RunTSNE(expr, dims = 1:npc)
expr <- Seurat::RunTSNE(expr, perplexity=60, max_iter=2000, stop_lying_iter=1000,dims = 1:npc) # tsne
expr <- Seurat::RunUMAP(expr, dims = 1:npc) # umap
saveRDS(expr,file = "/home/rlhua/project/cholangiocarcinoma/result/expr.RDS")

##### scrublet去除doublet
fs <- list.files("/home/rlhua/project/cholangiocarcinoma/result/scrublet","doublet.txt",full.names = T)
fs <- fs[-3:-4]
dt_scrublet <- lapply(1:8, function(x){
  dt <- fread(fs[x],sep = ",")
  dt$sample <- gsub(".*/|_doublet.txt","",fs[x])
  dt$barcode <- gsub("-1$",paste0("_",x),dt$barcode)
  return(dt)
}) %>% rbindlist()
dt <- merge(dt_scrublet[doublet_scores > threshold,.N,by=sample],dt_scrublet[,.N,by=sample],by="sample")
colnames(dt)[2:3] <- c("doublet","total")

dt_scrublet <- dt_scrublet[barcode %in% colnames(expr)]
dt <- merge(dt_scrublet[doublet_scores > threshold,.N,by=sample],dt_scrublet[,.N,by=sample],by="sample")
colnames(dt)[2:3] <- c("doublet","total")

dt_scrublet$predicted_doublets <- ifelse(dt_scrublet$doublet_scores >　dt_scrublet$threshold,T,F)

cellManifest <- fread("~/scRNApip/M1scCancer/stat/190065F_C6/cellManifest-all.txt",sep="\t")
cellManifest <- cellManifest[droplet.type=="cell"]
mito.thres <- min(grDevices::boxplot.stats(cellManifest$mito.percent[cellManifest$mito.percent < 0.75])$out)

fs <- list.files("/home/rlhua/project/cholangiocarcinoma/M2classification","RDS$",full.names = T)
for (i in fs) {
  cat(i,"\n")
  sr <- readRDS(i)
  cat(ncol(sr),"\n")
  sr <- sr[,colnames(sr) %in% dt_scrublet$barcode[!dt_scrublet$predicted_doublets]]
  cat(ncol(sr),"\n")
  sr <- sr[,!(sr$mito.percent > mito.thres & sr$orig.ident=="190065F_Ca6")]
  cat(ncol(sr),"\n")
  fn <- gsub(".RDS","1.RDS",i)
  saveRDS(sr,fn)
}

##### 各细胞类型亚类tsneplot
fs <- list.files("/home/rlhua/project/cholangiocarcinoma/M2classification","1.RDS$",full.names = T)
lst <- list()
for (i in fs) {
  cat(i,"\n")
  sr <- readRDS(i)
  colors <- scCancer::getDefaultColors(n =length(unique(sr$seurat_clusters)),type = 2)
  lst[[i]] <- DimPlot(sr, reduction = "tsne", group.by = "seurat_clusters",label = T,label.size=3,cols = colors)+
    ggtitle(gsub(".*/M2|1.RDS","",i))+
    theme(aspect.ratio = 1,plot.title = element_text(hjust = 0.5),legend.position = "none",text = element_text(size=7),axis.text = element_text(size=7))
}
p <- plot_grid(plotlist = lst,ncol = 4,nrow = 2,align = "hv",axis="tblr")
save_plot("/home/rlhua/project/cholangiocarcinoma/result/celltype_subcluster_tsneplot.pdf",p,ncol = 4,nrow =2,base_height = 3,base_width = 3,limitsize = FALSE)

##### 各细胞类型亚类markers生存分析
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
colnames(survival_info) <- c("sample","time","status","type")
survival_info$status <- ifelse(survival_info$status=="Alive",0,1)

## M3_clustermarker2
fs <- list.files("/home/rlhua/project/cholangiocarcinoma/M3markers","^M3_clustermarker2.RDS",all.files = T,full.names = T,recursive = T)
lst <- list()
for (i in 1:length(fs)) {
  df <- readRDS(fs[i])
  df[["celltype"]] <- gsub("/home/rlhua/project/cholangiocarcinoma/M3markers/M3markers_M2","",fs[i]) %>%
    gsub(".RDS/M3_clustermarker2.RDS","",.)
  lst[[i]] <- df
}
df <- rbindlist(lst)
df$cluster <- paste(df$celltype,df$cluster)
df$cluster <- gsub("_"," ",df$cluster)

## 绘图
clusters <- unique(df$cluster)
lst <- list()
lst1 <- list()
for (i in 1:length(clusters)) {
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
  
}
p1 <- cowplot::plot_grid(plotlist = lst,ncol = 8,nrow = 8,align = "v")
save_plot("/home/rlhua/project/cholangiocarcinoma/result/celltype_subcluster_markers2_boxplot.pdf",p1,ncol = 8,nrow = 8,base_width=2,base_height=2)
p2 <- arrange_ggsurvplots(lst1, print = T,ncol = 8, nrow = 8)
ggsave("/home/rlhua/project/cholangiocarcinoma/result/celltype_subcluster_markers2_survplot.pdf",p2,width = 24,height = 24)





