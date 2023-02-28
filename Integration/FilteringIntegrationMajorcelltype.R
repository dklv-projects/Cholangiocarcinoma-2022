library(scCancer)
library(ggplot2)
library(magrittr)
library(Seurat)
library(cowplot)
library(data.table)
library(survival)
library(survminer)


##### 读取每个样本sccancer seurat对象
single.savePaths <- list.files("/home/rlhua/scRNApip/M1scCancer/anno","190065",full.names=T)
expr.list <- list()
for(i in 1:length(single.savePaths)){
  sampleName <- gsub(".*/","",single.savePaths[i])
  print(sampleName)
  expr.list[[sampleName]] <- readRDS(paste0(single.savePaths[i], "/expr.RDS"))
} # 读取单样本数据

##### 去除190065F_C6样本线粒体基因比例高的细胞
cellManifest <- fread("~/scRNApip/M1scCancer/stat/190065F_C6/cellManifest-all.txt",sep="\t")
cellManifest <- cellManifest[droplet.type=="cell"]
mito.thres <- min(grDevices::boxplot.stats(cellManifest$mito.percent[cellManifest$mito.percent < 0.75])$out)
expr.list[["190065F_C6"]] <- expr.list[["190065F_C6"]][,expr.list[["190065F_C6"]]@meta.data$mito.percent < mito.thres]

##### scrublet去除doublet
fs <- list.files("/home/rlhua/project/cholangiocarcinoma/result/scrublet","doublet.txt",full.names = T)
fs <- fs[-3:-4]
for(i in 1:length(fs)){
  dt <- fread(fs[i],sep = ",")
  dt$barcode <- gsub("-1$","",dt$barcode)
  cells_singlet <- dt$barcode[dt$doublet_scores < dt$threshold]
  expr.list[[i]] <- expr.list[[i]][,colnames(expr.list[[i]]) %in% cells_singlet]
}

##### 多样本整合
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
genes <- unique(genes)
expr <- Seurat::IntegrateData(anchorset = expr.anchors,dims = 1:npc,features.to.integrate =genes) # 整合数据
expr <- Seurat::ScaleData(expr,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent")) # 归一化
expr <- Seurat::RunPCA(expr,npcs=100) # pca
expr <- Seurat::FindNeighbors(expr, reduction = "pca", dims = 1:npc)
expr <- Seurat::FindClusters(expr, resolution = 0.8)
expr <- Seurat::FindClusters(expr, resolution = 0.1)# 聚类
expr <- Seurat::RunTSNE(expr, perplexity=60, max_iter=2000, stop_lying_iter=1000,dims = 1:npc) # tsne
expr <- Seurat::RunUMAP(expr, dims = 1:npc) # umap

##### 细胞类型注释
expr$celltype <- ""
expr$celltype[expr$seurat_clusters %in% c(0,1,13)] <- "T cells"
expr$celltype[expr$seurat_clusters %in% c(2,10)] <- "B cells"
expr$celltype[expr$seurat_clusters %in% c(3)] <- "Myeloid"
expr$celltype[expr$seurat_clusters %in% c(4)] <- "Fibroblasts"
expr$celltype[expr$seurat_clusters %in% c(5)] <- "Endothelial"
expr$celltype[expr$seurat_clusters %in% c(6)] <- "Mast cells"
expr$celltype[expr$seurat_clusters %in% c(11)] <- "Schwann"
expr$celltype[expr$seurat_clusters %in% c(7,8,9,12)] <- "Epithelial"
expr$celltype[expr$integrated_snn_res.0.8==15] <- "Epithelial"

##### 保存整合后seurat对象
saveRDS(expr,file = "/home/rlhua/project/cholangiocarcinoma/result/new/expr.RDS")
