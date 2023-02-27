#manually curate clusters with the cell annotation and marker distribution
#根据细胞类型注释和常见细胞marker的表达分布手动确认细胞类型，保存到
library(Seurat)
library(paran)
expr=readRDS("/home/rlhua/project/cholangiocarcinoma/result/new/expr.RDS")
setwd("/home/rlhua/project/cholangiocarcinoma/result/new/")

#细胞类型已经保存在celltype中
for(cellkind in sort(unique(expr$celltype))){
  #取一大类所有细胞进行亚分类
  aclass <- expr[,expr$celltype == cellkind]
  #DefaultAssay(aclass)="RNA"
  aclass <- FindVariableFeatures(aclass, selection.method = "vst", nfeatures = 1000, verbose = F)
  aclass <- ScaleData(aclass,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"))
  mtx <- as.matrix(aclass@assays$integrated@scale.data)
  #library(paran)
  pr <- paran(t(mtx),iterations = 10)
  npca <- sum(pr$Ev/pr$RndEv > 1.5)
  aclass$npca=npca
  aclass <- Seurat::RunPCA(aclass,npcs=100) # pca
  #ElbowPlot(aclass,ndims = 100)
  aclass <- Seurat::FindNeighbors(aclass, reduction = "pca", dims = 1:npca)
  #删除在大类中的分类结果
  ccols=grep("integrated_snn_res",colnames(aclass@meta.data))
  print(ccols)
  aclass@meta.data=aclass@meta.data[,-ccols]
  #重新再分亚类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.02) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.05) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.1) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.2) # 聚类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.5) # 聚类
  aclass <- Seurat::RunTSNE(aclass, dims = 1:npca)
  #aclass <- Seurat::RunTSNE(aclass, perplexity=90, max_iter=2000, stop_lying_iter=1000,dims = 1:npca)
  aclass <- Seurat::RunUMAP(aclass, dims = 1:npca) # umap
  saveRDS(aclass,file = paste0(sub(" ","_",cellkind),".RDS"))
}