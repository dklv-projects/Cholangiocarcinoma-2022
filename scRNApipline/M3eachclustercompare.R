#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <seurat对象> <分类用meta.data列名>
#功能：计算所有类tumor和control的差异基因和signal2noise
#input：seurat对象
#input: 用于分类的列名，如integrated_snn_res.0.5
#output：M3_cluster_mean_expression_TumorOrNot.RDS
#output：M3_cluster_signal2noise_TumorOrNot.RDS
#output：M3_cluster_mean_expression_TumorOrNot.csv
#output：M3_diffgene_TumorOrNot.RDS
#output：M3_diffgene_TumorOrNot.csv

#保存在<包含seurat对象目录> 所在目录中的M3results目录中

#Rscript ../script/M3eachclustercompare.R ../M2classification/M2B_cells.RDS integrated_snn_res.0.5

Args<-commandArgs(T)
if(length(Args)!=2){stop("parameters number incorrect")}
workdir <- dirname(Args[1]) #项目目录
exprobj <- basename(Args[1])
setresolution <- Args[2] #项目目录

pwd=getwd()
setwd(workdir)
savePath=paste0(pwd,"/M3_TumorVsControl",exprobj)
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}

## 读取Seurat对象
library(Seurat)
expr_comb=readRDS(exprobj)
rewrite=FALSE
#注意：设置分类选项，需要根据需要选择！！！，建议跑完M2模块后修改M2模块
expr_comb@meta.data[,setresolution]=as.factor(expr_comb@meta.data[,setresolution])
if(!identical(expr_comb@meta.data$seurat_clusters,expr_comb@meta.data[,setresolution])){
  expr_comb@meta.data$seurat_clusters=expr_comb@meta.data[,setresolution]
  expr_comb@active.ident = expr_comb$seurat_clusters #设置默认分类
  rewrite=TRUE
}
if(rewrite==TRUE){
  saveRDS(expr_comb,file = exprobj)
}

## 第一步，计算每个细胞cluster的基因平均表达量和平均表达量最大的cluster
#input：seurat对象或表达谱矩阵，expr_comb
#output：基因平均表达量矩阵保存在cluster_mean_expression，平均表达量最大的类保存在maxcell,结果保存在M3_fig1cluster_mean_expression.RDS(csv)

#input:seurat对象，expr_comb
#input:细胞分类
cluster_mean_expression=data.frame()
cluster_signal2noise=data.frame()
diff_genes=data.frame()
meta=expr_comb@meta.data

matexpr=as.matrix(expr_comb@assays$integrated@data)
tissue="TUMOR"
allcells=rownames(meta)[meta$tissues == tissue]
allothercells=rownames(meta)[meta$tissues != tissue]


diff.expr.genes = tryCatch(
  {
    Seurat::FindMarkers(matexpr,cells.1=allcells,cells.2=allothercells,only.pos = T,logfc.threshold = 0.25)
    #diff.expr.genes = diff.expr.genes[diff.expr.genes$p_val_adj<0.05,]
    diff.expr.genes$cluster="all"
  },
    error=function(e){NA})
if(!is.na(diff.expr.genes)){
  diff_genes=diff.expr.genes
}

exprmat=matexpr[,allcells]
otherexprmat=matexpr[,allothercells]
rowmean=rowMeans(exprmat)
otherrowmean=rowMeans(otherexprmat)
sds=apply(exprmat,1,sd)
othersds=apply(otherexprmat,1,sd)
signal2noise=(rowmean-otherrowmean)/(sds+othersds)
signal2noisedf=data.frame(value=signal2noise)
colnames(signal2noisedf)="all"
rowmeandf=data.frame(value1=rowmean,value2=otherrowmean)
colnames(rowmeandf)=c("allTUMOR","allnotTUMOR")
cluster_signal2noise=signal2noisedf
cluster_mean_expression=rowmeandf

for(clust in sort(unique(meta$seurat_clusters))){
  print(clust)
  
    tissue_cells=rownames(meta)[meta$seurat_clusters==clust & meta$tissues == tissue]
    tissue_othercells=rownames(meta)[meta$seurat_clusters==clust & meta$tissues != tissue]

  print(length(tissue_cells))
  print(length(tissue_othercells))
  if(length(tissue_cells)<3 | length(tissue_othercells)<3){next}

  exprmat=matexpr[,tissue_cells]
  otherexprmat=matexpr[,tissue_othercells]
  rowmean=rowMeans(exprmat)
  otherrowmean=rowMeans(otherexprmat)
  sds=apply(exprmat,1,sd)
  othersds=apply(otherexprmat,1,sd)
  signal2noise=(rowmean-otherrowmean)/(sds+othersds)
  signal2noisedf=data.frame(value=signal2noise)
  
  colnames(signal2noisedf)=clust
  rowmeandf=data.frame(value1=rowmean,value2=otherrowmean)
  colnames(rowmeandf)=c(paste0(clust,c("TUMOR","notTUMOR")))
  if(ncol(cluster_signal2noise)==0){
    cluster_signal2noise=signal2noisedf
  }else{
    cluster_signal2noise=cbind(cluster_signal2noise,signal2noisedf)    
  }
  
  if(ncol(cluster_mean_expression)==0){
    cluster_mean_expression=rowmeandf
  }else{
    cluster_mean_expression=cbind(cluster_mean_expression,rowmeandf)    
  }

  diff.expr.genes <- Seurat::FindMarkers(matexpr,cells.1=tissue_cells,cells.2=tissue_othercells,only.pos = T,logfc.threshold = 0.25)
  #diff.expr.genes = diff.expr.genes[diff.expr.genes$p_val_adj<0.05,]
  if(nrow(diff.expr.genes)>0){
    diff.expr.genes$cluster=clust
    diff_genes=rbind(diff_genes,diff.expr.genes)  
  }
}

saveRDS(cluster_mean_expression, file = file.path(savePath,"M3_cluster_mean_expression_TumorOrNot.RDS"))
saveRDS(cluster_signal2noise, file = file.path(savePath,"M3_cluster_signal2noise_TumorOrNot.RDS"))
write.csv(cluster_mean_expression,file = file.path(savePath,"M3_cluster_mean_expression_TumorOrNot.csv"))

saveRDS(diff_genes, file = file.path(savePath,"M3_diffgene_TumorOrNot.RDS"))
write.csv(diff_genes, file = file.path(savePath,"M3_diffgene_TumorOrNot.csv"))













