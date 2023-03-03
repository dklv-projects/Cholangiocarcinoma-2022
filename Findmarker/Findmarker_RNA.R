#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <seurat对象> <分类用meta.data列名>
#功能：细胞类型降维展示
#input：seurat对象
#input: 用于分类的列名，如integrated_snn_res.0.2
#input: FindMarkers中min.pct设置，如0.1
#output：包括2个pdf文件和6个数据文件
#保存在<包含seurat对象目录> 所在目录中的M3results目录中

#Rscript ../script/M3findmarkers.R ../M2classification/M2B_cells.RDS integrated_snn_res.0.5 0.1
#Rscript /home/dklv/scPip/script/M3findmarkers.R M1scCancer/comb/expr_combined.RDS integrated_snn_res.0.2 0.25

Args<-commandArgs(T)
if(length(Args)!=3){stop("parameters number incorrect")}
workdir <- dirname(Args[1]) #项目目录
exprobj <- basename(Args[1])
setresolution <- Args[2]
setpct <- as.numeric(Args[3])
pwd=getwd()
setwd(workdir)
savePath=paste0(pwd,"/RNAmarkers_",gsub("\\..*","",exprobj))
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}

## 读取Seurat对象
library(Seurat)
expr_comb=readRDS(exprobj)
DefaultAssay(expr_comb)="RNA"
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
cluster_mean_expression_org=data.frame()
cluster_signal2noise=data.frame()
meta=expr_comb@meta.data
matexpr=as.matrix(expr_comb@assays$RNA@data)
for(clust in sort(unique(meta$seurat_clusters))){
  print(clust)
  cells=rownames(meta)[meta$seurat_clusters==clust]
  othercells=rownames(meta)[meta$seurat_clusters!=clust]
  print(length(cells))
  print(length(othercells))
  exprmat=matexpr[,cells]
  otherexprmat=matexpr[,othercells]
  rowmean=rowMeans(exprmat)
  otherrowmean=rowMeans(otherexprmat)
  sds=apply(exprmat,1,sd)
  othersds=apply(otherexprmat,1,sd)
  signal2noise=(rowmean-otherrowmean)/(sds+othersds)
  signal2noise[is.na(signal2noise)]=0
  signal2noisedf=data.frame(value=signal2noise)
  colnames(signal2noisedf)=clust
  rowmeandf=data.frame(value1=rowmean,value2=otherrowmean)
  colnames(rowmeandf)=c(clust,paste0("not_",clust))
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
  
  for(sample in unique(meta$orig.ident)){
    print(sample)
    cells.org=rownames(meta)[meta$seurat_clusters==clust & meta$orig.ident == sample]
    if(length(cells.org)>1){
      exprmat.org=matexpr[,cells.org]
      rowmeandf.org=data.frame(value=rowMeans(exprmat.org))
      colnames(rowmeandf.org)=paste(sample,clust,sep="_")
    }else if(length(cells.org)==1){
      exprmat.org=matexpr[,cells.org]
      rowmeandf.org=data.frame(value=exprmat.org)
      colnames(rowmeandf.org)=paste(sample,clust,sep="_")
    }else{next}
    if(ncol(cluster_mean_expression_org)==0){
      cluster_mean_expression_org=rowmeandf.org
    }else{
      cluster_mean_expression_org=cbind(cluster_mean_expression_org,rowmeandf.org)    
    }
  }
}

maxcell=apply(cluster_mean_expression,1,function(x){
  colnames(cluster_mean_expression)[order(x,decreasing = T)[1]]
})
#cluster_mean_expression$maxcell=maxcell
saveRDS(cluster_mean_expression, file = file.path(savePath,"cluster_mean_expression.RDS"))

saveRDS(cluster_mean_expression_org, file = file.path(savePath,"cluster_mean_expression_org.RDS"))
saveRDS(cluster_signal2noise, file = file.path(savePath,"cluster_signal2noise.RDS"))


write.csv(cluster_mean_expression,file = file.path(savePath,"cluster_mean_expression.csv"))




## 第二步，筛选每个cluster的marker基因
#input：seurat对象，expr_comb
#output：每个cluster的marker基因矩阵，diff.expr.genes和diff.expr.genes2，结果保存到clustermarker.RDS和clustermarker2.RDS文件

#筛选方法1A，宽松
expr_comb@active.ident = expr_comb$seurat_clusters #设置默认分类
# 找每个类的marker基因,only.pos是否只找阳性marker，logfc.threshold设置倍数阈值
diff.expr.genes <- Seurat::FindAllMarkers(expr_comb,only.pos = T,logfc.threshold = 0.25, min.pct = setpct) 
#进一步筛选：marker基因p_val_adj<0.05,所在类是其表达水平最高的类
diff.expr.genes = diff.expr.genes[diff.expr.genes$p_val_adj<0.05,]
#diff.expr.genes =diff.expr.genes[diff.expr.genes$avg_logFC>=0.58 &diff.expr.genes$p_val_adj<0.05,]
#进一步筛选
diff.expr.genes$maxcell=maxcell[diff.expr.genes$gene]
diff.expr.genes = diff.expr.genes[as.character(diff.expr.genes$cluster)==diff.expr.genes$maxcell,]
saveRDS(diff.expr.genes,file = file.path(savePath,"clustermarker.RDS"))
write.csv(diff.expr.genes,file = file.path(savePath,"clustermarker.csv"))

#筛选方法1B，严谨
diff.expr.genes2 <- Seurat::FindAllMarkers(expr_comb,only.pos = T,logfc.threshold = 0.25,test.use = "roc", min.pct = setpct)
diff.expr.genes2$maxcell=maxcell[diff.expr.genes2$gene]
diff.expr.genes2 = diff.expr.genes2[as.character(diff.expr.genes2$cluster)==diff.expr.genes2$maxcell,]
diff.expr.genes2 = diff.expr.genes2[diff.expr.genes2$myAUC>0.5,]
diff.expr.genes2=diff.expr.genes2[ order(diff.expr.genes2$myAUC,decreasing=T),]
diff.expr.genes2=diff.expr.genes2[ !duplicated(diff.expr.genes2$gene),]
diff.expr.genes2=diff.expr.genes2[order(diff.expr.genes2$cluster),]


saveRDS(diff.expr.genes2,file = file.path(savePath,"clustermarker2.RDS"))
write.csv(diff.expr.genes2,file = file.path(savePath,"clustermarker2.csv"))

## 绘图1，heatmap，保存到M3_top5mark1.pdf,M3_top5mark2.pdf
#input：marker矩阵diff.expr.genes，基因平均表达量矩阵cluster_mean_expression
topdf<-function(diff.expr.genes,retain=10,avg_logFC=T){
  top10df=data.frame()
  for(clust in unique(diff.expr.genes$cluster)){
    subdf=diff.expr.genes[diff.expr.genes$cluster==clust,]
    if(nrow(subdf)>0){
      number=ifelse(nrow(subdf)>retain,retain,nrow(subdf))
      if(avg_logFC){
        if("avg_logFC" %in% colnames(subdf)){
          subdf=subdf[order(subdf$avg_logFC,decreasing = T),]#按照avg_logFC重新排序，可选，如果是AUC筛选的矩阵，此句省略
        }
        #print(dim(subdf))
      }
      top10df=rbind(top10df,head(subdf,number))
    }
  }
  top10df
}

library(pheatmap)
library(viridis)
top10df=topdf(diff.expr.genes,retain = 5)
top10expr=cluster_mean_expression[top10df$gene,]
pheatmap(top10expr,
         color=viridis(100),
         border_color=NA,
         scale="row",
         cluster_rows=F,cluster_cols=F,
         width=4,
         height = 15,
         filename=file.path(savePath,"top5mark1.pdf"))

top10df=topdf(diff.expr.genes2,retain = 5)
top10expr=cluster_mean_expression[top10df$gene,]
pheatmap(top10expr,
         color=viridis(100),
         border_color=NA,
         scale="row",
         cluster_rows=F,cluster_cols=F,
         width=4,
         height = 15,
         filename=file.path(savePath,"top5mark2.pdf"))

#################################################################
#2绘制所有类marker的热图
rowcluster<-function(df,Ind=F){
  #对数据的行进行层次聚类，输出聚类后的数据
  hc<-hclust(dist(df),method = "ward.D2") #对行进行聚类
  rowInd<-hc$order
  if(Ind==T){
    return(rowInd)
  }else{
    df[rowInd,]
  }
}
colcluster<-function(df,Ind=F){
  #对数据的列进行层次聚类，输出聚类后的数据
  hc<-hclust(dist(t(df)),method = "ward.D2")#对列进行聚类
  colInd<-hc$order 
  if(Ind==T){
    return(colInd)
  }else{
    df[,colInd]
  }
}
mycolor=c('#352a86', '#342b89', '#342d8c', '#342e8f', '#343092', '#343196', '#34339a', 
          '#34349d', '#3436a0', '#3438a4', '#3439a7', '#343baa', '#333dad', '#2f3fb2', 
          '#2b42b6', '#2745ba', '#2349be', '#1f4cc2', '#1b4fc6', '#1752ca', '#1355ce', 
          '#0f58d3', '#0a5bd7', '#065edb', '#0361de', '#0364de', '#0568dd', '#066cdc', 
          '#076fdb', '#0972da', '#0a75d8', '#0b78d7', '#0d7cd6', '#0e7fd5', '#1082d4',
          '#1185d3', '#1289d2', '#148ccd', '#1790c9', '#1a94c6', '#1c98c2', '#1e9cbe', 
          '#209fba', '#22a3b6', '#24a7b3', '#26abae', '#28afaa', '#2ab3a6', '#2eb6a2', 
          '#37b79e', '#41b899', '#4bb894', '#54b990', '#5fba8b', '#69ba86', '#73bb81', 
          '#7dbb7c', '#87bc78', '#90bc73', '#9abd6f', '#a4bd6a', '#abbd66', '#b3bd63',
          '#b9bc60', '#c0bc5c', '#c7bc59', '#cebc56', '#d5bb53', '#dcbb50', '#e3bb4d', 
          '#e9ba49', '#f1ba45', '#f6ba43', '#f7bc40', '#f7bf3d', '#f7c23a', '#f7c438',
          '#f7c735', '#f7c933', '#f6cc30', '#f6cf2d', '#f6d22a', '#f6d428', '#f6d725',
          '#f6d923', '#f6dc21', '#f6df1f', '#f6e11d', '#f6e41b', '#f6e719', '#f7e917',
          '#f7ed15', '#f7ef14', '#f7f212', '#f7f410', '#f7f70e', '#f8fa0d')

cluster_mean_expression=cluster_mean_expression[,seq(1, ncol(cluster_mean_expression), by = 2)]
expr_mean_mark=cluster_mean_expression[diff.expr.genes$gene,]

tmpdf=data.frame()
ranknames=colnames(expr_mean_mark)[colcluster(expr_mean_mark,Ind = T)]
tabs=table(diff.expr.genes$cluster)[ranknames]
rowgap=cumsum(tabs[1:(length(tabs)-1)])
for(clust in ranknames){
  print(clust)
  clustdf=expr_mean_mark[diff.expr.genes$cluster==clust,]
  print(nrow(clustdf))
  if(nrow(clustdf)==0){next}
  if(nrow(clustdf)>1){
    clustdf=rowcluster(clustdf)
  }
  if(nrow(tmpdf)==0){
    tmpdf=clustdf
  }else{
    tmpdf=rbind(tmpdf,clustdf)
  }
}
tmpdf=colcluster(tmpdf)
annotation_col=data.frame(Cluster=factor(colnames(tmpdf)))
rownames(annotation_col)=annotation_col$Cluster
pheatmap(tmpdf,scale = "row",
         cluster_rows = F,cluster_cols = F,show_rownames=F,
         color = mycolor,
         gaps_col = 1:ncol(tmpdf),gaps_row = rowgap,
         annotation_col = annotation_col,
         filename = file.path(savePath,"mark.pdf"))

expr_mean_mark=cluster_mean_expression[diff.expr.genes2$gene,]

tmpdf=data.frame()
ranknames=colnames(expr_mean_mark)[colcluster(expr_mean_mark,Ind = T)]
tabs=table(diff.expr.genes2$cluster)[ranknames]
rowgap=cumsum(tabs[1:(length(tabs)-1)])
for(clust in ranknames){
  print(clust)
  clustdf=expr_mean_mark[diff.expr.genes2$cluster==clust,]
  print(nrow(clustdf))
  if(nrow(clustdf)==0){next}
  if(nrow(clustdf)>1){
    clustdf=rowcluster(clustdf)
  }
  if(nrow(tmpdf)==0){
    tmpdf=clustdf
  }else{
    tmpdf=rbind(tmpdf,clustdf)
  }
}
tmpdf=colcluster(tmpdf)
annotation_col=data.frame(Cluster=factor(colnames(tmpdf)))
rownames(annotation_col)=annotation_col$Cluster
pheatmap(tmpdf,scale = "row",
         cluster_rows = F,cluster_cols = F,show_rownames=F,
         color = mycolor,
         gaps_col = 1:ncol(tmpdf),gaps_row = rowgap,
         annotation_col = annotation_col,
         filename = file.path(savePath,"mark2.pdf"))


#USAGE:Rscript cmdfile <seurat对象> <分类用meta.data列名>
#功能：细胞类型降维展示
#input：seurat对象
#input: 用于分类的列名，如integrated_snn_res.0.2
#input: FindMarkers中min.pct设置，如0.1
#output：包括2个pdf文件和6个数据文件
#保存在<包含seurat对象目录> 所在目录中的M3results目录中

#Rscript ../script/M3findmarkers.R ../M2classification/M2B_cells.RDS integrated_snn_res.0.5 0.1
#Rscript /home/dklv/scPip/script/M3findmarkers.R M1scCancer/comb/expr_combined.RDS integrated_snn_res.0.2 0.25

Args<-commandArgs(T)
if(length(Args)!=3){stop("parameters number incorrect")}
workdir <- dirname(Args[1]) #项目目录
exprobj <- basename(Args[1])
setresolution <- Args[2]
setpct <- as.numeric(Args[3])
pwd=getwd()
setwd(workdir)
savePath=paste0(pwd,"/RNAmarkers_",gsub("\\..*","",exprobj))
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}

## 读取Seurat对象
library(Seurat)
expr_comb=readRDS(exprobj)
DefaultAssay(expr_comb)="RNA"
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
cluster_mean_expression_org=data.frame()
cluster_signal2noise=data.frame()
meta=expr_comb@meta.data
matexpr=as.matrix(expr_comb@assays$RNA@data)
for(clust in sort(unique(meta$seurat_clusters))){
  print(clust)
  cells=rownames(meta)[meta$seurat_clusters==clust]
  othercells=rownames(meta)[meta$seurat_clusters!=clust]
  print(length(cells))
  print(length(othercells))
  exprmat=matexpr[,cells]
  otherexprmat=matexpr[,othercells]
  rowmean=rowMeans(exprmat)
  otherrowmean=rowMeans(otherexprmat)
  sds=apply(exprmat,1,sd)
  othersds=apply(otherexprmat,1,sd)
  signal2noise=(rowmean-otherrowmean)/(sds+othersds)
  signal2noise[is.na(signal2noise)]=0
  signal2noisedf=data.frame(value=signal2noise)
  colnames(signal2noisedf)=clust
  rowmeandf=data.frame(value1=rowmean,value2=otherrowmean)
  colnames(rowmeandf)=c(clust,paste0("not_",clust))
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
  
  for(sample in unique(meta$orig.ident)){
    print(sample)
    cells.org=rownames(meta)[meta$seurat_clusters==clust & meta$orig.ident == sample]
    if(length(cells.org)>1){
      exprmat.org=matexpr[,cells.org]
      rowmeandf.org=data.frame(value=rowMeans(exprmat.org))
      colnames(rowmeandf.org)=paste(sample,clust,sep="_")
    }else if(length(cells.org)==1){
      exprmat.org=matexpr[,cells.org]
      rowmeandf.org=data.frame(value=exprmat.org)
      colnames(rowmeandf.org)=paste(sample,clust,sep="_")
    }else{next}
    if(ncol(cluster_mean_expression_org)==0){
      cluster_mean_expression_org=rowmeandf.org
    }else{
      cluster_mean_expression_org=cbind(cluster_mean_expression_org,rowmeandf.org)    
    }
  }
}

maxcell=apply(cluster_mean_expression,1,function(x){
  colnames(cluster_mean_expression)[order(x,decreasing = T)[1]]
})
#cluster_mean_expression$maxcell=maxcell
saveRDS(cluster_mean_expression, file = file.path(savePath,"cluster_mean_expression.RDS"))

saveRDS(cluster_mean_expression_org, file = file.path(savePath,"cluster_mean_expression_org.RDS"))
saveRDS(cluster_signal2noise, file = file.path(savePath,"cluster_signal2noise.RDS"))


write.csv(cluster_mean_expression,file = file.path(savePath,"cluster_mean_expression.csv"))




## 第二步，筛选每个cluster的marker基因
#input：seurat对象，expr_comb
#output：每个cluster的marker基因矩阵，diff.expr.genes和diff.expr.genes2，结果保存到clustermarker.RDS和clustermarker2.RDS文件

#筛选方法1A，宽松
expr_comb@active.ident = expr_comb$seurat_clusters #设置默认分类
# 找每个类的marker基因,only.pos是否只找阳性marker，logfc.threshold设置倍数阈值
diff.expr.genes <- Seurat::FindAllMarkers(expr_comb,only.pos = T,logfc.threshold = 0.25, min.pct = setpct) 
#进一步筛选：marker基因p_val_adj<0.05,所在类是其表达水平最高的类
diff.expr.genes = diff.expr.genes[diff.expr.genes$p_val_adj<0.05,]
#diff.expr.genes =diff.expr.genes[diff.expr.genes$avg_logFC>=0.58 &diff.expr.genes$p_val_adj<0.05,]
#进一步筛选
diff.expr.genes$maxcell=maxcell[diff.expr.genes$gene]
diff.expr.genes = diff.expr.genes[as.character(diff.expr.genes$cluster)==diff.expr.genes$maxcell,]
saveRDS(diff.expr.genes,file = file.path(savePath,"clustermarker.RDS"))
write.csv(diff.expr.genes,file = file.path(savePath,"clustermarker.csv"))

#筛选方法1B，严谨
diff.expr.genes2 <- Seurat::FindAllMarkers(expr_comb,only.pos = T,logfc.threshold = 0.25,test.use = "roc", min.pct = setpct)
diff.expr.genes2$maxcell=maxcell[diff.expr.genes2$gene]
diff.expr.genes2 = diff.expr.genes2[as.character(diff.expr.genes2$cluster)==diff.expr.genes2$maxcell,]
diff.expr.genes2 = diff.expr.genes2[diff.expr.genes2$myAUC>0.5,]
diff.expr.genes2=diff.expr.genes2[ order(diff.expr.genes2$myAUC,decreasing=T),]
diff.expr.genes2=diff.expr.genes2[ !duplicated(diff.expr.genes2$gene),]
diff.expr.genes2=diff.expr.genes2[order(diff.expr.genes2$cluster),]


saveRDS(diff.expr.genes2,file = file.path(savePath,"clustermarker2.RDS"))
write.csv(diff.expr.genes2,file = file.path(savePath,"clustermarker2.csv"))

## 绘图1，heatmap，保存到M3_top5mark1.pdf,M3_top5mark2.pdf
#input：marker矩阵diff.expr.genes，基因平均表达量矩阵cluster_mean_expression
topdf<-function(diff.expr.genes,retain=10,avg_logFC=T){
  top10df=data.frame()
  for(clust in unique(diff.expr.genes$cluster)){
    subdf=diff.expr.genes[diff.expr.genes$cluster==clust,]
    if(nrow(subdf)>0){
      number=ifelse(nrow(subdf)>retain,retain,nrow(subdf))
      if(avg_logFC){
        if("avg_logFC" %in% colnames(subdf)){
          subdf=subdf[order(subdf$avg_logFC,decreasing = T),]#按照avg_logFC重新排序，可选，如果是AUC筛选的矩阵，此句省略
        }
        #print(dim(subdf))
      }
      top10df=rbind(top10df,head(subdf,number))
    }
  }
  top10df
}

library(pheatmap)
library(viridis)
top10df=topdf(diff.expr.genes,retain = 5)
top10expr=cluster_mean_expression[top10df$gene,]
pheatmap(top10expr,
         color=viridis(100),
         border_color=NA,
         scale="row",
         cluster_rows=F,cluster_cols=F,
         width=4,
         height = 15,
         filename=file.path(savePath,"top5mark1.pdf"))

top10df=topdf(diff.expr.genes2,retain = 5)
top10expr=cluster_mean_expression[top10df$gene,]
pheatmap(top10expr,
         color=viridis(100),
         border_color=NA,
         scale="row",
         cluster_rows=F,cluster_cols=F,
         width=4,
         height = 15,
         filename=file.path(savePath,"top5mark2.pdf"))

#################################################################
#2绘制所有类marker的热图
rowcluster<-function(df,Ind=F){
  #对数据的行进行层次聚类，输出聚类后的数据
  hc<-hclust(dist(df),method = "ward.D2") #对行进行聚类
  rowInd<-hc$order
  if(Ind==T){
    return(rowInd)
  }else{
    df[rowInd,]
  }
}
colcluster<-function(df,Ind=F){
  #对数据的列进行层次聚类，输出聚类后的数据
  hc<-hclust(dist(t(df)),method = "ward.D2")#对列进行聚类
  colInd<-hc$order 
  if(Ind==T){
    return(colInd)
  }else{
    df[,colInd]
  }
}
mycolor=c('#352a86', '#342b89', '#342d8c', '#342e8f', '#343092', '#343196', '#34339a', 
          '#34349d', '#3436a0', '#3438a4', '#3439a7', '#343baa', '#333dad', '#2f3fb2', 
          '#2b42b6', '#2745ba', '#2349be', '#1f4cc2', '#1b4fc6', '#1752ca', '#1355ce', 
          '#0f58d3', '#0a5bd7', '#065edb', '#0361de', '#0364de', '#0568dd', '#066cdc', 
          '#076fdb', '#0972da', '#0a75d8', '#0b78d7', '#0d7cd6', '#0e7fd5', '#1082d4',
          '#1185d3', '#1289d2', '#148ccd', '#1790c9', '#1a94c6', '#1c98c2', '#1e9cbe', 
          '#209fba', '#22a3b6', '#24a7b3', '#26abae', '#28afaa', '#2ab3a6', '#2eb6a2', 
          '#37b79e', '#41b899', '#4bb894', '#54b990', '#5fba8b', '#69ba86', '#73bb81', 
          '#7dbb7c', '#87bc78', '#90bc73', '#9abd6f', '#a4bd6a', '#abbd66', '#b3bd63',
          '#b9bc60', '#c0bc5c', '#c7bc59', '#cebc56', '#d5bb53', '#dcbb50', '#e3bb4d', 
          '#e9ba49', '#f1ba45', '#f6ba43', '#f7bc40', '#f7bf3d', '#f7c23a', '#f7c438',
          '#f7c735', '#f7c933', '#f6cc30', '#f6cf2d', '#f6d22a', '#f6d428', '#f6d725',
          '#f6d923', '#f6dc21', '#f6df1f', '#f6e11d', '#f6e41b', '#f6e719', '#f7e917',
          '#f7ed15', '#f7ef14', '#f7f212', '#f7f410', '#f7f70e', '#f8fa0d')

cluster_mean_expression=cluster_mean_expression[,seq(1, ncol(cluster_mean_expression), by = 2)]
expr_mean_mark=cluster_mean_expression[diff.expr.genes$gene,]

tmpdf=data.frame()
ranknames=colnames(expr_mean_mark)[colcluster(expr_mean_mark,Ind = T)]
tabs=table(diff.expr.genes$cluster)[ranknames]
rowgap=cumsum(tabs[1:(length(tabs)-1)])
for(clust in ranknames){
  print(clust)
  clustdf=expr_mean_mark[diff.expr.genes$cluster==clust,]
  print(nrow(clustdf))
  if(nrow(clustdf)==0){next}
  if(nrow(clustdf)>1){
    clustdf=rowcluster(clustdf)
  }
  if(nrow(tmpdf)==0){
    tmpdf=clustdf
  }else{
    tmpdf=rbind(tmpdf,clustdf)
  }
}
tmpdf=colcluster(tmpdf)
annotation_col=data.frame(Cluster=factor(colnames(tmpdf)))
rownames(annotation_col)=annotation_col$Cluster
pheatmap(tmpdf,scale = "row",
         cluster_rows = F,cluster_cols = F,show_rownames=F,
         color = mycolor,
         gaps_col = 1:ncol(tmpdf),gaps_row = rowgap,
         annotation_col = annotation_col,
         filename = file.path(savePath,"mark.pdf"))

expr_mean_mark=cluster_mean_expression[diff.expr.genes2$gene,]

tmpdf=data.frame()
ranknames=colnames(expr_mean_mark)[colcluster(expr_mean_mark,Ind = T)]
tabs=table(diff.expr.genes2$cluster)[ranknames]
rowgap=cumsum(tabs[1:(length(tabs)-1)])
for(clust in ranknames){
  print(clust)
  clustdf=expr_mean_mark[diff.expr.genes2$cluster==clust,]
  print(nrow(clustdf))
  if(nrow(clustdf)==0){next}
  if(nrow(clustdf)>1){
    clustdf=rowcluster(clustdf)
  }
  if(nrow(tmpdf)==0){
    tmpdf=clustdf
  }else{
    tmpdf=rbind(tmpdf,clustdf)
  }
}
tmpdf=colcluster(tmpdf)
annotation_col=data.frame(Cluster=factor(colnames(tmpdf)))
rownames(annotation_col)=annotation_col$Cluster
pheatmap(tmpdf,scale = "row",
         cluster_rows = F,cluster_cols = F,show_rownames=F,
         color = mycolor,
         gaps_col = 1:ncol(tmpdf),gaps_row = rowgap,
         annotation_col = annotation_col,
         filename = file.path(savePath,"mark2.pdf"))


