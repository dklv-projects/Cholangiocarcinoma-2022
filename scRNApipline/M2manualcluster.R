#manually curate clusters with the cell annotation and marker distribution
#根据细胞类型注释和常见细胞marker的表达分布手动确认细胞类型，保存到
library(Seurat)
library(data.table)
setwd("/home/dklv/scPip/M1scCancer/comb/")
#可以自制一张表格，列名为meta.data中的列名，类名，新的自定义类名
#读取seurat对象并把手工分类信息写入到meta.data中，列名为manual_clusters1
expr_comb <- readRDS("expr_combined.RDS")
dt1 <- data.table(cell=colnames(expr_comb),clusters.res.0.1=expr_comb$integrated_snn_res.0.1)
dt1$myCellType <- 0
dt1$myCellType[dt1$clusters.res.0.1==0] <- "T cells"
dt1$myCellType[dt1$clusters.res.0.1==1] <- "T cells"
dt1$myCellType[dt1$clusters.res.0.1==2] <- "B cells"
dt1$myCellType[dt1$clusters.res.0.1==3] <- "Epithelial"
dt1$myCellType[dt1$clusters.res.0.1==4] <- "Fibroblasts"
dt1$myCellType[dt1$clusters.res.0.1==5] <- "Myeloid"
dt1$myCellType[dt1$clusters.res.0.1==6] <- "Endothelial"
dt1$myCellType[dt1$clusters.res.0.1==7] <- "Epithelial"
dt1$myCellType[dt1$clusters.res.0.1==8] <- "Myeloid"
dt1$myCellType[dt1$clusters.res.0.1==9] <- "Mast cells"
dt1$myCellType[dt1$clusters.res.0.1==10] <- "Schwann"
dt1$myCellType[dt1$clusters.res.0.1==11] <- "B cells"
dt1$myCellType[dt1$clusters.res.0.1==12] <- "Epithelial"

expr_comb$manual_clusters1 <- dt1$myCellType

##添加样本信息,这部分应该从本脚本中分离出
expr_comb$tissues="TUMOR"
expr_comb$tissues[expr_comb$orig.ident %in% c("190065A_P1","190065D_P4","190065E_P5")]="CONTROL"
saveRDS(expr_comb,"/home/dklv/scPip/M2classification/M2expr_combination.RDS")

expr_comb<-readRDS("/home/dklv/scPip/M2classification/M2expr_combination.RDS")

newparasdf=data.frame()
for(cellkind in sort(unique(expr_comb$manual_clusters1))){
  
  if(cellkind %in% c('Schwann','Mast cells')){
    pointsize=0.5
  }else{
    pointsize=0.2
  }
  
  #取一大类所有细胞进行亚分类
  aclass <- expr_comb[,expr_comb$manual_clusters1 == cellkind]
  aclass <- FindVariableFeatures(aclass, selection.method = "vst", nfeatures = 1000, verbose = F)
  aclass <- ScaleData(aclass,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"))
  mtx <- as.matrix(aclass@assays$integrated@scale.data)
  library(paran)
  pr <- paran(t(mtx),iterations = 10)
  npca <- sum(pr$Ev/pr$RndEv > 1.5)
  aclass$npca=npca
  aclass <- Seurat::RunPCA(aclass,npcs=100) # pca
  #ElbowPlot(aclass,ndims = 100)
  aclass <- Seurat::FindNeighbors(aclass, reduction = "pca", dims = 1:npca)
  #删除大类的分类结果
  ccols=grep("integrated_snn_res",colnames(aclass@meta.data))
  print(ccols)
  aclass@meta.data=aclass@meta.data[,-ccols]
  #从新分类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.02) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.05) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.1) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.2) # 聚类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.5) # 聚类
  aclass <- Seurat::RunTSNE(aclass, dims = 1:npca)
  #aclass <- Seurat::RunTSNE(aclass, perplexity=90, max_iter=2000, stop_lying_iter=1000,dims = 1:npca)
  aclass <- Seurat::RunUMAP(aclass, dims = 1:npca) # umap
  saveRDS(aclass,file = paste0("/home/rlhua/scRNApip/M2classification/M2",sub(" ","_",cellkind),".RDS"))
}
  
stromal=c("Fibroblasts","Endothelial","Schwann")
immune=c("T cells","B cells","Myeloid","Mast cells")

#提取stromal细胞
  aclass <- expr_comb[,expr_comb$manual_clusters1 %in% stromal]
  aclass <- FindVariableFeatures(aclass, selection.method = "vst", nfeatures = 2000, verbose = F)
  aclass <- ScaleData(aclass,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"))
  mtx <- as.matrix(aclass@assays$integrated@scale.data)
  library(paran)
  pr <- paran(t(mtx),iterations = 10)
  npca <- sum(pr$Ev/pr$RndEv > 1.5)
  aclass$npca=npca
  aclass <- Seurat::RunPCA(aclass,npcs=100) # pca
  #ElbowPlot(aclass,ndims = 100)
  aclass <- Seurat::FindNeighbors(aclass, reduction = "pca", dims = 1:npca)
  #删除大类的分类结果
  ccols=grep("integrated_snn_res",colnames(aclass@meta.data))
  print(ccols)
  aclass@meta.data=aclass@meta.data[,-ccols]
  #从新分类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.1) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.2) # 聚类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.5) # 聚类
  aclass <- Seurat::RunTSNE(aclass, dims = 1:npca)
  #aclass <- Seurat::RunTSNE(aclass, perplexity=90, max_iter=2000, stop_lying_iter=1000,dims = 1:npca)
  aclass <- Seurat::RunUMAP(aclass, dims = 1:npca) # umap
  saveRDS(aclass,file = "/home/dklv/scPip/M2classification/M2stromal.RDS")

  
  #提取immune细胞
  aclass <- expr_comb[,expr_comb$manual_clusters1 %in% immune]
  aclass <- FindVariableFeatures(aclass, selection.method = "vst", nfeatures = 2000, verbose = F)
  aclass <- ScaleData(aclass,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"))
  mtx <- as.matrix(aclass@assays$integrated@scale.data)
  library(paran)
  pr <- paran(t(mtx),iterations = 10)
  npca <- sum(pr$Ev/pr$RndEv > 1.5)
  aclass$npca=npca
  aclass <- Seurat::RunPCA(aclass,npcs=100) # pca
  #ElbowPlot(aclass,ndims = 100)
  aclass <- Seurat::FindNeighbors(aclass, reduction = "pca", dims = 1:npca)
  #删除大类的分类结果
  ccols=grep("integrated_snn_res",colnames(aclass@meta.data))
  print(ccols)
  aclass@meta.data=aclass@meta.data[,-ccols]
  #从新分类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.1) 
  aclass <- Seurat::FindClusters(aclass, resolution = 0.2) # 聚类
  aclass <- Seurat::FindClusters(aclass, resolution = 0.5) # 聚类
  aclass <- Seurat::RunTSNE(aclass, dims = 1:npca)
  #aclass <- Seurat::RunTSNE(aclass, perplexity=90, max_iter=2000, stop_lying_iter=1000,dims = 1:npca)
  aclass <- Seurat::RunUMAP(aclass, dims = 1:npca) # umap
  saveRDS(aclass,file = "/home/dklv/scPip/M2classification/M2immune.RDS")
  table(aclass$manual_clusters1,aclass$integrated_snn_res.0.2)

  #函数1.最简单的展示细胞分类的函数,没有注释
  embed_cluster_seurat<-function(obj,reduction="tSNE",legend.position="bottom"){
    library(ggplot2)
    emb.meta=cbind(obj@meta.data,obj@reductions$tsne@cell.embeddings,obj@reductions$umap@cell.embeddings)
    resolutions=colnames(emb.meta)[grep("integrated_snn_res|manual_clusters",colnames(emb.meta))]
    dim1=paste(reduction,1,sep="_")
    dim2=paste(reduction,2,sep="_")
    plotlist=list()
    for(resolut in resolutions){
      plotlist[[resolut]]=ggplot(emb.meta,aes_string(dim1,dim2,color=resolut))+
        geom_point(shape=16,size=0.1)+
        theme_bw(base_size = 10)+theme(aspect.ratio = 1,panel.grid = element_blank(),axis.text= element_text(color="black"),panel.border = element_rect(color = "black",size=1))+ 
        theme(legend.title=element_blank(),legend.text = element_text(size = 10),legend.position = legend.position,plot.title = element_text(hjust = 0.5)) + 
        guides(colour = guide_legend(override.aes = list(size=5)))+
        labs(title=resolut)
    }
    plotlist
  }
  
  plotlist_tsne=embed_cluster_seurat(aclass,reduction="tSNE")
  plotlist_umap=embed_cluster_seurat(aclass,reduction="UMAP")
  
  plotlist=c(plotlist_tsne,plotlist_umap)
  pg=plot_grid(plotlist=plotlist,ncol=length(plotlist)/2,nrow=2)
  ggsave(plot=pg,filename = "aclass_embed.pdf",width=length(plotlist),height=8)
  
  
  #评价分类效果
  silhouette_width<-function(seuratobj,
                             res=NULL,
                             npca=NULL,
                             res.low = .05, res.high=2, res.n = 20
  ){
    library(Seurat)
    require(cluster)
    
    srobj.tmp = seuratobj
    
    if(is.null(npca)){
      if("npca" %in% colnames(srobj.tmp@meta.data)){
        npca=unique(srobj.tmp$npca)
      }else{
        library(paran)
        pr <- paran(t(as.matrix(srobj.tmp@assays$integrated@scale.data)),iterations = 10)
        npca <- sum(pr$Ev/pr$RndEv > 1.5)
      }
    }
    PCAmatrix=srobj.tmp@reductions$pca@cell.embeddings[,1:npca]
    dis=dist(PCAmatrix)
    
    if(is.null(res)){
      set.res = round(exp(seq(log(res.low), log(res.high), length.out=res.n)), digits=2)
    }else{
      set.res=res
    }
    sillist=list()
    for(i in set.res){
      srobj.tmp = Seurat::FindClusters(srobj.tmp, resolution=i)
      clusters= as.numeric(srobj.tmp@meta.data[,paste0("integrated_snn_res.",i)])
      if(length(table(clusters))>1){
        sil = cluster::silhouette(clusters, dis)
        sillist[[paste0("resolution",i)]]=sil
      }
    }
    sillist
    
  }
  # colors <- scCancer::getDefaultColors(n = length(unique(aclass$seurat_clusters)),type = 2)
  # library(ggplot2)
  # p1=DimPlot(aclass, reduction = "tsne", group.by = "seurat_clusters",cols = colors)+
  #   theme(aspect.ratio=1,legend.position="bottom")
  # colors <- scCancer::getDefaultColors(n = length(unique(aclass$orig.ident)),type = 2)
  # p2=DimPlot(aclass, reduction = "tsne", group.by = "orig.ident",cols = colors)+
  #   theme(aspect.ratio=1,legend.position="bottom")
  # ggsave("aclassroblasts_subcluster_sample_tsne-1000.pdf")
  # p3=DimPlot(aclass, reduction = "tsne", group.by = "tissues",cols = colors)+
  #   theme(aspect.ratio=1,legend.position="bottom")
  # ggsave("aclassroblasts_subcluster_tissue_tsne-1000.pdf")
  
  #aclass$tissues="TUMOR"
  #aclass$tissues[aclass$orig.ident %in% c("190065A_P1","190065D_P4","190065E_PA5")]="CONTROL"
  diff.expr.genes <- Seurat::FindAllMarkers(aclass, only.pos = TRUE, min.pct = 0.25,logfc.threshold = 0.25) # 找每个类的marker基因
  diff.expr.genes <- diff.expr.genes[, c("cluster", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")]
  levels(diff.expr.genes$cluster)=1:length(levels(diff.expr.genes$cluster))
  library(openxlsx)
  write.xlsx(diff.expr.genes,file=paste0(cellkind,".xlsx"))
  library(cowplot)
  
  paras=cbind(aclass@reductions$tsne@cell.embeddings,aclass@reductions$umap@cell.embeddings,aclass@meta.data)
  
  library(stringr)
  library(ggpubr)
  paras$id=str_split(rownames(paras),'_',simplify = T)[,2]
  newparas=data.frame()
  for(i in unique(paras$id)){
    print(i)
    partcells=paras[paras$id==i,]
    cellnum=nrow(partcells)
    totalcount=sum(partcells$nCount_RNA)
    partcells$normalized_nCount=partcells$nCount_RNA*cellnum/totalcount
    newparas=rbind(newparas,partcells)
  }
  newparasdf=rbind(newparasdf,newparas)
  colpan=c("#DB3B32","#CED951","#82BC5E","#476AAE","#8E579B","#EB9B3F","#85CBDB","#DA408F")
  colpan2=c("#DB3B32","#CED951","#476AAE","#8E579B","#EB9B3F","#85CBDB")
  
  get_label_pos <- function(data, emb = "tSNE", group.by="seurat_clusters",add="tissues",labeltype='cluster') {
    cv=table(data[,group.by])
    names(cv)=1:length(names(cv))
    new.data <- data[, c(paste(emb, 1:2, sep = "_"), group.by,add)]
    colnames(new.data) <- c("x","y","cluster",add)
    clusters <- names(table(new.data$cluster))
    if(labeltype=='cluster'){
      new.pos <- lapply(clusters, function(i) {
        tmp.data = subset(new.data, cluster == i)
        data.frame(
          x = median(tmp.data$x),
          y = median(tmp.data$y),
          label = paste0((as.integer(i)+1),'\n','(n=',cv[as.integer(i)+1],')')
        )
      })
    }else{
      new.pos <- lapply(clusters, function(i) {
        tmp.data = subset(new.data, cluster == i)
        data.frame(
          x = median(tmp.data$x),
          y = median(tmp.data$y),
          label = paste0((as.integer(i)+1),'\n','(',paste(table(tmp.data[,add]),collapse = ","),')')
        )
      })
    }
    do.call(rbind, new.pos)
  }
  
  p1=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "seurat_clusters",size=pointsize)+theme(aspect.ratio = 1,legend.position = "None" )+geom_text(inherit.aes = F, data = get_label_pos(newparas, emb = "tSNE"), aes(x,y,label=label), size=3, lineheight=.75) 
  p2=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "orig.ident",size=pointsize)+xlab("")+ylab("")+theme(aspect.ratio = 1 ,legend.position = "None" )+scale_colour_manual(values = colpan)
  p3=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "tissues",size=pointsize)+xlab("")+ylab("")+theme(aspect.ratio = 1 ,legend.position = "None" )+scale_color_manual(values = rev(c("#84BC56","#47559F")))+geom_text(inherit.aes = F, data = get_label_pos(newparas, emb = "tSNE",labeltype="tissues"), aes(x,y,label=label), size=3, lineheight=.75,colour="#67000d")
  p4=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "nCount_RNA",size=pointsize)+scale_color_viridis_c()+xlab("")+ylab("")+theme(aspect.ratio = 1 ,legend.position = "None" )
  p5=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "normalized_nCount",size=pointsize)+scale_color_viridis_c()+xlab("")+ylab("")+theme(aspect.ratio = 1 ,legend.position = "None" )
  p11=ggscatter(newparas[newparas$id %in% c('1','2','4','5','6','7'),],x = "tSNE_1",y="tSNE_2",color = "seurat_clusters",size=pointsize)+theme(aspect.ratio = 1,legend.position = "None" )+geom_text(inherit.aes = F, data = get_label_pos(newparas[newparas$id %in% c('1','2','4','5','6','7'),], emb = "tSNE"), aes(x,y,label=label), size=3, lineheight=.75)
  p21=ggscatter(newparas[newparas$id %in% c('1','2','4','5','6','7'),],x = "tSNE_1",y="tSNE_2",color = "orig.ident",size=pointsize)+xlab("")+ylab("")+theme(aspect.ratio = 1,legend.position = "None" )+scale_colour_manual(values = colpan2)
  p31=ggscatter(newparas[newparas$id %in% c('1','2','4','5','6','7'),],x = "tSNE_1",y="tSNE_2",color = "tissues",size=pointsize)+xlab("")+ylab("")+theme(aspect.ratio = 1,legend.position = "None" )+scale_color_manual(values = rev(c("#84BC56","#47559F")))+geom_text(inherit.aes = F, data = get_label_pos(newparas[newparas$id %in% c('1','2','4','5','6','7'),], emb = "tSNE",labeltype="tissues"), aes(x,y,label=label), size=3, lineheight=.75,colour="#67000d")
  p41=ggscatter(newparas[newparas$id %in% c('1','2','4','5','6','7'),],x = "tSNE_1",y="tSNE_2",color = "nCount_RNA",size=pointsize)+scale_color_viridis_c()+xlab("")+ylab("")+theme(aspect.ratio = 1,legend.position = "None" )
  p51=ggscatter(newparas[newparas$id %in% c('1','2','4','5','6','7'),],x = "tSNE_1",y="tSNE_2",color = "normalized_nCount",size=pointsize)+scale_color_viridis_c()+xlab("")+ylab("")+theme(aspect.ratio = 1,legend.position = "None" )
  pg=plot_grid(p1,p3,p4,p5,p2,p11,p31,p41,p51,p21,NULL,rel_heights=c(3,3,1),ncol=5,nrow=3)
  #pg=plot_grid(p1,p3,p4,p5,p2,p11,p31,p41,p51,p21,ncol=5,nrow=2,axis="lrtb",align="hv")
  ggsave(plot = pg,paste0(cellkind,'.pdf'),width = 12,height = 7)
}
saveRDS(newparasdf,"newparasdf2.RDS")

l2=ggscatter(paras,x = "tSNE_1",y="tSNE_2",color = "orig.ident",size=pointsize)+xlab("")+ylab("")+theme(legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=10))+guides(color = guide_legend(override.aes = list(size=5)))+scale_colour_manual(values = colpan)

l3=ggscatter(paras,x = "tSNE_1",y="tSNE_2",color = "tissues",size=pointsize)+xlab("")+ylab("")+theme(aspect.ratio = 1 )+scale_color_manual(values = rev(c("#84BC56","#47559F")))+theme(legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=10))+guides(color = guide_legend(override.aes = list(size=5)))
l4=ggscatter(paras,x = "tSNE_1",y="tSNE_2",color = "nCount_RNA",size=pointsize)+scale_color_viridis_c()+xlab("")+ylab("")+theme(aspect.ratio = 1)+theme(legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=10))+guides(color = guide_colourbar(barheight  = 1.5, barwidth = 10))
l5=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "normalized_nCount",size=pointsize)+scale_color_viridis_c()+xlab("")+ylab("")+theme(aspect.ratio = 1)+theme(legend.position = "bottom",legend.title = element_blank(),legend.text = element_text(size=10))+guides(color = guide_colourbar(barheight  = 1.5, barwidth = 10))
pg=plot_grid(l3,l4,l5,l2,ncol=4,axis="tblr",align="hv")
ggsave(plot = pg,"legends.pdf",width = 15,height = 5)

cellkind="Fibroblasts"
tsne=read.table('result1/s3.tsne.txt',sep='\t',header=T,row.names=1)
aucell=read.table('result1/s3.AUCell.txt',row.names=1,header=T)
pointsize=0.5
for(cellkind in unique(newparasdf$manual_clusters1)){
  
  newparas=newparasdf[newparasdf$manual_clusters1==cellkind,]
  newparas$X=tsne[rownames(newparas),]$X  
  newparas$Y=tsne[rownames(newparas),]$Y 
  p1=ggscatter(newparas,x = "tSNE_1",y="tSNE_2",color = "seurat_clusters",size=pointsize)+theme(aspect.ratio = 1,legend.position = "None" )+geom_text(inherit.aes = F, data = get_label_pos(newparas, emb = "tSNE"), aes(x,y,label=label), size=3) 
  p2=ggscatter(newparas,x = "X",y="Y",color = "seurat_clusters",size=pointsize)+theme(aspect.ratio = 1,legend.position = "None" )
  newaucell=aucell[rownames(newparas),]
  tsneaucell=Rtsne(newaucell,pca=T)
  newparas$aucell_tSNE_1=tsneaucell$Y[,1]
  newparas$aucell_tSNE_2=tsneaucell$Y[,2]
  p3=ggscatter(newparas,x = "aucell_tSNE_1",y="aucell_tSNE_2",color = "seurat_clusters",size=pointsize)+theme(aspect.ratio = 1,legend.position = "None" )+geom_text(inherit.aes = F, data = get_label_pos(newparas, emb = "aucell_tSNE"), aes(x,y,label=label), size=3) 
  pg=plot_grid(p1,p2,p3,ncol=3,nrow=1)
  ggsave(plot = pg,paste0(cellkind,".pdf"),width = 9,height = 3)
}





