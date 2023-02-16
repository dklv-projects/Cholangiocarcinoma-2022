#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <包含seurat对象> <分类用分辨率>
#功能：细胞类型降维展示
#input：seurat对象
#input: 分类用meta.data列名,如integrated_snn_res.0.2
#output：包括3个pdf文件
#保存在<包含seurat对象目录> 所在目录中的M2results目录中
#M2_fig1cluterwithdifrevolutions.pdf展示不同分辨率下的分类结果
#M2_fig2samplEmbeding.pdf展示不同样本
#M2_fig3cellTypeEmbeding.pdf展示细胞注释结果,分类用meta.data列名

Args<-commandArgs(T)
if(length(Args)>2){stop("parameters number incorrect")}
workdir <- dirname(Args[1]) #项目目录
exprobj <- basename(Args[1])
if(length(Args)==2){
  setresolution <- Args[2] #分类用分辨率设置,如果不设置则使用所有的分辨率做注释展示
}
pwd=getwd()
setwd(workdir)
savePath=paste0(pwd,"/M2classification_",exprobj)
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}

#函数1.最简单的展示细胞分类的函数,没有注释
embed_cluster<-function(emb.meta,reduction="tSNE",colorVar,title="",legend.position="bottom"){
  dim1=paste(reduction,1,sep="_")
  dim2=paste(reduction,2,sep="_")
  pg=ggplot(emb.meta,aes_string(dim1,dim2,color=colorVar))+
    geom_point(shape=16,size=0.1)+
    theme_bw(base_size = 10)+theme(aspect.ratio = 1,panel.grid = element_blank(),axis.text= element_text(color="black"),panel.border = element_rect(color = "black",size=1))+ 
    theme(legend.title=element_blank(),legend.text = element_text(size = 10),legend.direction="vertical",legend.position = legend.position,plot.title = element_text(hjust = 0.5)) + 
    guides(colour = guide_legend(override.aes = list(size=5)))+
    labs(title=title)
  pg
}

#函数2.突出一类细胞的函数,也没有注释
highlight1cluster<-function(emb.meta,reduction="tSNE",colorVar){
  dim1=paste(reduction,1,sep="_")
  dim2=paste(reduction,2,sep="_")
  plot_list=list()
  for(clust in unique(emb.meta[,colorVar])){
    subsample=emb.meta[emb.meta[,colorVar]==clust,]
    subsample$tag="s"
    whole=emb.meta
    whole$tag="w"
    plot_list[[clust]]=ggplot(rbind(whole,subsample),aes_string(dim1,dim2,color="tag"))+
      geom_point(shape=16,size=0.1)+
      theme_bw(base_size = 10)+theme(aspect.ratio = 1,panel.grid = element_blank(),axis.text= element_text(color="black"),panel.border = element_rect(color = "black",size=1))+ 
      theme(legend.title=element_blank(),legend.text = element_text(size = 10),legend.direction="vertical",legend.position = "None",plot.title = element_text(hjust = 0.5)) + 
      guides(colour = guide_legend(override.aes = list(size=5)))+
      labs(title=clust)+scale_color_manual(values = c("#69AD58","grey"))
  }
  plot_list
}

## get_label_pos函数获得每个细胞分类的中心点和细胞数量
#获得类的中心位置和类细胞的数量
get_label_pos <- function(emb.meta, reduction = "tSNE", group.by="seurat_clusters",count=TRUE) {
  #emb.meta:
  #reduction: 降维方法，tSNE、UMAP、pca三选一，仅取前2维
  #group.by:细胞分类依据，来自seurat对象中的meta.data列的名称
  #count:是否显示不同类型细胞的数量
  cols=c(paste(reduction,1:2,sep="_"),group.by)
  new.data = emb.meta[,cols]
  colnames(new.data) = c("x","y","cluster")
  cv=table(new.data$cluster)
  clusters = names(table(new.data$cluster))
  new.pos = lapply(clusters, function(i){
    tmp.data = subset(new.data, cluster == i)
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = ifelse(count,paste0(i,"(",cv[i],")"),i)
    )
  })
  do.call(rbind, new.pos)
}

#函数3.绘制单图,每个类用不同颜色表示，图层上有注释信息
DimRedunc_type<-function(emb.meta,reduction="tSNE",colorVar,labelSize=3,count=T){
  dim1=paste(reduction,1,sep="_")
  dim2=paste(reduction,2,sep="_")
  plot=ggplot(emb.meta,aes_string(dim1,dim2,color=colorVar))+
    geom_point(shape=16,size=0.1)+
    theme_bw(base_size = 10)+
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          legend.position = "None" ,
          axis.text= element_text(color="black"),
          panel.border = element_rect(color = "black",size=1),
          plot.title = element_text(hjust = 0.5)
    )+labs(title=colorVar)
  if(count==TRUE){
    plot=plot+
      geom_text_repel(inherit.aes = F, 
                      data = get_label_pos(emb.meta,reduction = reduction,group.by=colorVar,count = count), 
                      aes(x,y,label=label), size=labelSize)
  }else{
    plot=plot+
      geom_text(inherit.aes = F, 
                data = get_label_pos(emb.meta,reduction = reduction,group.by=colorVar,count = count), 
                aes(x,y,label=label), size=labelSize)
  }
  plot
}

##细胞类型定性分析
get_label_pos_anno <- function(emb.meta, emb = "tSNE", group.by="seurat_clusters",anno="celltype_anno") {
  new.data = emb.meta[, c(paste(emb, 1:2, sep = "_"), group.by, anno)]
  colnames(new.data) = c("x","y","cluster","anno")
  cv=table(new.data$cluster)
  clusters = names(table(new.data$cluster))[table(new.data$cluster)>0]
  new.pos = lapply(clusters, function(i){
    tmp.data = subset(new.data, cluster == i)
    tab = sort(table(tmp.data$anno),decreasing = T)
    cum = cumsum(tab)
    j = which(cum > (0.7*nrow(tmp.data) ))[1]
    new.tab = tab[1:j]
    lab = paste(paste0(i,"(",cv[i],")"),
                paste(names(new.tab),collapse = "\n"),
                paste(round(100*new.tab/nrow(tmp.data)),collapse = ","),
                sep = "\n")
    data.frame(
      x = median(tmp.data$x),
      y = median(tmp.data$y),
      label = lab
    )
  })
  do.call(rbind, new.pos)
}

DimRedunc_type2<-function(emb.meta,reduction="tSNE",colorVar="seurat_clusters",labelSize=1,anno="celltype_anno",title=NULL){
  dim1=paste(reduction,1,sep="_")
  dim2=paste(reduction,2,sep="_")
  plot=ggplot(emb.meta,aes_string(dim1,dim2,color=colorVar))+
    geom_point(shape=16,size=0.1)+
    theme_bw(base_size = 10)+
    theme(aspect.ratio = 1,
          panel.grid = element_blank(),
          legend.position = "None" ,
          axis.text= element_text(color="black"),
          panel.border = element_rect(color = "black",size=1),
          plot.title = element_text(hjust = 0.5)
    ) +
    theme(legend.title=element_blank(),
          legend.text = element_text(size = 10),
          legend.direction="vertical",
          legend.position = "None") +
    geom_text_repel(inherit.aes = F, data = get_label_pos_anno(emb.meta, emb = reduction,group.by=colorVar,anno=anno), aes(x,y,label=label), size=labelSize)+
    labs(title=ifelse(is.null(title),colorVar,title))
  plot
}
##################################################################
## 读取Seurat对象,提取绘图对象，保存到emb.meta
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

emb.meta=cbind(expr_comb@meta.data,expr_comb@reductions$tsne@cell.embeddings,expr_comb@reductions$umap@cell.embeddings)

##开始使用ggplot2绘图
library(ggplot2)
library(ggrepel)
library(cowplot)
#1.下面把所有分辨率的分类降维结果绘第1张图，保存到M2_fig1cluterwithdifrevolutions.pdf文件
resolutions=colnames(emb.meta)[grep("integrated_snn_res|manual_clusters",colnames(emb.meta))]
plotlist_tsne=list()
plotlist2_tsne=list()
plotlist_umap=list()
plotlist2_umap=list()
for(resol in resolutions){
  plotlist_tsne[[resol]]<-DimRedunc_type(emb.meta,reduction="tSNE",colorVar = resol,count = F)
  plotlist2_tsne[[resol]]<-DimRedunc_type(emb.meta,reduction="tSNE",colorVar = resol,count = T)
  plotlist_umap[[resol]]<-DimRedunc_type(emb.meta,reduction="UMAP",colorVar = resol,count = F)
  plotlist2_umap[[resol]]<-DimRedunc_type(emb.meta,reduction="UMAP",colorVar = resol,count = T)
}
pg=plot_grid(plotlist=c(plotlist_tsne,plotlist2_tsne,plotlist_umap,plotlist2_umap),nrow = 4)
ggsave(pg,filename=file.path(savePath,paste0("M2_fig1cluterwithdifrevolutions_",exprobj,".pdf")),width=3*length(resolutions),height=12)

##2.下面把细胞样本来源降维结果绘第2张图，保存到M2_fig2samplEmbeding.pdf文件
##细胞样本来源分布

plot_tsne=embed_cluster(emb.meta,reduction="tSNE",colorVar="orig.ident",title="Sample origin",legend.position = "bottom")
plot_umap=embed_cluster(emb.meta,reduction="UMAP",colorVar="orig.ident",title="Sample origin",legend.position = "bottom")

plot.list.tsne=highlight1cluster(emb.meta,colorVar="orig.ident")
plot.list.tsne[["All"]]=plot_tsne
plot.list.umap=highlight1cluster(emb.meta,reduction="UMAP",colorVar="orig.ident")
plot.list.umap[["All"]]=plot_umap

pg=plot_grid(plotlist = c(plot.list.tsne,plot.list.umap),nrow=2,align="hv",axis = "tb")
ggsave(plot = pg,filename = file.path(savePath,paste0("M2_fig2samplEmbeding_",exprobj,".pdf")),width = 3*length(plot.list.tsne),height = 9)

##3.下面把注释结果按照选定的分类和样本来源绘出，保存到M2_fig3cellTypeEmbeding.pdf文件
plot.tsne=DimRedunc_type2(emb.meta,reduction="tSNE",colorVar="seurat_clusters",labelSize=1,anno="celltype_anno",title="All samples")
plot.umap=DimRedunc_type2(emb.meta,reduction="UMAP",colorVar="seurat_clusters",labelSize=1,anno="celltype_anno",title="All samples")
#pg=plot_grid(plot.tsne,plot.umap)
#ggsave(plot = pg,filename ="M2_meta_anno.pdf",width = 10,height = 5)

plot.tsne.list=list()
plot.umap.list=list()
for(orig in unique(emb.meta$orig.ident)){
  subsample=emb.meta[emb.meta$orig.ident==orig,]
  plot.tsne.list[[orig]]=DimRedunc_type2(subsample,reduction="tSNE",colorVar="seurat_clusters",labelSize=1,anno="celltype_anno",title=orig)
  plot.umap.list[[orig]]=DimRedunc_type2(subsample,reduction="UMAP",colorVar="seurat_clusters",labelSize=1,anno="celltype_anno",title=orig)
}
plot.tsne.list[["All"]]=plot.tsne
plot.umap.list[["All"]]=plot.umap
pg=plot_grid(plotlist = c(plot.tsne.list,plot.umap.list),nrow = 2)
ggsave(plot = pg,filename = file.path(savePath,paste0("M2_fig3cellTypeEmbeding_",exprobj,"_",setresolution,".pdf")),width = 5*length(plot.tsne.list),height = 10)


