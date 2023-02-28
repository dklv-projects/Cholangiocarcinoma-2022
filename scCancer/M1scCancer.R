#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <项目目录> <物种(human/mouse)>
#功能：数据预处理，包括细胞筛选、基因筛选、normaliation、sample integration、scale、降维、细胞分类和分类后细胞类型注释
#input：cellranger的结果
#output：在程序运行时的当前目录中新建M1scCancer作为输出目录，主要结果保存在comb目录中的expr_comb.RDS
#项目目录：包含所有cellranger结果的目录，其中的每个子目录下有outs目录
#
Args<-commandArgs(T)
if(length(Args)!=2){stop("parameters number incorrect")}
workdir <- Args[1] #项目目录
species <- Args[2]
if(!(species %in% c("human","mouse"))){stop("incorrect species")}
if(species=="human"){genome="hg38"}else{genome="mm10"}
print(genome)
#保存当前目录，进入项目目录
pwd=getwd()
setwd(workdir)
#project_dir=getwd()
#project_name=basename(project_dir)
#在程序运行时的当前目录中新建M1scCancer作为输出目录
savePath=paste0(pwd,"/M1scCancer")
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}

dir=list.dirs(recursive=F,full.names=F)
library(scCancer)
#1.质量控制、注释,结果写入_stat/和_anno/
for(fold in dir){
  #输入文件夹为project_dir
  #输出文件夹为project_name_stat
  stat.result <- runScStatistics(
    #这里的dataPath是cellrange跑出来的output文件夹，
    #将它重命名后，整理到一个文件夹内，包含所有文库对应的文件夹，
    #dataPath的上级目录是这个项目的分析结果文件夹
    dataPath = file.path(fold,"outs"),
    savePath = paste0(savePath,"/stat/",fold),
    sampleName = fold,
    species=species,
    mix.anno = c("human" = "hg38", "mouse" = "mm10")
    )
  
  #输入文件夹1为project_dir
  #输入文件夹2为project_name_stat
  #输出文件夹为project_name_anno
  anno.results <- runScAnnotation(
    dataPath = file.path(fold,"outs"),
    statPath = paste0(savePath,"/stat/",fold),
    savePath = paste0(savePath,"/anno/",fold),
    sampleName = fold,
    geneSet.method = "average",
    pc.use = 50,
    species = species,
    genome = genome
  )
}

#2.多样本整合,结果写入expr.comb
#combinedata函数的用途：features.to.integrate的参数是scCancer中runScCombination无法传递的
combinedata <- function(single.savePaths,
                        savePath,
                        genome="hg38",
                        npc=50,
                        vars.to.regress=c("nCount_RNA", "mito.percent", "ribo.percent")
){
  #输入文件夹为上一步anno的输出文件夹
  #输出为保存expr_comb的文件夹
  expr.list <- list()
  for(sp in single.savePaths){
    sampleName <- basename(sp)
    print(sampleName)
    expr.list[[sampleName]] <- readRDS(paste0(sp, "/expr.RDS"))
  } # 读取单样本数据
  expr.anchors <- Seurat::FindIntegrationAnchors(object.list = expr.list,dims = 1:npc) # 找anchor
  #下面对anchor进行筛选后更新
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
  #这一步features.to.integrate的参数是scCancer中runScCombination无法传递的，也是combinedata存在的唯一原因
  expr_comb <- Seurat::IntegrateData(anchorset = expr.anchors,dims = 1:npc,features.to.integrate =genes) # 整合数据
  expr_comb <- Seurat::ScaleData(expr_comb,vars.to.regress = vars.to.regress) # 中心化
  
  #函数用于计算有效PCA数
  number.PCA=function(matrix){
    library(paran)
    mtx=as.matrix(matrix)
    pr=paran(t(mtx),iterations = 10)
    npca=sum(pr$Ev/pr$RndEv > 1.5)
    return(npca)
  }
  npca=number.PCA(expr_comb@assays$integrated@scale.data)
  expr_comb$npca=npca
  #npca==38
  print(npca)
  #所有细胞聚类
  expr_comb <- Seurat::RunPCA(expr_comb,npcs=npca) # pca
  expr_comb <- Seurat::FindNeighbors(expr_comb, reduction = "pca", dims = 1:npca)
  expr_comb <- Seurat::FindClusters(expr_comb, resolution = 0.5) # 聚类
  expr_comb <- Seurat::FindClusters(expr_comb, resolution = 0.2) # 聚类
  expr_comb <- Seurat::FindClusters(expr_comb, resolution = 0.1) 
  expr_comb <- Seurat::FindClusters(expr_comb, resolution = 0.05) 
  expr_comb <- Seurat::FindClusters(expr_comb, resolution = 0.02) 
  #注意：设置分类选项，需要根据需要选择！！！，建议跑完M2模块后修改M2模块
  expr_comb$seurat_clusters=expr_comb$integrated_snn_res.0.2
  expr_comb@active.ident = expr_comb$seurat_clusters #设置默认分类
  #所有细胞降纬
  expr_comb <- Seurat::RunTSNE(expr_comb, dims = 1:npca)
  expr_comb <- Seurat::RunUMAP(expr_comb, dims = 1:npca) # umap
  #expr_comb <- Seurat::RunTSNE(expr_comb, perplexity=60, max_iter=2000, stop_lying_iter=1000,dims = 1:npc) # tsne
  # if(genome=="hg38"){
  #   #所有细胞类型注释
  #   library(scHCL)
  #   hcl_anno<-scHCL(scdata = as.matrix(expr_comb@assays$RNA@counts), numbers_plot = 1)
  #   expr_comb$hcl_anno<- hcl_anno$scHCL_probility[,2]
  #   expr_comb$celltype_anno<-expr_comb$hcl_anno
  # }
  # if(genome=="mm10"){
  #   library(scMCA)
  #   mca_anno<-scMCA(scdata = as.matrix(expr_comb@assays$RNA@counts), numbers_plot = 1)
  #   expr_comb$mca_anno<- mca_anno$scMCA_probility[,2]
  #   expr_comb$celltype_anno<-expr_comb$mca_anno
  # }
  
  #保存整合后数据
  if(!dir.exists(savePath)){
    dir.create(savePath,recursive = T)
  }
  saveRDS(expr_comb,file = file.path(savePath,"expr_combined.RDS"))
  return(expr_comb)
}

if(length(dir)>1){
  single.savePaths=list.dirs(paste0(savePath,"/anno"),full.names = T,recursive = F)
  expr_comb <- combinedata(single.savePaths,
                           savePath = paste0(savePath,"/comb"))
}


