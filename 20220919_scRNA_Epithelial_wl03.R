library(Seurat)
library(data.table)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(paran)
library(pheatmap)
library(limma)
library(survival)
library(survminer)

sr <- readRDS("/home/rlhua/project/cholangiocarcinoma/M2classification/M2Epithelial1.RDS")

sr <- FindVariableFeatures(sr, selection.method = "vst", nfeatures = 1000, verbose = F)
sr <- ScaleData(sr,vars.to.regress = c("nCount_RNA", "mito.percent", "ribo.percent"))
mtx <- as.matrix(sr@assays$integrated@scale.data)
pr <- paran(t(mtx),iterations = 10)
npca <- sum(pr$Ev/pr$RndEv > 1.5)
sr$npca=npca
sr <- Seurat::RunPCA(sr,npcs=100) # pca
ElbowPlot(sr,ndims = 100)
sr <- Seurat::FindNeighbors(sr, reduction = "pca", dims = 1:npca)
#删除大类的分类结果
ccols=grep("integrated_snn_res",colnames(sr@meta.data))
print(ccols)
sr@meta.data=sr@meta.data[,-ccols]
#从新分类
sr <- Seurat::FindClusters(sr, resolution = 0.02) 
sr <- Seurat::FindClusters(sr, resolution = 0.05) 
sr <- Seurat::FindClusters(sr, resolution = 0.1) 
sr <- Seurat::FindClusters(sr, resolution = 0.2) # 聚类
sr <- Seurat::FindClusters(sr, resolution = 0.5) # 聚类
sr <- Seurat::RunTSNE(sr, dims = 1:npca,seed.use=10)
# sr <- Seurat::RunTSNE(sr, perplexity=90, max_iter=2000, stop_lying_iter=1000,dims = 1:npca)
sr <- Seurat::RunUMAP(sr, dims = 1:npca) # umap

##### infercnv
endo <- readRDS("/home/rlhua/project/cholangiocarcinoma/M2classification/M2Endothelial1.RDS")
fib <- readRDS("/home/rlhua/project/cholangiocarcinoma/M2classification/M2Fibroblasts1.RDS")
ep <- sr
set.seed(1234)
endo_mtx <- endo@assays$RNA@counts[,sample(1:ncol(endo),500)]
set.seed(1234)
fib_mtx <- fib@assays$RNA@counts[,sample(1:ncol(fib),500)]
ep_mtx <- ep@assays$RNA@counts
cnt_mtx <- cbind(endo_mtx,fib_mtx,ep_mtx)
cellanno <- rbind(data.table(colnames(endo_mtx),"endothelial"),
                  data.table(colnames(fib_mtx),"fibroblast"),
                  data.table(colnames(ep_mtx),paste0("C",ep$seurat_clusters)))
fwrite(cellanno,"/home/rlhua/project/cholangiocarcinoma/result/infercnv/cellAnnotations.txt",sep = "\t",col.names = F,quote = F)

infercnv_obj = CreateInfercnvObject(raw_counts_matrix=cnt_mtx,
                                    annotations_file="/home/rlhua/project/cholangiocarcinoma/result/infercnv/cellAnnotations.txt",
                                    delim="\t",
                                    gene_order_file="/home/rlhua/project/cholangiocarcinoma/result/infercnv/gene_ordering_file.txt",
                                    ref_group_names=c("endothelial","fibroblast"))
infercnv_obj = infercnv::run(infercnv_obj,
                             cutoff=0.1,
                             out_dir="/home/rlhua/project/cholangiocarcinoma/result/infercnv/",
                             cluster_by_groups=T,
                             hclust_method="ward.D2",
                             num_threads=10,
                             denoise=T,
                             HMM=T,
                             analysis_mode="subclusters",
                             sd_amplifier=3,
                             noise_logistic=TRUE
)

##### 根据infercnv结果区分癌和上皮细胞
## 读入infercnv结果
df_cnv <- fread("/home/rlhua/project/cholangiocarcinoma/result/infercnv/infercnv.observations.txt",sep=" ")
df_hmm <- fread("/home/rlhua/project/cholangiocarcinoma/result/infercnv/infercnv.19_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt")
df_cnv <- as.data.frame(df_cnv)
rownames(df_cnv) <- df_cnv$V1
df_cnv$V1 <- NULL
df_hmm <- as.data.frame(df_hmm)
rownames(df_hmm) <- df_hmm$V1
df_hmm$V1 <- NULL

## 计算cnv打分和相关系数
sr$patient <- gsub("_.*","",sr$orig.ident)
df <- df_cnv
df[df_hmm==1] <- 1 ## 选择可信的CNV
df <- df -1 
cnv_score <- apply(df, 2, function(x){sum(x^2)}) # 每个细胞的CNV打分
cnv_cor <- lapply(unique(sr$patient), function(x){
  pt_score <- cnv_score[names(cnv_score) %in% colnames(sr)[sr$patient==x]] #取出同一病人样本的细胞
  topcells <- names(sort(pt_score,decreasing=T))[1:ceiling(length(pt_score)*0.05)] #CNV打分最高的前5%细胞
  cat(ceiling(length(pt_score)*0.05),"\n")
  toppattern <- apply(df[,topcells], 1, mean) # 前5%细胞的CNV模式
  cnv_cor <- apply(df[,colnames(df) %in% colnames(sr)[sr$patient==x]],2,cor,y=toppattern) #其他细胞和前5%细胞CNV模式的相关系数
  return(cnv_cor)
}) # 每个细胞的cnv模式和前5%细胞的相关系数
cnv_cor <- unlist(cnv_cor)
cnv_cor <- cnv_cor[colnames(sr)]
cnv_score <- cnv_score[colnames(sr)]

## cnv打分和相关系数的散点图
dt <- data.table(cnv_cor=cnv_cor,
                 cnv_score=cnv_score,
                 cell=names(cnv_cor),
                 patient=sr$patient,
                 group=sr$tissues)
dt$group <- factor(dt$group,levels=c("TUMOR","CONTROL"))
lst <- lapply(unique(dt$patient), function(x){
  p <- ggplot(dt[patient==x],aes(x=cnv_score,y=cnv_cor,color=group))+
    geom_point(size=0.1)+
    geom_hline(yintercept = 0.5,size=0.4)+
    geom_vline(xintercept = 1,size=0.4)+
    scale_color_manual(values = c('dodgerblue4','palegreen3'))+
    ggtitle(x)+
    xlab("CNV score")+ylab("CNV correlation")+
    xlim(0,7)+
    theme_bw()+
    theme(aspect.ratio = 0.8,
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_rect(color="black"),
          legend.position = "none",
          text = element_text(size = 7))
})
p <- plot_grid(plotlist = lst,ncol = 1,nrow = 5,align = "hv",axis = "tblr")
fn <- "/home/rlhua/project/cholangiocarcinoma/result/Epithelial_CNVscore_corr_tissue_scatterplot.pdf"
save_plot(fn,p,ncol = 1,nrow = 5,base_height = 2,base_width = 2,limitsize = FALSE)

dt$status <- "uncertain"
dt$status[dt$cnv_cor > 0.5 & dt$cnv_score > 1] <- "cancer" # 根据阈值筛选癌细胞和正常上皮
dt$status[dt$cnv_cor < 0.5 & dt$cnv_score < 1] <- "epithelial"
dt$status <- factor(dt$status,levels = c("cancer","epithelial","uncertain"))
dt <- dt[order(status,decreasing = T)]
lst <- lapply(unique(dt$patient), function(x){
  p <- ggplot(dt[patient==x],aes(x=cnv_score,y=cnv_cor,color=status))+
    geom_point(size=0.1)+
    geom_hline(yintercept = 0.5,size=0.4)+
    geom_vline(xintercept = 1,size=0.4)+
    scale_color_manual(values = c("red","blue","black"))+
    ggtitle(x)+
    xlab("CNV score")+ylab("CNV correlation")+
    xlim(0,7)+
    theme_bw()+
    theme(aspect.ratio = 0.8,
          plot.title = element_text(hjust = 0.5),
          panel.grid = element_blank(),
          panel.border = element_rect(color="black"),
          legend.position = "none",
          text = element_text(size = 7))
})
p <- plot_grid(plotlist = lst,ncol = 1,nrow = 5,align = "hv",axis = "tblr")
fn <- "/home/rlhua/project/cholangiocarcinoma/result/Epithelial_CNVscore_corr_scatterplot.pdf"
save_plot(fn,p,ncol = 1,nrow = 5,base_height = 2,base_width = 2,limitsize = FALSE)

## 根据阈值区分癌和上皮细胞
sr$status <- "uncertain"
sr$status[cnv_cor>0.5 & cnv_score>1] <- "cancer"
sr$status[cnv_cor<0.5 & cnv_score<1] <- "epithelial"

## 每个亚类确定细胞类型
sr$seurat_clusters <- as.character(sr$integrated_snn_res.0.2)
sr$seurat_clusters[sr$integrated_snn_res.0.5==9] <- "9" # 手动分开一个cluster
sr$seurat_clusters <- factor(sr$seurat_clusters,levels=names(sort(table(sr$seurat_clusters),decreasing=T)))
sr$seurat_clusters <- as.numeric(sr$seurat_clusters)-1
sr$seurat_clusters <- factor(sr$seurat_clusters)

sr$sample <- NA
for (i in unique(sr$seurat_clusters)) {
  status <- sr$status[sr$seurat_clusters==i]
  cat()
  aa <- sum(status=="cancer")/sum(status=="epithelial") ## 根据癌和上皮细胞比例定义亚类性质
  cat(i,aa,"\n")
  if(aa>1){sr$sample[sr$seurat_clusters==i] <- "cancer"}else{
    sr$sample[sr$seurat_clusters==i] <- "cholangiocyte"
  }
}
sr$sample[sr$seurat_clusters==9] <- "SMC"
sr$sample[sr$seurat_clusters==3] <- "HSC"

## 计算CNV矩阵的细胞距离
lst <- lapply(unique(sr$patient), function(x){
  cat(x,"\n")
  mtx <- as.matrix(df[,colnames(sr)[sr$patient==x]]) #按病人样本计算细胞距离
  mtx <- mtx[apply(mtx, 1, var)>0,] #去除在所有细胞无cnv的基因
  out.dist <- dist(t(mtx))
  return(out.dist)
})
names(lst) <- unique(sr$patient)
saveRDS(lst,"/home/rlhua/project/cholangiocarcinoma/result/Epithelial_CNV_dist.RDS") #保存数据

## CNV热图
out.dist <- readRDS("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_CNV_dist.RDS") #读入细胞cnv距离数据
lst <- lapply(unique(sr$patient), function(x){
  cat(x,"\n")
  hc <- hclust(out.dist[[x]],method = "ward.D2") #对每个病人样本细胞根据CNV聚类
  out.id <- cutree(hc,k=3) # 分3类
  dt1 <- dt[cell %in% names(out.id)] 
  dt1$cluster <- as.character(out.id[dt1$cell])
  p <- ggplot(dt1,aes(x=cnv_score,y=cnv_cor,color=cluster))+
    geom_point(size=0.2)+
    theme_bw()+
    theme(aspect.ratio = 1) ## 聚类后的散点图
  anno <- data.frame(cluster=as.character(out.id),
                     sample=sr$orig.ident[names(out.id)])
  rownames(anno) <- names(out.id)
  mtx <- as.matrix(df[,hc$labels[hc$order]])
  # 热图
  pheatmap(t(df[,hc$labels[hc$order]]),
           cluster_rows = F,cluster_cols = F,
           annotation_row = anno,
           filename = paste0("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_",x,"_infercnv_pheatmap.png"))
  return(p)
})
p <- plot_grid(plotlist = lst,ncol = 1,nrow = 5,align = "hv",axis = "tblr")
fn <- "/home/rlhua/project/cholangiocarcinoma/result/Epithelial_patient_CNVscore_corr_clustering_scatterplot.pdf"
save_plot(fn,p,ncol = 1,nrow = 5,base_height = 4,base_width = 4,limitsize = FALSE)

ca <- lapply(unique(dt$patient), function(x){
  cat(x,"\n")
  cells <- dt[patient==x & status=="cancer",cell]
  dist1 <- as.dist(as.matrix(out.dist[[x]])[cells,cells])
  hc <- hclust(dist1,method = "ward.D2")
  cells <- hc$labels[hc$order]
  return(cells)
}) %>% unlist() #癌细胞根据患者细胞聚类
ep <- lapply(unique(dt$patient), function(x){
  cat(x,"\n")
  cells <- dt[patient==x & status=="epithelial",cell]
  dist1 <- as.dist(as.matrix(out.dist[[x]])[cells,cells])
  hc <- hclust(dist1,method = "ward.D2")
  cells <- hc$labels[hc$order]
  return(cells)
}) %>% unlist() #正常上皮根据患者细胞聚类

# 癌细胞热图
anno <- data.frame(patient=sr$patient[ca],
                   sample=dt$group[match(ca,dt$cell)])
rownames(anno) <- ca
pheatmap(t(df[,ca]),cluster_rows = F,cluster_cols = F,
         annotation_row = anno,
         show_rownames = F,show_colnames = F,
         filename = "/home/rlhua/project/cholangiocarcinoma/result/Epithelial_infercnv_cancer_pheatmap.png")
# 正常上皮细胞热图
anno <- data.frame(patient=sr$patient[ep],
                   sample=dt$group[match(ep,dt$cell)])
rownames(anno) <- ep
pheatmap(t(df[,ep]),cluster_rows = F,cluster_cols = F,
         annotation_row = anno,
         show_rownames = F,show_colnames = F,
         filename = "/home/rlhua/project/cholangiocarcinoma/result/Epithelial_infercnv_normal_pheatmap.png")

##### SingleR annotation
mtx <- sr@assays$RNA@data[,sr$seurat_clusters==3]
pt <- sr$orig.ident[sr$seurat_clusters==3]

library(data.table)
library(magrittr)
library(Seurat)
library(SingleR)
library(scater)
library(pheatmap)

hpca.se <- celldex::HumanPrimaryCellAtlasData()
Blue.se <- celldex::BlueprintEncodeData()
Immune.se <- celldex::DatabaseImmuneCellExpressionData()
Nover.se <- celldex::NovershternHematopoieticData()
MonacoIm.se <- celldex::MonacoImmuneData()
ImmGen.se <- celldex::ImmGenData() #(鼠)
Mouse.se <- celldex::MouseRNAseqData() #(鼠)
save(hpca.se,Blue.se,Immune.se,Nover.se,MonacoIm.se,ImmGen.se,Mouse.se,file = "/home/rlhua/project/cholangiocarcinoma/result/SingleR/database.Rdata")

pred_hpca <- SingleR(test = mtx, ref = hpca.se, labels = hpca.se$label.fine,de.method="wilcox")
pred_Immune <- SingleR(test = mtx, ref = Immune.se, labels = Immune.se$label.fine,de.method="wilcox")
pred_Nover <- SingleR(test = mtx, ref = Nover.se, labels = Nover.se$label.fine,de.method="wilcox")
pred_Blue <- SingleR(test = mtx, ref = Blue.se, labels = Blue.se$label.fine,de.method="wilcox")
pred_MonacoIm <- SingleR(test = mtx, ref = MonacoIm.se, labels = MonacoIm.se$label.fine,de.method="wilcox")
save(pred_hpca,pred_Immune,pred_Nover,pred_Blue,pred_MonacoIm,file = "/home/rlhua/project/cholangiocarcinoma/result/SingleR/HSPC/HSPC.Rdata")

pdf("/home/rlhua/project/cholangiocarcinoma/result/SingleR/HSPC/HSPC.pdf",height = 14,width = 14)
plotScoreHeatmap(pred_hpca,cluster_cols=T,clusters=pt)
plotScoreHeatmap(pred_Immune,cluster_cols=T,clusters=pt)
plotScoreHeatmap(pred_Nover,cluster_cols=T,clusters=pt)
plotScoreHeatmap(pred_Blue,cluster_cols=T,clusters=pt)
plotScoreHeatmap(pred_MonacoIm,cluster_cols=T,clusters=pt)
dev.off()

##### celltypist annotation
## Immune_All_AddPIP
predicted_labels <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Immune_All_AddPIP/predicted_labels.csv",sep=",",header=T)
colnames(predicted_labels)[1] <- "cell"
predicted_labels$subcluster  <- sr$seurat_clusters[predicted_labels$cell]
probability_matrix <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Immune_All_AddPIP/probability_matrix.csv",sep=",",header=T)
probability_matrix <- as.data.frame(probability_matrix)
rownames(probability_matrix) <- probability_matrix$V1
probability_matrix$V1 <- NULL
probability_matrix <- as.matrix(probability_matrix)
probability_matrix <- t(probability_matrix)
probability_matrix <- probability_matrix[,predicted_labels[order(subcluster,predicted_labels),cell]]
anno_col <- data.frame(subcluster=sr$seurat_clusters[colnames(probability_matrix)])
pheatmap(probability_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=8,
         fontsize_row=8,
         scale="column",
         show_colnames=F,
         show_rownames = T,
         fontsize_col=10,
         cellheight = 10,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=anno_col,
         filename = "~/project/cholangiocarcinoma/result/celltypist/Epithelial/Immune_All_AddPIP/probability_matrix_pheatmap.png")
decision_matrix <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Immune_All_AddPIP/decision_matrix.csv",sep=",",header=T)
decision_matrix <- as.data.frame(decision_matrix)
rownames(decision_matrix) <- decision_matrix$V1
decision_matrix$V1 <- NULL
decision_matrix <- as.matrix(decision_matrix)
decision_matrix <- t(decision_matrix)
decision_matrix <- decision_matrix[,predicted_labels[order(subcluster,predicted_labels),cell]]
pheatmap(decision_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=8,
         fontsize_row=8,
         scale="none",
         show_colnames=F,
         show_rownames = T,
         fontsize_col=10,
         cellheight = 10,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=anno_col,
         filename = "~/project/cholangiocarcinoma/result/celltypist/Epithelial/Immune_All_AddPIP/decision_matrix_pheatmap.png")

## Cells_Intestinal_Tract
predicted_labels <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Cells_Intestinal_Tract/predicted_labels.csv",sep=",",header=T)
colnames(predicted_labels)[1] <- "cell"
predicted_labels$subcluster  <- sr$seurat_clusters[predicted_labels$cell]
probability_matrix <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Cells_Intestinal_Tract/probability_matrix.csv",sep=",",header=T)
probability_matrix <- as.data.frame(probability_matrix)
rownames(probability_matrix) <- probability_matrix$V1
probability_matrix$V1 <- NULL
probability_matrix <- as.matrix(probability_matrix)
probability_matrix <- t(probability_matrix)
probability_matrix <- probability_matrix[,predicted_labels[order(subcluster,predicted_labels),cell]]
anno_col <- data.frame(subcluster=sr$seurat_clusters[colnames(probability_matrix)])
pheatmap(probability_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=8,
         fontsize_row=8,
         scale="column",
         show_colnames=F,
         show_rownames = T,
         fontsize_col=10,
         cellheight = 10,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=anno_col,
         filename = "~/project/cholangiocarcinoma/result/celltypist/Epithelial/Cells_Intestinal_Tract/probability_matrix_pheatmap.png")
decision_matrix <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Cells_Intestinal_Tract/decision_matrix.csv",sep=",",header=T)
decision_matrix <- as.data.frame(decision_matrix)
rownames(decision_matrix) <- decision_matrix$V1
decision_matrix$V1 <- NULL
decision_matrix <- as.matrix(decision_matrix)
decision_matrix <- t(decision_matrix)
decision_matrix <- decision_matrix[,predicted_labels[order(subcluster,predicted_labels),cell]]
pheatmap(decision_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=8,
         fontsize_row=8,
         scale="none",
         show_colnames=F,
         show_rownames = T,
         fontsize_col=10,
         cellheight = 10,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=anno_col,
         filename = "~/project/cholangiocarcinoma/result/celltypist/Epithelial/Cells_Intestinal_Tract/decision_matrix_pheatmap.png")

## Healthy_COVID19_PBMC
predicted_labels <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Healthy_COVID19_PBMC/predicted_labels.csv",sep=",",header=T)
colnames(predicted_labels)[1] <- "cell"
predicted_labels$subcluster  <- sr$seurat_clusters[predicted_labels$cell]
probability_matrix <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Healthy_COVID19_PBMC/probability_matrix.csv",sep=",",header=T)
probability_matrix <- as.data.frame(probability_matrix)
rownames(probability_matrix) <- probability_matrix$V1
probability_matrix$V1 <- NULL
probability_matrix <- as.matrix(probability_matrix)
probability_matrix <- t(probability_matrix)
probability_matrix <- probability_matrix[,predicted_labels[order(subcluster,predicted_labels),cell]]
anno_col <- data.frame(subcluster=sr$seurat_clusters[colnames(probability_matrix)])
pheatmap(probability_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=8,
         fontsize_row=8,
         scale="column",
         show_colnames=F,
         show_rownames = T,
         fontsize_col=10,
         cellheight = 10,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=anno_col,
         filename = "~/project/cholangiocarcinoma/result/celltypist/Epithelial/Healthy_COVID19_PBMC/probability_matrix_pheatmap.png")
decision_matrix <- fread("~/project/cholangiocarcinoma/result/celltypist/Epithelial/Healthy_COVID19_PBMC/decision_matrix.csv",sep=",",header=T)
decision_matrix <- as.data.frame(decision_matrix)
rownames(decision_matrix) <- decision_matrix$V1
decision_matrix$V1 <- NULL
decision_matrix <- as.matrix(decision_matrix)
decision_matrix <- t(decision_matrix)
decision_matrix <- decision_matrix[,predicted_labels[order(subcluster,predicted_labels),cell]]
pheatmap(decision_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=8,
         fontsize_row=8,
         scale="none",
         show_colnames=F,
         show_rownames = T,
         fontsize_col=10,
         cellheight = 10,
         treeheight_row = 0,treeheight_col = 0,
         annotation_col=anno_col,
         filename = "~/project/cholangiocarcinoma/result/celltypist/Epithelial/Healthy_COVID19_PBMC/decision_matrix_pheatmap.png")


##### tsneplot
sr$subcluster <- ""
sr$subcluster[sr$seurat_clusters==0] <- "ANXA1-C0-CA"
sr$subcluster[sr$seurat_clusters==1] <- "MYO1A-C1-CHO"
sr$subcluster[sr$seurat_clusters==2] <- "KIT-C2-HSPC"
sr$subcluster[sr$seurat_clusters==3] <- "NFKB1-C3-CA"
sr$subcluster[sr$seurat_clusters==4] <- "FXYD2-C4-CHO"
sr$subcluster[sr$seurat_clusters==5] <- "PSCA-C5-CA"
sr$subcluster[sr$seurat_clusters==6] <- "XAF1-C6-CHO"
sr$subcluster[sr$seurat_clusters==7] <- "SAA1-C7-CHO"
sr$subcluster[sr$seurat_clusters==8] <- "KRT19-C8-CA"
sr$subcluster[sr$seurat_clusters==9] <- "no-C9-CHO"
sr$subcluster[sr$seurat_clusters==10] <- "no-C10-CA"
sr$subcluster[sr$seurat_clusters==11] <- "ACTA1-C11-SMC"

FeaturePlot(sr,reduction="tsne",pt.size=0.1,features=c("ANXA1","MYO1A","CD34","NFKB1","FXYD2","PSCA","XAF1","SAA1","KRT19","ACTA1"),ncol=5)
colors <- c("#21855F","#549024","#D19512","#8D601C","#525252","#BF4812","#C7146B","#9BCB73","#5C5496","#2C8629","#205B94")
colors <- c("#BF4812","#549024","#5C5496","#525252","#C7146B","#D19512","#21855F","#9BCB73","#8D601C","#205B94")
names(colors) <- 0:9
idx <- sample(1:ncol(sr))
ggplot()+
  geom_point(aes(x=sr@reductions$tsne@cell.embeddings[,"tSNE_1"][idx],
                 y=sr@reductions$tsne@cell.embeddings[,"tSNE_2"][idx],
                 color=sr$subcluster[idx]),shape=16,size=0.001)+
  scale_color_manual(values = colors)+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position="none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_subcluster_tsneplot.png",width = 1.7,height = 1.7,dpi=600)
ggplot()+
  geom_point(aes(x=sr@reductions$tsne@cell.embeddings[,"tSNE_1"][idx],
                 y=sr@reductions$tsne@cell.embeddings[,"tSNE_2"][idx],
                 color=sr$tissues[idx]),shape=16,size=0.001)+
  scale_color_manual(values = c('palegreen3','dodgerblue4'))+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_tissue_tsneplot.png",width = 1.7,height = 1.7,dpi=600)
ggplot()+
  geom_point(aes(x=sr@reductions$tsne@cell.embeddings[,"tSNE_1"][idx],
                 y=sr@reductions$tsne@cell.embeddings[,"tSNE_2"][idx],
                 color=sr$sample[idx]),shape=16,size=0.001)+
  scale_color_brewer(palette="Set1")+
  xlab("tSNE_1")+ylab("tSNE_2")+
  theme_classic()+
  theme(aspect.ratio = 1,legend.position = "none",axis.text = element_blank(),axis.line = element_blank(),axis.title = element_blank(),axis.ticks = element_blank())
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_celltype_tsneplot.png",width = 1.7,height = 1.7,dpi=600)

##### heatmap of subcluster markers
markers <- readRDS("/home/rlhua/project/cholangiocarcinoma/M3markers/M3markers_M2Epithelial.RDS/M3_clustermarker2.RDS")
ca <- sr[,sr$sample=="cancer"]
ca <- ScaleData(ca,features=rownames(ca@assays$RNA@data),assay="RNA")
ca$subcluster <- factor(ca$subcluster,levels = names(sort(table(ca$subcluster),decreasing=T)))
top_gene <- lapply(sort(unique(ca$seurat_clusters)), function(x){
  gene <- markers[markers$cluster==x,"gene"][1:10]
}) %>% unlist()

colors <- scCancer::getDefaultColors(n =length(unique(ca$subcluster)),type = 2)
p <- DoHeatmap(ca, features = top_gene,assay = "RNA",group.by="subcluster",group.colors=colors,size=2)+theme(text = element_text(size=7))
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_subcluster_markers_heatmap.pdf",height = 40,width = 10)

##### survival analysis of subcluster markers
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

## markers of each cluster
markers <- readRDS("/home/rlhua/project/cholangiocarcinoma/M3markers/M3markers_M2Epithelial.RDS/M3_clustermarker2.RDS")

## survival plot
lst <- lapply(sort(unique(markers$cluster)), function(x){
  genes <- markers[markers$cluster==x,"gene"] ## 取marker基因
  genes <- intersect(genes,rownames(mtx)) ## marker基因与表达矩阵基因取交集
  gsva_score <- GSVA::gsva(mtx,list(genes)) ## GSVA
  dt <- survival_info[type=="cancer"]
  dt$gsva_score <- gsva_score[,dt$sample]
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
  p <- ggsurvplot(fit,
                  conf.int = F, #不画置信区间，想画置信区间就把F改成T
                  #conf.int.style = "step",#置信区间的类型，还可改为ribbon
                  censor = F, #不显示观察值所在的位置
                  palette = c("#D95F02","#1B9E77"), #线的颜色对应高、低
                  xlab="Months",
                  legend=c(0.85,0.95),
                  title=paste("Cluster",x),
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
  return(p)
})
p <- arrange_ggsurvplots(lst, print = T,ncol = 5, nrow = 2)
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_subcluster_markers_GSVA_survivalplot.pdf",p,width = 15,height = 6,limitsize = FALSE)

## boxplot of markers of each cluster
lst <- lapply(sort(unique(markers$cluster)), function(x){
  genes <- markers$gene[markers$cluster==x] ## 取marker基因
  genes <- intersect(genes,rownames(mtx)) ## marker基因与表达矩阵基因取交集
  gsva_score <- GSVA::gsva(mtx,list(genes)) ## GSVA
  dt <- survival_info
  dt$gsva_score <- gsva_score[,dt$sample]
  dt$type <- factor(dt$type,levels = c("normal","paratumor","cancer"))
  dt_n <- dt[type=="normal"]
  dt_p <- dt[type=="paratumor"]
  dt_t <- dt[type=="cancer"]
  dt_n <- dt_n[!gsva_score %in% boxplot.stats(dt_n$gsva_score)$out] 
  dt_p <- dt_p[!gsva_score %in% boxplot.stats(dt_p$gsva_score)$out]
  dt_t <- dt_t[!gsva_score %in% boxplot.stats(dt_t$gsva_score)$out]
  ymax <- max(c(dt_n$gsva_score,dt_p$gsva_score,dt_t$gsva_score))*1.2
  ## 箱线图
  p <- ggboxplot(dt,x="type",y="gsva_score",fill = "type",xlab = F,ylab = "GSVA score",
                                  title = x,
                                  bxp.errorbar=T,outlier.shape=NA,ylim=c(min(dt$gsva_score)+0.0001,ymax))+
    theme_classic(base_size = 10)+theme(plot.title = element_text(hjust = 0.5),legend.position = "none")+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.x.npc = 0.5,label.y = ymax*0.95,size=3)
  return(p)
})
p <- cowplot::plot_grid(plotlist = lst,ncol = 4,nrow = 1,align = "v")
save_plot("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_celltype_markers_boxplot.pdf",p,ncol = 4,nrow = 1,base_width=2,base_height=2)

##### violinplot for mucins
mucins <- rownames(sr@assays$RNA@data)[rownames(sr@assays$RNA@data) %like% "^MUC"]
dt <- reshape2::melt(as.matrix(sr@assays$RNA@data[mucins,])) %>% as.data.table()
colnames(dt)[1:2] <- c("gene","cell")
dt$cluster <- sr$subcluster[dt$cell]

ggviolin(dt,'cluster','value',scale='width',fill='cluster',color='transparent')+
  facet_grid(gene~.)+
  xlab('')+ylab('')+
  scale_y_continuous(limit=c(-0.01,5),breaks=c(-0.01,5),labels = c("",5))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        strip.background = element_blank(),legend.position='None',
        strip.text.y = element_text(angle = 0,size=10,hjust=1))

##### violinplot for inhibitory receptors
ir <- c("BTLA","CD200R1","CD22","CD300A","CD300LF",
            "CD33","CD5","CD72","CEACAM1","CLEC12A",
            "CLEC4A","CTLA4","FCGR2B","KLRB1",
            "KLRC1","KLRG1","LAIR1","LILRB1","LILRB2",
            "LILRB4","LILRB5","NCR2","PDCD1","PECAM1",
            "PILRA","PVR","SIGLEC11","SIGLEC5","SIGLEC7",
            "SIGLEC8","SIGLEC9","SIRPA","TIGIT","VSTM1"
) ## from "Functional categories of immune inhibitory receptors"
ir <- ir[ir %in% rownames(sr@assays$RNA@counts)]
dt <- reshape2::melt(as.matrix(sr@assays$RNA@data[ir,])) %>% as.data.table()
colnames(dt)[1:2] <- c("gene","cell")
dt$cluster <- sr$subcluster[dt$cell]

p <- ggviolin(dt,'cluster','value',scale='width',fill='cluster',color='transparent')+
  facet_grid(gene~.)+
  xlab('')+ylab('')+
  scale_y_continuous(limit=c(-0.01,5),breaks=c(-0.01,5),labels = c("",5))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        strip.background = element_blank(),legend.position='None',
        strip.text.y = element_text(angle = 0,size=10,hjust=1))
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_subcluster_inhibitory_receptors_violinplot.pdf",height=13,width=5)

##### violinplot for L-R
lr <- fread("~/project/cholangiocarcinoma/result/t_cell_interaction.tsv",sep="\t",skip=1)
negative_ligand <- lr$ligand[1:24]
negative_ligand <- negative_ligand[negative_ligand %in% rownames(sr@assays$RNA@data)]
dt <- reshape2::melt(as.matrix(sr@assays$RNA@data[negative_ligand,])) %>% as.data.table()
colnames(dt)[1:2] <- c("gene","cell")
dt$cluster <- sr$seurat_clusters[dt$cell]
p <- ggviolin(dt,'cluster','value',scale='width',fill='cluster',color='transparent')+
  facet_grid(gene~.)+
  xlab('')+ylab('')+
  scale_y_continuous(limit=c(-0.01,6),breaks=c(-0.01,6),labels = c("",6))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        strip.background = element_blank(),legend.position='None',
        strip.text.y = element_text(angle = 0,size=10,hjust=1))
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_subcluster_negative_ligand_violinplot.pdf",height=13,width=5)

positive_ligand <- lr$ligand[26:44]
positive_ligand <- positive_ligand[positive_ligand %in% rownames(sr@assays$RNA@data)]
dt <- reshape2::melt(as.matrix(sr@assays$RNA@data[positive_ligand,])) %>% as.data.table()
colnames(dt)[1:2] <- c("gene","cell")
dt$cluster <- sr$seurat_clusters[dt$cell]
p <- ggviolin(dt,'cluster','value',scale='width',fill='cluster',color='transparent')+
  facet_grid(gene~.)+
  xlab('')+ylab('')+
  scale_y_continuous(limit=c(-0.01,5),breaks=c(-0.01,5),labels = c("",5))+
  theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),
        strip.background = element_blank(),legend.position='None',
        strip.text.y = element_text(angle = 0,size=10,hjust=1))
ggsave("/home/rlhua/project/cholangiocarcinoma/result/Epithelial_subcluster_positive_ligands_violinplot.pdf",height=13,width=5)

##### dCCA and pCCA
ca <- sr[,sr$sample=="cancer"]
ca$type <- ifelse(ca$seurat_clusters==4 ,"pCCA","dCCA")
Idents(ca) <- "type"

## marker genes of pCCA and dCCA
cluster_mean_expression=data.frame()
cluster_mean_expression_org=data.frame()
cluster_signal2noise=data.frame()
meta=ca@meta.data
matexpr=as.matrix(ca@assays$RNA@data)
for(clust in sort(unique(meta$type))){
  print(clust)
  cells=rownames(meta)[meta$type==clust]
  othercells=rownames(meta)[meta$type!=clust]
  print(length(cells))
  print(length(othercells))
  exprmat=matexpr[,cells]
  otherexprmat=matexpr[,othercells]
  rowmean=rowMeans(exprmat)
  otherrowmean=rowMeans(otherexprmat)
  sds=apply(exprmat,1,sd)
  othersds=apply(otherexprmat,1,sd)
  signal2noise=(rowmean-otherrowmean)/(sds+othersds)
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
}

cluster_mean_expression <- cluster_mean_expression[,unique(meta$type)]
maxcell=apply(cluster_mean_expression,1,function(x){
  colnames(cluster_mean_expression)[order(x,decreasing = T)[1]]
})

diff.expr.genes2 <- Seurat::FindAllMarkers(ca,only.pos = T,logfc.threshold = 0.25,test.use = "roc",assay="RNA")
diff.expr.genes2$maxcell=maxcell[diff.expr.genes2$gene]
diff.expr.genes2 = diff.expr.genes2[as.character(diff.expr.genes2$cluster)==diff.expr.genes2$maxcell,]
diff.expr.genes2 = diff.expr.genes2[diff.expr.genes2$myAUC>0.7,]
diff.expr.genes2=diff.expr.genes2[ order(diff.expr.genes2$myAUC,decreasing=T),]
diff.expr.genes2=diff.expr.genes2[ !duplicated(diff.expr.genes2$gene),]
diff.expr.genes2=diff.expr.genes2[order(diff.expr.genes2$cluster),]

gsva_score <- GSVA::gsva(mtx[,survival_info$sample[survival_info$type=="cancer" & survival_info$histology %in% c("pCCA","dCCA")]], 
                         list(dCCA=diff.expr.genes2$gene[diff.expr.genes2$cluster=="dCCA"],pCCA=diff.expr.genes2$gene[diff.expr.genes2$cluster=="pCCA"]))

lst <- list()
for (i in c("pCCA","dCCA")) {
  dt <- data.table(sample=colnames(gsva_score),
                   histology=survival_info$histology[survival_info$type=="cancer" & survival_info$histology %in% c("pCCA","dCCA")],
                   value=gsva_score[i,])
  ymax <- max(dt$value)*1.2
  lst[[i]] <- ggboxplot(dt,x="histology",y="value",fill = "histology",xlab = F,ylab = "GSVA score",
                        title = i,
                        bxp.errorbar=T,outlier.shape=NA,ylim=c(min(dt$value),ymax))+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.x.npc = 0.5,size=3,label.y = ymax*0.95)+
    theme(panel.border = element_rect(fill = NA),plot.title = element_text(hjust = 0.5),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),text = element_text(size = 7))
  
}
p <- cowplot::plot_grid(plotlist = lst,ncol = 2,nrow = 1,align = "v")
save_plot("/home/rlhua/project/cholangiocarcinoma/result/Cancer_pCCA2dCCA_GSVA_72samples_boxplot.pdf",p,ncol = 2,nrow = 1,base_width=1.5,base_height=2)

lst <- list()
for (i in mucins) {
  dt <- data.table(sample=survival_info$sample[survival_info$type=="cancer" & survival_info$histology %in% c("pCCA","dCCA")],
                   histology=survival_info$histology[survival_info$type=="cancer" & survival_info$histology %in% c("pCCA","dCCA")],
                   value=mtx[i,survival_info$sample[survival_info$type=="cancer" & survival_info$histology %in% c("pCCA","dCCA")]])
  ymax <- max(dt$value)*1.2
  lst[[i]] <- ggboxplot(dt,x="histology",y="value",fill = "histology",xlab = F,
                        title = i,
                        bxp.errorbar=T,outlier.shape=NA,ylim=c(min(dt$value),ymax))+
    stat_compare_means(aes(label = paste0("p = ", ..p.format..)),label.x.npc = 0.5,size=3,label.y = ymax*0.95)+
    theme(panel.border = element_rect(fill = NA),plot.title = element_text(hjust = 0.5),legend.position = "none",axis.title.x = element_blank(),axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1),text = element_text(size = 7))
}
p <- cowplot::plot_grid(plotlist = lst,ncol = 4,nrow = ceiling(length(lst)/4),align = "v")
save_plot("/home/rlhua/project/cholangiocarcinoma/result/Cancer_pCCA2dCCA_mucins_72samples_boxplot.pdf",p,ncol = 4,nrow = ceiling(length(lst)/4),base_width=1.5,base_height=2)

##### GSVA for hallmark pathways
## hallmark genesets
geneSets <- readLines(system.file("txt", "hallmark-pathways.txt", package = "scCancer")) # 从scCancer包获取hallmark基因集
geneSets <- strsplit(geneSets, "\t")
geneSets <- as.data.frame(geneSets, stringsAsFactors = F)
colnames(geneSets) <- geneSets[1, ]
geneSets <- as.list(geneSets[2, ])
geneSets <- sapply(geneSets, function(x) strsplit(x, ", ")) # 将hallmark基因集整理成列表格式
gene_cnt <- unlist(geneSets) %>% table() # 计算基因在50个基因集中出现次数
pathway_specific_genes <- names(gene_cnt)[gene_cnt==1] # 仅保留每个基因集独有的基因
geneSets <- lapply(geneSets, function(x){x[x %in% pathway_specific_genes]}) # 保留独有基因后的基因集

## GSVA
t.scores <- GSVA::gsva(as.matrix(sr@assays$RNA@data), geneSets)
saveRDS(t.scores,"/home/rlhua/project/cholangiocarcinoma/result/Epithelial_hallmark_GSVA.RDS")

## cancer cells
# limma统计检验
lst <- list()
for (i in unique(sr$seurat_clusters[sr$sample=="cancer"])) {
  group_list <- data.table(cell=colnames(sr)[sr$sample=="cancer"],
                           subcluster=paste0("C",sr$seurat_clusters[sr$sample=="cancer"]),
                           patient=gsub("_.*","",sr$orig.ident[sr$sample=="cancer"]))
  group_list$subcluster <- ifelse(group_list$subcluster==paste0("C",i),group_list$subcluster,"else")
  design <- model.matrix(~0+subcluster+patient,group_list)
  rownames(design) <- group_list$cell
  contrast.matrix <- makeContrasts(contrasts = paste0("subclusterC",i,"-subclusterelse"), levels = design)
  fit <- lmFit(t.scores[,group_list$cell], design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  toptable <- topTable(fit2, coef= 1 ,number=Inf,adjust.method = "BH", sort.by="P")
  toptable[["pathway"]] <- rownames(toptable)
  toptable[["subcluster"]] <- i
  lst[[i]] <- toptable[,c("pathway","t","P.Value","adj.P.Val","subcluster")]
}
dt <- rbindlist(lst)

# heatmap
dt <- dt[pathway %in% unique(dt$pathway[dt$adj.P.Val<0.05])]
dt <- as.data.frame(dt)
mtx <- reshape2::dcast(dt,pathway~subcluster,fun.aggregate = mean,value.var = "t")
rownames(mtx) <- mtx[,"pathway"]
mtx <- mtx[,-1]
mtx <- as.matrix(mtx)
pheatmap(mtx,
         color =colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=7,
         fontsize_row=7,
         scale="none",
         show_colnames=T,
         show_rownames = T,
         fontsize_col=10,
         main = "Cluster:",
         cellheight = 10,cellwidth = 18,
         treeheight_row = 0,treeheight_col = 0,
         clustering_method = "ward.D2",
         filename = "/home/rlhua/project/cholangiocarcinoma/result/cancer_subcluster_hallmark_GSVA_pheatmap.pdf")

## cholangiocyte
# limma统计检验
lst <- list()
for (i in unique(sr$seurat_clusters[sr$sample=="cholangiocyte"])) {
  group_list <- data.table(cell=colnames(sr)[sr$sample=="cholangiocyte"],
                           subcluster=paste0("C",sr$seurat_clusters[sr$sample=="cholangiocyte"]),
                           patient=gsub("_.*","",sr$orig.ident[sr$sample=="cholangiocyte"]))
  group_list$subcluster <- ifelse(group_list$subcluster==paste0("C",i),group_list$subcluster,"else")
  design <- model.matrix(~0+subcluster+patient,group_list)
  rownames(design) <- group_list$cell
  contrast.matrix <- makeContrasts(contrasts = paste0("subclusterC",i,"-subclusterelse"), levels = design)
  fit <- lmFit(t.scores[,group_list$cell], design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  toptable <- topTable(fit2, coef= 1 ,number=Inf,adjust.method = "BH", sort.by="P")
  toptable[["pathway"]] <- rownames(toptable)
  toptable[["subcluster"]] <- i
  lst[[i]] <- toptable[,c("pathway","t","P.Value","adj.P.Val","subcluster")]
}
dt <- rbindlist(lst)

# heatmap
dt <- dt[pathway %in% unique(dt$pathway[dt$adj.P.Val<0.05])]
dt <- as.data.frame(dt)
mtx <- reshape2::dcast(dt,pathway~subcluster,fun.aggregate = mean,value.var = "t")
rownames(mtx) <- mtx[,"pathway"]
mtx <- mtx[,-1]
mtx <- as.matrix(mtx)
pheatmap(mtx,
         color =colorRampPalette(c("blue", "white", "red"))(100),
         cluster_cols =F,
         cluster_rows = T,
         fontsize=7,
         fontsize_row=7,
         scale="none",
         show_colnames=T,
         show_rownames = T,
         fontsize_col=10,
         main = "Cluster:",
         cellheight = 10,cellwidth = 18,
         treeheight_row = 0,treeheight_col = 0,
         clustering_method = "ward.D2",
         filename = "/home/rlhua/project/cholangiocarcinoma/result/cholangiocyte_subcluster_hallmark_GSVA_pheatmap.pdf")

##### comparison of cancer and cholangiocyte
tbar<-function(df,cutoff=3.27,title=""){
  colnames(df)=c('ID','score')
  #cutoff <- 3.27
  df$group <- droplevels(cut(df$score, breaks = c(-Inf, -cutoff, cutoff, Inf),labels = c(1,2,3)))
  colvec=c('palegreen3', 'snow3', 'dodgerblue4')
  names(colvec)=1:3
  colpan=colvec[levels(df$group)]
  coltext=c("black","snow3","black")
  names(coltext)=1:3
  coltextpan=coltext[levels(df$group)]
  #按照score排序
  sortdf <- df[order(df$score),]
  sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
  #head(sortdf)
  
  ggplot(sortdf, aes(ID, score, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = colpan, guide = FALSE) + 
  
  #画2条虚线.
  geom_hline(yintercept = c(-cutoff,cutoff), 
             color="white",
             linetype = 2, #画虚线
             size = 0.3) + #线的粗细
    
  #写label
  geom_text(inherit.aes = F,data = subset(df, score < 0),
            aes(x=ID, y= 0.2, label= paste0(" ", ID), color = group),#bar跟坐标轴间留出间隙
            size = 1.5, #字的大小
            hjust = 0 ) +  #字的对齐方式
    geom_text(inherit.aes = F,data = subset(df, score > 0),
              aes(x=ID, y= -0.2, label=ID, color = group),
              size = 1.5, hjust = 1) +  
  scale_colour_manual(values = coltextpan, guide = FALSE) +
  xlab("") +ylab("t value of GSVA score, tumor \n versus non-malignant")+
  theme_bw() + #去除背景色
    theme(panel.grid =element_blank()) + #去除网格线
    theme(panel.border = element_rect(size = 0.6)) + #边框粗细
    theme(plot.title = element_text(hjust = .5),axis.title = element_text(size=7),text = element_text(size = 7),axis.line.y = element_blank(), axis.ticks.y = element_blank(), axis.text.y = element_blank()) + #去除y轴
    ggtitle(title)
} ## 绘图函数

cells <- colnames(sr)[sr$sample %in% c("cancer","cholangiocyte")]
group_list <- data.table(cell=cells,
                         tissue=sr$sample[cells],
                         patient=gsub("_.*","",sr$orig.ident[cells]))
design <- model.matrix(~0+tissue+patient,group_list)
rownames(design) <- group_list$cell
contrast.matrix <- makeContrasts(tissuecancer-tissuecholangiocyte, levels = design)
fit <- lmFit(t.scores[,cells], design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
toptable <- topTable(fit2, coef= 1 ,number=Inf,adjust.method = "BH", sort.by="P")
toptable[["pathway"]] <- rownames(toptable)
df <- toptable[,c("pathway","t")]
cutoff <- min(abs(toptable[toptable$adj.P.Val<0.05,"t"]))-0.00001
p <- tbar(df=df,cutoff = cutoff)+ylab("t value of GSVA score, cancer\nversus cholangiocyte")
ggsave("/home/rlhua/project/cholangiocarcinoma/result/cancerVScholangiocyte_hallmark_GSVA_barplot.pdf",width = 2.5,height = 4)

