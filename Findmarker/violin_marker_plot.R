##### 亚群markers
## markers by RNA assay 
fs <- list.files("/home/rlhua/project/cholangiocarcinoma/result/new","^clustermarker2.RDS",all.files = T,full.names = T,recursive = T)
lst <- list()
for (i in 1:length(fs)) {
  df <- readRDS(fs[i])
  df[["celltype"]] <- gsub("/home/rlhua/project/cholangiocarcinoma/result/new/RNAmarkers_","",fs[i]) %>% gsub("/clustermarker2.RDS","",.)
  lst[[i]] <- df
}
df_RNA <- rbindlist(lst)
df_RNA$cluster <- paste(df_RNA$celltype,df_RNA$cluster)
df_RNA$cluster <- gsub("_"," ",df_RNA$cluster)

## markers by integrated assay 
fs <- list.files("/home/rlhua/project/cholangiocarcinoma/result/new","^M3_clustermarker2.RDS",all.files = T,full.names = T,recursive = T)
lst <- list()
for (i in 1:length(fs)) {
  df <- readRDS(fs[i])
  df[["celltype"]] <- gsub("/home/rlhua/project/cholangiocarcinoma/result/new/M3markers_","",fs[i]) %>% gsub("/M3_clustermarker2.RDS","",.)
  lst[[i]] <- df
}
df_integrated <- rbindlist(lst)
df_integrated$cluster <- paste(df_integrated$celltype,df_integrated$cluster)
df_integrated$cluster <- gsub("_"," ",df_integrated$cluster)

## violinplot of markers
# integrated
for (i in sort(unique(df_integrated$cluster))) {
  markers <- df_integrated$gene[df_integrated$cluster==i]
  ct <- unique(df_integrated$celltype[df_integrated$cluster==i])
  sr <- readRDS(paste0("/home/rlhua/project/cholangiocarcinoma/result/new/seurat_RDSdata/",ct,".RDS"))
  VlnPlot(sr,markers,pt.size=0,group.by="seurat_clusters",ncol=8,assay="integrated")
  fn <- paste0("/home/rlhua/project/cholangiocarcinoma/result/new/clustermarkers/",i,"_integrated_violinplot.pdf")
  ggsave(fn,width=24,height=ceiling(length(markers)/8)*3,limitsize = FALSE)
  avgexp <- apply(sr@assays$integrated@data[markers,,drop=F], 1, function(x) mean(x[x>0]))
  markers <- markers[avgexp>1]
  if (length(markers)==0) {next}
  VlnPlot(sr,markers,pt.size=0,group.by="seurat_clusters",ncol=8,assay="integrated")
  fn <- paste0("/home/rlhua/project/cholangiocarcinoma/result/new/clustermarkers/",i,"_integrated_highexpression_violinplot.pdf")
  ggsave(fn,width=24,height=ceiling(length(markers)/8)*3,limitsize = FALSE)
}

# RNA
for (i in sort(unique(df_RNA$cluster))) {
  markers <- df_RNA$gene[df_RNA$cluster==i]
  ct <- unique(df_RNA$celltype[df_RNA$cluster==i])
  sr <- readRDS(paste0("/home/rlhua/project/cholangiocarcinoma/result/new/seurat_RDSdata/",ct,".RDS"))
  VlnPlot(sr,markers,pt.size=0,group.by="seurat_clusters",ncol=8,assay="RNA")
  fn <- paste0("/home/rlhua/project/cholangiocarcinoma/result/new/clustermarkers/",i,"_RNA_violinplot.pdf")
  ggsave(fn,width=24,height=ceiling(length(markers)/8)*3,limitsize = FALSE)
}
