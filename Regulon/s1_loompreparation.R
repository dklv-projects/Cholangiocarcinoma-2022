#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <SeuratObject.RDS> <number>
#功能：将seurat对象转换成loom，提取部分(默认5000)细胞做成loom
#input：seurat对象，来源物种必须是人
#output：loom文件query.loom
#Rscript s1_loompreparation.R ../../M2classification/M2Schwann.RDS 5000

#reference.loom是关键，第一种，使用所有细胞做，
#                      第二种，抽样，（此脚本用法）
#                      第三种，根据分群，1:1，
#                      第四种，分群后取一群做
Args<-commandArgs()
print(Args)
cmdfile=strsplit(Args[4], split="=" )[[1]][2]
cmddir=dirname(normalizePath(cmdfile))
print(cmddir)#命令文件所在目录
expr=readRDS(Args[6])
num=as.integer(Args[7])
print(num)
pwd=getwd()
#保存结果目录
savePath=paste0(pwd,"/M3regulon")
savefile=paste0(basename(Args[6]),".loom")
print(savefile)
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}
setwd(savePath)
## loom文件的生成
#input：Seurat对象或者单细胞表达谱
#output：loom
##第一步，读取seurat对象或单细胞表达谱，写入expr
options(stringsAsFactors = FALSE)
library(plyr)
library(permute)
library(data.table)
library(SCopeLoomR)
library(SCENIC)
library(Seurat)

#这里使用的是样本整合后的数据
expr=as.matrix(expr@assays$integrated@data)

##第二步，过滤基因后写入loom文件
#参数初始化
scenicOptions <- initializeScenic(
  org="hgnc", # 物种名 or hgnc, or dmel
  dbDir=file.path(cmddir,"cisTarget_databases_hg38"), # RcisTarget databases location
  dbs="hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather", # file na+ datasetTitle="SCENIC on Human Cells"
  nCores=1
)
#过滤基因
genesKept <- geneFiltering(
  exprMat=expr, 
  scenicOptions=scenicOptions,
  minCountsPerGene = 1,
  minSamples = 50
)
dge <- expr[genesKept, ]

#此外，抽取部分细胞预测regulon，这里保留5000个细胞
if(ncol(dge)>num){
  set.seed(666)
  colsample=sample(1:ncol(dge),num)
  dge.sample=dge[,colsample]
}else(dge.sample=dge)

loom <- build_loom(savefile, dgem=dge.sample)
close_loom(loom)


