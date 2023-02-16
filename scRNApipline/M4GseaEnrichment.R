#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <包含signal2noise的数据框RDS>
#功能：所有类的marker基因进行功能富集分析
#input：包含marker的数据框,至少包含gene和cluster两列
#output：


Args<-commandArgs()
print(Args)
signal2noise=readRDS(Args[6])
cmdfile=strsplit(Args[4], split="=" )[[1]][2]
cmddir=dirname(normalizePath(cmdfile))
print(cmddir)
pwd=getwd()
savePath=paste0(pwd,"/M4gseaEnrichment_",basename(Args[6]))
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}
#将结果保存文件夹设置为工作目录
setwd(savePath)

for(colume in colnames(signal2noise)){
  sig2noi=signal2noise[,colume]
  names(sig2noi)=rownames(signal2noise)
  filename=paste0("sigal2noise_",colume,".RDS")
  saveRDS(sig2noi,file = filename)
  GOenrichfile=paste0(cmddir,"/human_functional_class_scoring/GO_GSEA.R")
  system(paste("Rscript",GOenrichfile,filename))
  KEGGenrichfile=paste0(cmddir,"/human_functional_class_scoring/KEGG_GSEA.R")
  system(paste("Rscript",KEGGenrichfile,filename))
  Hallmarkenrichfile=paste0(cmddir,"/human_functional_class_scoring/Hallmark_GSEA.R")
  system(paste("Rscript",Hallmarkenrichfile,filename))
  #Reactomeenrichfile=paste0(cmddir,"/human_functional_class_scoring/Reactome_GSEA.R")
  #system(paste("Rscript",Reactomeenrichfile,filename))
  system(paste("rm",filename))
}
