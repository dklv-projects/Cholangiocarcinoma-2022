#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <包含表达谱的数据框RDS>
#功能：所有类的marker基因进行基因集活性分析
#input：包含表达谱的数据框,行名为gene symbol,列名为样本
#output：
#for fl in `ls `;do echo $fl;cd $fl;Rscript ../../script/M4Gsva.R M3_cluster_mean_expression.RDS;cd ../;done

Args<-commandArgs()
print(Args)
expr.file=normalizePath(Args[6])
cmdfile=strsplit(Args[4], split="=" )[[1]][2]
cmddir=dirname(normalizePath(cmdfile))
print(cmddir)
pwd=getwd()
savePath=paste0(pwd,"/M4gsva_",basename(Args[6]))
if(!dir.exists(savePath)){
  dir.create(savePath,recursive = T)
}
#将结果保存文件夹设置为工作目录
setwd(savePath)

gsvaKEGGfile=paste0(cmddir,"/human_geneset_activity/KEGG_GSVA.R")
system(paste("Rscript",gsvaKEGGfile,expr.file))
gsvaHallmarkfile=paste0(cmddir,"/human_geneset_activity/Hallmark_GSVA.R")
system(paste("Rscript",gsvaHallmarkfile,expr.file))
gsvaGOfile=paste0(cmddir,"/human_geneset_activity/GO_GSVA.R")
system(paste("Rscript",gsvaGOfile,expr.file))
gsvaReactomefile=paste0(cmddir,"/human_geneset_activity/Reactome_GSVA.R")
system(paste("Rscript",gsvaReactomefile,expr.file))
