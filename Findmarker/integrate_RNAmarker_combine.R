#!/usr/bin/env Rscript
options(warn=-1)
#USAGE:Rscript cmdfile <包含integrate marker的数据框RDS> <包含RNA marker的数据框RDS>
#功能：合并RDS
#input：包含marker的数据框,至少包含gene和cluster两列
#output：


Args<-commandArgs(T)
print(Args)
markerdf.integrate=readRDS(Args[1])
markerdf.integrate=markerdf.integrate[,1:9]
markerdf.RNA=readRDS(Args[2])
integrate.gene=unique(markerdf.integrate$gene)
markerdf.RNA=markerdf.RNA[!(markerdf.RNA$gene %in% integrate.gene),]
markerdf=rbind(markerdf.integrate,markerdf.RNA)
markerdf=markerdf[order(markerdf$cluster),]
saveRDS(markerdf,"marker2.RDS")
