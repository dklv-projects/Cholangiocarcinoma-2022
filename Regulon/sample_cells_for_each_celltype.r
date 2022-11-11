setwd("/home/rlhua/tmp/RDSdata")
colname_type=c()
for(fl in list.files(pattern="RDS")){
  print(fl)
  akind=readRDS(fl)
  #print(tail(colnames(akind)))
  type1=sub("M2","",fl)
  type2=sub("1\\.RDS","",type1)
  vec=rep(type2,ncol(akind))
  names(vec)=colnames(akind)
  print(type2)
  colname_type=c(colname_type,vec)
}

expr <- readRDS("/home/rlhua/project/cholangiocarcinoma/result/expr.RDS")
new_expr=expr[,names(colname_type)]
new_expr$celltype=colname_type
saveRDS(new_expr,file = "/home/rlhua/project/cholangiocarcinoma/result/expr20221111.RDS")

sample_types=c()
for(type in unique(colname_type)){
  print(type)
  col_type=colname_type[colname_type==type]
  sample_col=sample(names(col_type),700)
  sample_type=col_type[sample_col]
  sample_types=c(sample_types,sample_type)
}
sample_expr=new_expr[,names(sample_types)]
saveRDS(sample_expr,file = "/home/rlhua/project/cholangiocarcinoma/result/expr_sample20221111.RDS")