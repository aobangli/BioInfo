rm(list=ls())

############### 1.下载GEO数据  ############
library(GEOquery)
eSet <- getGEO('GSE14520', destdir=".", AnnotGPL = F, getGPL = F)

b = eSet[[1]]
c = eSet[[2]]
raw_exprSet = exprs(b) 
phe = pData(b)
exprSet = raw_exprSet

############### 2.ID转换  ################
##### GPL96 ####
library(hgu133a.db)
ids=toTable(hgu133aSYMBOL)

#### GPL570 ####
library(hgu133plus2.db) 
ids=toTable(hgu133plus2SYMBOL)

#### 非GPL96和GPL570 的id转换 ####
library(Biobase)
library(GEOquery)
gpl <- getGEO('GPL3921', destdir=".")

write.csv(Table(gpl), 'GPL3921.csv')

ids=Table(gpl)[,c(1,11)]
head(ids)     ## you need to check this , which column do you need


##### 转换  #####
colnames(ids)=c('probe_id','symbol')
exprSet=exprSet[rownames(exprSet) %in% ids$probe_id,]

ids=ids[match(rownames(exprSet),ids$probe_id),]

jimmy <- function(exprSet,ids){
  tmp = by(exprSet, ids$symbol,
           function(x) rownames(x)[which.max(rowMeans(x))] )
  probes = as.character(tmp)
  print(dim(exprSet))
  exprSet=exprSet[rownames(exprSet) %in% probes ,]
  
  print(dim(exprSet))
  rownames(exprSet)=ids[match(rownames(exprSet),ids$probe_id),2]
  return(exprSet)
}

geo_exprSet <- jimmy(exprSet,ids)

write.csv(geo_exprSet, 'gene_expr.csv')
write.csv(phe, 'phenotype.csv')

save(geo_exprSet, file = 'gene_expr.Rdata')
save(phe, file = 'phenotype.Rdata')
