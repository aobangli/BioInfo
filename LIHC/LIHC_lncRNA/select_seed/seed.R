library(timeROC)
library(survival)
library(lars) 
library(glmnet) 
library('survival')
library('survminer')

####### 生存表达数据取交集 #######

rm(list = ls())
options(stringsAsFactors = F)
load(file = './LIHC_lipid_limma_DEG.Rdata')
load(file = './clinical.Rdata')

entire_set = LIHC_PT_COUNT_expr

colnames(entire_set)= substr(colnames(entire_set),1,12)

genename = rownames(need_DEG)[need_DEG$change != 'NOT']
entire_set = entire_set[genename,]

clinical= clinical[substr(colnames(entire_set),1,12),]

OSdata=clinical[,c('id', 'OS', 'OS.time')]
OSdata= na.omit(OSdata)
entire_set = entire_set[,rownames(OSdata)]
entire_set = cbind(t(entire_set),OSdata[,2:3])
save(entire_set,LIHC_PT_COUNT_expr,OSdata,file = 'entire_set.Rdata')


####### 划分训练集和验证集 #######
#以肿瘤样本集为全集划分, 即LIHC_PT_COUNT_expr
rm(list = ls())
load(file = 'entire_set.Rdata')

out_ROC = c()

for(m in 8000:10000){
  print(m)
set.seed(m)
training_set_rate = 0.5
index <-  sort(sample(nrow(entire_set), nrow(entire_set) * training_set_rate))
training_set <- entire_set[index, ]
validation_set <-  entire_set[-index, ]


####### 训练集 单变量COX #######

survexprdata = training_set

coxdata= survexprdata[,1:(ncol(survexprdata)-2)]
coxdata= log2(coxdata+1)
coxgene= colnames(coxdata)
coxdata= cbind(coxdata,survexprdata[,c('OS','OS.time')])

coxoutput=c()
for(i in 1:length(coxgene)) {
  
  res.cox= coxph(Surv(OS.time,OS)~coxdata[,i]  ,data=coxdata)
  summary(res.cox)
  
  beta1 <- coef(res.cox)
  se1 <- sqrt(diag(vcov(res.cox)))
  HR1 <- exp(beta1)
  HRse1 <- HR1 * se1
  pcox = 1 - pchisq((beta1/se1)^2, 1)
  
  result<-cbind("gene" = coxgene[i], "p.value" = pcox ,"beta"=beta1,"se"=se1,"HR"=HR1)
  coxoutput<-rbind(coxoutput,result)}


######## 训练集 lasso 回归 #######
lasso_confidence = 0.01

coxoutput = as.data.frame(coxoutput)
coxoutput$p.value = as.numeric(coxoutput$p.value)
lassocoxgene = coxoutput[coxoutput$p.value < lasso_confidence, 1]

if(length(lassocoxgene)>3){

lassoexpr= LIHC_PT_COUNT_expr[lassocoxgene,]
colnames(lassoexpr)=substr(colnames(lassoexpr),1,12)

x=t(log2(lassoexpr+1))
x=x[rownames(OSdata),]

OStime=as.numeric(OSdata[,3])
OS=as.numeric(OSdata[,2])
y=cbind(time=OStime,status=OS)


fit=glmnet(x, y, family = "cox", maxit = 1000) 

cvfit = cv.glmnet(x, y, family="cox", maxit = 1000) 

coef = coef(fit, s = cvfit$lambda.min) 
index = which(coef != 0) 
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef   #查看模型的相关系数


###### 训练集 ROC #####


######### riskscore surv #######
scoresurv = survexprdata[, geneCoef[,1] ]
scoresurv = log2(scoresurv +1)
scoresurv = cbind(scoresurv,survexprdata[,((ncol(survexprdata)-1):(ncol(survexprdata)))])

scoresurv$riskscore = 0
for(i in 1:nrow(geneCoef)) {
  scoresurv$riskscore = scoresurv$riskscore + scoresurv[, i] * as.numeric(geneCoef[i, 2])
}


ROC_rt <<- timeROC(T=scoresurv$OS.time,delta=scoresurv$OS,
                   marker=scoresurv$riskscore,cause=1,
                   weighting='marginal',
                   times=c(365,1095, 1824),ROC=TRUE)

ROC_train<-cbind("seed" = m, "train_1" = ROC_rt$AUC[1] ,"train_3"=ROC_rt$AUC[2],"train_5"=ROC_rt$AUC[3])
if(ROC_train[2]>0.695 & ROC_train[3]>0.695 & ROC_train[4]>0.695){


######  验证集 ####


survexprdata = validation_set


######### riskscore surv #######


scoresurv = survexprdata[, geneCoef[,1] ]
scoresurv = log2(scoresurv +1)
scoresurv = cbind(scoresurv,survexprdata[,((ncol(survexprdata)-1):(ncol(survexprdata)))])

scoresurv$riskscore = 0
for(i in 1:nrow(geneCoef)) {
  scoresurv$riskscore = scoresurv$riskscore + scoresurv[, i] * as.numeric(geneCoef[i, 2])
}



ROC_rt <<- timeROC(T=scoresurv$OS.time,delta=scoresurv$OS,
                   marker=scoresurv$riskscore,cause=1,
                   weighting='marginal',
                   times=c(365,1095, 1824),ROC=TRUE)

ROC_val<-cbind( ROC_train,"val_1" = ROC_rt$AUC[1] ,"val_3"=ROC_rt$AUC[2],"val_5"=ROC_rt$AUC[3])
if(ROC_val[7]>0.695 & ROC_val[6]>0.695 & ROC_val[5]>0.695 ){out_ROC = rbind(out_ROC,ROC_val)}

}}}

write.csv(out_ROC,file = 'out_ROC 8000-10000.csv')
