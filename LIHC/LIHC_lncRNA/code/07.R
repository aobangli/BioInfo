# nomogram 绘制列线图函数
# f 上述已构建好的COX比例风险回归模型
# fun 将预测的生存概率转化为线性，用于绘制新的生存概率坐标轴，
# 若要绘制2条及以上生存概率坐标轴，则用list函数连接；
# 需要根据数据的时间单位来设置预测时间点的生存概率
# fun.at 设置生存概率坐标轴的范围及步长
# funlabel 设置生存概率坐标轴的名称


# 加载与匹配 临床数据和表达数据
rm(list = ls())
clinical= read.csv("./clinical_change.csv")
rownames(clinical) = clinical[,1]

load("./survexprdata_geneCoef.Rdata")

clinicalexpr = log2(survexprdata[, geneCoef[,1] ]+1)
clinicalgene = colnames(clinicalexpr)

rownames(clinicalexpr) = substring(rownames(clinicalexpr),1,12) 
clinicalexpr = clinicalexpr[rownames(clinical),]
clinicalexpr = na.omit(clinicalexpr)
clinical = clinical[rownames(clinicalexpr),]
clinicalexpr = cbind(clinicalexpr,clinical)    #临床特征与表达数据merge

# 计算risk score得分
clinicalexpr$Riskscore = 0
for(i in 1:nrow(geneCoef)) {
        clinicalexpr$Riskscore = clinicalexpr$Riskscore + clinicalexpr[, i] * as.numeric(geneCoef[i, 2])
}

clinicalexpr$OS = as.numeric(clinicalexpr$OS) + 1
clinicalexpr$OS.time = as.numeric(clinicalexpr$OS.time)

nomomatrix = clinicalexpr[,c('agegroup','gender','pM','pN1','pT','pT1','stage1','Riskscore','OS','OS.time')]
nomomatrix$stage1 = factor(nomomatrix$stage1)
nomomatrix$stage2 = gsub('1','Stage I',nomomatrix$stage1)
nomomatrix$stage2 = gsub('2','Stage II',nomomatrix$stage2)
nomomatrix$stage2 = gsub('3','Stage III',nomomatrix$stage2)
nomomatrix$stage2 = gsub('4','Stage IV',nomomatrix$stage2)

nomomatrix$stage2 = factor(nomomatrix$stage2)


library(survival)
library(rms)


nomomatrix = nomomatrix[,c('agegroup','stage2','Riskscore','OS','OS.time')]
colnames(nomomatrix) = c('Age','Stage','RiskScore','OS','OS.time')
dd=datadist(nomomatrix)
options(datadist="dd") 


# 构建cox模型
library(survival)
library(rms)
f = cph(Surv(OS.time,OS) ~ Age +Stage +RiskScore  ,data=nomomatrix,
        x=T,y=T,surv = T)

#画图
survival = Survival(f)
survival1 = function(x)survival(365,x)
survival2 = function(x)survival(1095,x)
survival3 = function(x)survival(1825,x)


nom = nomogram(f,fun = list(survival1,survival2,survival3),lp=F,
               fun.at = c(0.1,seq(0.1,0.9,by=0.1),0.95),
               funlabel = c('1 year survival' ,'3 year survival','5 year survival'))
pdf(file="nomo.pdf",width=45,height=25,pointsize = 60)
plot(nom, xfrac=.2)
dev.off()


###### 内部验证 #####
### 1. 区分度 (discrimination)
validate(f, method="boot", B=1000, dxy=T)

rcorrcens(Surv(OS.time,OS) ~ predict(f), data = nomomatrix)


library(timeROC)
library(survival)
#### 1year ROC  ####
nomoROC =  clinicalexpr[,c('agegroup','gender','pM','pN1','pT','pT1','stage1','Riskscore','OS','OS.time')]
nomoROC$predict = predict(f)
nomoROC$OS = as.numeric(nomoROC$OS) - 1  #将event改为0和1

ROC_rt_stage <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$stage1,
                         cause=1,
                         weighting='marginal',
                         times= 365,
                         ROC=TRUE)

ROC_rt_score <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$Riskscore,
                         cause=1,
                         weighting='marginal',
                         times= 365,
                         ROC=TRUE)

ROC_rt_predict <<- timeROC(T=nomoROC$OS.time,
                           delta=nomoROC$OS,   #时间 结局
                           marker=nomoROC$predict,
                           cause=1,
                           weighting='marginal',
                           times= 365,
                           ROC=TRUE)


pdf(file="Fig5D1_1year.pdf",width = 10,height = 11,pointsize = 25)

plot(ROC_rt_stage,time=365,col='black',title=FALSE,lwd=3)+box(lwd=4)
par(new=TRUE)
plot(ROC_rt_score,time=365,col='blue',title=FALSE,lwd=3)
par(new=TRUE)
plot(ROC_rt_predict,time=365,col='red',title=FALSE,lwd=3)


legend('bottomright',
       c(paste0('stage: ',round(ROC_rt_stage$AUC[2],3)),
         paste0('score: ',round(ROC_rt_score$AUC[2],3)),
         paste0('nomogram: ',round(ROC_rt_predict$AUC[2],3))
       ),
       col=c('black','blue','red'),lwd=4,bty = 'n')

dev.off()


#### 3year ROC  ####
nomoROC =  clinicalexpr[,c('agegroup','gender','pM','pN1','pT','pT1','stage1','Riskscore','OS','OS.time')]
nomoROC$predict = predict(f)
nomoROC$OS = as.numeric(nomoROC$OS) - 1  #将event改为0和1

ROC_rt_stage <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$stage1,
                         cause=1,
                         weighting='marginal',
                         times= 1095,
                         ROC=TRUE)

ROC_rt_score <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$Riskscore,
                         cause=1,
                         weighting='marginal',
                         times= 1095,
                         ROC=TRUE)

ROC_rt_predict <<- timeROC(T=nomoROC$OS.time,
                           delta=nomoROC$OS,   #时间 结局
                           marker=nomoROC$predict,
                           cause=1,
                           weighting='marginal',
                           times= 1095,
                           ROC=TRUE)


pdf(file="Fig5E1_3year.pdf",width = 10,height = 11,pointsize = 25)

plot(ROC_rt_stage,time=1095,col='black',title=FALSE,lwd=3)+box(lwd=4)
par(new=TRUE)
plot(ROC_rt_score,time=1095,col='blue',title=FALSE,lwd=3)
par(new=TRUE)

plot(ROC_rt_predict,time=1095,col='red',title=FALSE,lwd=3)


legend('bottomright',
       c(paste0('stage: ',round(ROC_rt_stage$AUC[2],3)),
         paste0('score: ',round(ROC_rt_score$AUC[2],3)),
         paste0('nomogram: ',round(ROC_rt_predict$AUC[2],3))
       ),
       col=c('black','blue','red'),lwd=4,bty = 'n')

dev.off()

#### 5year ROC  ####
nomoROC =  clinicalexpr[,c('agegroup','gender','pM','pN1','pT','pT1','stage1','Riskscore','OS','OS.time')]
nomoROC$predict = predict(f)
nomoROC$OS = as.numeric(nomoROC$OS) - 1  #将event改为0和1

ROC_rt_stage <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$stage1,
                         cause=1,
                         weighting='marginal',
                         times= 1825,
                         ROC=TRUE)

ROC_rt_score <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$Riskscore,
                         cause=1,
                         weighting='marginal',
                         times= 1825,
                         ROC=TRUE)

ROC_rt_predict <<- timeROC(T=nomoROC$OS.time,
                           delta=nomoROC$OS,   #时间 结局
                           marker=nomoROC$predict,
                           cause=1,
                           weighting='marginal',
                           times= 1825,
                           ROC=TRUE)


pdf(file="Fig5F1_5year.pdf",width = 10,height = 11,pointsize = 25)

plot(ROC_rt_stage,time=1825,col='black',title=FALSE,lwd=3)+box(lwd=4)
par(new=TRUE)
plot(ROC_rt_score,time=1825,col='blue',title=FALSE,lwd=3)
par(new=TRUE)

plot(ROC_rt_predict,time=1825,col='red',title=FALSE,lwd=3)


legend('bottomright',
       c(paste0('stage: ',round(ROC_rt_stage$AUC[2],3)),
         paste0('score: ',round(ROC_rt_score$AUC[2],3)),
         paste0('nomogram: ',round(ROC_rt_predict$AUC[2],3))
       ),
       col=c('black','blue','red'),lwd=4,bty = 'n')

dev.off()

### 2. 一致性 (calibration)
# 1-year survival calibration curve
f1 <- cph(Surv(OS.time,OS) ~ Age +Stage + RiskScore, x=T, y=T, surv=T, data=nomomatrix, time.inc=365)

cal1 <- calibrate(f1, cmethod="KM", method="boot", u=365, m=50, B=181)

pdf(file="1yearnomo.pdf",width=15,height=15,pointsize = 32)
plot(cal1,xlim=c(0,1),ylim=c(0,1),xlab = 'Nomogram-predicted probability of 1-year survival',
     ylab = 'Observed 1-year survival',lwd=6,
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=0.5)+box(lwd=4)
abline(a = 0, b = 1, col = "gray50",lwd=3,lty=2)
dev.off()

# 3-year survival calibration curve
f3 <- cph(Surv(OS.time,OS) ~ Age +Stage + RiskScore, x=T, y=T, surv=T, data=nomomatrix, time.inc=1095)

cal3 <- calibrate(f3, cmethod="KM", method="boot", u=1095, m=50, B=181)

pdf(file="3yearnomo.pdf",width=15,height=15,pointsize = 32)
plot(cal3,xlim=c(0,1),ylim=c(0,1),xlab = 'Nomogram-predicted probability of 3-year survival',
     ylab = 'Observed 3-year survival',lwd=6,
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=0.5)+box(lwd=4)
abline(a = 0, b = 1, col = "gray50",lwd=3,lty=2)
dev.off()

# 5-year survival calibration curve
f5 <- cph(Surv(OS.time,OS) ~ Age +Stage + RiskScore, x=T, y=T, surv=T, data=nomomatrix, time.inc=1825)

cal5 <- calibrate(f5, cmethod="KM", method="boot", u=1825, m=50, B=181)

pdf(file="5yearnomo.pdf",width=15,height=15,pointsize = 32)
plot(cal5,xlim=c(0,1),ylim=c(0,1),xlab = 'Nomogram-predicted probability of 5-year survival',
     ylab = 'Observed 5-year survival',lwd=6,
     cex.lab=1.2, cex.axis=1.2, cex.main=1.2, cex.sub=0.5)+box(lwd=4)
abline(a = 0, b = 1, col = "gray50",lwd=3,lty=2)
dev.off()

