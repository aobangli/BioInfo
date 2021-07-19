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
clinical = clinical[ , -1]

load("./survexprdata_lipid.Rdata")

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
survival1 = function(x)survival(1095,x)
survival2 = function(x)survival(1825,x)


nom = nomogram(f,fun = list(survival1,survival2),lp=F,
               fun.at = c(0.1,seq(0.1,0.9,by=0.1),0.95),
               funlabel = c('3 year survival','5 year survival'))
pdf(file="nomo.pdf",width=45,height=25,pointsize = 60)
plot(nom, xfrac=.2)
dev.off()


###### 内部验证 #####
### 1. 区分度 (discrimination)
validate(f, method="boot", B=1000, dxy=T)

rcorrcens(Surv(OS.time,OS) ~ predict(f), data = nomomatrix)


library(timeROC)
library(survival)
#### 5year ROC  ####
nomoROC =  clinicalexpr[,c('agegroup','gender','pM','pN1','pT','pT1','stage1','Riskscore','OS','OS.time')]
nomoROC$predict = predict(f)
nomoROC$OS = as.numeric(nomoROC$OS) - 1  #将event改为0和1

ROC_rt_stage <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$stage1,
                         cause=1,
                         weighting='marginal',
                         times= 1760,
                         ROC=TRUE)



ROC_rt_predict <<- timeROC(T=nomoROC$OS.time,
                           delta=nomoROC$OS,   #时间 结局
                           marker=nomoROC$predict,
                           cause=1,
                           weighting='marginal',
                           times= 1760,
                           ROC=TRUE)


png(file="Fig5E1.png",width=1800,height=2000,pointsize = 60)

plot(ROC_rt_stage,time=1760,col='black',title=FALSE,lwd=8)+box(lwd=6)
par(new=TRUE)

plot(ROC_rt_predict,time=1760,col='red',title=FALSE,lwd=8)


legend('bottomright',
       c(paste0('stage: ',round(ROC_rt_stage$AUC[2],2)),
         
         paste0('nomogram: ',round(ROC_rt_predict$AUC[2],2))
       ),
       col=c('black','red'),lwd=8,bty = 'n')

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
                         times= 1000,
                         ROC=TRUE)

ROC_rt_score <<- timeROC(T=nomoROC$OS.time,
                         delta=nomoROC$OS,   #时间 结局
                         marker=nomoROC$Riskscore,
                         cause=1,
                         weighting='marginal',
                         times= 1000,
                         ROC=TRUE)

ROC_rt_predict <<- timeROC(T=nomoROC$OS.time,
                           delta=nomoROC$OS,   #时间 结局
                           marker=nomoROC$predict,
                           cause=1,
                           weighting='marginal',
                           times= 1000,
                           ROC=TRUE)


png(file="Fig5D1.png",width=1800,height=2000,pointsize = 60)

plot(ROC_rt_stage,time=1000,col='black',title=FALSE,lwd=8)+box(lwd=6)
par(new=TRUE)

plot(ROC_rt_predict,time=1000,col='red',title=FALSE,lwd=8)


legend('bottomright',
       c(paste0('stage: ',round(ROC_rt_stage$AUC[2],2)),
         
         paste0('nomogram: ',round(ROC_rt_predict$AUC[2],2))
       ),
       col=c('black','red'),lwd=8,bty = 'n')

dev.off()



### 2. 一致性 (calibration)
# 3-year survival calibration curve
f3 <- cph(Surv(OS.time,OS) ~ Age +Stage + RiskScore, x=T, y=T, surv=T, data=nomomatrix, time.inc=1095)

cal3 <- calibrate(f3, cmethod="KM", method="boot", u=1095, m=75, B=375)

pdf(file="3yearnomo.pdf",width=15,height=15,pointsize = 32)
plot(cal3,xlim=c(0,1),ylim=c(0,1),xlab = 'Nomogram-predicted probability of 3-year survival',
     ylab = 'Observed 3-year survival',lwd=6)+box(lwd=4)
abline(a = 0, b = 1, col = "gray50",lwd=3,lty=2)
dev.off()

# 5-year survival calibration curve
f5 <- cph(Surv(OS.time,OS) ~ Age +Stage + RiskScore, x=T, y=T, surv=T, data=nomomatrix, time.inc=1825)

cal5 <- calibrate(f5, cmethod="KM", method="boot", u=1825, m=75, B=228)

pdf(file="5yearnomo.pdf",width=15,height=15,pointsize = 32)
plot(cal5,xlim=c(0,1),ylim=c(0,1),xlab = 'Nomogram-predicted probability of 5-year survival',
     ylab = 'Observed 5-year survival',lwd=6)+box(lwd=4)
abline(a = 0, b = 1, col = "gray50",lwd=3,lty=2)
dev.off()

