rm(list = ls())
options(stringsAsFactors = F)
clinical= read.csv("./clinical_change.csv")
rownames(clinical) = clinical[,1]
clinical = clinical[ , -1]

#删除MX, NX, TX的样本
clinical = clinical[clinical$pM != 'MX', ]
clinical = clinical[clinical$pN != 'NX', ]
clinical = clinical[clinical$pT != 'TX', ]

########## 计算 ##########

load("./survexprdata_lipid.Rdata")

clinicalexpr = log2(survexprdata[, geneCoef[,1] ]+1)
clinicalgene = colnames(clinicalexpr)

rownames(clinicalexpr) = substring(rownames(clinicalexpr),1,12) 
clinicalexpr = clinicalexpr[rownames(clinical),]
clinicalexpr = na.omit(clinicalexpr)
clinical = clinical[rownames(clinicalexpr),]
clinicalexpr = cbind(clinicalexpr,clinical)    #临床特征与表达数据merge

clinicalexpr$Riskscore = 0
for(i in 1:nrow(geneCoef)) {
  clinicalexpr$Riskscore = clinicalexpr$Riskscore + clinicalexpr[, i] * as.numeric(geneCoef[i, 2])
}


clinicalexpr$Risk = ifelse(clinicalexpr$Riskscore < median(clinicalexpr$Riskscore),'Low','High')
clinicalexpr$ppT = as.numeric(substring(clinicalexpr$pT,2,2))
clinicalexpr$ppM = as.numeric(substring(clinicalexpr$pM,2,2))
clinicalexpr$ppN = as.numeric(substring(clinicalexpr$pN,2,2))


####### 临床指标与Risk ROC曲线 ######
library(timeROC)
library(survival)

ROC_TIME = 365

ROC_rt_score <<- timeROC(T=clinicalexpr$OS.time,
                         delta=clinicalexpr$OS,   #时间 结局
                         marker=clinicalexpr$Riskscore,
                         cause=1,
                         weighting='marginal',
                         times= ROC_TIME,
                         ROC=TRUE)

ROC_rt_stage <<- timeROC(T=clinicalexpr$OS.time,
                         delta=clinicalexpr$OS,   #时间 结局
                         marker=clinicalexpr$stage1,
                         cause=1,
                         weighting='marginal',
                         times= ROC_TIME,
                         ROC=TRUE)

ROC_rt_age <<- timeROC(T=clinicalexpr$OS.time,
                       delta=clinicalexpr$OS,   #时间 结局
                       marker=clinicalexpr$age,
                       cause=1,
                       weighting='marginal',
                       times= ROC_TIME,
                       ROC=TRUE)

ROC_rt_T <<- timeROC(T=clinicalexpr$OS.time,
                     delta=clinicalexpr$OS,   #时间 结局
                     marker=clinicalexpr$ppT,
                     cause=1,
                     weighting='marginal',
                     times= ROC_TIME,
                     ROC=TRUE)

ROC_rt_N <<- timeROC(T=clinicalexpr$OS.time,
                     delta=clinicalexpr$OS,   #时间 结局
                     marker=clinicalexpr$ppN,
                     cause=1,
                     weighting='marginal',
                     times= ROC_TIME,
                     ROC=TRUE)

ROC_rt_M <<- timeROC(T=clinicalexpr$OS.time,
                     delta=clinicalexpr$OS,   #时间 结局
                     marker=clinicalexpr$ppM,
                     cause=1,
                     weighting='marginal',
                     times= ROC_TIME,
                     ROC=TRUE)

png(file="Fig4D.png",width=1800,height=2000,pointsize = 60)
plot(ROC_rt_score,time=ROC_TIME,col= '#FF1717',title=FALSE,lwd=8)+box(lwd=5)
par(new=TRUE)
plot(ROC_rt_stage,time=ROC_TIME,col='#50FFAF',title=FALSE,lwd=8)
par(new=TRUE)
plot(ROC_rt_age,time=ROC_TIME,col='#FFDC16',title=FALSE,lwd=8)
par(new=TRUE)
plot(ROC_rt_T,time=ROC_TIME,col='#0089FF',title=FALSE,lwd=8)
par(new=TRUE)
plot(ROC_rt_N,time=ROC_TIME,col='#391496',title=FALSE,lwd=8)
par(new=TRUE)
plot(ROC_rt_M,time=ROC_TIME,col='#FF05D9',title=FALSE,lwd=8)

legend('bottomright',
       c(paste0('Risk score: ',round(ROC_rt_score$AUC[2],3)),
         paste0('Age: ',round(ROC_rt_age$AUC[2],3)),
         paste0('Stage: ',round(ROC_rt_stage$AUC[2],3)),
         paste0('T: ',round(ROC_rt_T$AUC[2],3)),
         paste0('N: ',round(ROC_rt_N$AUC[2],3)),
         paste0('M: ',round(ROC_rt_M$AUC[2],3))
       ),
       col=c('#FF1717','#FFDC16','#50FFAF','#0089FF','#391496','#FF05D9'),lwd=4,bty = 'n')

dev.off()



####### 临床指标与表达热图 ######
clinicalexpr = clinicalexpr[order(clinicalexpr$Riskscore), ]
library(pheatmap)
annotation_col <- data.frame(Age = clinicalexpr$agegroup,Gender = clinicalexpr$gender,Stage = clinicalexpr$stage,
                             'T'=clinicalexpr$pT,N = clinicalexpr$pN,M = clinicalexpr$pM,Risk = clinicalexpr$Risk )
rownames(annotation_col)= rownames(clinicalexpr)

ann_colors = list(Risk = c(High="#FF7070", Low="#7070FF"),
                  Stage = c("Stage I" = "#E2F0D9", "Stage II" = "#C5E0B4","Stage III" = "#A9D18E","Stage IV" ="#548235"),
                  M= c(M0 = "#F8CBAD", M1 = "#C55A11"),
                  N= c(N0 = "#FFE699", N1 = "#FFD966", N2 = "#BF9000", N3 = '#7F6000'),
                  'T'= c(T1 = "#DEEBF7", T2 = "#BDD7EE",T3 = "#9DC3E6",T4 ="#2E75B6" ),
                  Gender= c(male = "#8E44BC", female = "#DFCAEC"),
                  Age= c('>=66' = "#7C7C7C", '<66' = "#DBDBDB"))

clinical_gene_matrix = t(scale(clinicalexpr[,1:nrow(geneCoef)]))

pheatmap(clinical_gene_matrix,show_colnames = F,annotation_col = annotation_col,
         cluster_cols = F,treeheight_row =10 ,annotation_colors = ann_colors,
         cellheight = 60,cellwidth = 1.3,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         filename = "Fig4C.tiff")

####### 单变量COX #######

library('survival')
library('survminer')
clinicalsig = c('Riskscore','agegroup','gender','stage1','pT1','pN1','pM')
clinical_expr_ucox = clinicalexpr[,c(clinicalsig,'OS','OS.time')]

clinical_univ_cox=c()
for(i in 1:length(clinicalsig)) {
  
  res.cox= coxph(Surv(OS.time,OS)~clinical_expr_ucox[,i]  ,data=clinical_expr_ucox)
  summary(res.cox)
  
  beta1 <- coef(res.cox)
  se1 <- sqrt(diag(vcov(res.cox)))
  HR1 <- exp(beta1)
  HRse1 <- HR1 * se1
  pcox = 1 - pchisq((beta1/se1)^2, 1)
  
  result<-cbind("Factor" = clinicalsig[i], "p.value" = pcox ,"beta"=beta1,"se"=se1,"HR"=HR1,
                as.data.frame(summary(res.cox)[8]))
  clinical_univ_cox<-rbind(clinical_univ_cox,result)
}

write.csv(clinical_univ_cox,'clinical_univ_cox.csv')


#### 单变量森林图 #####
coxhr = read.csv("./clinical_univ_cox.csv",header = T)

coxhr = coxhr[ , c('Factor', 'p.value', 'HR', 'conf.int.exp.coef.', 'conf.int.lower..95', 'conf.int.upper..95')]
names(coxhr) = c('Factor', 'p.value', 'Hazard.ratio', 'median', 'lower', 'upper')

coxhr$p.value.text = ''
coxhr[coxhr$p.value < 0.001 , ]$p.value.text = '<0.001'
coxhr[coxhr$p.value >= 0.001 , ]$p.value.text = as.character(sprintf("%0.4f", coxhr[coxhr$p.value >= 0.001 , ]$p.value))

coxhr$HR.text = paste(sprintf("%0.3f", coxhr$Hazard.ratio) , paste('(' , sprintf("%0.3f", coxhr$lower) , '-' , sprintf("%0.3f", coxhr$upper) , ')' , sep = ''))

coxhr$Factor = c('Riskscore', 'Age', 'Gender', 'Stage', 'T', 'N', 'M')

tabletext <- cbind(c("Factor",as.vector(coxhr$Factor)),
                   c("P-value",as.vector(coxhr$p.value.text)),
                   c("Hazard ratio",as.vector(coxhr$HR.text)))

library(forestplot)
png(file="clinical_expr_ucox_forest.png",width=1000,height=550,pointsize = 25)
forestplot(tabletext,  #显示的文本
           c(NA,coxhr$median), #误差条的均值(此处为差值的中值)
           c(NA,coxhr$lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,coxhr$upper), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=T, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.3, ##误差条中的圆心点大小
           col=fpColors(line = "blue", #误差条的线的颜色
                        box="red"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 3,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1), #文本大小
           lineheight = "auto", #线的高度 
           xlab="Hazard ratio" ,#x轴的标题
           xticks = c(0.5,1.0,5,30)
)
dev.off()


####### 多变量COX #######

library('survival')
library('survminer')
clinicalsig = c('Riskscore','agegroup','gender','stage1','pT1','pN1','pM')
clinical_expr_ucox = clinicalexpr[,c(clinicalsig,'OS','OS.time')]

res.cox= coxph(Surv(OS.time,OS)~agegroup+gender+stage1+pT1+pN1+pM+Riskscore ,
               data=clinical_expr_ucox)
clinical_m_cox = cbind(summary(res.cox)[["coefficients"]], summary(res.cox)[["conf.int"]])

clinical_m_cox = as.data.frame(clinical_m_cox)
clinical_m_cox$Factor = row.names(clinical_m_cox)
clinical_m_cox$median = clinical_m_cox$`exp(coef)`
clinical_m_cox = clinical_m_cox[, c('Factor', 'Pr(>|z|)','coef' ,'se(coef)' , 'exp(coef)', 'median', 'lower .95', 'upper .95')]
colnames(clinical_m_cox) = c('Factor', 'p.value', 'beta', 'se', 'Hazard.ratio', 'median', 'lower', 'upper')
clinical_m_cox = clinical_m_cox[c('Riskscore', 'agegroup>=66', 'gendermale', 'stage1', 'pT1T3-T4', 'pN1N1-N3', 'pMM1'), ]
rownames(clinical_m_cox)= c('Riskscore', 'Age', 'Gender', 'Stage', 'T', 'N', 'M')


write.csv(clinical_m_cox,'clinical_m_cox.csv')


#### 多变量森林图 #####
coxhr = read.csv("./clinical_m_cox.csv",header = T)

coxhr$p.value.text = ''
coxhr[coxhr$p.value < 0.001 , ]$p.value.text = '<0.001'
coxhr[coxhr$p.value >= 0.001 , ]$p.value.text = as.character(sprintf("%0.4f", coxhr[coxhr$p.value >= 0.001 , ]$p.value))

coxhr$HR.text = paste(sprintf("%0.3f", coxhr$Hazard.ratio) , paste('(' , sprintf("%0.3f", coxhr$lower) , '-' , sprintf("%0.3f", coxhr$upper) , ')' , sep = ''))

coxhr$Factor = coxhr[ , 1]

tabletext <- cbind(c("Factor",as.vector(coxhr$Factor)),
                   c("P-value",as.vector(coxhr$p.value.text)),
                   c("Hazard ratio",as.vector(coxhr$HR.text)))

library(forestplot)
png(file="clinical_expr_m_forest.png",width=1000,height=550,pointsize = 25)
forestplot(tabletext,  #显示的文本
           c(NA,coxhr$median), #误差条的均值(此处为差值的中值)
           c(NA,coxhr$lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,coxhr$upper), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=T, #x轴的坐标不取对数
           fn.ci_norm = fpDrawCircleCI, #误差条显示方式
           boxsize = 0.3, ##误差条中的圆心点大小
           col=fpColors(line = "blue", #误差条的线的颜色
                        box="red"), #误差条的圆心点的颜色
           lty.ci = 7,   # 误差条的线的线型
           lwd.ci = 3,   # 误差条的线的宽度
           ci.vertices.height = 0.15, # # 误差条末端的长度
           txt_gp = fpTxtGp(ticks = gpar(cex = 1), xlab = gpar(cex = 1), cex = 1), #文本大小
           lineheight = "auto", #线的高度 
           xlab="Hazard ratio" ,#x轴的标题
           xticks = c(0.1,1.0,5)
)
dev.off()