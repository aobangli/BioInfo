rm(list = ls())
options(stringsAsFactors = F)
clinical= read.csv("../clinical_change.csv")
rownames(clinical) = clinical[,1]
clinical = clinical[ , -1]

#删除MX, NX, TX的样本
clinical = clinical[clinical$pM != 'MX', ]
clinical = clinical[clinical$pN != 'NX', ]
clinical = clinical[clinical$pT != 'TX', ]

########## 计算 ##########

load("./survexprdata_validation.Rdata")

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
         cellheight = 20,cellwidth = 1.3,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         filename = "Fig4C.tiff")