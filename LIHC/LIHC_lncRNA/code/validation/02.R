library('survival')
library('survminer')
rm(list = ls())
load(file = '../LIHC_lipid_limma_DEG.Rdata')
load(file = '../clinical.Rdata')
load(file = '../survexprdata_lipid.Rdata')
rm('survexprdata')

#survexprdata = validation_set
survexprdata = LIHC_PT_COUNT_expr

colnames(survexprdata)= substr(colnames(survexprdata),1,12)

genename = rownames(need_DEG)[need_DEG$change != 'NOT']
survexprdata= survexprdata[genename,]

clinical= clinical[substr(colnames(survexprdata),1,12),]

OSdata=clinical[,c('id', 'OS', 'OS.time')]
OSdata= na.omit(OSdata)
survexprdata= survexprdata[,rownames(OSdata)]
survexprdata=cbind(t(survexprdata),OSdata[,2:3])
write.csv(survexprdata,'survexprdata_lipid.csv')

save(survexprdata, geneCoef, file='survexprdata_validation.Rdata')

rm(list = ls())
#setwd(dir = 'D:/LIHC/COUNT')
load(file = 'survexprdata_validation.Rdata')
######### riskscore surv #######
scoresurv = survexprdata[, geneCoef[,1] ]
scoresurv = log2(scoresurv +1)
scoresurv = cbind(scoresurv,survexprdata[,((ncol(survexprdata)-1):(ncol(survexprdata)))])

scoresurv$riskscore = 0
for(i in 1:nrow(geneCoef)) {
  scoresurv$riskscore = scoresurv$riskscore + scoresurv[, i] * as.numeric(geneCoef[i, 2])
}

save(scoresurv, file = 'scoresurv.Rdata')

###### 按riskscore中位数分组 ######

scoresurv$scoregroup = ifelse(scoresurv$riskscore<median(scoresurv$riskscore),'low','high')

fit <- survfit(Surv(OS.time, OS) ~ scoregroup, data = scoresurv)
surv_diff <- survdiff(Surv(OS.time, OS) ~ scoregroup, data = scoresurv)
pvalue <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
# pvalue = 6.012701e-06
write.csv(scoresurv,'scoresurv.csv')
gg<- ggsurvplot(
  fit,pval = TRUE,
  ggtheme = theme(panel.background=element_blank(),
                  axis.line = element_line(size=1),
                  plot.title = element_text(hjust = 0.5),
                  axis.text = element_text(size = 11,face = "bold")),
  risk.table = TRUE,
  legend.title="Risk score",palette = c("red","blue"),
  font.legend=c(11,"bold","black"),
  legend.labs=c("High","Low"),
  title="TCGA-LIHC",font.title=c(18,"bold","black"),
  xlab="Time (days)",ylab="Survival probability",
  font.xlab=c(14,"bold","black"), font.ylab=c(14,"bold","black"))

pdf('kmplot.pdf')
print(gg,newpage= FALSE)
dev.off()