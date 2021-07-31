library('survival')
library('survminer')
rm(list = ls())
setwd('validation set')
load(file = '../survexprdata_geneCoef.Rdata')
rm('survexprdata')

load(file = '../survexprdata.Rdata')

survexprdata = validation_set
#survexprdata = entire_set

save(survexprdata, geneCoef, file='survexprdata_validation.Rdata')



rm(list = ls())
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
scoresurv$OS.time1 = scoresurv$OS.time/(365)
fit <- survfit(Surv(OS.time1, OS) ~ scoregroup, data = scoresurv)
surv_diff <- survdiff(Surv(OS.time1, OS) ~ scoregroup, data = scoresurv)
pvalue <- 1 - pchisq(surv_diff$chisq, length(surv_diff$n) -1)
# pvalue = 6.012701e-06
write.csv(scoresurv,'scoresurv.csv')
gg<- ggsurvplot(
  fit,pval = TRUE,
  pval.size = 8,
  ggtheme = theme(panel.background=element_blank(),
                  axis.line = element_line(size=1),
                  plot.title = element_text(hjust = 0.5),
                  axis.text = element_text(size = 16,face = "bold")),
  risk.table = TRUE,
  legend.title="Risk score",palette = c("red","blue"),
  font.legend=c(16,"bold","black"),
  legend.labs=c("High","Low"),
  title="TCGA-LIHC",font.title=c(18,"bold","black"),
  xlab="Time (years)",ylab="Survival probability",
  font.xlab=c(18,"bold","black"), font.ylab=c(18,"bold","black"),
  break.x.by = 1)

pdf('kmplot.pdf')
print(gg,newpage= FALSE)
dev.off()

setwd('../')
