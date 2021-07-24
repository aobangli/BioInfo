rm(list = ls())
load(file="survexprdata_geneCoef.Rdata")
######### riskscore surv #######
scoresurv = survexprdata[, geneCoef[,1] ]
scoresurv = log2(scoresurv +1)
scoresurv = cbind(scoresurv,survexprdata[,((ncol(survexprdata)-1):(ncol(survexprdata)))])

scoresurv$riskscore = 0
for(i in 1:nrow(geneCoef)) {
  scoresurv$riskscore = scoresurv$riskscore + scoresurv[, i] * as.numeric(geneCoef[i, 2])
}

#按riskscore中位数分组

scoresurv$scoregroup = ifelse(scoresurv$riskscore<median(scoresurv$riskscore),'low','high')

library(timeROC)
library(survival)

ROC_rt <<- timeROC(T=scoresurv$OS.time,delta=scoresurv$OS,
                   marker=scoresurv$riskscore,cause=1,
                   weighting='marginal',
                   times=c(365,1095, 1824),ROC=TRUE)

pdf(file="ROC.pdf",width = 10,height = 11,pointsize = 25)

plot(ROC_rt,time=365,title=FALSE,lwd=3)+box(lwd=4)
plot(ROC_rt,time=1095,col='blue',add=TRUE,title=FALSE,lwd=3)
plot(ROC_rt,time=1824,col='black',add=TRUE,title=FALSE,lwd=3)

legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],3)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],3)),
         paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],3))),
       col=c('red','blue','black'),lwd=4,bty = 'n')
dev.off()