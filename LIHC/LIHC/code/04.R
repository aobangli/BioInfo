rm(list = ls())
load(file="survexprdata_lipid.Rdata")
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
                   times=c(365,1095,1852),ROC=TRUE)

png(file="ROC.png",width=2000,height=2000,pointsize = 60)

plot(ROC_rt,time=365,title=FALSE,lwd=8)+box(lwd=5)
plot(ROC_rt,time=1095,col='blue',add=TRUE,title=FALSE,lwd=8)
plot(ROC_rt,time=1852,col='black',add=TRUE,title=FALSE,lwd=8)

legend('bottomright',
       c(paste0('AUC at 1 years: ',round(ROC_rt$AUC[1],2)),
         paste0('AUC at 3 years: ',round(ROC_rt$AUC[2],2)),
         paste0('AUC at 5 years: ',round(ROC_rt$AUC[3],2))),
       col=c('red','blue','black'),lwd=4,bty = 'n')
dev.off()