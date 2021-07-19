rm(list = ls())

coxhr = read.csv("Single_variable_cox_output.csv",header = T)

confidence = 0.05

coxhr = coxhr[coxhr$p.value < confidence , ]

coxhr$median = coxhr$HR
coxhr$lower = exp(coxhr$beta - 1.96 * coxhr$se)
coxhr$upper = exp(coxhr$beta + 1.96 * coxhr$se)

coxhr = coxhr[order(coxhr$p.value) , ]

coxhr$p.value.text = ''
coxhr[coxhr$p.value < 0.001 , ]$p.value.text = '<0.001'
coxhr[coxhr$p.value >= 0.001 , ]$p.value.text = as.character(sprintf("%0.4f", coxhr[coxhr$p.value >= 0.001 , ]$p.value))

coxhr$HR.text = paste(sprintf("%0.3f", coxhr$HR) , paste('(' , sprintf("%0.3f", coxhr$lower) , '-' , sprintf("%0.3f", coxhr$upper) , ')' , sep = ''))

tabletext <- cbind(c("Gene",as.vector(coxhr$gene)),
                   c("P-value",as.vector(coxhr$p.value.text)),
                   c("Hazard ratio",as.vector(coxhr$HR.text)))

# tabletext <- cbind(c(NA,"Gene",out_multi$id),
#                    c(NA,"Coefficient",round(out_multi$coef,3)),
#                    c(NA,"P value",ifelse(out_multi$pvalue<0.001,"P < 0.001",round(out_multi$pvalue,3))),
#                    c(NA,"Hazard Ratio(95% CI)",hz))




library(forestplot)
pdf(file="Fig2A.pdf",width=10,height=8,pointsize = 20)
forestplot(tabletext,  #显示的文本
           c(NA,coxhr$median), #误差条的均值(此处为差值的中值)
           c(NA,coxhr$lower), #误差条的下界(此处为差值的25%分位数)
           c(NA,coxhr$upper), #误差条的上界(此处为差值的75%分位数)
           zero = 1, #显示y=0的垂直线
           xlog=FALSE, #x轴的坐标不取对数
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
           xticks = c(0.8,1.0,1.2,1.4)
)
dev.off()