rm(list = ls())

#取全集并计算riskscore
load(file="../data/survexprdata_geneCoef.Rdata")
rm('survexprdata')

load('../data/LIHC_Count_expMatrix.Rdata')
load(file="../data/survexprdata.Rdata")

genes = c("PDCD1", "CD274", "PDCD1LG2", "HAVCR2", "CTLA4", 
          "IDO1", "TNFRSF18", "SOAT1", "CDK1", "HDAC2", "MMP9")

######### 样本筛选 ##########
LIHC_Count_expMatrix = LIHC_Count_expMatrix[genes, ]

colnames(LIHC_Count_expMatrix) = substr(colnames(LIHC_Count_expMatrix), 1, 12)
LIHC_Count_expMatrix = LIHC_Count_expMatrix[, row.names(entire_set)]
LIHC_Count_expMatrix = t(LIHC_Count_expMatrix)

#计算riskscore
scoresurv = entire_set[, geneCoef[,1] ]
scoresurv = log2(scoresurv +1)
scoresurv = cbind(scoresurv,entire_set[,((ncol(entire_set)-1):(ncol(entire_set)))])

scoresurv$riskscore = 0
for(i in 1:nrow(geneCoef)) {
  scoresurv$riskscore = scoresurv$riskscore + scoresurv[, i] * as.numeric(geneCoef[i, 2])
}
#删除OS, OS.time
scoresurv = scoresurv[ , -which(colnames(scoresurv) %in% c("OS", "OS.time"))]

#只取riskscore
riskscore = c(scoresurv[ , 'riskscore'])

gene_score_Matrix = cbind(LIHC_Count_expMatrix, riskscore)
colnames(gene_score_Matrix) = c("PD-1", "PD-L1", "PD-L2", "TIM3", "CTLA4", 
                                  "IDO1", "GITR", "SOAT1", "CDK1", "HDAC2", "MMP9", "RiskScore")

##### 计算gene_score_Matrix(immune check gene加上risk score)自身的相关性 #####
library(psych)
corr = corr.test(gene_score_Matrix, use = "pairwise", method = "pearson", adjust = "none")

p_matrix = corr$p
r_matrix = corr$r

save(p_matrix, r_matrix, file = 'immune_check_correlation_matrix_1.Rdata')
write.csv(p_matrix, file = 'immune_check_correlation_P_matrix_1.csv')
write.csv(r_matrix, file = 'immune_check_correlation_R_matrix_1.csv')

##### 相关性可视化 #####
library('corrplot')
color <- colorRampPalette(c("blue", "white", "red"))

pdf(file="immune_check_correlation_1.pdf",width = 10,height = 10,pointsize = 15)
corrplot(r_matrix, p.mat = p_matrix, insig = "label_sig", sig.level = c(.001, .01, .05),  # 设置显著性阈值为0.05
         pch.cex = .9, pch.col = "black",
         method = 'circle',
         type = "lower",
         diag = FALSE,
         tl.srt = 45,
         tl.col = 'black',
         tl.cex = 1.3,
         col.lim=c(-1, 1),
         col = color(100))
dev.off()

##### 计算模型内lin及score 和immune check gene的相关性 #####
colnames(LIHC_Count_expMatrix) = c("PD-1", "PD-L1", "PD-L2", "TIM3", "CTLA4", 
                                   "IDO1", "GITR", "SOAT1", "CDK1", "HDAC2", "MMP9")
corr = corr.test(scoresurv, LIHC_Count_expMatrix, use = "pairwise", method = "pearson", adjust = "none")
p_matrix = corr$p
r_matrix = corr$r

save(p_matrix, r_matrix, file = 'immune_check_correlation_matrix_2.Rdata')
write.csv(p_matrix, file = 'immune_check_correlation_P_matrix_2.csv')
write.csv(r_matrix, file = 'immune_check_correlation_R_matrix_2.csv')

##### 相关性可视化 #####
bk <- c(seq(-0.6, 0.6,by=0.01))
# 相关性热图：
library(pheatmap)
pheatmap(r_matrix, display_numbers = matrix(ifelse(p_matrix <= 0.001, "***", ifelse(p_matrix<= 0.01 ,"**",ifelse(p_matrix<= 0.05,"*",""))),nrow(p_matrix)),
         fontsize=16, treeheight_row =0,treeheight_col = 0,angle_col = "45",cluster_cols  = F,cluster_rows = F,
         scale = "none",cexRow = 5,cellheight = 30,cellwidth = 30,
         color = c(colorRampPalette(colors = c("blue","white"))((length(bk)-1)/2),
                   colorRampPalette(colors = c("white","red"))((length(bk)-1)/2)),
         legend_breaks=seq(-0.6,0.6,0.3),
         breaks=bk,width = 9,height = 10, filename = 'immune_check_correlation_2.pdf')

