#计算免疫细胞和基因的相关性
rm(list = ls())

#可选量化方法
# quantiseq
# xcell
# cibersort
# cibersort_abs
# mcp_counter
# epic
# timer
immunedeconv_method = 'xcell'

load(paste0('./', immunedeconv_method, '_res.Rdata'))
immunedeconv_res = as.data.frame(immunedeconv_res)
immunedeconv_res = t(immunedeconv_res)
colnames(immunedeconv_res) = immunedeconv_res['cell_type', ]
immunedeconv_res = immunedeconv_res[-1, ]

immunedeconv_res_row_names = row.names(immunedeconv_res)
immunedeconv_res = apply(immunedeconv_res, 2, as.numeric)
row.names(immunedeconv_res) = immunedeconv_res_row_names

# #计算M1/M2
# if('Macrophage M1' %in% colnames(immunedeconv_res) & 'Macrophage M2'  %in% colnames(immunedeconv_res)) {
#   #删除M2为0的样本, 避免后续出现M1/M2为无穷大的情况
#   immunedeconv_res = immunedeconv_res[immunedeconv_res[, 'Macrophage M2'] != 0 , ]
#   
#   #额外计算M1/M2的值作为新的一列
#   immunedeconv_colnames = colnames(immunedeconv_res)
#   immunedeconv_res = cbind(immunedeconv_res, immunedeconv_res[ , c('Macrophage M1')] / immunedeconv_res[ , c('Macrophage M2')])
#   colnames(immunedeconv_res) = c(immunedeconv_colnames, 'M1/M2')
# }
# 
# #计算T cell CD4+ Th1 / T cell CD4+ Th2
# if('T cell CD4+ Th1' %in% colnames(immunedeconv_res) & 'T cell CD4+ Th2'  %in% colnames(immunedeconv_res)) {
#   #删除Th2为0的样本, 避免后续出现Th1/Th2为无穷大的情况
#   immunedeconv_res = immunedeconv_res[immunedeconv_res[, 'T cell CD4+ Th2'] != 0 , ]
#   
#   #额外计算Th1/Th2的值作为新的一列
#   immunedeconv_colnames = colnames(immunedeconv_res)
#   immunedeconv_res = cbind(immunedeconv_res, immunedeconv_res[ , c('T cell CD4+ Th1')] / immunedeconv_res[ , c('T cell CD4+ Th2')])
#   colnames(immunedeconv_res) = c(immunedeconv_colnames, 'Th1/Th2')
# }
  

load(file = '../data/scoresurv.Rdata')
scoresurv = scoresurv[ , -which(colnames(scoresurv) %in% c("OS", "OS.time"))]

#保持两个数据集样本完全一致
ids = intersect(row.names(scoresurv), row.names(immunedeconv_res))
immunedeconv_res = immunedeconv_res[ids, ]
scoresurv = scoresurv[ids, ]

library(psych)

p_matrix = c()
r_matrix = c()
for(i in 1 : length(colnames(immunedeconv_res))) {
  corr = corr.test(immunedeconv_res[ , i], scoresurv, use = "pairwise", method = "pearson", adjust = "none")
  
  p_result = cbind("cell" = colnames(immunedeconv_res)[i], "p" = corr$p)
  r_result = cbind("cell" = colnames(immunedeconv_res)[i], "r" = corr$r)
  
  p_matrix = rbind(p_matrix, p_result)
  r_matrix = rbind(r_matrix, r_result)
}

row.names(p_matrix) = p_matrix[ , c('cell')]
p_matrix = p_matrix[ , -which('cell' %in% colnames(p_matrix))]
row.names(r_matrix) = r_matrix[ , c('cell')]
r_matrix = r_matrix[ , -which('cell' %in% colnames(r_matrix))]

row_names = row.names(p_matrix)
p_matrix = apply(p_matrix, 2, as.numeric)
r_matrix = apply(r_matrix, 2, as.numeric)
row.names(p_matrix) = row_names
row.names(r_matrix) = row_names

save(p_matrix, r_matrix, file = paste0('./immune_correlation_matrix_', immunedeconv_method, '.Rdata'))
write.csv(p_matrix, file = paste0('./immune_correlation_P_matrix_', immunedeconv_method, '.csv'))
write.csv(r_matrix, file = paste0('./immune_correlation_R_matrix_', immunedeconv_method, '.csv'))


#breaks
#bk <- c(seq(-0.6,-0.001,by=0.01),seq(0,0.6,by=0.01))
bk <- c(seq(-0.6, 0.6,by=0.01))
cell_num = ncol(immunedeconv_res)
# 相关性热图：
library(pheatmap)
pheatmap(r_matrix, display_numbers = matrix(ifelse(p_matrix <= 0.001, "***", ifelse(p_matrix<= 0.01 ,"**",ifelse(p_matrix<= 0.05,"*",""))),nrow(p_matrix)),
         fontsize=16, treeheight_row =0,treeheight_col = 0,angle_col = "45",cluster_cols  = F,cluster_rows = F,
         scale = "none",cexRow = 5,cellheight = 30,cellwidth = 30,
         color = c(colorRampPalette(colors = c("blue","white"))((length(bk)-1)/2),
                   colorRampPalette(colors = c("white","red"))((length(bk)-1)/2)),
         legend_breaks=seq(-0.6,0.6,0.3),
         breaks=bk,width = 9,height = cell_num * 0.6 + 2, filename = paste0('./immune_correlation_heatmap_', immunedeconv_method, '.pdf'))
