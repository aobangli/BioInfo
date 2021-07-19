rm(list = ls())
load(file = 'survexprdata_lipid.Rdata')

scoresurv = survexprdata[, geneCoef[,1] ]
scoresurv = log2(scoresurv +1)
scoresurv = cbind(scoresurv,survexprdata[,((ncol(survexprdata)-1):(ncol(survexprdata)))])

scoresurv$riskscore = 0
for(i in 1:nrow(geneCoef)) {
  scoresurv$riskscore = scoresurv$riskscore + scoresurv[, i] * as.numeric(geneCoef[i, 2])
}


distribution= scoresurv[order(scoresurv$riskscore), ]
distribution$Patients = 1:nrow(distribution)
distribution$scoregroup = ifelse(distribution$riskscore<median(distribution$riskscore),
                                 'Low risk','High risk')

library(ggplot2)
P<- ggplot(distribution,aes(x=Patients, y=riskscore,colour =scoregroup ,))+
  geom_point() +scale_colour_manual(values=c("#cc0000","#000099"))+
  geom_hline(yintercept=c(median(distribution$riskscore)),linetype=2)+
  geom_vline(xintercept=c(median(distribution$Patients)),linetype=2)+
  theme(legend.position = c(0.07, 0.9),legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = 'transparent',colour = 'black'),
        axis.text.x=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 12,face = 'bold', vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 12,face = 'bold', vjust = 0.5, hjust = 0.5))+
  xlab("Patients (increasing risk score)") + ylab("Risk score") 

ggsave(filename = "risk_distribution.pdf", 
       width = 25, height = 10, units="cm",dpi = 500)


######## score与生存 散点图 ######
distribution$event = ifelse((distribution$OS==0),"Alive","Dead")
distribution$survivaltime = distribution$OS.time/30.5

library(ggplot2)
Q<- ggplot(distribution,aes(x=Patients, y=survivaltime,colour =event,))+
  geom_point() +scale_colour_manual(values=c("#000099","#cc0000"))+
  geom_vline(xintercept=c(median(distribution$Patients)),linetype=2)+
  theme(legend.position = c(0.07, 0.9),legend.text = element_text(size = 12),
        legend.title = element_text(size = 0),
        panel.background = element_rect(fill = 'transparent',colour = 'black'),
        axis.text.x=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 12,face = 'bold', vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 12,face = 'bold', vjust = 0.5, hjust = 0.5))+
  xlab("Patients (increasing risk score)") + ylab("Survival time (months)") 

ggsave(filename = "time_distribution.pdf", 
       width = 25, height = 10, units="cm",dpi = 500)



######## sig_gene热图 ######

library(pheatmap)
scoregroup = ifelse(distribution$riskscore<median(distribution$riskscore),
                    'Low','High')

annotation_col = data.frame(scoregroup)
rownames(annotation_col)= rownames(distribution)
sig_gene_matrix = t(scale(distribution[, 1:nrow(geneCoef)]))
colnames(annotation_col)='Risk score'

# ann_colors = list(Time = c("white", "firebrick"), #连续数值型分组可设置成渐变
#                   CellType = c(CT1 = "#1B9E77", CT2 = "#D95F02"),
#                   GeneClass = c(Path1 = "#7570B3", Path2 = "#E7298A", Path3 = "#66A61E"))

ann_colors = list('Risk score' =c(High="#cc0000", Low="#000099" ))
tsig_gene_matrix = t(sig_gene_matrix)

#pheatmap(tsig_gene_matrix,show_rownames = F,annotation_row  = annotation_col,
#         cluster_rows = F,treeheight_col =12 ,annotation_colors = ann_colors,
#         cellheight = 1.3,cellwidth = 30,angle_col = 45,
#         color = colorRampPalette(colors = c("blue","white","red"))(100),
#         filename = "Fig2D.pdf",width = 8)

pheatmap(sig_gene_matrix,show_colnames = F,annotation_col = annotation_col,cluster_cols = F,annotation_colors = ann_colors, 
         cellheight = 16,cellwidth = 0.3,treeheight_row =5,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         filename = 'Fig2D.pdf')