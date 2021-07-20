rm(list=ls())
options(stringsAsFactors = F)
load('.\\lncRNA_matrix_to_cal.Rdata')

######### 样本筛选 ##########
library(stringr)
sampletype= substr(colnames(lncRNA_matrix_to_cal),14,15)
LIHC_PT_COUNT_expr= lncRNA_matrix_to_cal[,sampletype=='01']
LIHC_CN_COUNT_expr= lncRNA_matrix_to_cal[,sampletype=='11']

tmp1= substr(colnames(LIHC_PT_COUNT_expr),16,16)
LIHC_PT_COUNT_expr= LIHC_PT_COUNT_expr[,tmp1=='A']

tmp2= substr(colnames(LIHC_CN_COUNT_expr),16,16)
LIHC_CN_COUNT_expr= LIHC_CN_COUNT_expr[,tmp2=='A']

sampleID1= substr(colnames(LIHC_PT_COUNT_expr),1,12)
LIHC_PT_COUNT_expr<- LIHC_PT_COUNT_expr[,!duplicated(sampleID1)]  #删除表达矩阵中重复数据

sampleID2= substr(colnames(LIHC_CN_COUNT_expr),1,12)
LIHC_CN_COUNT_expr<- LIHC_CN_COUNT_expr[,!duplicated(sampleID2)]  #删除表达矩阵中重复数据

LIHC_TN_COUNT_expr<- cbind(LIHC_PT_COUNT_expr,LIHC_CN_COUNT_expr)

group_list<- c(rep('tumor',each=length(colnames(LIHC_PT_COUNT_expr))),
               rep('normal',each=length(colnames(LIHC_CN_COUNT_expr))))

########limma 差异表达分析 ######
library(edgeR)

exprSet= LIHC_TN_COUNT_expr
suppressMessages(library(limma))

design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(exprSet)


dge <- DGEList(counts=exprSet)
dge <- calcNormFactors(dge)
logCPM <- cpm(dge, log=TRUE, prior.count=3)

v <- voom(dge,design,plot=F, normalize="quantile")   #不画图
fit <- lmFit(v, design)

cont.matrix=makeContrasts(contrasts=c('tumor-normal'),levels = design)
fit2=contrasts.fit(fit,cont.matrix)
fit2=eBayes(fit2)

tempOutput = topTable(fit2, coef='tumor-normal', n=Inf)
DEG_limma_voom = na.omit(tempOutput)

nrDEG=DEG_limma_voom[,c(1,4)]
colnames(nrDEG)=c('log2FoldChange','pvalue') 

exprSet=LIHC_TN_COUNT_expr
need_DEG=nrDEG


####### 火山图 #######
logFC_cutoff=1.5
confidenceDEG = 0.01
need_DEG$change = as.factor(ifelse(need_DEG$pvalue < confidenceDEG & abs(need_DEG$log2FoldChange) > logFC_cutoff,
                                   ifelse(need_DEG$log2FoldChange > logFC_cutoff ,'UP','DOWN'),'NOT'))
this_tile <- paste0('Cutoff for logFC is ',round(logFC_cutoff,3),
                    '\nThe number of up gene is ',nrow(need_DEG[need_DEG$change =='UP',]) ,
                    '\nThe number of down gene is ',nrow(need_DEG[need_DEG$change =='DOWN',]))
library(ggplot2)
g = ggplot(data=need_DEG, 
           aes(x=log2FoldChange, y=-log10(pvalue), 
               color=change)) +
  geom_hline(yintercept=-log10(confidenceDEG),linetype=2)+  #增加阈值线
  geom_vline(xintercept=c(-logFC_cutoff , logFC_cutoff),linetype=2)+
  geom_point(alpha=1, size=2) +
  theme(panel.background = element_rect(fill = 'transparent',colour = 'black'),
        axis.text.x=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.text.y=element_text(size = 12, vjust = 0.5, hjust = 0.5),
        axis.title.x = element_text(size = 15, vjust = 0.5, hjust = 0.5),
        axis.title.y = element_text(size = 15, vjust = 0.5, hjust = 0.5))+
  xlab("log2(Fold change)") + ylab("-log10(Adjust P-value)") +
  ggtitle( this_tile ) + theme(plot.title = element_text(size=15,hjust = 0.5))+
  scale_colour_manual(values = c('blue','grey','red')) 

ggsave(g,filename = 'Fig1A.pdf',width = 15, height = 15, units="cm",dpi = 180)


###### DEG 热图 #####
library(pheatmap)
choose_DEG_gene=rownames(need_DEG)[need_DEG$change!='NOT']  #找出DEG
choose_DEG_matrix=exprSet[choose_DEG_gene,]

choose_DEG_matrix=t(scale(t(log2(choose_DEG_matrix+1)))) #画热图要先归一化
## http://www.bio-info-trainee.com/1980.html

annotation_coll = data.frame( group_list=group_list)#变成一个表
annotation_coll$sample = colnames(exprSet)
annotation_coll= annotation_coll[order(annotation_coll[,1]),]
annotation_coll$group_list = gsub('normal','Normal',annotation_coll$group_list)
annotation_coll$group_list = gsub('tumor','Tumor',annotation_coll$group_list)

choose_DEG_matrix= choose_DEG_matrix[,c(annotation_coll$sample)]
group_list_heatmap=annotation_coll[,1]
annotation_col = data.frame( group_list_heatmap=group_list_heatmap)#变成一个表
rownames(annotation_col)= colnames(choose_DEG_matrix)

colnames(annotation_col)=' '
pheatmap(choose_DEG_matrix,show_colnames = F,annotation_col = annotation_col,cluster_cols = F,
         cellheight = 9,cellwidth = 1.9,
         color = colorRampPalette(colors = c("blue","white","red"))(100),
         filename = 'Fig1B.pdf')

##### DEG_PCA #######
library(ggfortify)
df=as.data.frame(t(choose_DEG_matrix))#行列转换
df$group=group_list_heatmap
df$group=gsub('Normal','2 Normal',df$group)
df$group=gsub('Tumor','1 Tumor',df$group)


autoplot(prcomp( df[,1:(ncol(df)-1)] ), data=df,colour = "group")+theme_bw()#画图
ggsave('Fig1C.pdf',width = 13, height = 10, units="cm",dpi = 180)

####### GO 富集分析 #####
#source('D:/RData/r final/functions_kegg.R')
# library(ggplot2)
# library(clusterProfiler)
# library(org.Hs.eg.db)
# 
# deg=nrDEG
# colnames(deg)=c('logFC','P.Value')
# 
# ## 不同的阈值，筛选到的差异基因数量就不一样，后面的超几何分布检验结果就大相径庭。
# logFC_t = logFC_cutoff
# deg$g=ifelse(deg$P.Value>0.05,'stable',
#              ifelse( deg$logFC > logFC_t,'UP',
#                      ifelse( deg$logFC < -logFC_t,'DOWN','stable')))
# # table(deg$g)
# # head(deg)
# deg$symbol=rownames(deg)
# 
# df <- bitr(unique(deg$symbol), fromType = "SYMBOL",
#            toType = c( "ENTREZID"),
#            OrgDb = org.Hs.eg.db)
# 
# DEG=deg
# DEG=merge(DEG,df,by.y='SYMBOL',by.x='symbol')
# 
# DEG= DEG[DEG[,4]!='stable',]
# 
# 
# ####### ALL ########
# ALL <- enrichGO(gene=DEG$ENTREZID,
#                 OrgDb=org.Hs.eg.db,
#                 keyType = "ENTREZID",
#                 ont = 'ALL',
#                 pvalueCutoff = 0.3,
#                 pAdjustMethod = "BH",
#                 qvalueCutoff = 0.3,
#                 readable=T)
# 
# #泡泡图
# dotplot(ALL, split="ONTOLOGY",showCategory=5)+ facet_grid(ONTOLOGY~.,scale="free")
# ggsave('Fig1D.pdf',width = 28, height = 14, units="cm",dpi = 180)
# 
# 
# ####### kegg #######
# 
# gene_all=as.character(DEG[ ,'ENTREZID'] )
# ekk <- enrichKEGG(gene= gene_all,organism  = 'hsa', qvalueCutoff = 0.05)	 #KEGG富集分析
# dotplot(ekk,font.size=12,showCategory=9)	# 画气泡图
# ggsave('Fig1E.pdf',width = 18, height = 9.5, units="cm",dpi = 180)
# 


####### 划分训练集和验证集 #######
#以肿瘤样本集为全集划分, 即LIHC_PT_COUNT_expr
set.seed(20210717)
index <-  sort(sample(ncol(LIHC_PT_COUNT_expr), ncol(LIHC_PT_COUNT_expr)*.7))
training_set <- LIHC_PT_COUNT_expr[, index]
validation_set <-  LIHC_PT_COUNT_expr[, -index]

save(need_DEG, nrDEG, DEG_limma_voom, LIHC_TN_COUNT_expr, LIHC_PT_COUNT_expr,
     LIHC_CN_COUNT_expr, training_set, validation_set, group_list, file='LIHC_lipid_limma_DEG.Rdata')


####### 读取临床信息文件 #######
rm(list = ls())

# library('readr')
# 
# clinical_tsv = read_tsv('./clinical.cart.2021-07-09/clinical.tsv')
# alive_group = clinical_tsv[clinical_tsv$vital_status == 'Alive' & clinical_tsv$days_to_last_follow_up > 0, c('case_submitter_id', 'days_to_last_follow_up')]
# alive_group$OS = 0
# names(alive_group) = c('id', 'OS.time', 'OS')
# 
# dead_group = clinical_tsv[clinical_tsv$vital_status == 'Dead' & clinical_tsv$days_to_death > 0, c('case_submitter_id', 'days_to_death')]
# dead_group$OS = 1
# names(dead_group) = c('id', 'OS.time', 'OS')
# 
# clinical = rbind(alive_group, dead_group)
# clinical = clinical[!duplicated(clinical$id) , ]
# 
# clinical$OS.time = as.numeric(clinical$OS.time)
# clinical = as.data.frame(clinical)
# rownames(clinical) = clinical$id

clinical= read.csv("./clinical_change.csv")
rownames(clinical) = clinical[,1]
clinical = clinical[ , -1]
save(clinical, file = './clinical.Rdata')


####### 单变量COX #######

rm(list = ls())
load(file = './LIHC_lipid_limma_DEG.Rdata')
load(file = './clinical.Rdata')

survexprdata = training_set
#survexprdata = LIHC_PT_COUNT_expr

colnames(survexprdata)= substr(colnames(survexprdata),1,12)

genename = rownames(need_DEG)[need_DEG$change != 'NOT']
survexprdata= survexprdata[genename,]

clinical= clinical[substr(colnames(survexprdata),1,12),]

OSdata=clinical[,c('id', 'OS', 'OS.time')]
OSdata= na.omit(OSdata)
survexprdata= survexprdata[,rownames(OSdata)]
survexprdata=cbind(t(survexprdata),OSdata[,2:3])
write.csv(survexprdata,'survexprdata_lipid.csv')

library('survival')
library('survminer')

coxdata= survexprdata[,1:(ncol(survexprdata)-2)]
coxdata= log2(coxdata+1)
coxgene= colnames(coxdata)
coxdata= cbind(coxdata,survexprdata[,c('OS','OS.time')])

coxoutput=c()
for(i in 1:length(coxgene)) {
  
  res.cox= coxph(Surv(OS.time,OS)~coxdata[,i]  ,data=coxdata)
  summary(res.cox)
  
  beta1 <- coef(res.cox)
  se1 <- sqrt(diag(vcov(res.cox)))
  HR1 <- exp(beta1)
  HRse1 <- HR1 * se1
  pcox = 1 - pchisq((beta1/se1)^2, 1)
  
  result<-cbind("gene" = coxgene[i], "p.value" = pcox ,"beta"=beta1,"se"=se1,"HR"=HR1)
  coxoutput<-rbind(coxoutput,result)}
write.csv(coxoutput,'Single_variable_cox_output.csv')

######## lasso 回归 #######
lasso_confidence = 0.01

coxoutput = as.data.frame(coxoutput)
coxoutput$p.value = as.numeric(coxoutput$p.value)
lassocoxgene = coxoutput[coxoutput$p.value < lasso_confidence, 1]

lassoexpr= LIHC_PT_COUNT_expr[lassocoxgene,]
colnames(lassoexpr)=substr(colnames(lassoexpr),1,12)

library(lars) 
library(glmnet) 
x=t(log2(lassoexpr+1))
x=x[rownames(OSdata),]

OStime=as.numeric(OSdata[,3])
OS=as.numeric(OSdata[,2])
y=cbind(time=OStime,status=OS)


fit=glmnet(x, y, family = "cox", maxit = 1000) 
pdf(file="Fig2B.pdf",width=12,height=12,pointsize = 24)
plot(fit, xvar = "lambda", label = TRUE,lwd=6,cex.axis=1.2,cex.lab = 1.5)+box(lwd=4)
dev.off()


cvfit = cv.glmnet(x, y, family="cox", maxit = 1000) 
pdf(file="Fig2C.pdf",width=12,height=12,pointsize = 24)
plot(cvfit,lwd=6,cex.axis=1.2,cex.lab = 1.5)+box(lwd=4) 
#其中两条虚线分别指示了两个特殊的λ值 
abline(v = log(c(cvfit$lambda.min,cvfit$lambda.1se)),lty="dashed")
dev.off()

coef = coef(fit, s = cvfit$lambda.min) 
index = which(coef != 0) 
actCoef = coef[index] 
lassoGene = row.names(coef)[index] 
geneCoef = cbind(Gene=lassoGene,Coef=actCoef) 
geneCoef   #查看模型的相关系数

save(survexprdata, geneCoef, file='survexprdata_lipid.Rdata')

rm(list = ls())
#setwd(dir = 'D:/LIHC/COUNT')
load(file = 'survexprdata_lipid.Rdata')
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