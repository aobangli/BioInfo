# install.packages("remotes")
# remotes::install_github("icbi-lab/immunedeconv")
rm(list = ls())

########## 筛选出临床特征符合条件的样本的基因表达数据 ##########

clinical= read.csv("../data/clinical_change.csv")
rownames(clinical) = clinical[,1]
clinical = clinical[ , -1]
sample_ids = clinical$id

load('.\\LIHC_TPM_expMatrix.Rdata')

########### 筛选特定肿瘤样本 ###########
library(stringr)
sampletype= substr(colnames(LIHC_TPM_expMatrix),14,15)
LIHC_PT_TPM_expr= LIHC_TPM_expMatrix[,sampletype=='01']

tmp1= substr(colnames(LIHC_PT_TPM_expr),16,16)
LIHC_PT_TPM_expr= LIHC_PT_TPM_expr[,tmp1=='A']

colnames(LIHC_PT_TPM_expr) = substr(colnames(LIHC_PT_TPM_expr), 1, 12)

sample_ids = intersect(sample_ids, colnames(LIHC_PT_TPM_expr))

LIHC_PT_TPM_expr = LIHC_PT_TPM_expr[ , sample_ids]

#样本名去重
#LIHC_PT_TPM_expr = LIHC_PT_TPM_expr[ , !duplicated(colnames(LIHC_PT_TPM_expr))]

library('immunedeconv')

########### quantiseq ###########

immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "quantiseq")

save(immunedeconv_res, file = './quantiseq_res.Rdata')

########### xcell ###########
immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "xcell")

save(immunedeconv_res, file = './xcell_res.Rdata')


########### cibersort ###########
source("../code/CIBERSORT.R")
#write.table(LIHC_PT_TPM_expr, file = 'LIHC_PT_TPM_expr.txt')
res = CIBERSORT('LIHC_PT_TPM_expr.txt', "LM22.txt", perm = 10, QN = F)


# set_cibersort_binary("../code/CIBERSORT.R")
# set_cibersort_mat("./LM22.txt")
# 
# immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "cibersort")

save(immunedeconv_res, file = './cibersort_res.Rdata')



########### mcp_counter ###########
immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "mcp_counter")

save(immunedeconv_res, file = './mcp_counter_res.Rdata')


########### epic ###########
immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "epic")

save(immunedeconv_res, file = './epic_res.Rdata')

########### TIMER ###########
immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "timer", indications = rep('LIHC', ncol(LIHC_PT_TPM_expr)))

save(immunedeconv_res, file = './timer_res.Rdata')