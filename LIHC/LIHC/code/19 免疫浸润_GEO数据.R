rm(list = ls())

load('gene_expr.Rdata')
load('phenotype.Rdata')

######### 筛选肿瘤样本 #########
tumor_samples = phe[phe$`Tissue:ch1` == 'Liver Tumor Tissue', ]
tumor_samples = subset(tumor_samples, tumor_samples$`Tissue:ch1`!="NA")
sample_ids = row.names(tumor_samples)

tumor_exp_matrix = geo_exprSet[ , sample_ids]


library('immunedeconv')

########### quantiseq ###########

immunedeconv_res = deconvolute(tumor_exp_matrix, "quantiseq")

save(immunedeconv_res, file = './quantiseq_res.Rdata')

########### xcell ###########
immunedeconv_res = deconvolute(tumor_exp_matrix, "xcell")

save(immunedeconv_res, file = './xcell_res.Rdata')


########### cibersort ###########
source("../code/CIBERSORT.R")
#write.table(tumor_exp_matrix, file = 'tumor_exp_matrix.txt')
res = CIBERSORT('tumor_exp_matrix.txt', "LM22.txt", perm = 10, QN = F)


# set_cibersort_binary("../code/CIBERSORT.R")
# set_cibersort_mat("./LM22.txt")
# 
# immunedeconv_res = deconvolute(LIHC_PT_TPM_expr, "cibersort")

save(immunedeconv_res, file = './cibersort_res.Rdata')



########### mcp_counter ###########
immunedeconv_res = deconvolute(tumor_exp_matrix, "mcp_counter")

save(immunedeconv_res, file = './mcp_counter_res.Rdata')


########### epic ###########
immunedeconv_res = deconvolute(tumor_exp_matrix, "epic")

save(immunedeconv_res, file = './epic_res.Rdata')

########### TIMER ###########
immunedeconv_res = deconvolute(tumor_exp_matrix, "timer", indications = rep('LIHC', ncol(tumor_exp_matrix)))

save(immunedeconv_res, file = './timer_res.Rdata')
