options(stringsAsFactors = F)

#############  将所有文件移动到同一个文件夹SampleFiles下  #############
dir.create("SampleFiles")
filepath <- dir(path =".\\Gene_Expression_Quantification",full.names = TRUE)
for(wd in filepath){
  files <- dir(path =wd,pattern = "gz$") #查看满足条件的文件
  fromfilepath <- paste(wd,"\\",files,sep ="")
  tofilepath <- paste(".\\SampleFiles\\",files,sep="")
  file.copy(fromfilepath,tofilepath)
}

#######################  解压所有文件并删除原文件  ####################
setwd(".\\SampleFiles")
countsFiles <- dir(path = ".\\",pattern="gz$")
library(R.utils)
sapply(countsFiles,gunzip)


##########################  处理json文件  #############################
library(rjson)
metadata_json_File <- fromJSON(file="..\\metadata.cart.2021-07-09.json")
json_File_Info <- data.frame(filesName =c(),TCGA_Baecode =c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <-metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info,data.frame(filesName = file_name,TCGA_Barcode =TCGA_Barcode))
}
rownames(json_File_Info) <- json_File_Info[,1]
write.csv(json_File_Info,file = "..\\json_File_Info.csv")


##########################  获取Counts矩阵  ###########################
filesName_TO_TCGA_BarcodeFile <- json_File_Info[-1]
countsFileNames <- dir(pattern="counts$")

allSampleRawCounts <- data.frame()
for(txtFile in countsFileNames){
  #每一个循环读取一个文件
  SampleCounts <- read.table(txtFile,header = FALSE)
  rownames(SampleCounts) <- SampleCounts[,1]
  SampleCounts <- SampleCounts[-1]
  #根据fileName_TO_TCGA_BarcodeFile文件中文件名称与barcode对应关系，命名列名
  colnames(SampleCounts) <- filesName_TO_TCGA_BarcodeFile[paste(txtFile,".gz",sep = ""),]
  if(dim(allSampleRawCounts)[1]==0){
    allSampleRawCounts <- SampleCounts
  }
  else
  {allSampleRawCounts <- cbind(allSampleRawCounts,SampleCounts)}
}
#write.csv(allSampleRawCounts,file = "..\\allSampleRawCounts.csv")
ensembl_id <- substr(row.names(allSampleRawCounts),1,15)
rownames(allSampleRawCounts) <- ensembl_id
#RawCounts.cs文件与allSampleRawCounts.csv文件的区别在于行名的ensembl去掉了版本号
write.csv(allSampleRawCounts,file = "..\\RawCounts.csv")


##############################  ID转换  ###############################
#添加一列Ensemble_ID到RawCounts数据框中
RawCounts <- allSampleRawCounts
Ensembl_ID <- data.frame(Ensembl_ID = row.names(RawCounts))
rownames(Ensembl_ID) <- Ensembl_ID[,1]
RawCounts <- cbind(Ensembl_ID,RawCounts)

#一个函数，通过gtf文件获取Ensemble_ID与基因名称的对应关系
get_map = function(input){
  if(is.character(input)){
    if(!file.exists(input))stop("Bad input file.")
    message("Treat input as file")
    input = data.table::fread(input,header = FALSE)
  }else{
    data.table::setDT(input)
  }
  
  input = input[input[[3]] == "gene",]
  
  pattern_id = ".*gene_id \"([^;]+)\";.*"
  pattern_name = ".*gene_name \"([^;]+)\";.*"
  
  gene_id = sub(pattern_id, "\\1", input[[9]])
  gene_name = sub(pattern_name, "\\1", input[[9]])
  
  Ensembl_ID_TO_Genename <- data.frame(gene_id = gene_id,
                                       gene_name = gene_name,
                                       stringsAsFactors = FALSE)
  return(Ensembl_ID_TO_Genename)
  
}

setwd('../')
Ensembl_ID_TO_Genename <- get_map("../../../gtf_data/gencode.v34lift37.annotation.gtf")

gtf_Ensembl_ID <- substr(Ensembl_ID_TO_Genename[,1],1,15)
Ensembl_ID_TO_Genename <- data.frame(gtf_Ensembl_ID,Ensembl_ID_TO_Genename[,2])
colnames(Ensembl_ID_TO_Genename) <- c("Ensembl_ID","gene_id")
write.csv(Ensembl_ID_TO_Genename,file = ".\\Ensemble_ID_TO_Genename.csv")

#融合数据
mergeRawCounts <- merge(Ensembl_ID_TO_Genename,RawCounts,by="Ensembl_ID")

#按照gene_id列进行排序
mergeRawCounts <- mergeRawCounts[order(mergeRawCounts[,"gene_id"]),]
#根据gene_id列进行索引
index<-duplicated(mergeRawCounts$gene_id)
mergeRawCounts <- mergeRawCounts[!index,]
#利用基因名称作为行名
rownames(mergeRawCounts) <- mergeRawCounts[,"gene_id"]
#删除前两列保存文件
LIHC_Count_expMatrix <- mergeRawCounts[,-c(1:2)]
write.csv(LIHC_Count_expMatrix,file = ".\\LIHC_Count_expMatrix.csv")

save(LIHC_Count_expMatrix,file=".\\LIHC_Count_expMatrix.Rdata")


lncRNAs = get_map('../../../gtf_data/gencode.v38.long_noncoding_RNAs.gtf')

lncRNAs_in_expMatrix = intersect(rownames(LIHC_Count_expMatrix), lncRNAs$gene_name)
lncRNA_matrix = LIHC_Count_expMatrix[rownames(LIHC_Count_expMatrix) %in% lncRNAs_in_expMatrix, ]

notLncRNA_matrix = LIHC_Count_expMatrix[! rownames(LIHC_Count_expMatrix) %in% lncRNAs_in_expMatrix, ]

lncRNA_matrix = t(lncRNA_matrix)
notLncRNA_matrix = t(notLncRNA_matrix)

#筛出有25%样本表达量大于10的基因
lncRNA_matrix = lncRNA_matrix[, which(apply(lncRNA_matrix,2,function(x){return(sum(x>10))} ) > nrow(lncRNA_matrix)*0.25)]

save(lncRNA_matrix, notLncRNA_matrix, file = './lncRNA_matrix.Rdata')

####### 选择花生四烯酸基因 #####

arachidonic_gene = c("PLA2G10","PLA2G2D","PLA2G2E","PLA2G3",
               "PLA2G2F","PLA2G12A","PLA2G12B","PLA2G1B",
               "PLA2G5","PLA2G2A","PLA2G2C","PLA2G4E","PLA2G4A",
               "JMJD7-PLA2G4B","PLA2G4B","PLA2G4C","PLA2G4D",
               "PLA2G4F","PLA2G6","PLB1","PTGS1","PLAAT3",
               "PTGS2","PTGES","PTGES2","PTGES3","CBR1","CBR3",
               "PRXL2B","TBXAS1","PTGDS","HPGDS","AKR1C3",
               "PTGIS","ALOX5","LTA4H","CYP4F2","CYP4F3","LTC4S",
               "GGT1","GGT5","GPX6","GPX7","GPX2","GPX3","GPX1",
               "GPX5","GPX8","CYP2E1","CYP2J2","CYP2U1","CYP4A11",
               "CYP2C19","CYP4F8","ALOX12","ALOX12B","ALOX15B","CYP2B6",
               "CYP2C8","CYP2C9","EPHX2","ALOX15")
arachidonicRNA_matrix = notLncRNA_matrix[ , arachidonic_gene]


library(psych)
#corr = corr.test(lncRNA_matrix, notLncRNA_matrix, use = "pairwise", method = "pearson", adjust = "none")

########### 将相关系数与置信度满足条件的lncRNA筛选出来 ###########
p_matrix = c()
r_matrix = c()
rm(corr)
for(i in 1 : length(colnames(lncRNA_matrix))) {
  print(i)
  corr = corr.test(lncRNA_matrix[ , i], arachidonicRNA_matrix, use = "pairwise", method = "pearson", adjust = "none")

  p_result = cbind("gene" = colnames(lncRNA_matrix)[i], "p" = corr$p)
  r_result = cbind("gene" = colnames(lncRNA_matrix)[i], "r" = corr$r)
  
  p_matrix = rbind(p_matrix, p_result)
  r_matrix = rbind(r_matrix, r_result)
}

row.names(p_matrix) = p_matrix[ , c('gene')]
p_matrix = p_matrix[ , -which('gene' %in% colnames(p_matrix))]
row.names(r_matrix) = r_matrix[ , c('gene')]
r_matrix = r_matrix[ , -which('gene' %in% colnames(r_matrix))]

row_names = row.names(p_matrix)

p_matrix = apply(p_matrix, 2, as.numeric)
r_matrix = apply(r_matrix, 2, as.numeric)

row.names(p_matrix) = row_names
row.names(r_matrix) = row_names

save(p_matrix, r_matrix, file = './correlation_coefficient_matrix.Rdata')

p_threshold = 0.001
r_threshold = 0.5
relative_lncRNAs = c()
for(i in 1 : length(row.names(p_matrix))) {
  p_genes = colnames(p_matrix)[p_matrix[i, ] < p_threshold]
  r_genes = colnames(r_matrix)[abs(r_matrix[i, ]) > r_threshold]
  if(length(intersect(p_genes, r_genes)) > 0) {
    relative_lncRNAs = c(relative_lncRNAs, row.names(p_matrix)[i])
  }
}

lncRNA_matrix_to_cal = lncRNA_matrix[ , relative_lncRNAs]

lncRNA_matrix_to_cal = t(lncRNA_matrix_to_cal)

save(lncRNA_matrix_to_cal, file = './lncRNA_matrix_to_cal.Rdata')


########## 读取临床TSV文件 ##########

rm(list = ls())
options(stringsAsFactors = F)
library('readr')
clinical_tsv = read_tsv('./clinical.cart.2021-07-09/clinical.tsv')
clinical_tsv = as.data.frame(clinical_tsv)

########## 清洗数据 ##########
clinical_tsv = clinical_tsv[!duplicated(clinical_tsv$case_submitter_id) , ]
clinical_tsv = clinical_tsv[ , c('case_submitter_id', 'age_at_index', 'gender', 'ajcc_pathologic_m', 'ajcc_pathologic_n', 'ajcc_pathologic_stage', 'ajcc_pathologic_t', 'vital_status', 'days_to_last_follow_up', 'days_to_death')]
clinical_tsv = clinical_tsv[clinical_tsv$age_at_index > 0, ]
clinical_tsv = clinical_tsv[clinical_tsv$gender == 'male' | clinical_tsv$gender == 'female', ]
clinical_tsv = clinical_tsv[clinical_tsv$vital_status == 'Alive'| clinical_tsv$vital_status == 'Dead', ]
clinical_tsv = clinical_tsv[(clinical_tsv$vital_status == 'Alive' & clinical_tsv$days_to_last_follow_up > 0) | (clinical_tsv$vital_status == 'Dead' & clinical_tsv$days_to_death > 0), ]

#grepl('^M\\d\\w?$', 'M0')
clinical_tsv = clinical_tsv[grepl('^M(\\d|X)\\w?$', clinical_tsv$ajcc_pathologic_m), ]
clinical_tsv = clinical_tsv[grepl('^N(\\d|X)\\w?$', clinical_tsv$ajcc_pathologic_n), ]
clinical_tsv = clinical_tsv[grepl('^T(\\d|X)\\w?$', clinical_tsv$ajcc_pathologic_t), ]
clinical_tsv = clinical_tsv[grepl('^Stage\\s\\w+$', clinical_tsv$ajcc_pathologic_stage), ]

########## 从clinical_tsv装载数据到clinical ##########
clinical = as.data.frame(matrix(nrow = nrow(clinical_tsv), ncol = 13))
colnames(clinical) = c('id', 'age', 'agegroup', 'gender', 'pM', 'pN', 'pN1', 'stage', 'stage1', 'pT', 'pT1', 'OS.time', 'OS')
clinical$id = clinical_tsv$case_submitter_id
clinical$age = as.numeric(clinical_tsv$age_at_index)
clinical$agegroup = ifelse(clinical$age >= 66, '>=66', '<66')
clinical$gender = clinical_tsv$gender

clinical$pM = regmatches(clinical_tsv$ajcc_pathologic_m, regexpr('M(\\d|X)', clinical_tsv$ajcc_pathologic_m))
clinical$pN = regmatches(clinical_tsv$ajcc_pathologic_n, regexpr('N(\\d|X)', clinical_tsv$ajcc_pathologic_n))
clinical$pN1 = ifelse(clinical$pN == 'N1' | clinical$pN == 'N2' | clinical$pN == 'N3' | clinical$pN == 'NX', 'N1-N3', clinical$pN)
clinical$pT = regmatches(clinical_tsv$ajcc_pathologic_t, regexpr('T(\\d|X)', clinical_tsv$ajcc_pathologic_t))
clinical[clinical$pT == 'T1' | clinical$pT == 'T2', ]$pT1 = 'T1-T2'
clinical[clinical$pT == 'T3' | clinical$pT == 'T4' | clinical$pT == 'TX', ]$pT1 = 'T3-T4'

clinical$stage = regmatches(clinical_tsv$ajcc_pathologic_stage, regexpr('Stage\\s[IV]+', clinical_tsv$ajcc_pathologic_stage))
stage_list = c('Stage I', 'Stage II', 'Stage III', 'Stage IV')
for(i in 1:nrow(clinical)) {
  clinical[i, ]$stage1 = which(stage_list == clinical[i, ]$stage)
}

clinical$OS = clinical_tsv$vital_status
clinical[clinical$OS == 'Alive', ]$OS = 0
clinical[clinical$OS == 'Dead', ]$OS = 1
clinical[clinical$OS == 0, ]$OS.time = clinical_tsv[clinical_tsv$vital_status == 'Alive', ]$days_to_last_follow_up
clinical[clinical$OS == 1, ]$OS.time = clinical_tsv[clinical_tsv$vital_status == 'Dead', ]$days_to_death
clinical$OS.time = as.numeric(clinical$OS.time)
clinical$OS = as.numeric(clinical$OS)

rownames(clinical) = clinical$id
write.csv(clinical,'clinical_change.csv')
