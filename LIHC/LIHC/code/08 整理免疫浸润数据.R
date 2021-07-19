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
FPKMFiles <- dir(path = ".\\",pattern="gz$")
library(R.utils)
sapply(FPKMFiles,gunzip)


##########################  处理json文件  #############################
library(rjson)
metadata_json_File <- fromJSON(file="..\\metadata.cart.2021-07-14.json")
json_File_Info <- data.frame(filesName =c(),TCGA_Baecode =c())
for(i in 1:length(metadata_json_File)){
  TCGA_Barcode <- metadata_json_File[[i]][["associated_entities"]][[1]][["entity_submitter_id"]]
  file_name <-metadata_json_File[[i]][["file_name"]]
  json_File_Info <- rbind(json_File_Info,data.frame(filesName = file_name,TCGA_Barcode =TCGA_Barcode))
}
rownames(json_File_Info) <- json_File_Info[,1]
#write.csv(json_File_Info,file = "..\\json_File_Info.csv")


##########################  获取FPKM矩阵  ###########################
filesName_TO_TCGA_BarcodeFile <- json_File_Info[-1]
FPKMFileNames <- dir(pattern="FPKM.txt$")

allSampleRawFPKM <- data.frame()
for(txtFile in FPKMFileNames){
  #每一个循环读取一个文件
  SampleFPKM <- read.table(txtFile,header = FALSE)
  rownames(SampleFPKM) <- SampleFPKM[,1]
  SampleFPKM <- SampleFPKM[-1]
  #根据fileName_TO_TCGA_BarcodeFile文件中文件名称与barcode对应关系，命名列名
  colnames(SampleFPKM) <- filesName_TO_TCGA_BarcodeFile[paste(txtFile,".gz",sep = ""),]
  if(dim(allSampleRawFPKM)[1]==0){
    allSampleRawFPKM <- SampleFPKM
  }
  else
  {allSampleRawFPKM <- cbind(allSampleRawFPKM,SampleFPKM)}
}
#write.csv(allSampleRawFPKM,file = "..\\allSampleRawFPKM.csv")
ensembl_id <- substr(row.names(allSampleRawFPKM),1,15)
rownames(allSampleRawFPKM) <- ensembl_id
#RawFPKM.cs文件与allSampleRawFPKM.csv文件的区别在于行名的ensembl去掉了版本号
#write.csv(allSampleRawFPKM,file = "..\\RawFPKM.csv")


##############################  ID转换  ###############################
#添加一列Ensemble_ID到RawFPKM数据框中
RawFPKM <- allSampleRawFPKM
Ensembl_ID <- data.frame(Ensembl_ID = row.names(RawFPKM))
rownames(Ensembl_ID) <- Ensembl_ID[,1]
RawFPKM <- cbind(Ensembl_ID,RawFPKM)

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
Ensembl_ID_TO_Genename <- get_map("../data/gencode.v34lift37.annotation.gtf")

gtf_Ensembl_ID <- substr(Ensembl_ID_TO_Genename[,1],1,15)
Ensembl_ID_TO_Genename <- data.frame(gtf_Ensembl_ID,Ensembl_ID_TO_Genename[,2])
colnames(Ensembl_ID_TO_Genename) <- c("Ensembl_ID","gene_id")
#write.csv(Ensembl_ID_TO_Genename,file = ".\\Ensemble_ID_TO_Genename.csv")

#融合数据
mergeRawFPKM <- merge(Ensembl_ID_TO_Genename,RawFPKM,by="Ensembl_ID")

#按照gene_id列进行排序
mergeRawFPKM <- mergeRawFPKM[order(mergeRawFPKM[,"gene_id"]),]
#根据gene_id列进行索引
index<-duplicated(mergeRawFPKM$gene_id)
mergeRawFPKM <- mergeRawFPKM[!index,]
#利用基因名称作为行名
rownames(mergeRawFPKM) <- mergeRawFPKM[,"gene_id"]
#删除前两列保存文件
LIHC_FPKM_expMatrix <- mergeRawFPKM[,-c(1:2)]
write.csv(LIHC_FPKM_expMatrix,file = ".\\LIHC_FPKM_expMatrix.csv")

save(LIHC_FPKM_expMatrix,file=".\\LIHC_FPKM_expMatrix.Rdata")

fpkmToTpm <- function(fpkm) {
  exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
}


LIHC_TPM_expMatrix <- apply(LIHC_FPKM_expMatrix, 2, fpkmToTpm)

write.csv(LIHC_TPM_expMatrix,file = ".\\LIHC_TPM_expMatrix.csv")

save(LIHC_TPM_expMatrix,file=".\\LIHC_TPM_expMatrix.Rdata")

