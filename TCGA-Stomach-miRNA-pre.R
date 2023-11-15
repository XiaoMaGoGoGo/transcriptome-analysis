#***-1--TCGA下载的miRNA数据整合与分组--------------------------------------------------------

#nstall.packages("stringr")
library("stringr")
#setwd('Transcriptome analysis/example/data/')
#install.packages("rjson")
library("rjson")
json <- jsonlite::fromJSON("metadata.cart.2023-10-30.json")

#id <- json$associated_entities[[1]][,1]
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)  
#得到file_name和sample_id的对应矩阵，然后将sample_id添加到对应文件名的矩阵中。

count_file <- list.files('gdc_download_20231030_063423.619995',pattern = '*mirnas.quantification.txt',recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})

# 先预定义一个矩阵，用于后面遍历文件时合并
path = paste0('gdc_download_20231030_063423.619995//',count_file[1])
data0<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)# 读取文件
colnames(data0)<-data0[1,]#将第一行作为列名
matrix<-data0[,"read_count",drop=FALSE]#提取count列
colnames(matrix) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[1])]#添加样本名称

# 遍历文件，读取数据后表示在一个矩阵中。所要提取的数据预预定义矩阵中的提取一致
for (i in 1:length(count_file_name)){
  path = paste0('gdc_download_20231030_063423.619995//',count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[1,]
  data<-data[,"read_count",drop=FALSE]# 提取count列，一般而言，修改此处
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]#添加样本名称
  matrix <- cbind(matrix,data)
}
matrix <- matrix[-1,]
matrix <- matrix[,-1]#miRNA在样本中的表达情况
miRNA_expression <- apply(matrix, 2, function(x) as.numeric(as.character(x)))

#----------------根据样本名称，分组：tumor和normal------------------------

#table(str_sub(colnames(matrix0)))
group_list <- ifelse(as.numeric(str_sub(colnames(matrix),14,15))<10,'tumor','normal')#样本的分组情况
table(group_list)
matrix_miRNA_grouped <- rbind(matrix,group_list)
write.csv(matrix_miRNA_grouped,'Group_miRNA_stomach_matrix.csv',row.names = TRUE)

#***-1----------------TCGA数据预处理与tumor和normal分组完成！-----------------


#***-2--------------DESeq2表达数据标准化-------------------------------------

#BiocManager::install("DESeq2")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DESeq2")
library(DESeq2)

#2-1获取整合好的miRNA表达矩阵miRNA_expression
#miRNA_names <- rownames(matrix)
#miRNA_expression <- cbind(miRNA_names,matrix)#miRNA表达数据
miRNA_expression <- apply(matrix, 2, function(x) as.numeric(as.character(x)))
rownames(miRNA_expression) <- rownames(matrix)

#2-2获取样本类别
sample_names <- colnames(matrix)
#table(substring(sample_names,14,15)) #看下样本类型
group_list <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')# 样本分组
sample_Class <- data.frame(Sample = sample_names, Class = group_list)

#2-3 标准化
miRNA_expression <- miRNA_expression[rowSums(miRNA_expression)>0,]#预处理，去除在所有样本中表达量为0的值

dds <- DESeqDataSetFromMatrix(countData=miRNA_expression,
                             colData=sample_Class,
                             design=~Class,
                             tidy=FALSE)

#vsd <- vst(dds, blind = TRUE)
vsd <- varianceStabilizingTransformation(dds)
plotPCA(vsd, "Class")#PCA降维按class（tumor和normal）着色
miRNA_exprSet_vst <- as.data.frame(assay(vsd))
write.csv(miRNA_exprSet_vst,file="miRNA_exprSet_vst.csv",fileEncoding = "UTF-8")
#***-2-----------标准化完成，包括去除在所有样本中0的miRNA和方差平稳转换------------


