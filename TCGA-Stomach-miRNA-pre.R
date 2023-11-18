### TCGA下载的胃癌miRNA数据预处理：表达矩阵整合、去除空值、重复样本取均值、标准化

library("stringr")
library("rjson")

setwd('example/')

# miRNA表达矩阵

#1 读取meta下的json文件信息（文件名和样本名）
json <- jsonlite::fromJSON("metadata.cart.2023-10-30.json")
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
file_sample <- data.frame(sample_id,file_name=json$file_name)
#得到file_name和sample_id的对应矩阵

#2 合并文件
count_file <- list.files('gdc_download_20231030_063423.619995',pattern = '*mirnas.quantification.txt',recursive = TRUE)
count_file_name <- strsplit(count_file,split='/')
count_file_name <- sapply(count_file_name,function(x){x[2]})#返回列表每个元素的第二位（即文件名）用于与json文件比较

#3 先预定义一个矩阵，用于后面遍历文件时合并
path = paste0('gdc_download_20231030_063423.619995//',count_file[1])
data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)# 读取文件
colnames(data)<-data[1,]#将第一行作为列名
miRNA_expression<-data[,"read_count",drop=FALSE]#提取count列
colnames(miRNA_expression) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[1])]#添加样本名称

#4 遍历文件，读取数据后表示在一个矩阵中。所要提取的数据预预定义矩阵中的提取一致
for (i in 1:length(count_file_name)){
  path = paste0('gdc_download_20231030_063423.619995//',count_file[i])
  data<- read.delim(path,fill = TRUE,header = FALSE,row.names = 1)
  colnames(data)<-data[1,]
  data<-data[,"read_count",drop=FALSE]# 提取count列，一般而言，修改此处
  colnames(data) <- file_sample$sample_id[which(file_sample$file_name==count_file_name[i])]#添加样本名称
  miRNA_expression <- cbind(miRNA_expression,data)
}
miRNA_expression <- miRNA_expression[-1,]
miRNA_expression <- miRNA_expression[,-1]
miRNA_names <- rownames(miRNA_expression)
miRNA_expression <- apply(miRNA_expression, 2, function(x) as.numeric(as.character(x)))#数据类型转换
rownames(miRNA_expression) <- miRNA_names

#5 补充：根据样本名称，分组：tumor和normal

group_list <- ifelse(as.numeric(str_sub(colnames(miRNA_expression),14,15))<10,'tumor','normal')#样本的分组情况
table(group_list)
print(miRNA_expression[,group_list=='normal'])

#数据预处理

#1 数据清洗和缺失值处理
miRNA_expression <- miRNA_expression[rowSums(miRNA_expression)>0,] #保留表达量大于0的

#2 重复样本判断和取平均值
duplicate_samples <- colnames(miRNA_expression)[duplicated(colnames(miRNA_expression))]
  #若不为空，则需要对重复样本取平均值

#3 DESeq2表达数据标准化

#if (!requireNamespace("BiocManager", quietly = TRUE))
  #install.packages("BiocManager")
#BiocManager::install("DESeq2")
library(DESeq2)

sample_names <- colnames(miRNA_expression)
group_list <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')# 样本分组
sample_Class <- data.frame(Sample = sample_names, Class = group_list)
dds <- DESeqDataSetFromMatrix(countData=miRNA_expression,
                              colData=sample_Class,
                              design=~Class,
                              tidy=FALSE)
#构建DESeq对象
#vsd <- vst(dds, blind = TRUE)
#标准化
vsd <- varianceStabilizingTransformation(dds)

#结果保存
plotPCA(vsd, "Class")#PCA降维按class（tumor和normal）着色
miRNA_exprSet_vst <- as.data.frame(assay(vsd))
write.csv(miRNA_exprSet_vst,file="miRNA_exprSet_vst.csv",fileEncoding = "UTF-8")
write.csv(miRNA_exprSet_vst,file="miRNA_exprSet_vst.csv",fileEncoding = "UTF-8")
#***-2-----------标准化完成，包括去除在所有样本中0的miRNA和方差平稳转换------------


