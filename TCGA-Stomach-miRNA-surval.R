#*对DESeq得到的miRNA与样本做生存分析
#*表达矩阵整理-仅取tumor样本
#*临床样本整理
#*生存分析（gender、age、high-expression和low-expression
#*单因素Cox+多因素Cox
library("rjson")
library(tidyverse) #可读取tsv文件

setwd('D:/Bioinfo/Transcriptome analysis/TCGA-Stomach/')

# 表达矩阵
miRNA_exprSet_vst<-read.csv('miRNA_exprSet_vst.csv',header =FALSE)
colnames(miRNA_exprSet_vst) <- miRNA_exprSet_vst[1,]
miRNA_exprSet_vst <- miRNA_exprSet_vst[-1,]
rownames(miRNA_exprSet_vst) <- miRNA_exprSet_vst[,1]
miRNA_exprSet_vst <- miRNA_exprSet_vst[,-1]
group_list=ifelse(as.numeric(substr(colnames(miRNA_exprSet_vst),14,15)) < 10,'tumor','normal')
table(group_list)
#仅取tumor样本，446个
miRNA_exprSet_vst=miRNA_exprSet_vst[,group_list=='tumor']

# 临床信息
clinical <- read_tsv('clinical.cart.2023-10-30/clinical.tsv') #读取clinical文件
clinical <- clinical[!duplicated(clinical$case_submitter_id), ]
clinical <- column_to_rownames(clinical,var = "case_submitter_id")
clinical=clinical[,colnames(clinical) %in% c("vital_status",
                                 "days_to_last_follow_up",
                                 "days_to_death",
                                 "race",
                                 "gender",
                                 "age_at_index",
                                 "tumor_stage")]
clinical <- subset(clinical, vital_status == "Dead" | vital_status == "Alive")
dim(clinical) # 435个
clinical=clinical[match(substr(colnames(miRNA_exprSet_vst),1,12),rownames(clinical)),]
clinical <- clinical[!grepl("\\.", rownames(clinical)) & !rownames(clinical) == 'NA', ]
#临床样本435个

#生存分析
library(survival)
library(survminer)

# 生存分析
#1、生存时间
clinical$days <- as.numeric(ifelse(clinical$vital_status=='Alive',clinical$days_to_last_follow_up,clinical$days_to_death))
#时间以月份记，保留两位小数
clinical$time=round(clinical$days/30,2)
clinical <- na.omit(clinical) # 删除空值

#2、根据生死定义活着是0，死的是1
clinical$event=ifelse(clinical$vital_status=='Alive',0,1)
table(clinical$event)

#3 年龄分组(部分样本缺失，考虑可能的影响应该不大)
clinical$age_at_index[is.na(clinical$age_at_index)] <- 0
clinical$age_at_index=as.numeric(clinical$age_at_index)
clinical$age_group=ifelse(clinical$age_at_index>median(clinical$age_at_index),'older','younger')
table(clinical$age_group)

#4 癌症阶段 缺失
table(clinical$tumor_stage)

#5 race 人种
table(clinical$race)

#6 性别 gender
table(clinical$gender)

# gender
sfit <- survfit(Surv(time, event)~race, data=clinical)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)

#基因表达的生存分析
#1 指定感兴趣的基因——DESeq计算的上调基因
miRNA_up <- as.vector(read.table('DESeq_up_miRNA.txt',header = TRUE)$miRNA_names)
clinical_sample <- as.vector(rownames(clinical))
miRNA_up_exp <- miRNA_exprSet_vst[miRNA_up,grep(paste0("^", clinical_sample, collapse = "|"), names(miRNA_exprSet_vst))]
colnames(miRNA_up_exp) <- substr(colnames(miRNA_up_exp), 1, 12)
write.table(miRNA_up_exp,'miRNA_up_exp.txt',sep='\t')

#166个上调基因在427个样本中的表达矩阵
miRNA=miRNA_up[match(substr(colnames(miRNA_exprSet_vst),1,12),clinical_sample),]
miRNA_up_exp <- miRNA_exprSet_vst[miRNA_up,grep(paste0("^", clinical_sample, collapse = "|"), names(miRNA_exprSet_vst))]

# 找出重复的列名
duplicate_cols <- duplicated(names(miRNA_up_exp))
duplicate_col_names <- names(miRNA_up_exp)[duplicate_cols]

# 输出重复的列名
print(duplicate_col_names)

# 对重复的列进行平均值计算
for (col_name in duplicate_col_names) {
  duplicate_col <- df[[col_name]]
  avg_col <- rowMeans(duplicate_col, na.rm = TRUE)
  
  # 在数据框中添加新的列，命名为原始列名 + "_avg"
  new_col_name <- paste0(col_name, "_avg")
  df[[new_col_name]] <- avg_col
}

# 删除重复的列
df <- df[, !duplicate_cols]


#test
exprSet= miRNA_exprSet_vst #套流程
for (col in names(exprSet)) {
  exprSet[[col]] <- as.numeric(exprSet[[col]])
}
g = rownames(exprSet)[150] # 随便选一个
clinical$gene = ifelse(exprSet[1,]> median(exprSet[1,]),'high','low')
sfit1=survfit(Surv(time, event)~gene, data=clinical)
ggsurvplot(sfit1,pval =TRUE, data = clinical, risk.table = TRUE)




# 差异基因
miRNA_up <- read.table('DESeq_up_miRNA.txt',header = TRUE)
miRNA_down <- read.table('DESeq_down_miRNA.txt',header = TRUE)
miRNA_diff <- rbind(miRNA_up,miRNA_down)
rownames(miRNA_diff) <- miRNA_diff$miRNA_names#上调和下调的基因，共240个
miRNA_diff <- miRNA_diff[,!colnames(miRNA_diff) %in% "miRNA_names"]

#提取240个基因在clinical中的表达矩阵
#先提取基因
miRNA_diff_exp <- miRNA_exprSet_vst[match(rownames(miRNA_exprSet_vst),rownames(miRNA_diff)),]
miRNA_diff_exp <- miRNA_diff_exp[!grepl("\\.", rownames(miRNA_diff_exp)) & !rownames(miRNA_diff_exp) == 'NA', ]
#再保留样本
miRNA_diff_exp <- miRNA_diff_exp[match(substr(colnames(miRNA_diff_exp),1,12),rownames(clinical)),]





#差异基因

clinical$os_status <- ifelse(clinical$vital_status == "Alive", 0, 1) # 1 表示终点事件
clinical$os_time <- ifelse(clinical$vital_status=='Alive',clinical$days_to_last_follow_up,clinical$days_to_death)
clinical$os_time <- as.numeric(clinical$os_time)
#去除生存时间为0、NA和重复的数据
clinical <- subset(clinical, !duplicated(clinical))
clinical <- na.omit(clinical)

#读取对应的case_id和sample_id
json <- jsonlite::fromJSON("metadata.cart.2023-10-30.json")
sample_id <- sapply(json$associated_entities,function(x){x[,1]})
case_id <- sapply(json$associated_entities,function(x){x[,3]})
case_simple <- cbind(case_id,sample_id)
clinical <- merge(clinical, case_simple, by = "case_id", all.x = TRUE)

#读取miRNA
miRNA_exprSet_vst<-read.csv('miRNA_exprSet_vst.csv',header =FALSE)
colnames(miRNA_exprSet_vst) <- miRNA_exprSet_vst[1,]
miRNA_exprSet_vst <- miRNA_exprSet_vst[-1,]
rownames(miRNA_exprSet_vst) <- miRNA_exprSet_vst[,1]
miRNA_exprSet_vst <- miRNA_exprSet_vst[,-1]

#提取上调和下调miRNA（标准化后）
miRNA_down <- as.vector(read.csv('DESeq_down_miRNA.txt',sep='\t')$miRNA_names)
miRNA_up <- as.vector(read.csv('DESeq_up_miRNA.txt',sep='\t')$miRNA_names)
selected_miRNAs <- miRNA_exprSet_vst[miRNA_exprSet_vst$V1 %in%c(miRNA_down,miRNA_up), ]#根据行名提取上调和下调的miRNA
selected_miRNAs <- t(selected_miRNAs)
sample_id <- rownames(selected_miRNAs)
selected_miRNAs <- cbind(sample_id,selected_miRNAs)
data_for_Cox <- merge(clinical, selected_miRNAs, by = "sample_id", all.x = TRUE)
rownames(data_for_Cox) <- data_for_Cox$sample_id
cols_to_convert <- 10:ncol(data_for_Cox)

# 循环遍历需要转换的列
for (i in cols_to_convert) {
  data_for_Cox[, i] <- as.numeric(data_for_Cox[, i])
}
data_for_Cox$os_status <- as.integer(data_for_Cox$os_status)
data_for_Cox$os_time <- as.integer(data_for_Cox$os_time)

#2--单因素COX回归分析----
library('survival')
df <- data_for_Cox
#设置p值的阈值
pfilter <- 0.05   
#创建数据框，存放单因素Cox结果
Single_Cox_Allresult <- data.frame()  
#使用for循环对输入数据中的100个基因依次进行单因素COX分析
#单因素COX回归分析中p值＜0.05的基因，其分析结果输入到之前新建的空白数据框uniresult中
for(i in colnames(df[,10:ncol(df)])){   
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ df[,i], data = df)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  Single_Cox_Allresult <- rbind(Single_Cox_Allresult,
                      cbind(gene=i,
                            HR=unisum$coefficients[,2],
                            L95CI=unisum$conf.int[,3],
                            H95CI=unisum$conf.int[,4],
                            pvalue=unisum$coefficients[,5]
                       ))

}
Single_Cox_Selected <- Single_Cox_Allresult[Single_Cox_Allresult$pvalue <0.05,]
#得到20个显著的miRNA

#提取各样本的生存状态、生存时间、
#以及单因素COX回归分析中p值＜0.05的基因在各样本中的表达情况，作为多因素COX回归分析的输入数据
Single_Cox_Selected_unigene <- subset(df,select = c("os_status","os_time",Single_Cox_Selected$gene))

#3-多因素COX回归分析----  
multicox <- coxph(Surv(time = os_time,event = os_status) ~ ., data = Single_Cox_Selected_unigene) 
multisum <- summary(multicox)
coef_muti<-coef(multicox)

#提取20个miRNA的多因素COX回归分析结果至multiresult对象中
gene <- colnames(Single_Cox_Selected_unigene)[3:ncol(Single_Cox_Selected_unigene)]
HR <- multisum$coefficients[,2]
L95CI <- multisum$conf.int[,3]
H95CI <- multisum$conf.int[,4]
pvalue <- multisum$coefficients[,5]
multiresult <- data.frame(gene=gene,
                          HR=HR,
                          L95CI=L95CI,
                          H95CI=H95CI,
                          pvalue=pvalue,
                          coef = coef_muti)
multi_Cox_result <- multiresult[multiresult$pvalue<pfilter,]
# 得到6个miRNA，hsa-mir-137、hsa-mir-378a、hsa-mir-4777、hsa-mir-548v、hsa-mir-549a、hsa-mir-7641-1
# 分析6个miRNA在肿瘤/非肿瘤样本中的表达情况
miRNA_KM <- multi_Cox_result$gene










##确定6个miRNA在样本中的高/低表达量

##绘制6个miRNA下的k-M曲线
library('survminer')
miRNA_KM <- multiresult$gene
data_for_KM <- data_for_Cox[,c('os_time','os_status',miRNA_KM)]#生存时间和表达量标准化后的数据
#colnames(data_for_KM)
#str(data_for_KM)
#data_for_KM$time <- as.numeric(data_for_KM$os_time)
#data_for_KM$status <- as.numeric(data_for_KM$os_status)
#data_for_KM$sex <- as.numeric(data_for_KM$`hsa-mir-137`)
#Surv(data_for_KM$time,data_for_KM$status)
fit <- survfit(Surv(data_for_KM$os_time,data_for_KM$os_status)~data_for_KM$`hsa-mir-137`, data = data_for_KM)
#注意fit里的time
ggsurvplot(fit,
           conf.int = TRUE,# 显示置信区间
           linetype = "strata", # 根据性别分组自动设置曲线类型
           surv.median.line = "hv", # 设置中位生存期显示
           ggtheme = theme_bw(), # 设置ggplot2主题
           #palette = c("#E7B800", "#2E9FDF")
           )

#保存单因素COX回归分析结果
write.csv(uniresult,file = "result/单因素COX分析结果.csv",row.names = F)



#生存曲线绘制实例
d<- lung
str(lung)
fit <- survfit(Surv(time,status) ~ sex, data = d)
ggsurvplot(fit,
           conf.int = TRUE,# 显示置信区间
           linetype = "strata", # 根据性别分组自动设置曲线类型
           surv.median.line = "hv", # 设置中位生存期显示
           ggtheme = theme_bw(), # 设置ggplot2主题
           palette = c("#E7B800", "#2E9FDF"))