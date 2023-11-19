###针对肿瘤样本的生存分析
###clinical信息整理
###gender、race等的K-M曲线绘制
###单基因/多基因K-M曲线绘制
###基因与生存率的p值计算

library("rjson")
library(tidyverse) #可读取tsv文件
library(survival)
library(survminer)
setwd('D:/Bioinfo/Transcriptome analysis/TCGA-Stomach/')

# 标准化后的表达矩阵
miRNA_exprSet_vst<-read.csv('miRNA_exprSet_vst.csv',row.names=1,
                            as.is=TRUE,check.names=FALSE)
group_list=ifelse(as.numeric(substr(colnames(miRNA_exprSet_vst),14,15)) < 10,'tumor','normal')
table(group_list)
# 仅取tumor样本，446个
miRNA_exprSet_tumor <- miRNA_exprSet_vst[,group_list=='tumor']

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
clinical=clinical[match(substr(colnames(miRNA_exprSet_tumor),1,12),rownames(clinical)),]
clinical <- clinical[!grepl("\\.", rownames(clinical)) & !rownames(clinical) == 'NA', ]
dim(clinical)
#临床样本427个

#生存分析
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

# gender下的K-M曲线
sfit <- survfit(Surv(time, event)~gender, data=clinical)
print(sfit)
ggsurvplot(sfit, conf.int=F, pval=TRUE)


# 基因表达生存分析
#1 差异基因
miRNA_up <- read.table('DESeq_up_miRNA.txt',header = TRUE)$miRNA_names
miRNA_down <- read.table('DESeq_down_miRNA.txt',header = TRUE)$miRNA_names

#2 差异基因在clinical样本中的表达情况:列是样本，行是miRNA
DE_miRNA_exp <- miRNA_exprSet_tumor
colnames(DE_miRNA_exp) <- substr(colnames(miRNA_exprSet_tumor),1,12)
clinical_sample <- rownames(clinical)
DE_miRNA_exp <- DE_miRNA_exp[c(miRNA_up,miRNA_down),clinical_sample]

#3 单基因表达生存分析
exprSet=DE_miRNA_exp  #套流程
g = rownames(exprSet)[150] # 随便选一个
clinical$gene = t(ifelse(exprSet[g,]>median(as.numeric(exprSet[g,])),'high','low'))
sfit1=survfit(Surv(time, event)~gene, data=clinical)
ggsurvplot(sfit1,pval =TRUE, data = clinical, risk.table = TRUE)

#4 指定多基因绘制在一张图上
gs=rownames(exprSet)[1:4]
splots <- lapply(gs, function(g){
  clinical$gene=t(ifelse(exprSet[g,]>median(as.numeric(exprSet[g,])),'high','low'))
  sfit1=survfit(Surv(time, event)~gene, data=clinical)
  ggsurvplot(sfit1,pval =TRUE, data = clinical, risk.table = TRUE)
}) 
arrange_ggsurvplots(splots, print = TRUE,  
                    ncol = 2, nrow = 2, risk.table.height = 0.4)


#批量计算所有表达矩阵的生存分析p值，从而挑选出“显著基因”
mySurv=with(clinical,Surv(time, event))
log_rank_p <- apply(exprSet , 1 , function(gene){
  # gene=exprSet[1,]
  clinical$group=ifelse(gene>median(gene),'high','low')  
  data.survdiff=survdiff(mySurv~group,data=clinical)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
  return(p.val)
})
#上述如果返回报错说只有一组，需要将一些低表达数据删除点
#expr = expr[apply(expr, 1, function(x) {
#sum(x > 1) > 9  #过滤掉低表达基因
#}), ]  
log_rank_p=data.frame(sort(log_rank_p))#miRNA及其与生存分析的p值
