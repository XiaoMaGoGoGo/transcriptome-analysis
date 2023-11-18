###差异分析

library(tidyverse)
library("rstatix")
library(dplyr)
library('ggrepel')
library('ggplot2')
library('ggsci')
# 读取表达矩阵及样本信息
miRNA_expression <- read.csv('miRNA_exprSet_vst.csv', 
                             row.names=1,as.is=TRUE,check.names=FALSE)#首行首列做列名和行名，不改变数据格式且不修改列名格式
miRNA_names <- rownames(miRNA_expression)
miRNA_expression <- cbind(miRNA_names,miRNA_expression)
sample_names <- colnames(miRNA_expression)
group_list <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')# 样本分组
sample_Class <- data.frame(Sample = sample_names, Class = group_list)

# miRNA按照tumor和normal分组整理
df = miRNA_expression %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "Sample",values_to = "value") %>%  # 转换成长数据
  left_join(sample_Class,by=c("Sample" = "Sample")) # 与分组数据合并

#方法一 基于FC的差异分析
#1 计算每个miRNA的FC log2FC
dfFC = df %>%
  group_by(miRNA_names,Class) %>%                         #按照样本类别分组
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算两类平均表达水平
  pivot_wider(names_from = Class,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = tumor/normal)                           # 实验组/对照组 计算差异倍数FC
dfFC$log2FC <- log2(dfFC$FC)

#2 根据阈值，确定差异基因
DE_miRNA_FC <- dfFC[dfFC$log2FC  > 1 | dfFC$log2FC  < -1, ]
#共133个


# 方法二、T检验下的差异分析
#1 对每一个miRNA 进行T检验并校正p值
dfT = df %>%
  group_by(miRNA_names) %>%
  t_test(value ~ Class,var.equal=T) # t_test 方法源于rstatix包; var.equal=T等方差
dfT= dfT %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  # p.adjust对P值FDR校正,算法选择BH
#2 根据p值，确定差异基因
DE_miRNA_T <- dfP[dfP$FDR<=0.05,]

#方法三、limma差异分析

library(edgeR) 
library(limma)
library(dplyr)
library("statmod")

#1 导入标准化后的表达数据，整理成limma需要的格式
miRNA_exp_limma <- miRNA_expression[,-1]#为计算每个miRNA，在之前的表达矩阵中增加了miRNA列，此处应删除

#2 把design中的样本信息转换成行名为样本名 ，列名为样本类型0/1
design <- model.matrix(~0+factor(group_list))
colnames(design)=levels(factor(group_list))
rownames(design)=colnames(miRNA_exp_limma)

#3 数据拟合上去
fit <- lmFit(miRNA_exp_limma,design) #将数据拟合上去
contrast.matrix <- makeContrasts(tumor - normal,levels = design) 

#4 非线性最小二乘法
fit <- contrasts.fit(fit, contrast.matrix)
fit <- eBayes(fit)#用经验贝叶斯调整t-test中方差的部分

#5 提取差异分析结果
dflimma <- topTable(fit, coef = 1,n = Inf,sort.by="logFC")
dflimma  <- na.omit(dflimma)
dflimma<- data.frame(dflimma, stringsAsFactors = FALSE, check.names = FALSE)

#6 根据p值筛选差异基因
DE_miRNA_limma <- dflimma[dflimma$adj.P.Val <=0.05 & !is.na(dflimma$adj.P.Val),]
#495个

#7 上调基因和下调基因
DE_miRNA_limma$change <-as.factor(ifelse(DE_miRNA_limma$logFC > 1 ,'Up', 
                                         ifelse(DE_miRNA_limma$logFC < -1,'Down','No change')))

table(DE_miRNA_limma$change)
#55下调和65上调

#8 绘制火山图这里的logFC就是log2FC
DE_miRNA_limma$miRNA_names <- rownames(DE_miRNA_limma)#添加miRNA_names列，用于绘图添加标签
ggplot(DE_miRNA_limma,aes(logFC, -log10(adj.P.Val)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  # 添加标签：
  geom_text_repel(data = filter(DE_miRNA_limma, abs(logFC) > 1 & -log10(adj.P.Val) > 38),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = miRNA_names, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(adj.P.val)"+
         xlim(-3, 3) + ylim(0, 35))


# 方法四、基于DESeq的差异分析

library(DESeq2)

#1 构建DESeq对象(用未标准化的表达数据)
miRNA_exp_DESeq <- read.csv('miRNA_expression.csv', 
                            row.names=1,as.is=TRUE,check.names=FALSE)#首行首列做列名和行名，不改变数据格式且不修改列名格式
sample_names <- colnames(miRNA_exp_DESeq)
group_list <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')# 样本分组
sample_Class <- data.frame(Sample = sample_names, Class = group_list)

dds <- DESeqDataSetFromMatrix(countData=miRNA_exp_DESeq,
                              colData=sample_Class,
                              design=~Class,
                              tidy=FALSE)
#2 执行差异分析
dds <- DESeq(dds)

#3 检查差异分析结果
dfDESeq <- results(dds, contrast = c('Class', 'tumor', 'normal'))# tumor在前，normal 在后，意为 tumor相较于 normal 中哪些基因上调/下调

dfDESeq  <- data.frame(dfDESeq , stringsAsFactors = FALSE, check.names = FALSE)
#4 根据p值筛选差异基因
DE_miRNA_DESeq <- dfDESeq[dfDESeq$padj <=0.05 & !is.na(dfDESeq$padj),]
#452个

#5 上调基因和下调基因
DE_miRNA_DESeq$change <-as.factor(ifelse(DE_miRNA_DESeq$log2FoldChange > 1 ,'Up', 
                                         ifelse(DE_miRNA_DESeq$log2FoldChange < -1,'Down','No change')))
table(DE_miRNA_DESeq$change)
#74下调和166个上调

#6 绘制火山图这里的logFC就是log2FC
DE_miRNA_DESeq$miRNA_names <- rownames(DE_miRNA_DESeq)#添加miRNA_names列，用于绘图添加标签
ggplot(DE_miRNA_DESeq,aes(log2FoldChange, -log10(padj)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  # 添加标签：
  geom_text_repel(data = filter(DE_miRNA_DESeq, abs(log2FoldChange) > 1 & -log10(padj) > 38),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = miRNA_names, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(padj)"+
         xlim(-3, 3) + ylim(0, 35))


# 方法五 edge
library(edgeR)

#1 数据整理
miRNA_exp_edge <- read.csv('miRNA_expression.csv', 
                            row.names=1,as.is=TRUE,check.names=FALSE)#首行首列做列名和行名，不改变数据格式且不修改列名格式
sample_names <- colnames(miRNA_exp_edge)
group_list <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')

#2 构建edge对象
dgelist <- DGEList(counts =miRNA_exp_edge, group=group_list)
dgelist <- calcNormFactors(dgelist,method='TMM')#TMM标准化
design <- model.matrix(~group_list)

#3 筛选差异基因
#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist, design, robust = TRUE)

#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型筛选差异基因
fit <- glmFit(dge, design, robust = TRUE)
dfedge <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
dfedge<- data.frame(dfedge, stringsAsFactors = FALSE, check.names = FALSE)

#4 根据p值筛选差异基因
DE_miRNA_edge <- dfedge[dfedge$FDR <=0.05 & !is.na(dfedge$FDR),]
#509个

#5 上调基因和下调基因
DE_miRNA_edge$change <-as.factor(ifelse(DE_miRNA_edge$logFC > 1 ,'Up', 
                                         ifelse(DE_miRNA_edge$logFC < -1,'Down','No change')))
table(DE_miRNA_edge$change)
#96下调和184个上调

#6 绘制火山图这里的logFC就是log2FC
DE_miRNA_edge$miRNA_names <- rownames(DE_miRNA_edge)#添加miRNA_names列，用于绘图添加标签
ggplot(DE_miRNA_edge,aes(logFC, -log10(FDR)))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1.2,1.2), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 0.2, 
             alpha = 0.5) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  # 添加标签：
  geom_text_repel(data = filter(DE_miRNA_edge, abs(logFC) > 1 & -log10(FDR) > 38),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = miRNA_names, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(FDR)"+
         xlim(-3, 3) + ylim(0, 35))


#五种方法下的差异分析结果保存 ---
write.table(DA_FC_p, 'DA_FC_p.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(DA_DEseq, 'DA_DEseq.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(DA_limma, 'DA_limma.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(lrt_glmQLFit, 'lrt_glmQLFit.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(lrt_glmLRT, 'lrt_glmLRT.txt', sep = '\t', col.names = NA, quote = FALSE)


#4-基于校正后的p分析DESeq结果，提取上调和下调miRNA提取
#基于校正后的p分析DESeq结果
col_A_index <- which(colnames(DA_DEseq) == "padj")
DESeq_p <- DA_DEseq[!is.na(DA_DEseq[, col_A_index]), ]
DESeq_p <- DESeq_p[DESeq_p$padj<0.05,]
DEseq_miRNA_sta <- 
  DESeq_p %>% 
  mutate(change = as.factor(ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                                   ifelse(log2FoldChange > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('miRNA_names')

#差异分析结果分别在DE_miRNA_FC/T/limma/DESeq/edge中
