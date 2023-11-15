#***-3---差异分析--------------

library(tidyverse)
library("rstatix")
library(dplyr)
library('ggrepel')
library('ggplot2')
library('ggsci')

#1.整理样本
# 读取样本信息:表达量+分组
miRNA_names <- rownames(miRNA_exprSet_vst)
miRNA_exp_ <- cbind(miRNA_names,miRNA_exprSet_vst)
sample_names <- colnames(miRNA_exprSet_vst)
group_list <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')# 样本分组
sample_Class <- data.frame(Sample = sample_names, Class = group_list)

#2 基于FC+T检验的 差异分析

# 表达量+分组 合并，整理
df = miRNA_exp_ %>%
  as_tibble() %>%
  pivot_longer(-1,names_to = "Sample",values_to = "value") %>%  # 转换成长数据
  left_join(sample_Class,by=c("Sample" = "Sample")) # 与分组数据合并
#将value转换成数值型
df$value <- as.numeric(df$value)
dfFC = df %>%
  group_by(miRNA_names,Class) %>%  
  summarise(mean = mean(value,na.rm=T)) %>%               # 计算平均值
  pivot_wider(names_from = Class,values_from = mean) %>%  # 转换成宽数据
  summarise(FC = tumor/normal)                           # 实验组/对照组 计算差异倍数FC
#基于T检验，计算p值
dfP = df %>%
  group_by(miRNA_names) %>%
  t_test(value ~ Class,var.equal=T) # t_test 方法源于rstatix包; var.equal=T等方差
#对p进行校正
dfP_FDR = dfP %>%
  select(1,last_col()) %>%
  mutate(FDR = p.adjust(.$p,method = "BH"))  # p.adjust对P值FDR校正,算法选择BH
DA_FC_p <- cbind(dfFC,dfP_FDR)
DA_FC_p[, "Log2_FC"] <- log2(DA_FC_p[, "FC"])
DA_FC_p <- DA_FC_p[, !(names(DA_FC_p) %in% "miRNA_names.1")]
# 绘制火山图
data_FC_p<- 
  DA_FC_p %>% 
  mutate(change = as.factor(ifelse(p < 0.05 & abs(Log2_FC) > 1,
                                   ifelse(Log2_FC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('number')
table(data_FC_p$change)

ggplot(data_FC_p,aes(Log2_FC, -log10(FDR)))+
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
  geom_text_repel(data = filter(data_FC_p, abs(Log2_FC) > 1 & -log10(FDR) > 10),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = miRNA_names, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)"+
         xlim(-3, 3) + ylim(0, 35))

#基于DESeq的差异分析

dds <- DESeq(dds)
sizeFactors(dds)

res <- results(dds, contrast = c('Class', 'tumor', 'normal'))# tumor在前，normal 在后，意为 tumor相较于 normal 中哪些基因上调/下调
DA_DEseq <- data.frame(res, stringsAsFactors = FALSE, check.names = FALSE)

data_DESeq <- 
  DA_DEseq %>% 
  mutate(change = as.factor(ifelse(pvalue < 0.05 & abs(log2FoldChange) > 1,
                                   ifelse(log2FoldChange > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('miRNA_names')
table(data_DESeq$change)

ggplot(data_DESeq,aes(log2FoldChange, -log10(padj)))+
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
  geom_text_repel(data = filter(data_DESeq, abs(log2FoldChange) > 1 & -log10(padj) > 38),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = miRNA_names, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(padj)"+
         xlim(-3, 3) + ylim(0, 35))


#基于limma的差异分析
#if (!require("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("edgeR")
if (!requireNamespace("airway", quietly = TRUE))
  BiocManager::install("airway")
library(edgeR) 
library(limma)
library(dplyr)
install.packages("statmod")
library("statmod")

#导入miRNA标准化后的矩阵（去不表达+方差平稳转换）和分组信息
miRNA_exp_limma <- miRNA_exprSet_vst
sample_names <- colnames(miRNA_exprSet_vst)
group <- ifelse(as.numeric(substring(sample_names,14,15))<10,'tumor','normal')# 样本分组
#suppressMessages(library(limma))
#把design中的样本信息转换成行名为样本名 ，列名为样本类型0/1
design <- model.matrix(~0+factor(group))
colnames(design)=levels(factor(group))
rownames(design)=colnames(miRNA_exp_limma)
fit <- lmFit(miRNA_exp_limma,design) #将数据拟合上去
contrast.matrix <- makeContrasts(tumor - normal,levels = design) 
#非线性最小二乘法
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)#用经验贝叶斯调整t-test中方差的部分
DA_limma <- topTable(fit2, coef = 1,n = Inf,sort.by="logFC")
DA_limma  <- na.omit(DA_limma )
head(DA_limma)
DA_limma_<- data.frame(DA_limma, stringsAsFactors = FALSE, check.names = FALSE)
data_limma_ <- 
  DA_limma_ %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('miRNA_names')
table(data_limma_$change)
#这里的logFC就是log2FC
ggplot(data_limma_,aes(logFC, -log10(adj.P.Val)))+
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
  geom_text_repel(data = filter(data_limma_, abs(logFC) > 1 & -log10(adj.P.Val) > 38),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = miRNA_names, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(adj.P.val)"+
         xlim(-3, 3) + ylim(0, 35))

#3-4edgeR筛选差异基因
library(edgeR)
dgelist <- DGEList(counts =miRNA_exprSet_vst, group=group)
dgelist_norm <- calcNormFactors(dgelist,method='TMM')#TMM标准化
design <- model.matrix(~group)

#（1）估算基因表达值的离散度
dge <- estimateDisp(dgelist_norm, design, robust = TRUE)

#（2）模型拟合，edgeR 提供了多种拟合算法
#负二项广义对数线性模型
fit <- glmFit(dge, design, robust = TRUE)
lrt <- topTags(glmLRT(fit), n = nrow(dgelist$counts))
lrt_glmLRT<- data.frame(lrt, stringsAsFactors = FALSE, check.names = FALSE)
data_edgeR_LRT_ <- 
  lrt_glmLRT %>% 
  mutate(change = as.factor(ifelse(PValue < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('miRNA_names')
table(data_limma_$change)
#拟似然负二项广义对数线性模型
fit <- glmQLFit(dge, design, robust = TRUE)
lrt <- topTags(glmQLFTest(fit), n = nrow(dgelist$counts))
write.table(lrt, 'control_treat.glmQLFit.txt', sep = '\t', col.names = NA, quote = FALSE)
lrt_glmQLFit <-data.frame(lrt, stringsAsFactors = FALSE, check.names = FALSE)
data_edgeR_QLFit_ <- 
  lrt_glmQLFit %>% 
  mutate(change = as.factor(ifelse(PValue < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('miRNA_names')
table(data_limma_$change)

#4-差异分析结果保存 ---
write.table(DA_FC_p, 'DA_FC_p.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(DA_DEseq, 'DA_DEseq.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(DA_limma, 'DA_limma.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(lrt_glmQLFit, 'lrt_glmQLFit.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(lrt_glmLRT, 'lrt_glmLRT.txt', sep = '\t', col.names = NA, quote = FALSE)


#5-基于校正后的p分析DESeq结果，提取上调和下调miRNA提取
#基于校正后的p分析DESeq结果
col_A_index <- which(colnames(DA_DEseq) == "padj")
DESeq_p <- DA_DEseq[!is.na(DA_DEseq[, col_A_index]), ]
DESeq_p <- DESeq_p[DESeq_p$padj<0.05,]
DEseq_miRNA_sta <- 
  DESeq_p %>% 
  mutate(change = as.factor(ifelse(padj < 0.05 & abs(log2FoldChange) > 1,
                                   ifelse(log2FoldChange > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('miRNA_names')
table(DEseq_miRNA_sta$change)

#按照上调、下调、不变对miRNA分组保存
DESeq_up_miRNA <- subset(DEseq_miRNA_sta, change == "Up")#change是因子列，基于银子列分组
DESeq_down_miRNA <- subset(DEseq_miRNA_sta, change == "Down")
DESeq_nochange_miRNA <- subset(DEseq_miRNA_sta, change == "No change")
write.table(DESeq_up_miRNA, 'DESeq_up_miRNA.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(DESeq_down_miRNA, 'DESeq_down_miRNA.txt', sep = '\t', col.names = NA, quote = FALSE)
write.table(DESeq_nochange_miRNA, 'DESeq_nochange_miRNA.txt', sep = '\t', col.names = NA, quote = FALSE)
