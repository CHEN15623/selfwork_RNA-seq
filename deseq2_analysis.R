# 加载必要的包
library(DESeq2)
library(tidyverse)

# 设置输出目录
out_dir <- "result/deseq2/"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 读取表达矩阵
counts <- read.delim("result/featureCounts/gene_counts.txt", comment.char = "#", check.names = FALSE)
counts <- counts[, -c(2:6)]  # 移除Chr, Start, End, Strand, Length列
rownames(counts) <- counts$Geneid
counts <- counts[, -1]

# 简化列名为样本名
colnames(counts) <- gsub("temp/STAR/|_dataAligned.sortedByCoord.out.bam", "", colnames(counts))

# 读取样本信息
coldata <- read.delim("result/mydata.txt", row.names = 1)

# 确保样本顺序匹配
counts <- counts[, rownames(coldata)]

# 构建 DESeq2 数据对象
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ Group)

# 预处理
dds <- dds[rowSums(counts(dds)) > 1, ]  # 过滤低表达基因

# 差异表达分析
dds <- DESeq(dds)
res <- results(dds, contrast = c("Group", "PA", "BC"))  # 比较PA vs BC

# 排序并保存结果
res_ordered <- res[order(res$padj), ]
res_ordered$FoldChange <- 2^res_ordered$log2FoldChange
write.csv(as.data.frame(res_ordered), file = file.path(out_dir, "DEG_results.csv"))

# 提取显著差异基因（padj < 0.05 且 |log2FC| > 1）
res_sig <- subset(res_ordered, pvalue < 0.05 & abs(log2FoldChange) > 0.8)#根据需求设置阈值
res_sig$FoldChange <- 2^res_sig$log2FoldChange
write.csv(as.data.frame(res_sig), file = file.path(out_dir, "DEG_significant.csv"))

# 保存归一化后的counts
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(as.data.frame(normalized_counts), file = file.path(out_dir, "normalized_counts.csv"))

# 火山图
res_df <- as.data.frame(res)
res_df$gene <- rownames(res_df)
res_df <- res_df %>%
  mutate(sig = ifelse(pvalue < 0.05 & abs(log2FoldChange) > 0.8, "Significant", "Not Significant"))#根据需求设置,pvalue,log2FC或者FC

ggplot(res_df, aes(x = log2FoldChange, y = -log10(pvalue), color = sig)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("gray", "red")) +
  theme_minimal() +
  ggtitle("Volcano Plot: PA vs BC") +
  ggsave(filename = file.path(out_dir, "volcano_plot.pdf"), width = 6, height = 5)
