#只用安装一次
#install.packages("BiocManager")
library(BiocManager)
#BiocManager::install(c("clusterProfiler", "org.Mm.eg.db", "enrichplot", "DOSE"))
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)
library(DOSE)
library(tidyverse)
setwd("C:/Users/86187/Desktop")#设置工作路径
deg <- read.csv("DEG_significant.csv", row.names = 1)#输入数据为Deseq2分析中的输出数据
deg_up <- deg %>% filter(log2FoldChange > 0.8 & pvalue < 0.05)#可以自己设置对应的上下调条件
deg_down <- deg %>% filter(log2FoldChange < -0.8 & pvalue < 0.05)#可以自己设置对应的上下调条件

# 提取 Ensembl ID
genes_up <- rownames(deg_up)
genes_down <- rownames(deg_down)
# 转换上调基因
genes_up_entrez <- bitr(genes_up, fromType = "ENSEMBL",
                        toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# 转换下调基因
genes_down_entrez <- bitr(genes_down, fromType = "ENSEMBL",
                          toType = "ENTREZID", OrgDb = org.Mm.eg.db)
head(genes_down_entrez)
# 上调基因 GO 富集
ego_up <- enrichGO(gene = genes_up_entrez$ENTREZID,
                   OrgDb = org.Mm.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05)

# 下调基因 GO 富集
ego_down <- enrichGO(gene = genes_down_entrez$ENTREZID,
                     OrgDb = org.Mm.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,
                     qvalueCutoff = 0.05)
# 上调基因 KEGG
ekegg_up <- enrichKEGG(gene = genes_up_entrez$ENTREZID,
                       organism = 'mmu', pvalueCutoff = 0.05)

# 下调基因 KEGG
ekegg_down <- enrichKEGG(gene = genes_down_entrez$ENTREZID,
                         organism = 'mmu', pvalueCutoff = 0.05)
# GO 富集 dotplot
dotplot(ego_up, showCategory = 10, title = "GO Upregulated") 
dotplot(ego_down, showCategory = 10, title = "GO Downregulated")

# KEGG dotplot
dotplot(ekegg_up, showCategory = 10, title = "KEGG Upregulated")
dotplot(ekegg_down, showCategory = 10, title = "KEGG Downregulated")
结果导出
write.csv(as.data.frame(ego_up), "GO_Up.csv")
write.csv(as.data.frame(ego_down), "GO_Down.csv")
write.csv(as.data.frame(ekegg_up), "KEGG_Up.csv")
write.csv(as.data.frame(ekegg_down), "KEGG_Down.csv")
