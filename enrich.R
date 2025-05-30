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
# 结果导出
# 函数：将 geneID 转换为 geneSymbol 并加入新列
add_gene_symbols <- function(df) {
  entrez_all <- df$geneID %>%
    str_split("/") %>%
    unlist() %>%
    unique()
  
  # ENTREZID -> SYMBOL 注释
  mapping <- bitr(entrez_all,
                  fromType = "ENTREZID",
                  toType = "SYMBOL",
                  OrgDb = org.Mm.eg.db)
  
  convert_geneid_to_symbol <- function(geneid_string) {
    ids <- str_split(geneid_string, "/")[[1]]
    symbols <- mapping %>%
      filter(ENTREZID %in% ids) %>%
      arrange(match(ENTREZID, ids)) %>%
      pull(SYMBOL)
    paste(symbols, collapse = "/")
  }
  
  df$geneSymbol <- sapply(df$geneID, convert_geneid_to_symbol)
  return(df)
}

# 将结果转换为 data.frame
go_up_df <- as.data.frame(ego_up)
go_down_df <- as.data.frame(ego_down)
kegg_up_df <- as.data.frame(ekegg_up)
kegg_down_df <- as.data.frame(ekegg_down)

# 添加 geneSymbol 列
go_up_df <- add_gene_symbols(go_up_df)
go_down_df <- add_gene_symbols(go_down_df)
kegg_up_df <- add_gene_symbols(kegg_up_df)
kegg_down_df <- add_gene_symbols(kegg_down_df)

# 保存为 CSV
write.csv(go_up_df, "GO_Up.csv", row.names = FALSE)
write.csv(go_down_df, "GO_Down.csv", row.names = FALSE)
write.csv(kegg_up_df, "KEGG_Up.csv", row.names = FALSE)
write.csv(kegg_down_df, "KEGG_Down.csv", row.names = FALSE)

