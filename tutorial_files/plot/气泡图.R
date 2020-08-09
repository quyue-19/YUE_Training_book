## 模式动物基因功能KEGG/GO注释(ClusterProfiler)

**模式动物有现成的OrgDb包，安装加载目标物种的R package就可以了。**
  
  ```R
### 模式动物(人or小鼠为例)使用clusterProfiler进行GO/KEGG分析
# Bioconductor 安装 clusterProfiler 
BioManager::install("clusterProfiler")

# 下载人or小鼠注释数据库
BioManager::install("org.Hs.eg.db")
install.packages("org.Hs.eg.db")

BioManager::install("org.Mm.eg.db")
install.packages("org.Mm.eg.db")

# org.Mm.eg.db因为网络问题下载不下来，手动下载以后直接install，实在不行直接下载以后library
install.packages("/Users/quyue/Desktop/BIOINFO/R_packages/org.Mm.eg.db_3.10.0.tar.gz", repos = NULL, type = "source")

# 加载clusterProfiler、绘图ggplot和基因注释包
library(clusterProfiler)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

# 读入基因数据
gene <- read.csv("/Users/quyue/Desktop/pericyte_LPS.csv", header = T, sep = "\t")
head(gene)

dim(gene)

# 筛选|FC|>= 2且FDR <= 0.05的为显著差异基因
gene_ensembl <- na.omit(gene[(abs(gene$log2FC) >= 2)&(gene$P.adj <= 0.05),1])

head(gene_ensembl)
length(gene_ensembl)

# 通过ENSEMBL，获取ENTREZID，GENENAME
genelist<- bitr(gene_ensembl,fromType="ENSEMBL",toType=c("ENTREZID","GENENAME"),OrgDb="org.Mm.eg.db")

head(genelist)

# 差异基因功能GO分析
go.all<- enrichGO(gene=genelist$ENTREZID, 
                  OrgDb = org.Mm.eg.db, 
                  ont='ALL', #ont='BP'
                  pAdjustMethod = 'BH',
                  pvalueCutoff = 0.05, 
                  qvalueCutoff = 0.2,
                  keyType = 'ENTREZID')
go.all

#随后对富集结果进行总览，查看BP，CC，MF的个数
dim(go.all[go.all$ONTOLOGY=='BP',]);dim(go.all[go.all$ONTOLOGY=='CC',]);dim(go.all[go.all$ONTOLOGY=='MF',])
#保存结果
write.csv(go.all@result,'/Users/quyue/Desktop/pericyte_DEG_go.all.result.csv',row.names=F)

#画图
dotplot(go.all,showCategory = 10, size = NULL,font.size = 10, title = "GO enrichment", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")
barplot(go.all,showCategory = 10, size = NULL,font.size = 10, title = "GO enrichment", split = "ONTOLOGY") + facet_grid(ONTOLOGY ~ ., scales = "free")

#kegg
enrich.kegg <- enrichKEGG(gene =genelist$ENTREZID,
                          organism ="mmu",
                          keyType = "kegg",
                          pvalueCutoff = 1,
                          pAdjustMethod = "BH",
                          minGSSize = 10,
                          maxGSSize = 500,
                          qvalueCutoff = 1,
                          use_internal_data =FALSE)
dim(enrich.kegg)

# 需要注意pvalue或p.adjust设置(filter未成功，直接跳过的)
sig.kegg<-filter(enrich.kegg,pvalue<.05)
dim(sig.kegg)

# 用未经filter的kegg画的图
kegg.bar<-barplot(enrich.kegg,showCategory=20,color = "pvalue")
kegg.dot<-dotplot(enrich.kegg,showCategory=20,color = "pvalue")

kegg.bar
kegg.dot







