library(clusterProfiler)
library(topGO)
library(Rgraphviz)
library(pathview)
library(org.Hs.eg.db)


adicar <- read.csv("/Users/quyue/Desktop/bonemarrowdata/geneID/adicar.txt",header = T, stringsAsFactors = F, sep = '\t')
head(adicar)

colnames(adicar)

# 用entrez-id
genelist <- adicar$To
class(genelist)

DEG.entrez_id = as.character(genelist)
DEG.entrez_id

enrich.go.BP = enrichGO(gene = DEG.entrez_id, 
                       OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", 
                       ont = "BP",
                       pvalueCutoff = 0.5, 
                       qvalueCutoff = 0.5) 
dotplot(enrich.go.BP)


# 用symbol
genecard <- adicar$From
class(genecard)
DEG.id = as.character(genecard)
DEG.id

enrich.go.BP = enrichGO(gene = DEG.id, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "SYMBOL", 
                        ont = "BP",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.1) 
dotplot(enrich.go.BP)
head(enrich.go.BP)

# top136 of adicar
adi <- read.csv("/Users/quyue/Desktop/adi.txt",header = T, stringsAsFactors = F, sep = '\t')
head(adi)
card <- adi$To
class(card)
adi.id = as.character(card)
adi.id

enrich.go.BP = enrichGO(gene = adi.id, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "SYMBOL", 
                        ont = "CC",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.1) 
barplot(enrich.go.BP)
head(enrich.go.BP)
plotGOgraph(enrich.go.BP)

# 所有adi-上调的gene
aadicar <- read.csv("/Users/quyue/Desktop/bonemarrowData/aadi.txt",header = T, stringsAsFactors = F, sep = '\t')
head(aadicar)

agenecard <- aadicar$From
class(agenecard)
aDEG.id = as.character(agenecard)
aDEG.id

aenrich.go.BP = enrichGO(gene = aDEG.id, 
                        OrgDb = org.Hs.eg.db, 
                        keyType = "SYMBOL", 
                        ont = "BP",
                        pvalueCutoff = 0.05, 
                        qvalueCutoff = 0.1) 
dotplot(aenrich.go.BP)
head(aenrich.go.BP)

# test
go <- enrichGO(aDEG.id, OrgDb = org.Hs.eg.db, ont='MF',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
               qvalueCutoff = 0.2,keyType = 'SYMBOL')

head(go)

dim(go[go$ONTOLOGY=='BP',])
dim(go[go$ONTOLOGY=='CC',])
dim(go[go$ONTOLOGY=='MF',])
barplot(go,showCategory=20,drop=T)
dotplot(go,showCategory=50)

