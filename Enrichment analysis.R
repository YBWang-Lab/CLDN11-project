library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(monocle)
library(org.Hs.eg.db)

dge.cluster <- FindMarkers(OL,ident.1 = "DAO1") 
sig_dge.celltype <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>0.25)
write.csv(sig_dge.celltype,file="DAO1_DEG.csv")


ego_ALL <- enrichGO(gene          = row.names(sig_dge.celltype),OrgDb = 'org.Hs.eg.db',keyType = 'SYMBOL',ont           = "ALL",pAdjustMethod = "BH",pvalueCutoff  = 0.01,qvalueCutoff  = 0.05)
ego_ALL <- data.frame(ego_ALL)
write.csv(ego_ALL,'enrichGO.csv')

genelist = bitr(row.names(sig_dge.celltype), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
genelist <- pull(genelist,ENTREZID)
ekegg <- enrichKEGG(gene = genelist, organism = 'hsa',pvalueCutoff = 0.01,qvalueCutoff = 0.01,minGSSize = 1)
write.csv(as.data.frame(ekegg@result), file="kegg.csv") 
genelist = bitr(genelist, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") 
write.csv(genelist,'KEGG_genelist.csv')

library(Seurat)
library(clusterProfiler)
library(tidyverse)
library(patchwork)
library(org.Hs.eg.db) 
library(msigdbr)
library(dplyr)
library(tibble)
library(fgsea)
library(ggplot2)
library(enrichplot)
geneset <- read.gmt("c2.cp.kegg.v7.5.1.entrez.gmt") 

dge.cluster <- FindMarkers(OL,ident.1 = "DAO1") 
markers <- subset(dge.cluster, p_val_adj<0.01&abs(avg_log2FC)>0.25)

gs <-bitr(rownames(markers), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") 
markers1<-cbind(markers[gs[,1],],gs) 
geneList = markers1$avg_log2FC 
names(geneList) = markers1$ENTREZID
geneList = sort(geneList,decreasing = T) 
egmt <- GSEA(geneList, TERM2GENE=geneset,verbose=F)
egmt1<- setReadable(egmt,OrgDb=org.Hs.eg.db, keyType = "ENTREZID")
y=data.frame(egmt1)
write.csv(y,file = "DAO1_GSEA.csv") 

for(i in seq_along(egmt@result$ID)){
  p <- gseaplot2(egmt, geneSetID = i, title = egmt@result$ID[i])
  filename <- paste0('GSEA', egmt@result$ID[i], '.png')
  ggsave(filename = filename, p, width = 4, height = 4)} 
