library(CellChat)
library(ggplot2)
library(ggalluvial)
library(svglite)
library(Seurat)
library(NMF)
options(stringsAsFactors = FALSE)

Ctrl=subset(HS,Group%in%c('Ctrl'))
HS=subset(HS,Group%in%c('HS'))

data.input  <- Ctrl@assays$RNA@data
identity = data.frame(group =Ctrl$cluster, row.names = names(Ctrl$cluster)) 
cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
groupSize
CellChatDB <- CellChatDB.human 
str(CellChatDB) 
showDatabaseCategory(CellChatDB)
cellchat@DB<-CellChatDB

cellchat<-subsetData(cellchat)
future::plan("multiprocess",workers=1)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.human)

cellchat<-computeCommunProb(cellchat,raw.use=FALSE)
cellchat<-filterCommunication(cellchat)
df.net<-subsetCommunication(cellchat)
write.csv(df.net,"net_lr.csv")

cellchat<-computeCommunProbPathway(cellchat)
df.netp<-subsetCommunication(cellchat,slot.name='netP')
write.csv(df.netp,'net_pathway.csv')

cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Number of interactions')
netVisual_circle(cellchat@net$weight,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Interaction weights/strength')

cellchat@netP$pathways
pathways.show.all<-cellchat@netP$pathways
levels(cellchat@idents)
dir.create("all_pathways_com_circle")
setwd("all_pathways_com_circle")
for ( i in 1:length(pathways.show.all)) {
  netVisual(cellchat,signaling=pathways.show.all[i],out.format=c('pdf'),
            vertex.receiver=vertex.reciver)
  gg<-netAnalysis_contribution(cellchat,signaling=pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width=5, height=2.5, units='in',dpi=300)
}
setwd('../')
Ctrl_cellchat=cellchat

data.input  <- HS@assays$RNA@data
identity = data.frame(group =HS$cluster, row.names = names(HS$cluster)) 
cellchat <- createCellChat(data.input)
cellchat <- addMeta(cellchat, meta = identity, meta.name = "labels")
cellchat <- setIdent(cellchat, ident.use = "labels")
groupSize <- as.numeric(table(cellchat@idents))
groupSize
CellChatDB <- CellChatDB.human 
str(CellChatDB) 
showDatabaseCategory(CellChatDB)
cellchat@DB<-CellChatDB

cellchat<-subsetData(cellchat)
future::plan("multiprocess",workers=1)
cellchat<-identifyOverExpressedGenes(cellchat)
cellchat<-identifyOverExpressedInteractions(cellchat)
cellchat<-projectData(cellchat,PPI.human)

cellchat<-computeCommunProb(cellchat,raw.use=FALSE)
cellchat<-filterCommunication(cellchat)
df.net<-subsetCommunication(cellchat)
write.csv(df.net,"net_lr.csv")

cellchat<-computeCommunProbPathway(cellchat)
df.netp<-subsetCommunication(cellchat,slot.name='netP')
write.csv(df.netp,'net_pathway.csv')

cellchat<-aggregateNet(cellchat)
groupSize<-as.numeric(table(cellchat@idents))
par(mfrow=c(1,1),xpd=TRUE)
netVisual_circle(cellchat@net$count,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Number of interactions')
netVisual_circle(cellchat@net$weight,vertex.weight=groupSize,weight.scale=T,label.edge=F,title.name='Interaction weights/strength')

cellchat@netP$pathways
pathways.show.all<-cellchat@netP$pathways
levels(cellchat@idents)
dir.create("all_pathways_com_circle")
setwd("all_pathways_com_circle")
for ( i in 1:length(pathways.show.all)) {
  netVisual(cellchat,signaling=pathways.show.all[i],out.format=c('pdf'),
            vertex.receiver=vertex.reciver)
  gg<-netAnalysis_contribution(cellchat,signaling=pathways.show.all[i])
  ggsave(filename=paste0(pathways.show.all[i], "_L-R_contribution.pdf"),
         plot=gg, width=5, height=2.5, units='in',dpi=300)
}
setwd('../')
HS_cellchat=cellchat

cco.list<-list(HS=HS_cellchat,Ctrl=Ctrl_cellchat)
cellchat<-mergeCellChat(cco.list,add.names=names(cco.list),cell.prefix=TRUE)

gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
ggsave("Overview_number_strength.pdf", p, width = 6, height = 4)

par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

par(mfrow=c(1,1))
h1<-netVisual_heatmap(cellchat)
h2<-netVisual_heatmap(cellchat,measure='weight')
h1+h2

par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F,
                   edge.weight.max = weight.max[2], edge.width.max = 12,
                   title.name = paste0("Number of interactions - ", names(cco.list)[i])) }

gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
ggsave("Compare_pathway_strengh.pdf", p, width = 10, height = 6)
