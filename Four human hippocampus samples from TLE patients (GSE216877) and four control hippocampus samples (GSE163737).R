library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)

TLE_1<-Read10X(data <- "~/TLE_1/")
TLE_2<-Read10X(data <- "~/TLE_2/")
TLE_3<-Read10X(data <- "~/TLE_3/")
TLE_4<-Read10X(data <- "~/TLE_4/")
Ctrl_1<-Read10X(data <- "~/Ctrl_1/")
Ctrl_2<-Read10X(data <- "~/Ctrl_2/")
Ctrl_3<-Read10X(data <- "~/Ctrl_3/")
Ctrl_4<-Read10X(data <- "~/Ctrl_4/")

TLE_1<-CreateSeuratObject(counts =TLE_1,project = "TLE_1",min.cells = 3,min.features = 200)
TLE_2<-CreateSeuratObject(counts =TLE_2,project = "TLE_2",min.cells = 3,min.features = 200)
TLE_3<-CreateSeuratObject(counts =TLE_3,project = "TLE_3",min.cells = 3,min.features = 200)
TLE_4<-CreateSeuratObject(counts =TLE_4,project = "TLE_4",min.cells = 3,min.features = 200)
Ctrl_1<-CreateSeuratObject(counts =Ctrl_1,project = "Ctrl_1",min.cells = 3,min.features = 200)
Ctrl_2<-CreateSeuratObject(counts =Ctrl_2,project = "Ctrl_2",min.cells = 3,min.features = 200)
Ctrl_3<-CreateSeuratObject(counts =Ctrl_3,project = "Ctrl_3",min.cells = 3,min.features = 200)
Ctrl_4<-CreateSeuratObject(counts =Ctrl_4,project = "Ctrl_4",min.cells = 3,min.features = 200)

TLE_1[["percent.mt"]]<-PercentageFeatureSet(TLE_1,pattern = "MT-")
TLE_2[["percent.mt"]]<-PercentageFeatureSet(TLE_2,pattern = "MT-")
TLE_3[["percent.mt"]]<-PercentageFeatureSet(TLE_3,pattern = "MT-")
TLE_4[["percent.mt"]]<-PercentageFeatureSet(TLE_4,pattern = "MT-")
Ctrl_1[["percent.mt"]]<-PercentageFeatureSet(Ctrl_1,pattern = "MT-")
Ctrl_2[["percent.mt"]]<-PercentageFeatureSet(Ctrl_2,pattern = "MT-")
Ctrl_3[["percent.mt"]]<-PercentageFeatureSet(Ctrl_3,pattern = "MT-")
Ctrl_4[["percent.mt"]]<-PercentageFeatureSet(Ctrl_4,pattern = "MT-")

HB.genes_total<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG2","HBM","HBQ1","HB2")

HB_m<-match(HB.genes_total,rownames(TLE_1@assays$RNA))
HB.genes<-rownames(TLE_1@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
TLE_1[["percent.HB"]]<-PercentageFeatureSet(TLE_1,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(TLE_2@assays$RNA))
HB.genes<-rownames(TLE_2@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
TLE_2[["percent.HB"]]<-PercentageFeatureSet(TLE_2,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(TLE_3@assays$RNA))
HB.genes<-rownames(TLE_3@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
TLE_3[["percent.HB"]]<-PercentageFeatureSet(TLE_3,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(TLE_4@assays$RNA))
HB.genes<-rownames(TLE_4@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
TLE_4[["percent.HB"]]<-PercentageFeatureSet(TLE_4,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(Ctrl_1@assays$RNA))
HB.genes<-rownames(Ctrl_1@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
Ctrl_1[["percent.HB"]]<-PercentageFeatureSet(Ctrl_1,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(Ctrl_2@assays$RNA))
HB.genes<-rownames(Ctrl_2@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
Ctrl_2[["percent.HB"]]<-PercentageFeatureSet(Ctrl_2,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(Ctrl_3@assays$RNA))
HB.genes<-rownames(Ctrl_3@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
Ctrl_3[["percent.HB"]]<-PercentageFeatureSet(Ctrl_3,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(Ctrl_4@assays$RNA))
HB.genes<-rownames(Ctrl_4@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
Ctrl_4[["percent.HB"]]<-PercentageFeatureSet(Ctrl_4,features = HB.genes)

VlnPlot(TLE_1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(TLE_2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(TLE_3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(TLE_4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
TLE_1 <- subset(TLE_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
TLE_2 <- subset(TLE_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
TLE_3 <- subset(TLE_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
TLE_4 <- subset(TLE_4, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
Ctrl_1 <- subset(Ctrl_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
Ctrl_2 <- subset(Ctrl_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
Ctrl_3 <- subset(Ctrl_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
Ctrl_4 <- subset(Ctrl_4, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)

TLE_1<-NormalizeData(TLE_1,normalization.method = "LogNormalize",scale.factor=10000)
TLE_2<-NormalizeData(TLE_2,normalization.method = "LogNormalize",scale.factor=10000)
TLE_3<-NormalizeData(TLE_3,normalization.method = "LogNormalize",scale.factor=10000)
TLE_4<-NormalizeData(TLE_4,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_1<-NormalizeData(Ctrl_1,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_2<-NormalizeData(Ctrl_2,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_3<-NormalizeData(Ctrl_3,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_4<-NormalizeData(Ctrl_4,normalization.method = "LogNormalize",scale.factor=10000)

TLE_1 <- FindVariableFeatures(TLE_1, selection.method = "vst", nfeatures = 2000)
TLE_2 <- FindVariableFeatures(TLE_2, selection.method = "vst", nfeatures = 2000)
TLE_3 <- FindVariableFeatures(TLE_3, selection.method = "vst", nfeatures = 2000)
TLE_4 <- FindVariableFeatures(TLE_4, selection.method = "vst", nfeatures = 2000)
Ctrl_1 <- FindVariableFeatures(Ctrl_1, selection.method = "vst", nfeatures = 2000)
Ctrl_2 <- FindVariableFeatures(Ctrl_2, selection.method = "vst", nfeatures = 2000)
Ctrl_3 <- FindVariableFeatures(Ctrl_3, selection.method = "vst", nfeatures = 2000)
Ctrl_4 <- FindVariableFeatures(Ctrl_4, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(TLE_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(TLE_2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(TLE_3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(TLE_4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(Ctrl_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(Ctrl_2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(Ctrl_3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(Ctrl_4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")

TLE<-FindIntegrationAnchors(object.list = list(TLE_1,TLE_2,TLE_3,TLE_4,Ctrl_1,Ctrl_2,Ctrl_3,Ctrl_4), dims = 1:20)
TLE<-IntegrateData(anchorset=TLE,dims=1:20)
DefaultAssay(TLE)<-"integrated"

all.genes<-rownames(TLE)
TLE<-ScaleData(TLE,features = all.genes,vars.to.regress = "percent.mt")

TLE<-RunPCA(TLE,npcs=30,verbose=FALSE) #PCA
print(TLE[["pca"]],dims = 1:5,nfeatures = 5) 
VizDimLoadings(TLE, dims = 1:2, reduction = "pca")
DimPlot(TLE,reduction = "pca")
DimHeatmap(TLE, dims = 1:6, cells = 500, balanced = TRUE)
TLE <- JackStraw(TLE, num.replicate = 100)
TLE <- ScoreJackStraw(TLE, dims = 1:20)
JackStrawPlot(TLE, dims = 1:15)
ElbowPlot(TLE)

TLE <- FindNeighbors(TLE, dims = 1:10)
TLE <- FindClusters(TLE, resolution = 0.15)

TLE<-RunUMAP(TLE,dims = 1:10)
DimPlot(TLE,reduction= "umap")
DimPlot(TLE,reduction= "umap",split.by='orig.ident')
TLE$Group<-ifelse(TLE$orig.ident=="Ctrl_1","Ctrl",ifelse(TLE$orig.ident=="Ctrl_2","Ctrl",ifelse(TLE$orig.ident=="Ctrl_3","Ctrl",ifelse(TLE$orig.ident=="Ctrl_4","Ctrl",'TLE'))))
DimPlot(TLE,reduction= "umap",split.by='Group')
TLE.markers <- FindAllMarkers(TLE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VlnPlot(TLE,features=c('PLP1','MBP'),pt.size=0)

#Subset oligodendrocytes in the dataset
OL=subset(TLE,cluster=="OL")
all.genes<-rownames(OL)
OL<-ScaleData(OL)

OL<-RunPCA(OL,npcs=30,verbose=FALSE) #PCA
print(OL[["pca"]],dims = 1:5,nfeatures = 5)
VizDimLoadings(OL, dims = 1:2, reduction = "pca")
DimPlot(OL,reduction = "pca")
DimHeatmap(OL, dims = 1:6, cells = 500, balanced = TRUE)
OL <- JackStraw(OL, num.replicate = 100)
OL <- ScoreJackStraw(OL, dims = 1:20)
JackStrawPlot(OL, dims = 1:15)
ElbowPlot(OL)

OL<- FindNeighbors(OL, dims = 1:10)
OL <- FindClusters(OL, resolution = 0.1)
OL<-RunUMAP(OL,dims = 1:10)
DimPlot(OL,reduction= "umap")
OL.markers <- FindAllMarkers(HS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#Define OL subcluster: OL3; DAO2
DotPlot(OL,features=c('PLP1','MBP','CCL3','CCL4'),group.by="cluster")