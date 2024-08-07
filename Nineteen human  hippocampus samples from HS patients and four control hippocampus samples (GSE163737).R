library(Seurat)
library(cowplot)
library(patchwork)
library(dplyr)

HS_1<-Read10X(data <- "~/HS_1/")
HS_2<-Read10X(data <- "~/HS_2/")
HS_3<-Read10X(data <- "~/HS_3/")
HS_4<-Read10X(data <- "~/HS_4/")
HS_5<-Read10X(data <- "~/HS_5/")
HS_6<-Read10X(data <- "~/HS_6/")
HS_7<-Read10X(data <- "~/HS_7/")
HS_8<-Read10X(data <- "~/HS_8/")
HS_9<-Read10X(data <- "~/HS_9/")
HS_10<-Read10X(data <- "~/HS_10/")
HS_11<-Read10X(data <- "~/HS_11/")
HS_12<-Read10X(data <- "~/HS_12/")
HS_13<-Read10X(data <- "~/HS_13/")
HS_14<-Read10X(data <- "~/HS_14/")
HS_15<-Read10X(data <- "~/HS_15/")
HS_16<-Read10X(data <- "~/HS_16/")
HS_17<-Read10X(data <- "~/HS_17/")
HS_18<-Read10X(data <- "~/HS_18/")
HS_19<-Read10X(data <- "~/HS_19/")
Ctrl_1<-Read10X(data <- "~/Ctrl_1/")
Ctrl_2<-Read10X(data <- "~/Ctrl_2/")
Ctrl_3<-Read10X(data <- "~/Ctrl_3/")
Ctrl_4<-Read10X(data <- "~/Ctrl_4/")

HS_1<-CreateSeuratObject(counts =HS_1,project = "HS_1",min.cells = 3,min.features = 200)
HS_2<-CreateSeuratObject(counts =HS_2,project = "HS_2",min.cells = 3,min.features = 200)
HS_3<-CreateSeuratObject(counts =HS_3,project = "HS_3",min.cells = 3,min.features = 200)
HS_4<-CreateSeuratObject(counts =HS_4,project = "HS_4",min.cells = 3,min.features = 200)
HS_5<-CreateSeuratObject(counts =HS_5,project = "HS_5",min.cells = 3,min.features = 200)
HS_6<-CreateSeuratObject(counts =HS_6,project = "HS_6",min.cells = 3,min.features = 200)
HS_7<-CreateSeuratObject(counts =HS_7,project = "HS_7",min.cells = 3,min.features = 200)
HS_8<-CreateSeuratObject(counts =HS_8,project = "HS_8",min.cells = 3,min.features = 200)
HS_9<-CreateSeuratObject(counts =HS_9,project = "HS_9",min.cells = 3,min.features = 200)
HS_10<-CreateSeuratObject(counts =HS_10,project = "HS_10",min.cells = 3,min.features = 200)
HS_11<-CreateSeuratObject(counts =HS_11,project = "HS_11",min.cells = 3,min.features = 200)
HS_12<-CreateSeuratObject(counts =HS_12,project = "HS_12",min.cells = 3,min.features = 200)
HS_13<-CreateSeuratObject(counts =HS_13,project = "HS_13",min.cells = 3,min.features = 200)
HS_14<-CreateSeuratObject(counts =HS_14,project = "HS_14",min.cells = 3,min.features = 200)
HS_15<-CreateSeuratObject(counts =HS_15,project = "HS_15",min.cells = 3,min.features = 200)
HS_16<-CreateSeuratObject(counts =HS_16,project = "HS_16",min.cells = 3,min.features = 200)
HS_17<-CreateSeuratObject(counts =HS_17,project = "HS_17",min.cells = 3,min.features = 200)
HS_18<-CreateSeuratObject(counts =HS_18,project = "HS_18",min.cells = 3,min.features = 200)
HS_19<-CreateSeuratObject(counts =HS_19,project = "HS_19",min.cells = 3,min.features = 200)
Ctrl_1<-CreateSeuratObject(counts =Ctrl_1,project = "Ctrl_1",min.cells = 3,min.features = 200)
Ctrl_2<-CreateSeuratObject(counts =Ctrl_2,project = "Ctrl_2",min.cells = 3,min.features = 200)
Ctrl_3<-CreateSeuratObject(counts =Ctrl_3,project = "Ctrl_3",min.cells = 3,min.features = 200)
Ctrl_4<-CreateSeuratObject(counts =Ctrl_4,project = "Ctrl_4",min.cells = 3,min.features = 200)

HS_1[["percent.mt"]]<-PercentageFeatureSet(HS_1,pattern = "MT-")
HS_2[["percent.mt"]]<-PercentageFeatureSet(HS_2,pattern = "MT-")
HS_3[["percent.mt"]]<-PercentageFeatureSet(HS_3,pattern = "MT-")
HS_4[["percent.mt"]]<-PercentageFeatureSet(HS_4,pattern = "MT-")
HS_5[["percent.mt"]]<-PercentageFeatureSet(HS_5,pattern = "MT-")
HS_6[["percent.mt"]]<-PercentageFeatureSet(HS_6,pattern = "MT-")
HS_7[["percent.mt"]]<-PercentageFeatureSet(HS_7,pattern = "MT-")
HS_8[["percent.mt"]]<-PercentageFeatureSet(HS_8,pattern = "MT-")
HS_9[["percent.mt"]]<-PercentageFeatureSet(HS_9,pattern = "MT-")
HS_10[["percent.mt"]]<-PercentageFeatureSet(HS_10,pattern = "MT-")
HS_11[["percent.mt"]]<-PercentageFeatureSet(HS_11,pattern = "MT-")
HS_12[["percent.mt"]]<-PercentageFeatureSet(HS_12,pattern = "MT-")
HS_13[["percent.mt"]]<-PercentageFeatureSet(HS_13,pattern = "MT-")
HS_14[["percent.mt"]]<-PercentageFeatureSet(HS_14,pattern = "MT-")
HS_15[["percent.mt"]]<-PercentageFeatureSet(HS_15,pattern = "MT-")
HS_16[["percent.mt"]]<-PercentageFeatureSet(HS_16,pattern = "MT-")
HS_17[["percent.mt"]]<-PercentageFeatureSet(HS_17,pattern = "MT-")
HS_18[["percent.mt"]]<-PercentageFeatureSet(HS_18,pattern = "MT-")
HS_19[["percent.mt"]]<-PercentageFeatureSet(HS_19,pattern = "MT-")
Ctrl_1[["percent.mt"]]<-PercentageFeatureSet(Ctrl_1,pattern = "MT-")
Ctrl_2[["percent.mt"]]<-PercentageFeatureSet(Ctrl_2,pattern = "MT-")
Ctrl_3[["percent.mt"]]<-PercentageFeatureSet(Ctrl_3,pattern = "MT-")
Ctrl_4[["percent.mt"]]<-PercentageFeatureSet(Ctrl_4,pattern = "MT-")

HB.genes_total<-c("HBA1","HBA2","HBB","HBD","HBE1","HBG2","HBM","HBQ1","HB2")

HB_m<-match(HB.genes_total,rownames(HS_1@assays$RNA))
HB.genes<-rownames(HS_1@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_1[["percent.HB"]]<-PercentageFeatureSet(HS_1,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_2@assays$RNA))
HB.genes<-rownames(HS_2@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_2[["percent.HB"]]<-PercentageFeatureSet(HS_2,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_3@assays$RNA))
HB.genes<-rownames(HS_3@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_3[["percent.HB"]]<-PercentageFeatureSet(HS_3,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_4@assays$RNA))
HB.genes<-rownames(HS_4@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_4[["percent.HB"]]<-PercentageFeatureSet(HS_4,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_5@assays$RNA))
HB.genes<-rownames(HS_5@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_5[["percent.HB"]]<-PercentageFeatureSet(HS_5,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_6@assays$RNA))
HB.genes<-rownames(HS_6@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_6[["percent.HB"]]<-PercentageFeatureSet(HS_6,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_7@assays$RNA))
HB.genes<-rownames(HS_7@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_7[["percent.HB"]]<-PercentageFeatureSet(HS_7,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_8@assays$RNA))
HB.genes<-rownames(HS_8@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_8[["percent.HB"]]<-PercentageFeatureSet(HS_8,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_9@assays$RNA))
HB.genes<-rownames(HS_9@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_9[["percent.HB"]]<-PercentageFeatureSet(HS_9,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_10@assays$RNA))
HB.genes<-rownames(HS_10@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_10[["percent.HB"]]<-PercentageFeatureSet(HS_10,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_11@assays$RNA))
HB.genes<-rownames(HS_11@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_11[["percent.HB"]]<-PercentageFeatureSet(HS_11,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_12@assays$RNA))
HB.genes<-rownames(HS_12@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_12[["percent.HB"]]<-PercentageFeatureSet(HS_12,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_13@assays$RNA))
HB.genes<-rownames(HS_13@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_13[["percent.HB"]]<-PercentageFeatureSet(HS_13,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_14@assays$RNA))
HB.genes<-rownames(HS_14@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_14[["percent.HB"]]<-PercentageFeatureSet(HS_14,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_15@assays$RNA))
HB.genes<-rownames(HS_15@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_15[["percent.HB"]]<-PercentageFeatureSet(HS_15,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_16@assays$RNA))
HB.genes<-rownames(HS_16@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_16[["percent.HB"]]<-PercentageFeatureSet(HS_16,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_17@assays$RNA))
HB.genes<-rownames(HS_17@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_17[["percent.HB"]]<-PercentageFeatureSet(HS_17,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_18@assays$RNA))
HB.genes<-rownames(HS_18@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_18[["percent.HB"]]<-PercentageFeatureSet(HS_18,features = HB.genes)
HB_m<-match(HB.genes_total,rownames(HS_19@assays$RNA))
HB.genes<-rownames(HS_19@assays$RNA)[HB_m]
HB.genes<-HB.genes[!is.na(HB.genes)]
HS_19[["percent.HB"]]<-PercentageFeatureSet(HS_19,features = HB.genes)
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

VlnPlot(HS_1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_5,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_6,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_7,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_8,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_9,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_10,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_11,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_12,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_13,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_14,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_15,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_16,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_17,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_18,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(HS_19,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_1,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_2,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_3,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
VlnPlot(Ctrl_4,features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.HB"))
HS_1 <- subset(HS_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 20)
HS_2 <- subset(HS_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_3 <- subset(HS_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_4 <- subset(HS_4, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_5 <- subset(HS_5, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_6 <- subset(HS_6, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_7 <- subset(HS_7, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_8 <- subset(HS_8, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_9 <- subset(HS_9, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_10 <- subset(HS_10, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_11 <- subset(HS_11, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_12 <- subset(HS_12, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_13 <- subset(HS_13, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_14 <- subset(HS_14, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_15 <- subset(HS_15, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_16 <- subset(HS_16, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_17 <- subset(HS_17, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_18 <- subset(HS_18, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
HS_19 <- subset(HS_19, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
Ctrl_1 <- subset(Ctrl_1, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
Ctrl_2 <- subset(Ctrl_2, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
Ctrl_3 <- subset(Ctrl_3, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)
Ctrl_4 <- subset(Ctrl_4, subset = nFeature_RNA > 200 & nFeature_RNA < 8000& percent.mt < 5)

HS_1<-NormalizeData(HS_1,normalization.method = "LogNormalize",scale.factor=10000)
HS_2<-NormalizeData(HS_2,normalization.method = "LogNormalize",scale.factor=10000)
HS_3<-NormalizeData(HS_3,normalization.method = "LogNormalize",scale.factor=10000)
HS_4<-NormalizeData(HS_4,normalization.method = "LogNormalize",scale.factor=10000)
HS_5<-NormalizeData(HS_5,normalization.method = "LogNormalize",scale.factor=10000)
HS_6<-NormalizeData(HS_6,normalization.method = "LogNormalize",scale.factor=10000)
HS_7<-NormalizeData(HS_7,normalization.method = "LogNormalize",scale.factor=10000)
HS_8<-NormalizeData(HS_8,normalization.method = "LogNormalize",scale.factor=10000)
HS_9<-NormalizeData(HS_9,normalization.method = "LogNormalize",scale.factor=10000)
HS_10<-NormalizeData(HS_10,normalization.method = "LogNormalize",scale.factor=10000)
HS_11<-NormalizeData(HS_11,normalization.method = "LogNormalize",scale.factor=10000)
HS_12<-NormalizeData(HS_12,normalization.method = "LogNormalize",scale.factor=10000)
HS_13<-NormalizeData(HS_13,normalization.method = "LogNormalize",scale.factor=10000)
HS_14<-NormalizeData(HS_14,normalization.method = "LogNormalize",scale.factor=10000)
HS_15<-NormalizeData(HS_15,normalization.method = "LogNormalize",scale.factor=10000)
HS_16<-NormalizeData(HS_16,normalization.method = "LogNormalize",scale.factor=10000)
HS_17<-NormalizeData(HS_17,normalization.method = "LogNormalize",scale.factor=10000)
HS_18<-NormalizeData(HS_18,normalization.method = "LogNormalize",scale.factor=10000)
HS_19<-NormalizeData(HS_19,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_1<-NormalizeData(Ctrl_1,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_2<-NormalizeData(Ctrl_2,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_3<-NormalizeData(Ctrl_3,normalization.method = "LogNormalize",scale.factor=10000)
Ctrl_4<-NormalizeData(Ctrl_4,normalization.method = "LogNormalize",scale.factor=10000)

HS_1 <- FindVariableFeatures(HS_1, selection.method = "vst", nfeatures = 2000)
HS_2 <- FindVariableFeatures(HS_2, selection.method = "vst", nfeatures = 2000)
HS_3 <- FindVariableFeatures(HS_3, selection.method = "vst", nfeatures = 2000)
HS_4 <- FindVariableFeatures(HS_4, selection.method = "vst", nfeatures = 2000)
HS_5 <- FindVariableFeatures(HS_5, selection.method = "vst", nfeatures = 2000)
HS_6 <- FindVariableFeatures(HS_6, selection.method = "vst", nfeatures = 2000)
HS_7 <- FindVariableFeatures(HS_7, selection.method = "vst", nfeatures = 2000)
HS_8 <- FindVariableFeatures(HS_8, selection.method = "vst", nfeatures = 2000)
HS_9 <- FindVariableFeatures(HS_9, selection.method = "vst", nfeatures = 2000)
HS_10 <- FindVariableFeatures(HS_10, selection.method = "vst", nfeatures = 2000)
HS_11 <- FindVariableFeatures(HS_11, selection.method = "vst", nfeatures = 2000)
HS_12 <- FindVariableFeatures(HS_12, selection.method = "vst", nfeatures = 2000)
HS_13 <- FindVariableFeatures(HS_13, selection.method = "vst", nfeatures = 2000)
HS_14 <- FindVariableFeatures(HS_14, selection.method = "vst", nfeatures = 2000)
HS_15 <- FindVariableFeatures(HS_15, selection.method = "vst", nfeatures = 2000)
HS_16 <- FindVariableFeatures(HS_16, selection.method = "vst", nfeatures = 2000)
HS_17 <- FindVariableFeatures(HS_17, selection.method = "vst", nfeatures = 2000)
HS_18 <- FindVariableFeatures(HS_18, selection.method = "vst", nfeatures = 2000)
HS_19 <- FindVariableFeatures(HS_19, selection.method = "vst", nfeatures = 2000)
Ctrl_1 <- FindVariableFeatures(Ctrl_1, selection.method = "vst", nfeatures = 2000)
Ctrl_2 <- FindVariableFeatures(Ctrl_2, selection.method = "vst", nfeatures = 2000)
Ctrl_3 <- FindVariableFeatures(Ctrl_3, selection.method = "vst", nfeatures = 2000)
Ctrl_4 <- FindVariableFeatures(Ctrl_4, selection.method = "vst", nfeatures = 2000)

plot1 <- VariableFeaturePlot(HS_1)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_2)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_3)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_4)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_6)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_7)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_8)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_9)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_10)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_11)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_12)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_13)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_14)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_15)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_16)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_17)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_18)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size=2.5)
CombinePlots(plots = list(plot1, plot2),legend="bottom")
plot1 <- VariableFeaturePlot(HS_19)
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

HS<-FindIntegrationAnchors(object.list = list(HS_1,HS_2,HS_3,HS_4,HS_5,HS_6,HS_7,HS_8,HS_9,HS_10,HS_11,HS_12,HS_13,HS_14,HS_15,HS_16,HS_17,HS_18,HS_19,Ctrl_1,Ctrl_2,Ctrl_3,Ctrl_4), dims = 1:20)
HS<-IntegrateData(anchorset=HS,dims=1:20)
DefaultAssay(HS)<-"integrated"

HS<-ScaleData(HS)

HS<-RunPCA(HS,npcs=30,verbose=FALSE) #PCA
print(HS[["pca"]],dims = 1:5,nfeatures = 5) 
VizDimLoadings(HS, dims = 1:2, reduction = "pca")
DimPlot(HS,reduction = "pca")
DimHeatmap(HS, dims = 1:6, cells = 500, balanced = TRUE)
HS <- JackStraw(HS, num.replicate = 100)
HS <- ScoreJackStraw(HS, dims = 1:20)
JackStrawPlot(HS, dims = 1:15)
ElbowPlot(HS)

HS <- FindNeighbors(HS, dims = 1:10)
HS <- FindClusters(HS, resolution = 0.15)

HS<-RunUMAP(HS,dims = 1:10)
DimPlot(HS,reduction= "umap")
DimPlot(HS,reduction= "umap",split.by='orig.ident')
HS$Group<-ifelse(HS$orig.ident=="Ctrl_1","Ctrl",ifelse(HS$orig.ident=="Ctrl_2","Ctrl",ifelse(HS$orig.ident=="Ctrl_3","Ctrl",ifelse(HS$orig.ident=="Ctrl_4","Ctrl",'HS'))))
DimPlot(HS,reduction= "umap",split.by='Group')
HS.markers <- FindAllMarkers(HS, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
VlnPlot(HS,features=c('SLC17A7','PLP1','GAD1','P2RY12','PDGFRA','VCAN','AQP4','VWF','LUM'),pt.size=0)
FeaturePlot(HS,features=c('SLC17A7','PLP1','GAD1','P2RY12','PDGFRA','VCAN','AQP4','VWF','LUM'))
# Define the major cell types: Excitatory neuron (EN); oligodendrocyte (OL); Inhibitory neuron (IN); Microglia; Oligodendrocyte precursor cell (OPC); Astrocyte; vascular cell (VC).

#Subset oligodendrocytes in the dataset
OL=subset(HS,cluster=="OL")
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
##Define OL subcluster: OL1; OL2; DAO1
DotPlot(OL,features=c('CCL3','CCL4',group.by="cluster")