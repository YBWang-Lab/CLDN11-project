library(SCENIC)
library(Seurat)
library(cowplot)
library(dplyr)

cells.use<-subset(OL, idents = c('OL1','OL2','DAO1'))
expr <- GetAssayData(object = cells.use, assay= "RNA", slot = "data")

expr <- as(Class = 'matrix', object = expr)
save(expr,file='expr.RData')
cellInfo <- data.frame(seuratCluster=Idents(cells.use))
cellInfo <- data.frame(cellInfo)

cellTypeColumn <- "Class"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "seuratCluster"

dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

colVars <- list(seuratCluster=c("OL1"="forestgreen", 
                                "OL2"="darkorange", 
                                "DAO1"="magenta4"))
colVars$seuratCluster <- colVars$seuratCluster[intersect(names(colVars$seuratCluster), cellInfo$seuratCluster)]
saveRDS(colVars, file="int/colVars.Rds")

org="hgnc" 
dbDir="./SCENIC_REF/" 

myDatasetTitle="SCENIC 0323"    # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=70) 

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds")

exprMat <- expr
genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=3*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]
rm(exprMat)

source('runSCENIC_2_createRegulons.R')

require(data.table)
library(data.table)

runCorrelation(exprMat_filtered, scenicOptions)
exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)# Run GENIE3

library(SCENIC)
scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"] # For toy run

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions, coexMethod=c("top5perTarget")) #** Only for toy run!!
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
runSCENIC_4_aucell_binarize(scenicOptions)
tsneAUC(scenicOptions, aucType="AUC") # choose settings

regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
rss <- calcRSS(AUC=getAUC(regulonAUC), 
               cellAnnotation=cellInfo[colnames(regulonAUC), 
               ]) 
write.csv(rss,file="rss.csv")


