library(Seurat)
library(SeuratData)
library(SingleR)
library(harmony)
library(tidyverse)
library(tibble)
library(ggplot2)



setwd("D:/R Windows/RStudio program/program/FENG LAB/machine/GEO_data")

{untar("GSE106118_RAW.tar", exdir = "GSE106118_RAW")


library(Seurat)
library(data.table)

## 获取文件列表
file_list <- list.files("./GSE106118_RAW/", pattern = "\\.txt\\.gz$")


seurat_list <- list()
for (file in file_list) {
  data.path <- paste0("./GSE106118_RAW/", file)
  data <- read.delim(gzfile(data.path), header = TRUE,row.names = 1)
  data <- as.matrix(data)
  sample_name <- substr(basename(file), 1, 10)
  seurat_obj <- CreateSeuratObject(counts = data,
                                   project = sample_name,
                                   min.features = 200,
                                   min.cells = 3)
  seurat_list <- append(seurat_list, seurat_obj)
}

file_list <- file_list[1:98]

sample_names <- substr(file_list, 1, 10)
sample_names <- gsub("-", "", sample_names)
scRNA <- merge(seurat_list[[1]],
               y = seurat_list[-1],
               add.cell.ids = sample_names)


scRNA <- JoinLayers(scRNA)
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^MT-")
scRNA[["percent.rp"]] <- PercentageFeatureSet(scRNA, pattern = "^RP[SL]")

p1 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
              ncol = 3, pt.size = 0, group.by = "orig.ident")
p1



scRNA <- subset(scRNA, percent.mt < 25 &  nFeature_RNA < 6000)


Idents(scRNA) <- scRNA$orig.ident


f <- "obj.Rdata"
library(harmony)
future::plan(sequential)

if (!file.exists(f)) {
  scRNA <- scRNA %>%
    NormalizeData() %>%
    FindVariableFeatures() %>%
    ScaleData(features = rownames(.)) %>%
    RunPCA(pc.genes = VariableFeatures(.)) %>%
    RunHarmony("orig.ident") %>%
    FindNeighbors(dims = 1:15, reduction = "harmony") %>%
    FindClusters(resolution = seq(from = 0.1, to = 1.0, by = 0.1)) %>%
    RunUMAP(dims = 1:15, reduction = "harmony") %>%
    RunTSNE(dims = 1:15, reduction = "harmony")
  
  save(scRNA, file = f)
}
load(f)

p1 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
              ncol = 3, pt.size = 0, group.by = "orig.ident")
p1


#-------------------------------------------------------------------------------
#去双细胞
#ncount表达太高了，九十五万，去双细胞试一下
library(DoubletFinder)

sweep.res.list <- paramSweep(scRNA,PCs = 1:30,sct = FALSE)
sweep.stats <- summarizeSweep(sweep.res.list,GT = FALSE)
bcmvn <- find.pK(sweep.stats)
pk_best = bcmvn%>%
  arrange(desc(BCmetric)) %>% 
  pull(pK) %>% 
  .[1] %>% as.character() %>% as.numeric()

#估算双细胞群中homo doublet比例
annotation <- scRNA$seurat_clusters
homo.prop <- modelHomotypic(annotation)
print(homo.prop)

#占据6%
nexp_poi <- round(0.06*nrow(scRNA@meta.data))
nexp_poi_adj <- round(nexp_poi*(1-homo.prop))

#找
scRNA <- doubletFinder(scRNA,
                       PCs = 1:30,
                       pN = 0.25,
                       pK =  pk_best,
                       nExp = nexp_poi_adj,
                       sct = FALSE)

colnames(scRNA@meta.data)[17] <- "Double_score"
colnames(scRNA@meta.data)[18] <- "Is_Double"



DimPlot(scRNA,reduction = "tsne",group.by = "Is_Double")


VlnPlot(scRNA,group.by = "Is_Double",
        features = c("nCount_RNA","nFeature_RNA"),
          ncol = 2)


scRNA <- subset(scRNA,Is_Double=="Singlet")

p1 <- VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "percent.rp"),
              ncol = 3, pt.size = 0, group.by = "orig.ident")
p1





#-------------------------------------------------------------------------------
ElbowPlot(scRNA, ndims = 50)

library(clustree)
clustree(scRNA)

scRNA <- FindClusters(scRNA, resolution = 0.3)
scRNA$seurat_clusters <- scRNA@active.ident

plot <- DimPlot(scRNA, label = TRUE)
plot
plot <- DimPlot(scRNA, reduction = "tsne",label = TRUE)
plot


g <- "doubletFinder.RData"
save(scRNA,file = g)
load(g)



#-------------------------------------------------------------------------------
#注释
library(SingleR)


SC_marker1 <- FindAllMarkers(scRNA,only.pos = TRUE,logfc.threshold = 1,)

SC_marker_top <- SC_marker1 %>% group_by(cluster)%>%top_n(n=2,wt=avg_log2FC)


#创建marker集合
markers <- c("DDAH1","MEF2C","HAND2","NKX2-5","TBX5","BMP4","MYH6","TNNI3","TNNT2","ACTN2",
             "POSTN","VIMENTIN","WT1","TBX18","CD31","PECAM1", "CD79A", "CD79B")

markers2 <- c("TTN",
  "MYH7",
  "MYH6",
  "TNNT2",
  "VWF",
  "IFI27",
  "PECAM1",
  "DCN",
  "C7",
  "LUM",
  "FBLN1",
  "COL1A2",
  "PTPRC",
  "CD163",
  "CCL4",
  "CXCL8",
  "ACTA2",
  "CALD1",
  "MYH11")

#绘图
##marker好的注释来比对下
p <- FeaturePlot(scRNA,features = markers2,ncol = 3)
p


#MYH6,TNNI3,TNNT2,POSTN,PECAM1,DCN,MYH7,C7,COL1A2,DCN,LUM,FBLN1


new.cluster.ids <- c("MYH6","TNNI3","TNNT2","POSTN","PECAM1","DCN"
                     ,"MYH7","COL1A2","DCN","LUM","C7","FBLN1","TTN","CALD1")
names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)
DimPlot(scRNA, reduction = "tsne", label = TRUE, pt.size = 0.5) + NoLegend()

}

library(data.table)
data1<- GetAssayData(scRNA, assay = "RNA", layer = "counts")
data1 <- data1 %>% data.frame(.)

data2 <- read.csv(choose.files(), header = T, sep = ",", row.names = 1, stringsAsFactors = FALSE)

expset <- merge(data1, data2, by = "row.names", all = TRUE)
rownames(expset) <- expset$Row.names
expset$Row.names <- NULL
expset[is.na(expset)] <- 0




group_list <- c(rep("embryo",4382),rep("aldult",4933))

metadata <- as.data.frame(group_list)
colnames(metadata) <- c("class")
metadata$sample <- colnames(expset)
rownames(metadata) <- metadata$sample
expr_mat <- t(expset)


group <- "class"
metadata[[group]] <- as.factor(metadata[[group]])
library(randomForest)
set.seed(304)
rf <- randomForest(expr_mat, metadata[[group]])
rf

setwd("D:/R Windows/RStudio program/program/FENG LAB/machine/随机森林/output")

save(rf,file="rf.RData")
#结果
#Call:
#randomForest(x = expr_mat, y = metadata[[group]]) 
#Type of random forest: classification
#Number of trees: 500
#No. of variables tried at each split: 172

#OOB estimate of  error rate: 0%
#Confusion matrix:
#  aldult embryo class.error
#aldult   4933      0           0
#embryo      0   4382           0



#划分训练集和测试集进行重复验证
library(caret)
seed <- 1
set.seed(seed)
train_index <- createDataPartition(metadata[[group]], p = 0.75, list = F)
train_data <- expr_mat[train_index, ]
train_data_group <- metadata[[group]][train_index]
test_data <- expr_mat[-train_index, ]
test_data_group <- metadata[[group]][-train_index]
dim(train_data)#[1]  6987 29891
dim(test_data)#[1]  2328 29891


#boruta选择鉴定关键变量
library(Boruta)
set.seed(1)
boruta <- Boruta(x = train_data, y = train_data_group, pValue = 0.01, mcAdj = T, maxRuns = 1000)
boruta



























