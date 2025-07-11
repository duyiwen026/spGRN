setwd('/data/xukun/-duyw/xin/net/allnet/scenic/network4')
LR <- read.csv('/data/network3/LR_scenic.csv')
RT <- read.csv('/data/network3/RT_scenic.csv')
TT <- read.csv('/data/network3/TT_scenic.csv')
DEGs <- read.csv('/data/myCAF.csv')
DEGs <- subset(DEGs,p_val_adj<0.05)

#way1
LR <- subset(LR,receptor%in%c(DEGs$X))
RT <- subset(RT,tf%in%c(DEGs$X))
TT <- subset(TT,dest%in%c(DEGs$X))

RT <- subset(RT,receptor%in%c(LR$receptor))
TT <- subset(TT,src%in%c(RT$tf))

write.csv(LR,'LR_DEG.csv')
write.csv(RT,'RT_DEG.csv')
write.csv(TT,'TT_DEG.csv')


LR <- subset(LR,receptor%in%c(DEGs$X))
RT <- subset(RT,tf%in%c(DEGs$X))
TT <- subset(TT,dest%in%c(DEGs$X))

RT <- subset(RT,receptor%in%c(LR$receptor))
TT <- subset(TT,src%in%c(RT$tf))

write.csv(LR,'LR_DEG.csv')
write.csv(RT,'RT_DEG.csv')
write.csv(TT,'TT_DEG.csv')




#spatial--------------

setwd('d:/data/spatial/')
brain <- readRDS('st.rds')

allen_reference <- readRDS("sc.rds")
table(allen_reference$cell_type)
#设置ncells=3000会标准化整个数据集，
allen_reference <- SCTransform(allen_reference, ncells = 3000, verbose = FALSE) %>%
  RunPCA(verbose = FALSE) %>%
  RunUMAP(dims = 1:30)

#亚组化后，我们重新规范化皮质
cortex <- brain
cortex <- SCTransform(cortex, assay = "Spatial", verbose = FALSE) %>%
  RunPCA(verbose = FALSE)
#注释存储在对象元数据的“subclass”列中
DimPlot(allen_reference, group.by = "cell_type", label = TRUE)

anchors <- FindTransferAnchors(reference = allen_reference, query = cortex, normalization.method = "SCT")
predictions.assay <- TransferData(anchorset = anchors, refdata = allen_reference$cell_type, prediction.assay = TRUE,
                                  weight.reduction = cortex[["pca"]], dims = 1:30)
cortex[["predictions"]] <- predictions.assay

#获取每个班级的每个点的预测分数。额叶皮层区域特别感兴趣的是层流兴奋性神经元。
#在这里，我们可以区分这些神经元亚型的不同顺序层，例如：
DefaultAssay(cortex) <- "predictions"
#基于这些预测分数，我们还可以预测其位置受空间限制的细胞类型。
#我们使用基于标记点过程的相同方法来定义空间可变特征，但使用细胞类型预测分数作为"标记"而不是基因表达。
cortex <- FindSpatiallyVariableFeatures(cortex, assay = "predictions", selection.method = "markvariogram",
                                        features = rownames(cortex), r.metric = 5, slot = "data")
top.clusters <- head(SpatiallyVariableFeatures(cortex), 4)

library(spacexr)
#setwd('/data/xukun/RCTD')
# set up reference
ref <- allen_reference
ref <- UpdateSeuratObject(ref)
Idents(ref) <- "cell_type"

# extract information to pass to the RCTD Reference function
counts <- ref[["RNA"]]@counts
cluster <- as.factor(ref$cell_type)
names(cluster) <- colnames(ref)
nUMI <- ref$nCount_RNA
names(nUMI) <- colnames(ref)
reference <- Reference(counts, cluster, nUMI)

# set up query with the RCTD function SpatialRNA
#slide.seq <- SeuratData::LoadData("ssHippo")

slide.seq <- cortex
counts <- slide.seq[["Spatial"]]@counts
coords <- GetTissueCoordinates(slide.seq)
colnames(coords) <- c("x", "y")
coords[is.na(colnames(coords))] <- NULL
query <- SpatialRNA(coords, counts, colSums(counts))


RCTD <- create.RCTD(query, reference, max_cores = 8)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
slide.seq <- AddMetaData(slide.seq, metadata = RCTD@results$results_df)

p1 <- SpatialDimPlot(slide.seq, group.by = "first_type")
p2 <- SpatialDimPlot(slide.seq, group.by = "second_type")
pdf('cell_type.pdf')
p1 | p2
dev.off()

write.csv(slide.seq$second_type,'labe.csv')
x <- as.data.frame(slide.seq$second_type)
colnames(x) <- 'predicted.id'
label <- read.csv('label_transfer_bc.csv',row.names = 1)
rownames(label) <- gsub("_", "-", rownames(label))
all(rownames(x)%in%rownames(label))
label_transfer_bc <- cbind(x,label)
write.csv(label_transfer_bc,'label_transfer.csv')

#靶基因功能分析------
library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("org.Hs.eg.db") 
setwd('d:/Ryy/spatalk-myBRCA/network4/GO_kegg/')
TT <- read.csv('/data/network4/way1/TT_DEG.csv')
DEGs <- read.csv('/data/scenic/myCAF.csv')

rt <- subset(DEGs,X%in%c(TT$dest)&avg_log2FC>0)
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)    #找出基因对应的id
entrezIDs <- as.character(entrezIDs)
out=cbind(rt,entrezID=entrezIDs)
write.table(out,file="id.txt",sep="\t",quote=F,row.names=F) 
rt=read.table("id.txt",sep="\t",header=T,check.names=F)           #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                                 #去除基因id为NA的基因
gene=rt$entrezID
#GO富集分析
kk <- enrichGO(gene = gene,
               OrgDb = org.Hs.eg.db, 
               pvalueCutoff =0.05, 
               qvalueCutoff = 0.05,
               ont="all",
               readable =T)
write.table(kk,file="GO.txt",sep="\t",quote=F,row.names = F)                 #保存富集结果
#柱状图
pdf(file="GO_barplot.pdf",width = 6,height = 8)
barplot(kk, drop = TRUE, showCategory =6,split="ONTOLOGY") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
#气泡图
pdf(file="GO_bubble.pdf",width = 6,height = 8)
dotplot(kk,showCategory = 6,split="ONTOLOGY",orderBy = "GeneRatio") + facet_grid(ONTOLOGY~., scale='free')
dev.off()
#kegg富集分析
rt=read.table("id.txt",sep="\t",header=T,check.names=F)       #读取id.txt文件
rt=rt[is.na(rt[,"entrezID"])==F,]                             #去除基因id为NA的基因
gene=rt$entrezID

options(clusterProfiler.download.method = "wget")
R.utils::setOption("clusterProfiler.download.method",'auto')
kegg <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =0.05, qvalueCutoff =2)      #富集分析
write.table(kegg,file="pvl_KEGGId.txt",sep="\t",quote=F,row.names = F)                          #保存富集结果
#柱状图
pdf(file="pvl_kegg_barplot.pdf",width = 10,height = 14)
barplot(kegg, drop = TRUE, showCategory = 30)
dev.off()