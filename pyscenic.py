setwd('/data/scenic')

write.csv(t(as.matrix(myCAF@assays$RNA@counts)),file = "myCAF.csv")
cd /data/xukun/fib/scenic
#ipython 
import os, sys 
os.getcwd()
os.listdir(os.getcwd())
import loompy as lp
import numpy as np
import scanpy as sc
x=sc.read_csv("myCAF.csv")
####https://zhuanlan.zhihu.com/p/537069794
row_attrs = {"Gene": np.array(x.var_names),}
col_attrs = {"CellID": np.array(x.obs_names)}
lp.create("sample.loom",x.X.transpose(),row_attrs,col_attrs)
exit   #退出ipython包

#
nano run_pyscenic.sh
#need：hs_hgnc_tfs.txt。
#hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
#reg.csv
#pyscenic--- GRN
pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
hs_hgnc_tfs.txt    
#PYSCENIC--- CISTARGET
pyscenic ctx \
adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname sample.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 3 \
--mask_dropouts
#PYSCENIC--- AUCELL
pyscenic aucell \
sample.loom \
reg.csv \
--output sample_SCENIC.loom \
--num_workers 3
# Ctrl + O(save)
# Ctrl + X 
#bash
chmod +x run_pyscenic.sh
./run_pyscenic.sh




#SCENIC---RT、TT
library(foreach)
library(Seurat) 
library(SCENIC)
packageVersion("SCENIC") 
library(doParallel)
library(pheatmap)
library(SCopeLoomR)
library(AUCell)
library(hdf5r)
library(tidydr)
library(dplyr)
library(data.table)
library(reshape)
library(reshape2)
setwd("D:/data/network3/")
scenicLoomPath='D:/data/sample_SCENIC.loom'

LR <- read.csv('/data/scenic/network2/LR_stlearn.csv')
RT <- read.csv('/data/scenic/network2/RT_stlearn.csv')
TT <- read.csv('/data/scenic/network2/TT_stlearn.csv')

loom <- open_loom('sample_SCENIC.loom')                 #导入sample_SCENIC.loom文件
#loom
#1
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons_incidMat[1:4,1:4]  
regulons <- regulonsToGeneLists(regulons_incidMat)      
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAUC[1:4,1:4] 
regulonAucThresholds <- get_regulon_thresholds(loom)     
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))]) 
x <- regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))]
x <- as.data.frame(x)
x$weight <- rownames(x)
top100_TF <- x %>% 
  arrange(desc(weight)) %>%  
  slice(1:100)  

xx1 <- list()
for(i in names(regulons)){
  xx1[[i]] <- data.frame(i,regulons[[i]])
}
scenic_TT <- do.call(rbind, xx1)
colnames(scenic_TT) <- c("x","y")
scenic_TT <- tidyr::unite(scenic_TT,"x_y",x,y,sep = "_",remove = F)
scenic_TT <- subset(scenic_TT,x%in%c(top100_TF$x))

#cellInfo <- read.csv("m1_meta.csv",row.names = 1)
#colnames(cellInfo)[9] <- "CellType"

RT$sign <- "(+)"
RT <- tidyr::unite(RT,"to_sign",tf,sign,sep = "",remove = F)
TT$sign <- "(+)"
TT <- tidyr::unite(TT,"from_sign",src,sign,sep = "",remove = F)
TT <- tidyr::unite(TT,"from_sign_to",from_sign,src,sep = "_",remove = F)#挑选不重复的TF

TT_scenic <- TT[TT$from_sign_to %in% scenic_TT$x_y,]

write.csv(TT_scenic,"TT_scenic.csv")
RT_scenic <- RT[RT$to_sign %in% TT_scenic$from_sign,]
write.csv(RT_scenic,"RT_scenic.csv")
