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

#编辑脚本
nano run_pyscenic.sh
#路径下需要：hs_hgnc_tfs.txt。
#hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather
#reg.csv
#pyscenic 的3个步骤之 GRN
pyscenic grn \
--num_workers 20 \
--output adj.sample.tsv \
--method grnboost2 \
sample.loom \
hs_hgnc_tfs.txt    #转录因子文件，1839  个基因的名字列表
#PYSCENIC 的3个步骤之 CISTARGET
pyscenic ctx \
adj.sample.tsv \
hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather \
--annotations_fname motifs-v9-nr.hgnc-m0.001-o0.0.tbl \
--expression_mtx_fname sample.loom \
--mode "dask_multiprocessing" \
--output reg.csv \
--num_workers 3 \
--mask_dropouts
#PYSCENIC 的3个步骤之 AUCELL
pyscenic aucell \
sample.loom \
reg.csv \
--output sample_SCENIC.loom \
--num_workers 3
#完成编辑后，按下 Ctrl + O(字母 O，不是数字 0)，这会提示你保存文件
#按 Ctrl + X 退出编辑器
#基于bash
chmod +x run_pyscenic.sh

#运行脚本
./run_pyscenic.sh




#利用SCENIC结果去验证RT、TT
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
#首先我们要把导入的loom处理成R中的数据
#获取regulon        regulon定义：TF与作用的genes
#1.提取每一个TF与每一个gene作用系数
regulons_incidMat <- get_regulons(loom, column.attr.name="Regulons")   
regulons_incidMat[1:4,1:4]    #在这里就可以看出 每一个TF与每一个gene的作用数值
regulons <- regulonsToGeneLists(regulons_incidMat)      #做成一个list  TF与其作用gene的list  TF+genes  个人感觉这里假如后面想分析这个TF，则这里可以画这个TF与其作用的gene的网络图                             
#2.获得regulon的AUC  即TF在每一个细胞的激活程度
regulonAUC <- get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAUC[1:4,1:4]  #regulonAUC这个文件含有每一个TF在各个细胞中的表达量  列名为细胞名   行名为TF
#3.找出在这单细胞数据中 高表达的TF
regulonAucThresholds <- get_regulon_thresholds(loom)     
tail(regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))])   #这里可以看出哪一些TF是在这个单细胞数据中高表达的
x <- regulonAucThresholds[order(as.numeric(names(regulonAucThresholds)))]
x <- as.data.frame(x)
x$weight <- rownames(x)
top100_TF <- x %>% 
  arrange(desc(weight)) %>%   # 根据指定列降序排序
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
#利用已有的靶基因筛选SCENIC中的转录因子
RT$sign <- "(+)"
RT <- tidyr::unite(RT,"to_sign",tf,sign,sep = "",remove = F)
TT$sign <- "(+)"
TT <- tidyr::unite(TT,"from_sign",src,sign,sep = "",remove = F)
TT <- tidyr::unite(TT,"from_sign_to",from_sign,src,sep = "_",remove = F)#挑选不重复的TF

TT_scenic <- TT[TT$from_sign_to %in% scenic_TT$x_y,]

write.csv(TT_scenic,"TT_scenic.csv")
RT_scenic <- RT[RT$to_sign %in% TT_scenic$from_sign,]
write.csv(RT_scenic,"RT_scenic.csv")
