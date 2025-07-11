#spatalk----------
.generate_ggi_res <- function(ggi_tf, cell_pair, receptor_name, st_data, max_hop, co_exp_ratio) {
  .co_exp <- function(x) {
    x_1 <- x[1:(length(x)/2)]
    x_2 <- x[(length(x)/2 + 1):length(x)]
    x_12 <- x_1 * x_2
    x_12_ratio <- length(x_12[x_12 > 0])/length(x_12)
    return(x_12_ratio)
  }
  .co_exp_batch <- function(st_data, ggi_res, cell_pair) {
    ggi_res_temp <- unique(ggi_res[, c("src", "dest")])
    cell_receiver <- unique(cell_pair$cell_receiver)
    m <- floor(nrow(ggi_res_temp)/5000)
    i <- 1
    res <- NULL
    while (i <= (m + 1)) {
      m_int <- 5000 * i
      if (m_int < nrow(ggi_res_temp)) {
        ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):(5000 * i), ]
      } else {
        if (m_int == nrow(ggi_res_temp)) {
          ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):(5000 * i), ]
          i <- i + 1
        } else {
          ggi_res_temp1 <- ggi_res_temp[((i - 1) * 5000 + 1):nrow(ggi_res_temp), ]
        }
      }
      ndata_src <- st_data[ggi_res_temp1$src, cell_receiver]
      ndata_dest <- st_data[ggi_res_temp1$dest, cell_receiver]
      ndata_gg <- cbind(ndata_src, ndata_dest)
      # calculate co-expression
      ggi_res_temp1$co_ratio <- NA
      ggi_res_temp1$co_ratio <- apply(ndata_gg, 1, .co_exp)
      res <- rbind(res, ggi_res_temp1)
      i <- i + 1
    }
    res$merge_key <- paste0(res$src, res$dest)
    ggi_res$merge_key <- paste0(ggi_res$src, ggi_res$dest)
    ggi_res <- merge(ggi_res, res, by = "merge_key", all.x = T, sort = F)
    ggi_res <- ggi_res[, c(2:6, 9)]
    colnames(ggi_res) <- c("src", "dest", "src_tf", "dest_tf", "hop", "co_ratio")
    return(ggi_res)
  }
  # generate ggi_res
  ggi_res <- NULL
  ggi_tf1 <- ggi_tf[ggi_tf$src == receptor_name, ]
  ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
  n <- 0
  ggi_tf1$hop <- n + 1
  while (n <= max_hop) {
    ggi_res <- rbind(ggi_res, ggi_tf1)
    ggi_tf1 <- ggi_tf[ggi_tf$src %in% ggi_tf1$dest, ]
    ggi_tf1 <- unique(ggi_tf1[ggi_tf1$dest %in% rownames(st_data), ])
    if (nrow(ggi_tf1) == 0) {
      break
    }
    ggi_tf1$hop <- n + 2
    n <- n + 1
  }
  ggi_res <- unique(ggi_res)
  # ndata_src and ndata_dest
  ggi_res_temp <- unique(ggi_res[, c("src", "dest")])
  if (nrow(ggi_res_temp) >= 5000) {
    ggi_res <- .co_exp_batch(st_data, ggi_res, cell_pair)
  } else {
    ndata_src <- st_data[ggi_res$src, cell_pair$cell_receiver]
    ndata_dest <- st_data[ggi_res$dest, cell_pair$cell_receiver]
    ndata_gg <- cbind(ndata_src, ndata_dest)
    # calculate co-expression
    ggi_res$co_ratio <- NA
    ggi_res$co_ratio <- apply(ndata_gg, 1, .co_exp)
  }
  ggi_res <- ggi_res[ggi_res$co_ratio > co_exp_ratio, ]
  return(ggi_res)
}

.generate_tf_gene_all <- function(ggi_res, max_hop) {
  tf_gene_all <- NULL
  ggi_hop <- ggi_res[ggi_res$hop == 1, ]
  for (k in 1:max_hop) {
    ggi_hop_yes <- ggi_hop[ggi_hop$dest_tf == "YES", ]
    if (nrow(ggi_hop_yes) > 0) {
      ggi_hop_tf <- ggi_res[ggi_res$hop == k + 1, ]
      if (nrow(ggi_hop_tf) > 0) {
        ggi_hop_yes <- ggi_hop_yes[ggi_hop_yes$dest %in% ggi_hop_tf$src, ]
        if (nrow(ggi_hop_yes) > 0) {
          tf_gene <- ggi_hop_yes$hop
          names(tf_gene) <- ggi_hop_yes$dest
          tf_gene_all <- c(tf_gene_all, tf_gene)
        }
      }
    }
    ggi_hop_no <- ggi_hop[ggi_hop$dest_tf == "NO", ]
    ggi_hop <- ggi_res[ggi_res$hop == k + 1, ]
    ggi_hop <- ggi_hop[ggi_hop$src %in% ggi_hop_no$dest, ]
  }
  return(tf_gene_all)
}


.generate_tf_res <- function(tf_gene_all, celltype_sender, celltype_receiver, receptor_name, ggi_res) {
  receptor_tf_temp <- data.frame(celltype_sender = celltype_sender, celltype_receiver = celltype_receiver,
                                 receptor = receptor_name, tf = names(tf_gene_all), n_hop = as.numeric(tf_gene_all), n_target = 0, stringsAsFactors = F)
  tf_names <- names(tf_gene_all)
  tf_n_hop <- as.numeric(tf_gene_all)
  for (i in 1:length(tf_names)) {
    ggi_res_tf <- ggi_res[ggi_res$src == tf_names[i] & ggi_res$hop == tf_n_hop[i] + 1, ]
    receptor_tf_temp$n_target[i] <- length(unique(ggi_res_tf$dest))
  }
  return(receptor_tf_temp)
}

.get_tf_path <- function(ggi_res, tf_gene, tf_hop, receptor) {
  tf_path <- NULL
  ggi_res1 <- ggi_res[ggi_res$dest == tf_gene & ggi_res$hop == tf_hop, ]
  if (tf_hop > 1) {
    tf_hop_new <- tf_hop - 1
    for (i in tf_hop_new:1) {
      ggi_res2 <- ggi_res[ggi_res$dest %in% ggi_res1$src & ggi_res$hop == i, ]
      ggi_res1 <- ggi_res1[ggi_res1$src %in% ggi_res2$dest, ]
      ggi_res2 <- ggi_res2[ggi_res2$dest %in% ggi_res1$src, ]
      if (i == tf_hop_new) {
        tf_path <- rbind(tf_path, ggi_res1, ggi_res2)
      } else {
        tf_path <- rbind(tf_path, ggi_res2)
      }
      ggi_res1 <- ggi_res2
    }
  } else {
    tf_path <- ggi_res1
  }
  tf_path_new <- NULL
  ggi_res1 <- tf_path[tf_path$src == receptor & tf_path$hop == 1, ]
  if (tf_hop > 1) {
    for (i in 2:tf_hop) {
      ggi_res2 <- tf_path[tf_path$src %in% ggi_res1$dest & tf_path$hop == i, ]
      ggi_res2 <- ggi_res2[ggi_res2$src %in% ggi_res1$dest, ]
      if (i == 2) {
        tf_path_new <- rbind(tf_path_new, ggi_res1, ggi_res2)
      } else {
        tf_path_new <- rbind(tf_path_new, ggi_res2)
      }
      ggi_res1 <- ggi_res2
    }
  } else {
    tf_path_new <- ggi_res1
  }
  ggi_res1 <- ggi_res[ggi_res$src == tf_gene & ggi_res$hop == (tf_hop + 1), ]
  tf_path_new <- rbind(tf_path_new, ggi_res1)
  tf_path_new$tf <- tf_gene
  return(tf_path_new)
}


save(.generate_ggi_res,.generate_tf_gene_all,.generate_tf_res,.get_tf_path,file = 'function.rda')
load('function.rda')

load('function.rda')
library(SpaTalk)
library(Seurat)

object <- readRDS('obj.rds')

LRT_path <- function(obj, celltype_sender, celltype_receiver, ligand, receptor, color = NULL, size = 5, arrow_length = 0.1) {
  pathways <- object@lr_path$pathways
  ggi_tf <- unique(pathways[, c("src", "dest", "src_tf", "dest_tf")])
  # get result from dec_celltype
  st_type <- object@para$st_type
  if_skip_dec_celltype <- object@para$if_skip_dec_celltype
  if (st_type == "single-cell") {
    st_meta <- object@meta$rawmeta
    if (if_skip_dec_celltype) {
      st_data <- object@data$rawdata
    } else {
      st_data <- object@data$rawndata
    }
  } else {
    if (if_skip_dec_celltype) {
      st_meta <- object@meta$rawmeta
      colnames(st_meta)[1] <- "cell"
      st_data <- object@data$rawdata
    } else {
      st_meta <- object@meta$newmeta
      st_data <- object@data$newdata
    }
  }
  if (!celltype_sender %in% st_meta$celltype) {
    stop("Please provide the correct name of celltype_sender!")
  }
  if (!celltype_receiver %in% st_meta$celltype) {
    stop("Please provide the correct name of celltype_receiver!")
  }
  if (!ligand %in% rownames(st_data)) {
    stop("Please provide the correct name of ligand!")
  }
  if (!receptor %in% rownames(st_data)) {
    stop("Please provide the correct name of receptor!")
  }
  max_hop <- object@para$max_hop
  cell_pair <- object@cellpair
  cell_pair <- cell_pair[[paste0(celltype_sender, " -- ", celltype_receiver)]]
  if (is.null(cell_pair)) {
    stop("No LR pairs found from the celltype_sender to celltype_receiver!")
  }
  co_exp_ratio <- object@para$co_exp_ratio
  
  
  ggi_res <- .generate_ggi_res(ggi_tf, cell_pair, receptor, st_data, max_hop, co_exp_ratio)
  tf_gene_all <- .generate_tf_gene_all(ggi_res, max_hop)
  tf_gene_all <- data.frame(gene = names(tf_gene_all), hop = tf_gene_all, stringsAsFactors = F)
  tf_gene_all_new <- unique(tf_gene_all)
  tf_gene_all <- tf_gene_all_new$hop
  names(tf_gene_all) <- tf_gene_all_new$gene
  tf_path_all <- NULL
  for (i in 1:length(tf_gene_all)) {
    tf_path <- .get_tf_path(ggi_res, names(tf_gene_all)[i], as.numeric(tf_gene_all[i]), receptor)
    tf_path_all <- rbind(tf_path_all, tf_path)
  }
  tf_path_all$hop <- tf_path_all$hop + 2
  node_x <- unique(tf_path_all[, c("dest", "hop")])
  plot_node <- data.frame(gene = c(node_x$dest, ligand, receptor), x = c(node_x$hop, 1, 2), stringsAsFactors = F)
  node_y <- as.data.frame(table(plot_node$x), stringsAsFactors = F)
  node_y_max <- max(node_y$Freq)
  plot_node$y <- 0
  plot_node_new <- NULL
  for (i in 1:max(plot_node$x)) {
    plot_node_temp <- plot_node[plot_node$x == i, ]
    if (nrow(plot_node_temp) == node_y_max) {
      plot_node_temp$y <- 1:nrow(plot_node_temp)
    } else {
      node_y_inter <- (node_y_max - 1)/(nrow(plot_node_temp) + 1)
      node_y_new <- 1 + node_y_inter
      for (j in 1:nrow(plot_node_temp)) {
        plot_node_temp$y[j] <- node_y_new
        node_y_new <- node_y_new + node_y_inter
      }
    }
    plot_node_new <- rbind(plot_node_new, plot_node_temp)
  }
  plot_node_new$Expression <- 0
  plot_node_new$Celltype <- celltype_receiver
  st_data_gene <- st_data[plot_node_new$gene, st_meta[st_meta$celltype == celltype_receiver, ]$cell]
  plot_node_new$Expression <- as.numeric(rowMeans(as.matrix(st_data_gene)))
  st_data_ligand <- st_data[ligand, st_meta[st_meta$celltype == celltype_sender, ]$cell]
  plot_node_new[plot_node_new$gene == ligand, ]$Expression <- mean(st_data_ligand)
  plot_node_new[plot_node_new$gene == ligand, ]$Celltype <- celltype_sender
  plot_node_new$Celltype <- factor(plot_node_new$Celltype, levels = c(celltype_sender, celltype_receiver))
  if (is.null(color)) {
    cellname <- unique(st_meta$celltype)
    cellname <- cellname[order(cellname)]
    if ("unsure" %in% cellname) {
      cellname <- cellname[-which(cellname == "unsure")]
    }
    col_manual <- ggpubr::get_palette(palette = "lancet", k = length(cellname))
    col_manual <- c(col_manual[which(cellname == celltype_sender)], col_manual[which(cellname == celltype_receiver)])
  } else {
    col_manual <- color
  }
  tf_path_all <- tf_path_all[, c("src", "dest", "hop")]
  colnames(tf_path_all)[3] <- "dest_x"
  tf_path_all$src_x <- tf_path_all$dest_x - 1
  tf_path_all$src_y <- 0
  tf_path_all$dest_y <- 0
  for (i in 1:nrow(tf_path_all)) {
    gene <- tf_path_all$src[i]
    gene_x <- tf_path_all$src_x[i]
    plot_node_temp <- plot_node_new[plot_node_new$gene == gene & plot_node_new$x == gene_x, ]
    tf_path_all$src_y[i] <- plot_node_temp$y
    gene <- tf_path_all$dest[i]
    gene_x <- tf_path_all$dest_x[i]
    plot_node_temp <- plot_node_new[plot_node_new$gene == gene & plot_node_new$x == gene_x, ]
    tf_path_all$dest_y[i] <- plot_node_temp$y
  }
  tf_path_all <- tf_path_all[, c("src", "src_x", "src_y", "dest", "dest_x", "dest_y")]
  ligand_xy <- plot_node_new[plot_node_new$gene == ligand, ]
  ligand_xy <- ligand_xy[order(ligand_xy$x), ]
  ligand_x <- ligand_xy$x[1]
  ligand_y <- ligand_xy$y[1]
  receptor_xy <- plot_node_new[plot_node_new$gene == receptor, ]
  receptor_xy <- receptor_xy[order(receptor_xy$x), ]
  receptor_x <- receptor_xy$x[1]
  receptor_y <- receptor_xy$y[1]
  tf_path_temp <- data.frame(src = ligand, src_x = ligand_x, src_y = ligand_y, dest = receptor, dest_x = receptor_x, dest_y = receptor_y, stringsAsFactors = F)
  tf_path_all <- rbind(tf_path_all, tf_path_temp)
  tf_target_all <<- tf_path_all
}

#ligands <- list('TGFB1','GRN')
#receptors <- list('TGFBR2', 'TNFRSF1A')
# 
Rec_TFs_Target <- list()
#
if (length(ligands) == length(receptors)) {
  for (i in seq_along(ligands)) {
    #  ligand and receptor
    current_ligand <- ligands[[i]]
    current_receptor <- receptors[[i]]
    # 
    result <- LRT_path(obj = object,
                       celltype_sender = 'Malignant_Epithelial',
                       celltype_receiver = 'cluster2',
                       ligand = current_ligand,
                       receptor = current_receptor)
    
    result_label <- paste0("result_", i)
    
    Rec_TFs_Target[[result_label]] <- result
  }
  
} else {
  print("The number of ligands and receptors should be equal.")
}



save(Rec_TFs_Target,file = 'Rec_TFs_Target.Rdata')

LRT_path_out <- function(obj, ligands, receptors, sender_celltype, receiver_celltype) {
  Rec_TFs_Target <- list()
  
  if (length(ligands) == length(receptors)) {
    for (i in seq_along(ligands)) {
      current_ligand <- ligands[[i]]
      current_receptor <- receptors[[i]]
      
      result <- LRT_path(obj = obj,
                         celltype_sender = sender_celltype,
                         celltype_receiver = receiver_celltype,
                         ligand = current_ligand,
                         receptor = current_receptor)
      
      result_label <- paste0("result_", i)
      Rec_TFs_Target[[result_label]] <- result
    }
  } else {
    print("The number of ligands and receptors should be equal.")
  }
  
  return(Rec_TFs_Target)
}


ligands <- list('TGFB1', 'GRN')
receptors <- list('TGFBR2', 'TNFRSF1A')
sender_celltype <- 'Malignant_Epithelial'
receiver_celltype <- 'CAFs_myCAF_like'

LRT_out <- LRT_path_out(object, ligands, receptors, sender_celltype, receiver_celltype)





#SpaTalk### 
# st_data: A matrix containing counts of st data
# st_meta: A data.frame containing x and y
# sc_data: A matrix containing counts of scRNA-seq data as the reference
# sc_celltype:  A character containing the cell types for scRNA-seq data
library(SpaTalk)
library(Seurat)
setwd('/data/spatalk-myBRCA')

#export OMP_NUM_THREADS=1

brain <- readRDS('/data /TumorST12_23.rds')
brain <- subset(brain,Location%in%c('Bdy'))

sc <- readRDS('/data /scRNA.rds')
# 10X mouse kidney spatial data
st_data <- brain@assays$Spatial@data
st_data <- rev_gene(data = st_data,data_type = "count",species = "Human",geneinfo = geneinfo)
st_meta <- brain@images[["image"]]@coordinates
st_meta <- st_meta[,c("tissue","imagerow","imagecol")]
colnames(st_meta) <- c("spot", "x", "y")


st_meta$spot <- rownames(st_meta)
rownames(st_meta) <- 1:nrow(st_meta)
obj <- createSpaTalk(st_data = st_data, st_meta = st_meta,species = "Human",if_st_is_sc = F,spot_max_cell = 30)
# sc_data: scRNA-seq data
# sc_celltype: cell type for each cell
sc_data <- sc@assays$RNA@counts
sc_meta <- sc@meta.data
sc_meta1 <- as.data.frame(sc_meta[,13])#celltype
rownames(sc_meta1) <- rownames(sc_meta)
colnames(sc_meta1) <- c('celltype')
sc_celltype <- as.character(sc_meta1$celltype)

#load("lrpairs.rda")
#load("pathways.rda")
load('/data/lrpairs.rda')
load('/data/pathways.rda')

#export OMP_NUM_THREADS=1

obj <- dec_celltype(object = obj,sc_data = sc_data,sc_celltype = sc_celltype)
obj <- find_lr_path(object = obj,lrpairs = lrpairs,pathways = pathways)
obj <- dec_cci_all(object = obj)
saveRDS(obj,'obj.rds')

#obj <- readRDS('data/obj.rds')

write.csv(obj@lrpair,'lrpair.csv')
write.csv(obj@lr_path$lrpairs,'lr_path_lrpairs.csv')
write.csv(obj@lr_path$pathways,'lr_path_pathways.csv')
write.csv(obj@tf,'tf.csv')

LRT_out <- LRT_path_out(object, ligands, receptors, sender_celltype, receiver_celltype)


LR <- read.csv('lrpair.csv',row.names = 1)
RT <- read.csv('tf.csv',row.names = 1)
LR <- read.csv('/data /lrpair.csv',row.names = 1)
RT <- read.csv('/data/tf.csv',row.names = 1)

LR <- subset(LR,celltype_sender%in%c('Malignant_Epithelial')&celltype_receiver%in%c('CAFs_myCAF_like'))
RT<- subset(RT,celltype_sender%in%c('Malignant_Epithelial')&celltype_receiver%in%c('CAFs_myCAF_like'))

# column_to_sort
sorted_LR <- LR[order(LR$score,decreasing = TRUE), ]

top_30_LR <- (head(sorted_LR$ligand, 30))
LR <- subset(LR,ligand%in%c(top_30_LR))
RT <- subset(RT,receptor%in%c(LR$receptor))

LRT_path_out <- function(object, ligands, receptors, sender_celltype, receiver_celltype) {
  Rec_TFs_Target <- list()
  
  if (length(ligands) == length(receptors)) {
    for (i in seq_along(ligands)) {
      current_ligand <- ligands[[i]]
      current_receptor <- receptors[[i]]
      
      result <- LRT_path(obj = obj,
                         celltype_sender = sender_celltype,
                         celltype_receiver = receiver_celltype,
                         ligand = current_ligand,
                         receptor = current_receptor)
      
      result_label <- paste0("result_", i)
      Rec_TFs_Target[[result_label]] <- result
    }
  } else {
    print("The number of ligands and receptors should be equal.")
  }
  
  return(Rec_TFs_Target)
}
LRT_out <- LRT_path_out(object, ligands, receptors, sender_celltype, receiver_celltype)

########

LR<- subset(LR,celltype_sender%in%c('Malignant_Epithelial')&celltype_receiver%in%c('CAFs_myCAF_like'))
RT<- subset(RT,celltype_sender%in%c('Malignant_Epithelial')&celltype_receiver%in%c('CAFs_myCAF_like'))

ligands <- as.list(LR$ligand)
receptors <- as.list(LR$receptor)
sender_celltype <- 'Malignant_Epithelial'
receiver_celltype <- 'CAFs_myCAF_like'

LRT_path_out <- function(object, ligands, receptors, sender_celltype, receiver_celltype) {
  Rec_TFs_Target <- list()
  
  if (length(ligands) == length(receptors)) {
    for (i in seq_along(ligands)) {
      current_ligand <- ligands[[i]]
      current_receptor <- receptors[[i]]
      
      result <- LRT_path(obj = obj,
                         celltype_sender = sender_celltype,
                         celltype_receiver = receiver_celltype,
                         ligand = current_ligand,
                         receptor = current_receptor)
      
      result_label <- paste0("result_", i)
      Rec_TFs_Target[[result_label]] <- result
    }
  } else {
    print("The number of ligands and receptors should be equal.")
  }
  
  return(Rec_TFs_Target)
}
LRT_out <- LRT_path_out(object, ligands, receptors, sender_celltype, receiver_celltype)
save(LRT_out,file = 'LRT_out.Rdata')
LRT_out <- do.call(rbind, LRT_out)
LRT_out <- subset(LRT_out,src %in%c(RT$tf))
unique(LRT_out$src)
#setwd('network/')
write.csv(LR,'LR_spatalk.csv')
