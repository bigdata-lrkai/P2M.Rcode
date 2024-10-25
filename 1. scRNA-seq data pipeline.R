###########################################################################################################################################################################################

### Set the working directory to the location where each dataset is stored
rm(list=ls())
GSE_number <- "GSE220946"
setwd(file.path("E:/Data/",GSE_number))

###########################################################################################################################################################################################

library(Seurat)   # 4.3.0.1
library(Matrix)   # 1.6-4
library(rtracklayer)   # 1.60.0
library(harmony)   # 1.2.0
library(SingleR)   # 2.2.0
library(ggplot2)   # 3.4.4
library(CellChat)   # 1.6.1
library(monocle3)   # 1.3.1
library(igraph)   # 1.6.0
library(ggpubr)   # 0.6.0
library(clusterProfiler)   # 4.8.1
library(org.Hs.eg.db)   # 3.17.0
library(GEOquery)  # 2.70.0

###########################################################################################################################################################################################

samples <- list.files(path=".",pattern=NULL)   # all sample names

### Create seurat object
seurat.obj <- sapply(samples,function(sample) {
  CreateSeuratObject(counts=Read10X(data.dir=paste0("./",sample)),project =sample,min.cells = 3,min.features = 200)})

###########################################################################################################################################################################################

### Quality control for each sample
seurat.obj <- lapply(seurat.obj,function(obj){
  obj[["percent_mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj <- subset(obj, subset = nCount_RNA > 200 & percent.mt < 20)
})

###########################################################################################################################################################################################

# Export the file format required for doublet detection in Python using DoubletDetection, with each file corresponding to an individual sample.
lapply(seurat.obj,function(obj){
  GSM_mtx<-Matrix(obj@assays$RNA@counts, sparse = T)
  writeMM(obj = GSM_mtx, file=paste0("./",levels(obj@active.ident),"/",levels(obj@active.ident),".mtx"))
  system(paste0("gzip ",levels(obj@active.ident),"/",levels(obj@active.ident),".mtx"),intern=TRUE)
})

###########################################################################################################################################################################################

# run 2.DoubletDetection.py

###########################################################################################################################################################################################

### Read the python doubletdetection results by sample and remove the identified doublets.
seurat.obj <- lapply(seurat.obj,function(obj){
  sample_name <- levels(obj@active.ident)
  sample_double_cells <- read.csv(paste0(sample_name,"/",sample_name,"_doubletdetection.csv"))
  obj<-obj[,which(sample_double_cells$X0=="False")]
}) 

###########################################################################################################################################################################################

### Extract the expression profiles of shared genes from each sample and integrate them into the Seurat object.
gene_common <- Reduce(intersect,sapply(seurat.obj,rownames))
seurat.obj <- lapply(seurat.obj,function(obj){
  obj<-obj[gene_common,]   # subset(obj, features = gene_common)
}) 
### save(seurat.obj,file="../seurat.obj.Rdata")
seurat_obj_merge <- merge(x = seurat.obj[[1]], y = seurat.obj[2:length(seurat.obj)])   # 13344 24089
View(seurat_obj_merge@meta.data)

seurat_obj_merge@meta.data$sample_type <- ifelse(seurat_obj_merge@meta.data$orig.ident=="GSM4959965","Pancreas(P)","Liver(M)")
seurat_obj_merge@meta.data$sample_class <- ifelse(seurat_obj_merge@meta.data$orig.ident=="GSM4959965","Primary","Metastasis")

###########################################################################################################################################################################################

### data standardization and dimensionality reduction
seurat_obj_merge <- NormalizeData(seurat_obj_merge, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj_merge <- FindVariableFeatures(seurat_obj_merge, selection.method = "vst", nfeatures = 2000)
seurat_obj_merge <- ScaleData(seurat_obj_merge, features = VariableFeatures(object = seurat_obj_merge))
seurat_obj_merge <- RunPCA(seurat_obj_merge, features = VariableFeatures(object = seurat_obj_merge))

### Batch effect removal
seurat_obj_merge <- RunHarmony(object=seurat_obj_merge, group.by.vars="orig.ident", plot_convergence = FALSE,lambda=1)

ElbowPlot(seurat_obj_merge, ndims = 25, reduction = "harmony")   # select the appropriate number of principal components
seurat_obj_merge <- FindNeighbors(seurat_obj_merge, dims = 1:20, reduction = "harmony")
seurat_obj_merge <- FindClusters(seurat_obj_merge, resolution = 0.8)
View(seurat_obj_merge@meta.data)
seurat_obj_merge <- RunUMAP(seurat_obj_merge, reduction = "harmony", dims = 1:20)

# Check the effectiveness of batch effect removal and manually adjust the lambda parameter.
DimPlot(seurat_obj_merge, reduction = "umap", group.by = c("orig.ident","sample_type","seurat_clusters")) 

###########################################################################################################################################################################################

### cell annotation (SingleR)
load("E:/Data/ref_Human_all.RData")
data <- GetAssayData(seurat_obj_merge, slot="data")
clusters <- seurat_obj_merge@meta.data$seurat_clusters
cell_singleR <- SingleR(test = data, ref = ref_Human_all, labels = ref_Human_all$label.main, clusters = clusters)
celltype = data.frame(ClusterID=rownames(cell_singleR), celltype=cell_singleR$labels, stringsAsFactors = FALSE)
seurat_obj_merge@meta.data$celltype_singler = celltype[match(seurat_obj_merge@meta.data$seurat_clusters,celltype$ClusterID),2]
seurat_obj_merge <- SetIdent(seurat_obj_merge, value = "celltype_singler")
View(seurat_obj_merge@meta.data)
table(seurat_obj_merge$celltype_singler)

# Review the cell annotation results and make appropriate adjustments based on the cancer context and cell markers from the Cell Marker 2.0 database.
DimPlot(seurat_obj_merge, reduction = "umap", group.by = c("celltype_singler"),label=T,repel=T)

###########################################################################################################################################################################################

### cell annotation (scType)
library(dplyr) # 1.1.4
library(HGNChelper) # 0.8.1
library(openxlsx) # 4.2.5.2

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";

tissue = "Immune system" 
gs_list = gene_sets_prepare(db_, tissue)
es.max = sctype_score(scRNAseqData = seurat_obj_merge[["RNA"]]@scale.data, scaled = TRUE, 
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative) 

cL_resutls = do.call("rbind", lapply(unique(seurat_obj_merge@meta.data$seurat_clusters), function(cl){
  es.max.cl = sort(rowSums(es.max[ ,rownames(seurat_obj_merge@meta.data[seurat_obj_merge@meta.data$seurat_clusters==cl, ])]), decreasing = !0)
  head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj_merge@meta.data$seurat_clusters==cl)), 10)
}))
sctype_scores = cL_resutls %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  

# set low-confident (low ScType score) clusters to "unknown"
sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells/4] = "Unknown"

seurat_obj_merge@meta.data$celltype_sctype = ""
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,]; 
  seurat_obj_merge@meta.data$celltype_sctype[seurat_obj_merge@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(seurat_obj_merge, reduction = "umap", label = TRUE, repel = TRUE, group.by = 'celltype_sctype')   

###########################################################################################################################################################################################

dir.create("./output")
seurat_obj_merge@meta.data <- cbind(seurat_obj_merge@meta.data,seurat_obj_merge@reductions[["umap"]]@cell.embeddings)
View(seurat_obj_merge@meta.data)
write.csv(seurat_obj_merge@meta.data, paste0("output/",GSE_number,"_dccc.csv") , row.names = T) 

# save(seurat_obj_merge,file="./seurat_obj_merge.Rdata")

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

### Identify the differentially expressed genes between metastatic and primary samples

### result of singleR
Met_site <- "Liver"

dir.create(file.path(".",GSE_number))
de_sr_dir <- file.path(GSE_number,paste0(GSE_number,"_",Met_site,"_sr"))
dir.create(de_sr_dir)

a<-unique(seurat_obj_merge$celltype_singler)

for (i in 1:length(a)) {
  
  cat("calculating: ",a[i],"\ndifferentially expressed genes\n")
  df.celltype0<-seurat_obj_merge[,seurat_obj_merge$celltype_singler%in%a[i]]
  dge.cluster <- FindMarkers(df.celltype0,ident.1 = "Metastasis",ident.2 = "Primary", group.by = 'sample_class')
  write.csv(dge.cluster, file.path(de_sr_dir,paste0(a[i],"_diffgenes.csv")), row.names = T)
  
}

#################################################################

DEG_Metastasis <- FindMarkers(seurat_obj_merge,   
                              slot = "data",  
                              ident.1 ="Metastasis", ident.2 ="Primary", group.by ='sample_class', 
                              test.use = "wilcox" )  
							  
write.csv(DEG_Metastasis,file.path(de_sr_dir,"Allcells_diffgenes.csv"), row.names = T)


###########################################################################################################################################################################################

### result of scType

de_st_dir <- file.path(GSE_number,paste0(GSE_number,"_",Met_site,"_st"))
dir.create(de_st_dir)

b<-unique(seurat_obj_merge$celltype_sctype)

for (i in 1:length(b)) {
  
  cat("calculating: ",b[i],"\ndifferentially expressed genes\n")
  df.celltype0<-seurat_obj_merge[,seurat_obj_merge$celltype_sctype%in%b[i]]
  dge.cluster <- FindMarkers(df.celltype0,ident.1 = "Metastasis",ident.2 = "Primary", group.by = 'sample_class') 
  write.csv(dge.cluster, file.path(de_st_dir,paste0(b[i],"_diffgenes.csv")), row.names = T) 
  
}

###########################################################################################################################################################################################

# run 3. DEG_formatting.R      4. FC_formatting.R

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

### Cell interaction analysis(singler)
View(seurat_obj_merge@meta.data)

### primary
seurat_obj_merge <- SetIdent(seurat_obj_merge, value = "sample_class")
seurat_primary <- subset(x = seurat_obj_merge, idents = 'Primary')

cellchat <- createCellChat(seurat_primary@assays$RNA@data, meta = seurat_primary@meta.data, group.by = "celltype_singler")  
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
# Count the number of communications between cells (how many ligand-receptor pairs) and the intensity (probability) of these interactions.
cellchat <- aggregateNet(cellchat)

interaction_mtx <- cellchat@net$weight
nonzero_rows <- rowSums(interaction_mtx != 0) > 0
nonzero_cols <- colSums(interaction_mtx != 0) > 0
interaction_mtx <- interaction_mtx[nonzero_rows,nonzero_cols]
interaction_mtx <- round(interaction_mtx,4)

interaction_edge <- data.frame(
  source = c(1:(nrow(interaction_mtx))^2),
  target = c(1:(nrow(interaction_mtx))^2),
  prob =   c(1:(nrow(interaction_mtx))^2)
)

interaction_edge$source <- rep(rownames(interaction_mtx), time = 1,  each = (nrow(interaction_mtx)))
interaction_edge$target <- rep(rownames(interaction_mtx), time = (nrow(interaction_mtx)),  each = 1)
interaction_edge$prob <- as.vector(t(interaction_mtx))

write.csv(interaction_mtx, file.path("output",paste0(GSE_number,"_",Met_site,"_sr_p.csv")), row.names = T)
write.csv(interaction_edge, file.path("output",paste0(GSE_number,"_",Met_site,"_sr_p_edge.csv")), row.names = F)

###########################################################################################################################################################################################

### metastatic
seurat_metastasis <- subset(x = seurat_obj_merge, idents = "Metastasis")

cellchat <- createCellChat(seurat_metastasis@assays$RNA@data, meta = seurat_metastasis@meta.data, group.by = "celltype_singler")  
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
# Count the number of communications between cells (how many ligand-receptor pairs) and the intensity (probability) of these interactions.
cellchat <- aggregateNet(cellchat)

interaction_mtx <- cellchat@net$weight
nonzero_rows <- rowSums(interaction_mtx != 0) > 0
nonzero_cols <- colSums(interaction_mtx != 0) > 0
interaction_mtx <- interaction_mtx[nonzero_rows,nonzero_cols]
interaction_mtx <- round(interaction_mtx,4)

interaction_edge <- data.frame(
  source = c(1:(nrow(interaction_mtx))^2),
  target = c(1:(nrow(interaction_mtx))^2),
  prob =   c(1:(nrow(interaction_mtx))^2)
)

interaction_edge$source <- rep(rownames(interaction_mtx), time = 1,  each = (nrow(interaction_mtx)))
interaction_edge$target <- rep(rownames(interaction_mtx), time = (nrow(interaction_mtx)),  each = 1)
interaction_edge$prob <- as.vector(t(interaction_mtx))

write.csv(interaction_mtx, file.path("output",paste0(GSE_number,"_",Met_site,"_sr_m.csv")), row.names = T)
write.csv(interaction_edge, file.path("output",paste0(GSE_number,"_",Met_site,"_sr_m_edge.csv")), row.names = F)

###########################################################################################################################################################################################
###########################################################################################################################################################################################

# Cell interaction analysis (sctype)

### primary

cellchat <- createCellChat(seurat_primary@assays$RNA@data, meta = seurat_primary@meta.data, group.by = "celltype_sctype")  
CellChatDB <- CellChatDB.human
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
# Count the number of communications between cells (how many ligand-receptor pairs) and the intensity (probability) of these interactions.
cellchat <- aggregateNet(cellchat)

interaction_mtx <- cellchat@net$weight
nonzero_rows <- rowSums(interaction_mtx != 0) > 0
nonzero_cols <- colSums(interaction_mtx != 0) > 0
interaction_mtx <- interaction_mtx[nonzero_rows,nonzero_cols]
interaction_mtx <- round(interaction_mtx,4)

interaction_edge <- data.frame(
  source = c(1:(nrow(interaction_mtx))^2),
  target = c(1:(nrow(interaction_mtx))^2),
  prob =   c(1:(nrow(interaction_mtx))^2)
)

interaction_edge$source <- rep(rownames(interaction_mtx), time = 1,  each = (nrow(interaction_mtx)))
interaction_edge$target <- rep(rownames(interaction_mtx), time = (nrow(interaction_mtx)),  each = 1)
interaction_edge$prob <- as.vector(t(interaction_mtx))

write.csv(interaction_mtx, file.path("output",paste0(GSE_number,"_",Met_site,"_st_p.csv")), row.names = T)
write.csv(interaction_edge, file.path("output",paste0(GSE_number,"_",Met_site,"_st_p_edge.csv")), row.names = F)

###########################################################################################################################################################################################

### metastatic

cellchat <- createCellChat(seurat_metastasis@assays$RNA@data, meta = seurat_metastasis@meta.data, group.by = "celltype_sctype")  
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multicore", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat,raw.use = FALSE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
# Count the number of communications between cells (how many ligand-receptor pairs) and the intensity (probability) of these interactions.
cellchat <- aggregateNet(cellchat)

interaction_mtx <- cellchat@net$weight
nonzero_rows <- rowSums(interaction_mtx != 0) > 0
nonzero_cols <- colSums(interaction_mtx != 0) > 0
interaction_mtx <- interaction_mtx[nonzero_rows,nonzero_cols]
interaction_mtx <- round(interaction_mtx,4)

interaction_edge <- data.frame(
  source = c(1:(nrow(interaction_mtx))^2),
  target = c(1:(nrow(interaction_mtx))^2),
  prob =   c(1:(nrow(interaction_mtx))^2)
)

interaction_edge$source <- rep(rownames(interaction_mtx), time = 1,  each = (nrow(interaction_mtx)))
interaction_edge$target <- rep(rownames(interaction_mtx), time = (nrow(interaction_mtx)),  each = 1)
interaction_edge$prob <- as.vector(t(interaction_mtx))

write.csv(interaction_mtx, file.path("output",paste0(GSE_number,"_",Met_site,"_st_m.csv")), row.names = T)
write.csv(interaction_edge, file.path("output",paste0(GSE_number,"_",Met_site,"_st_m_edge.csv")), row.names = F)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

### Calculate cell state scores.
View(seurat_obj_merge@meta.data)
CellState.markers <- read.csv("E:/Data/CellState.csv",check.names = F) %>% as.list()
seurat_obj_merge <- AddModuleScore( object = seurat_obj_merge, features = CellState.markers)

cs_mtx <- cbind(seurat_obj_merge@meta.data[,match(c("sample_class","celltype_singler","celltype_sctype"),colnames(seurat_obj_merge@meta.data))],
                seurat_obj_merge@meta.data[,c((ncol(seurat_obj_merge@meta.data)-17):ncol(seurat_obj_merge@meta.data))])

colnames(cs_mtx)[1] <- "sample_type"
colnames(cs_mtx)[4:21] <- names(CellState.markers)

write.csv(cs_mtx,file.path("output",paste0(GSE_number,"_cs_",Met_site,".csv")), row.names = T)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

### Trajectory Analysis
# 1. Merge samples from primary and all metastatic sites into a single Seurat object for trajectory analysis using Monocle3 
# 2. Use the cell with the highest Cytotrace score (indicating the lowest differentiation level) as the starting point for the trajectory.
# 3. Since the trajectory analysis includes all cell types, nearly all genes show statistically significant changes over time. Therefore, this function still utilizes genes with high variability for querying (which can be considered the intersection of differential genes and highly variable genes).
# 4. Output only two result files for each metastatic site.

### Function Section
Get_root_nodes <- function(cds,    # Monocle3 object.
                           rank0   # CytoTRACE scores should be in a data frame where the row names are the cell names, and there must be a column containing the CytoTRACE scores, labeled as "CytoTRACE."
){
  
  # View(closest_vertex)
  closest_vertex <-cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex %>% as.data.frame()
  closest_vertex$sample <- as.data.frame(colData(cds))$CellType1
  closest_vertex$cellname <- rownames(closest_vertex)
  closest_vertex <- closest_vertex[which(closest_vertex[2]=="Primary"),c(1,2)] %>% as.data.frame()
  
  # Group the primary cells that are adjacent to the same key points into a list.
  closest_vertex_cells <- split(rownames(closest_vertex), closest_vertex$V1)
  
  # Calculate the average CytoTRACE score of the cells near each key point.
  mean_scores <- lapply(closest_vertex_cells, function(cell_names) {
    mean(rank0$CytoTRACE[rownames(rank0) %in% cell_names])
  })
  # The key point with the highest CytoTRACE score is defined as the root.
  root_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names(which.max(mean_scores)))]
  root_nodes
}
##########################################

seurat_obj_merge <- SetIdent(seurat_obj_merge, value = "sample_class")
seurat_obj_pm <- subset(x = seurat_obj_merge, idents = c('Primary', 'Metastasis'))
# rm(seurat_obj_merge)

### monocle3

data <- GetAssayData(seurat_obj_pm, assay = 'RNA', slot = 'counts')
cell_metadata <- seurat_obj_pm@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

cds@int_colData$reducedDims$UMAP <- Embeddings(seurat_obj_pm, reduction = "umap")  
cds <- cluster_cells(cds)
cds <- learn_graph(cds,verbose = TRUE, use_partition = F)

root_node_select <- Get_root_nodes(cds,cell_state)
cds = order_cells(cds, root_pr_nodes = root_node_select)

Var_top1000_gene <- VariableFeatures(seurat_obj_pm)[1:1000]
genes_mtx <- seurat_obj_pm@assays$RNA@data[Var_top1000_gene,] %>% as.data.frame() %>% round(.,3) %>% t()

meta.data <- as.data.frame(colData(cds))
meta.data$pseudotime <- pseudotime(cds,reduction_method = 'UMAP')

key_point_umap <- t(as.matrix(cds@principal_graph_aux@listData[["UMAP"]][["dp_mst"]]))

##################################################

### Organizing the Output Result Format

meta0 <- cbind(meta.data[,match(c("celltype_singler","sample_class","celltype_sctype","UMAP_1","UMAP_2"),colnames(meta.data))],round(pseudotime(cds),5)) # 整合metadata
# View(meta0)
colnames(meta0) <- c("singler","site","sctype","umap_1","umap_2","pseudotime")
meta0$umap_1 <- round(meta0$umap_1,3)
meta0$umap_2 <- round(meta0$umap_2,3)

#################################################

# Combine the expression values of every 1,000 genes into a single column.
num_new_dfs <- 1
new_dfs <- list()
for (i in 1:num_new_dfs) {
  start_col <- (i - 1) * 100 + 1
  end_col <- i * 100
  new_df <- genes_mtx[, start_col:end_col]
  new_df_combined <- apply(new_df, 1, function(x) paste(x, collapse = ";"))
  new_dfs[[i]] <- data.frame(Combined = new_df_combined)
}
combined_df_list <- do.call(cbind, new_dfs)
colnames(combined_df_list) <- paste0("Column", 1:num_new_dfs)

################################################

tc_mtx <- cbind(meta0,combined_df_list)
tc_mtx_KP0 <- t(as.matrix(cds@principal_graph_aux@listData[["UMAP"]][["dp_mst"]]))

################################################

### Obtain the trajectory of the key points.

p2p <- principal_graph(cds)[["UMAP"]]
library(igraph)
# Convert the igraph graph object to a data frame.
edge_df <- get.data.frame(p2p, what = "edges")
edge_df <- edge_df[, -3]  

colnames(edge_df) <- c("Node1", "Node2")
edge_df <- edge_df[rev(rownames(edge_df)), ]
rownames(edge_df) <- 1:nrow(edge_df)

tc_mtx_KP <- cbind(key_point_umap[edge_df[,1],],key_point_umap[edge_df[,2],])
colnames(tc_mtx_KP) <- c("UMAP_1","UMAP_2","UMAP_3","UMAP_4")
tc_mtx_KP <- round(tc_mtx_KP,3)

write.csv(tc_mtx, file.path("output",paste0(GSE_number,"_TC_",Met_site,".csv")), row.names = T)
write.csv(tc_mtx_KP, file.path("output",paste0(GSE_number,"_TC_",Met_site,"_KP.csv")), row.names = T)

###########################################################################################################################################################################################
###########################################################################################################################################################################################
###########################################################################################################################################################################################

# run 5. GC_formatting_sr.R  6. GC_formatting_st.R

###########################################################################################################################################################################################

### end















