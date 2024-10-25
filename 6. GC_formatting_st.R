### 
# gene characteristics
# View(seurat_obj_merge@meta.data)
# 
# # This step is to extract the Seurat objects for primary and metastatic sites separately. [Note!] If there are multiple sites, they need to be processed separately; the idents='' here should be changed to the corresponding names of your data's metastatic sites.
# seurat_obj_merge <- SetIdent(seurat_obj_merge, value = "sample_class") # The "sample_class" is the column name in the Seurat object's meta.data that records the "distinction between primary and metastatic sites."
# seurat_primary <- subset(x = seurat_obj_merge, idents = 'Primary')
# seurat_metastasis <- subset(x = seurat_obj_merge, idents = "Metastasis")

P_mtx <- seurat_primary@assays$RNA@data
M_mtx <- seurat_metastasis@assays$RNA@data

Var_top1000_gene <- VariableFeatures(seurat_obj_merge)[1:1000]
cellType_all <- unique(seurat_obj_merge@meta.data$celltype_sctype)

#################################################################################################################################################

# The table corresponding to the middle graph (a CSV file needs to be output for each metastatic site and each annotation method).

#################################################################################################################################################
### 1.The expression ratio and mean expression of the 1,000 highly variable genes in various cell types within the primary samples. #############

library(progress)  

pb <- progress_bar$new(
  format = " Progress [:bar] :percent Time Remaining: :eta",
  total = length(Var_top1000_gene), clear = FALSE, width= 60)

exp_num_celltype_p <- NULL
exp_ratio_celltype_p <- NULL
exp_mean_celltype_p <- NULL

for(i in 1:length(Var_top1000_gene)){
  pb$tick()
  num_temp <- t(as.matrix(tapply(1:ncol(P_mtx), seurat_primary$celltype_sctype, function(x){
    sum(P_mtx[Var_top1000_gene[i],x]!=0)
  })))
  
  mean_temp <- t(as.matrix(tapply(1:ncol(P_mtx), seurat_primary$celltype_sctype, function(x){
    round (mean(P_mtx[Var_top1000_gene[i],x]),3)
  })))
  
  exp_mean_celltype_p <- rbind(exp_mean_celltype_p,mean_temp)
  exp_num_celltype_p <- rbind(exp_num_celltype_p,num_temp)
  ratio_temp <- round (num_temp /sum(num_temp),3)
  exp_ratio_celltype_p <- rbind(exp_ratio_celltype_p,ratio_temp)
  
}

rownames(exp_num_celltype_p) <- Var_top1000_gene
rownames(exp_ratio_celltype_p) <- Var_top1000_gene
rownames(exp_mean_celltype_p) <- Var_top1000_gene

#################################################################################################################################################
## 2.The expression ratio and mean expression of the 1,000 highly variable genes in various cell types within the lymphatic metastatic samples.

pb2 <- progress_bar$new(
  format = " Progress [:bar] :percent Time Remaining: :eta",
  total = length(Var_top1000_gene), clear = FALSE, width= 60) 

exp_num_celltype_m <- NULL
exp_ratio_celltype_m <- NULL
exp_mean_celltype_m <- NULL

for(i in 1:length(Var_top1000_gene)){
  
  pb2$tick()
  num_temp <- t(as.matrix(tapply(1:ncol(M_mtx), seurat_metastasis$celltype_sctype, function(x){
    sum(M_mtx[Var_top1000_gene[i],x]!=0)
  })))
  
  mean_temp <- t(as.matrix(tapply(1:ncol(M_mtx), seurat_metastasis$celltype_sctype, function(x){
    round (mean(M_mtx[Var_top1000_gene[i],x]),3)
  })))
  
  exp_mean_celltype_m <- rbind(exp_mean_celltype_m,mean_temp)
  exp_num_celltype_m <- rbind(exp_num_celltype_m,num_temp)
  ratio_temp <- round (num_temp /sum(num_temp),3)
  exp_ratio_celltype_m <- rbind(exp_ratio_celltype_m,ratio_temp)
  
}

rownames(exp_num_celltype_m) <- Var_top1000_gene
rownames(exp_ratio_celltype_m) <- Var_top1000_gene
rownames(exp_mean_celltype_m) <- Var_top1000_gene

#################################################################################################################################################
### Check whether the two ratio matrices have the same column names; if one column is missing, add that column and fill it with zeros. #####

# Save the original matrix name
exp_ratio_celltype_p_name <- "exp_ratio_celltype_p"
exp_ratio_celltype_m_name <- "exp_ratio_celltype_m"

# Determine which matrix is missing a column.
if (ncol(exp_ratio_celltype_p) < ncol(exp_ratio_celltype_m)) {
  missing_matrix <- exp_ratio_celltype_p
  complete_matrix <- exp_ratio_celltype_m
  missing_matrix_name <- exp_ratio_celltype_p_name
  complete_matrix_name <- exp_ratio_celltype_m_name
  cat("missing_matrix is",missing_matrix_name,"\n")
} else if (ncol(exp_ratio_celltype_m) < ncol(exp_ratio_celltype_p)) {
  missing_matrix <- exp_ratio_celltype_m
  complete_matrix <- exp_ratio_celltype_p
  missing_matrix_name <- exp_ratio_celltype_m_name
  complete_matrix_name <- exp_ratio_celltype_p_name
  cat("missing_matrix is",missing_matrix_name,"\n")
} else {
  missing_matrix <- NULL
  complete_matrix <- NULL
  missing_matrix_name <- ""
  complete_matrix_name <- ""
  cat("Both matrices have the same number of columns","\n")
}

# If there are missing matrices, retrieve the names of the missing columns.
if (!is.null(missing_matrix)) {
  missing_col <- setdiff(colnames(complete_matrix), colnames(missing_matrix))
  cat("missing_col is",missing_col,"\n")
} else {
  missing_col <- NULL
}

# Add the missing columns and set their values to 0.
if (!is.null(missing_col)) {
  missing_col_index <- match(missing_col, colnames(complete_matrix))
  missing_matrix <- cbind(missing_matrix, matrix(0, nrow(missing_matrix), length(missing_col)))
  colnames(missing_matrix)[(ncol(missing_matrix)-length(missing_col)+1):ncol(missing_matrix)] <- missing_col
  
  # Rearrange the column names to ensure that the two matrices have consistent ordering
  complete_matrix <- complete_matrix[, colnames(missing_matrix)]
  
  if(missing_matrix_name=="exp_ratio_celltype_p"){
    assign(exp_ratio_celltype_p_name, missing_matrix)
    assign(exp_ratio_celltype_m_name, complete_matrix)
  }
  if(missing_matrix_name=="exp_ratio_celltype_m"){
    assign(exp_ratio_celltype_m_name, missing_matrix)
    assign(exp_ratio_celltype_p_name, complete_matrix)
  }
  
}

#################################################################################################################################################
## Check whether the two mean matrices have the same column names; if a column is missing, add it and fill it with zeros. #######################

exp_mean_celltype_p_name <- "exp_mean_celltype_p"
exp_mean_celltype_m_name <- "exp_mean_celltype_m"

if (ncol(exp_mean_celltype_p) < ncol(exp_mean_celltype_m)) {
  missing_matrix <- exp_mean_celltype_p
  complete_matrix <- exp_mean_celltype_m
  missing_matrix_name <- exp_mean_celltype_p_name
  complete_matrix_name <- exp_mean_celltype_m_name
  cat("missing_matrix is",missing_matrix_name,"\n")
} else if (ncol(exp_mean_celltype_m) < ncol(exp_mean_celltype_p)) {
  missing_matrix <- exp_mean_celltype_m
  complete_matrix <- exp_mean_celltype_p
  missing_matrix_name <- exp_mean_celltype_m_name
  complete_matrix_name <- exp_mean_celltype_p_name
  cat("missing_matrix is",missing_matrix_name,"\n")
} else {
  missing_matrix <- NULL
  complete_matrix <- NULL
  missing_matrix_name <- ""
  complete_matrix_name <- ""
  cat("Both matrices have the same number of columns","\n")
}

if (!is.null(missing_matrix)) {
  missing_col <- setdiff(colnames(complete_matrix), colnames(missing_matrix))
  cat("missing_col is",missing_col,"\n")
} else {
  missing_col <- NULL
}

if (!is.null(missing_col)) {
  missing_col_index <- match(missing_col, colnames(complete_matrix))
  missing_matrix <- cbind(missing_matrix, matrix(0, nrow(missing_matrix), length(missing_col)))
  colnames(missing_matrix)[(ncol(missing_matrix)-length(missing_col)+1):ncol(missing_matrix)] <- missing_col
  
  complete_matrix <- complete_matrix[, colnames(missing_matrix)]
  
  if(missing_matrix_name=="exp_mean_celltype_p"){
    assign(exp_mean_celltype_p_name, missing_matrix)
    assign(exp_mean_celltype_m_name, complete_matrix)
  }
  if(missing_matrix_name=="exp_mean_celltype_m"){
    assign(exp_mean_celltype_m_name, missing_matrix)
    assign(exp_mean_celltype_p_name, complete_matrix)
  }
  
}

####################################################### Merge four matrices ################################################################################

# Replace the zero values with 0.0001.
exp_ratio_celltype_p[exp_ratio_celltype_p == 0] <- "0.0001"
exp_ratio_celltype_m[exp_ratio_celltype_m == 0] <- "0.0001"
exp_mean_celltype_p[exp_mean_celltype_p == 0] <- "0.0001"
exp_mean_celltype_m[exp_mean_celltype_m == 0] <- "0.0001"

gc_df <- matrix(paste(exp_ratio_celltype_p, exp_ratio_celltype_m, exp_mean_celltype_p, exp_mean_celltype_m, sep = ";"), nrow = nrow(exp_ratio_celltype_p))
# View(gc_df)
rownames(gc_df) <- rownames(exp_ratio_celltype_p)
colnames(gc_df) <- colnames(exp_ratio_celltype_p)

# The output format for the first table should be: gse197778_LymphNode_df.csv [Note]
# sr represents SingleR annotation, st represents scType annotation, at represents article original annotation

write.csv(gc_df,file.path("output",paste0(GSE_number,"_",Met_site,"_st_df.csv")), row.names = T)

#################################################################################################################################################

# The tables corresponding to the graphs on both sides.

#################################################################################################################################################
## 1. The number of cells expressing the gene / total number of cells in primary or metastatic samples (calculate separately for primary and metastatic samples and then merge).

# primary
Exp_in_P <- apply(P_mtx[Var_top1000_gene,],1,function(a){sum(a[]!=0)}) 
cell_num_P <- apply(P_mtx[Var_top1000_gene,],1,length) 
ratio_exp_p <- as.matrix(apply(P_mtx[Var_top1000_gene,],1,function(a){sum(a[]!=0)/length(a)}))  # 1000x2

# metastatic
Exp_in_M <- apply(M_mtx[Var_top1000_gene,],1,function(a){sum(a[]!=0)}) 
cell_num_M <- apply(M_mtx[Var_top1000_gene,],1,length) 
ratio_exp_M <- as.matrix(apply(M_mtx[Var_top1000_gene,],1,function(a){sum(a[]!=0)/length(a)}))  # 1000x2

ratio_exp_all <- cbind(ratio_exp_p,ratio_exp_M)
colnames(ratio_exp_all) <- c("ratio_exp_p","ratio_exp_m")

### This file does not need to differentiate between cell annotation methods; run it once for each metastatic site. [Note!]
write.csv(ratio_exp_all, file.path("output",paste0(GSE_number,"_",Met_site,".csv")), row.names = T)

## 2. The proportion of each cell type in the total number of cells (calculated separately for primary and metastatic samples before merging). 

num_celltype_P <- t(as.matrix(table(seurat_primary@meta.data$celltype_sctype)))
num_celltype_M <- t(as.matrix(table(seurat_metastasis@meta.data$celltype_sctype)))

# library(magrittr)
percent_celltype_P <- t(round(as.matrix(table(seurat_primary@meta.data$celltype_sctype))/ncol(seurat_primary),3)) %>% as.data.frame()  #1x7
percent_celltype_M <- t(round(as.matrix(table(seurat_metastasis@meta.data$celltype_sctype))/ncol(seurat_metastasis),3)) %>% as.data.frame()  #1x8

# Get the union of column names
all_columns <- union(colnames(percent_celltype_P), colnames(percent_celltype_M))

# Fill in the missing columns
percent_celltype_P <- cbind(percent_celltype_P, setNames(data.frame(matrix(0, nrow(percent_celltype_P), length(all_columns) - ncol(percent_celltype_P))), setdiff(all_columns, colnames(percent_celltype_P))))
percent_celltype_M <- cbind(percent_celltype_M, setNames(data.frame(matrix(0, nrow(percent_celltype_M), length(all_columns) - ncol(percent_celltype_M))), setdiff(all_columns, colnames(percent_celltype_M))))

# Merge the two data frames by columns
percent_celltype <- t(rbind(percent_celltype_P, percent_celltype_M))
colnames(percent_celltype) <- c("percent_celltype_p","percent_celltype_m") # 8x2

write.csv(percent_celltype, file.path("output",paste0(GSE_number,"_",Met_site,"_st_2.csv")), row.names = T)

##3. The number of cells expressing genes in each cell type divided by the total number of cells of that cell type (calculated separately for primary and metastatic samples, then merged).

# primary
ratio_ownType_P <- t(apply(exp_num_celltype_p,1,function(a){
  round(a/num_celltype_P,3)
}))

# View(ratio_ownType_P)
colnames(ratio_ownType_P) <- colnames(exp_num_celltype_p)  # 1000x7

# metastatic
ratio_ownType_M <- t(apply(exp_num_celltype_m,1,function(a){
  round(a/num_celltype_M,3)
}))
# View(ratio_ownType_M)
colnames(ratio_ownType_M) <- colnames(exp_num_celltype_m)  # 1000x8

# merge
#Check if the two matrices have the same column names. If one matrix is missing a column, add that column and fill it with 0

ratio_ownType_P_name <- "ratio_ownType_P"
ratio_ownType_M_name <- "ratio_ownType_M"

if (ncol(ratio_ownType_P) < ncol(ratio_ownType_M)) {
  missing_matrix <- ratio_ownType_P
  complete_matrix <- ratio_ownType_M
  missing_matrix_name <- ratio_ownType_P_name
  complete_matrix_name <- ratio_ownType_M_name
  cat("missing_matrix is",missing_matrix_name,"\n")
} else if (ncol(ratio_ownType_M) < ncol(ratio_ownType_P)) {
  missing_matrix <- ratio_ownType_M
  complete_matrix <- ratio_ownType_P
  missing_matrix_name <- ratio_ownType_M_name
  complete_matrix_name <- ratio_ownType_P_name
  cat("missing_matrix is",missing_matrix_name,"\n")
} else {
  missing_matrix <- NULL
  complete_matrix <- NULL
  missing_matrix_name <- ""
  complete_matrix_name <- ""
  cat("Both matrices have the same number of columns","\n")
}

if (!is.null(missing_matrix)) {
  missing_col <- setdiff(colnames(complete_matrix), colnames(missing_matrix))
  cat("missing_col is",missing_col,"\n")
} else {
  missing_col <- NULL
}

if (!is.null(missing_col)) {
  missing_col_index <- match(missing_col, colnames(complete_matrix))
  missing_matrix <- cbind(missing_matrix, matrix(0, nrow(missing_matrix), length(missing_col)))
  colnames(missing_matrix)[(ncol(missing_matrix)-length(missing_col)+1):ncol(missing_matrix)] <- missing_col
  
  complete_matrix <- complete_matrix[, colnames(missing_matrix)]
  
  if(missing_matrix_name=="ratio_ownType_P"){
    assign(ratio_ownType_P_name, missing_matrix)
    assign(ratio_ownType_M_name, complete_matrix)
  }
  if(missing_matrix_name=="ratio_ownType_M"){
    assign(ratio_ownType_M_name, missing_matrix)
    assign(ratio_ownType_P_name, complete_matrix)
  }
  
}

ratio_ownType_merge <- matrix(paste(ratio_ownType_P, ratio_ownType_M, sep = ";"), nrow = nrow(ratio_ownType_P))

rownames(ratio_ownType_merge) <- rownames(ratio_ownType_P)
colnames(ratio_ownType_merge) <- colnames(ratio_ownType_P)

write.csv(ratio_ownType_merge, file.path("output",paste0(GSE_number,"_",Met_site,"_st_3.csv")), row.names = T)

####################################################################################################################################################

### This file does not need to distinguish between metastatic sites and cell annotation methods; one set of data is sufficient for a single run.
Variance_top1000_gene <- seurat_obj_merge@assays$RNA@meta.features[Var_top1000_gene,]
Variance_top1000_gene <- cbind(rownames(Variance_top1000_gene),Variance_top1000_gene[,4])
colnames(Variance_top1000_gene) <- c("gene","variance")

write.csv(Variance_top1000_gene, file.path("output",paste0(GSE_number,"_var.csv")), row.names = F)

####################################################################################################################################################
