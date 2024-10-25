###############################################################################################################################
###############################################################################################################################

rm(list=ls()) 
GSE_Number <- 'GSE161097'
options(stringsAsFactors = F)

setwd(file.path("E:/Data2/",GSE_Number))
dir.create("output")

library(GEOquery) # 2.70.0
library(stringr) # 1.5.1
library(ggplot2) # 3.4.4
library(reshape2) # 1.4.4
library(limma) # 3.58.1

###############################################################################

### Data download and grouping
gset = getGEO( GSE_Number, destdir=".", AnnotGPL = F, getGPL = F)
exp<-exprs(gset[[1]]) %>% as.data.frame()
gene.names <- read.csv("E:/Data2/GSE161097/GSE161097_family.xml/GPL29376-tbl-1.txt",header = F,row.names = 1,sep = "\t")
rownames(exp) <- gene.names[,2]

# Remove genes that are expressed as 0 in more than 30% of the samples.
exp_persent<-apply(exp,1,function(i){length(which(i==0))/length(i)})
exp <- exp[-which(exp_persent>0.3),]


# log-transform
summary(as.vector(as.matrix(exp))) 
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# -1.858   5.353   6.844   6.761   8.169  15.127 
### Conclusion: This data has already undergone a log-like normalization.
# exp <- log2(exp+1)

# extract clinical information
pdata<-pData(gset[[1]]) 

# Set the sample grouping, which needs to be modified based on each dataset!!!
group_list<-pdata$source_name_ch1
group_list <- gsub("peritoneal metastases", "metastasis", group_list)
group_list <- gsub("primary colon cancer", "primary", group_list)

group_list = factor(group_list,levels = c("primary","metastasis"))
table(group_list)

###############################################################################################################################
# Remove batch effects in bulk data while retaining differences between primary and metastatic samples.

## Check for the presence of batch effects.
### 1
exp_L=melt(exp)      
colnames(exp_L)=c('probe','sample','value')  
exp_L$group=rep(group_list,each=nrow(exp))   
head(exp_L)
p=ggplot(exp_L,aes(x=sample,y=value,fill=group))+  
  geom_boxplot()+
  ggtitle("校正前数据分布")+  
  theme(plot.title = element_text(hjust = 0.5)) +   
  #theme(axis.text.x = element_blank()  
  theme(axis.text.x=element_text(angle=90,hjust = 1,size=7)   
  ) 
print(p)

### 2
library(FactoMineR) # 2.10
library(factoextra) # 1.0.7
dat=as.data.frame(t(exp))
dat.pca <- PCA(dat, graph = FALSE)
fviz_pca_ind(dat.pca,
             geom.ind = "point",
             col.ind = group_list,
             addEllipses = TRUE,
             palette = "jco", 
             legend.title = "Patients"
)

###############################################################################
## Manually assess for any significant batch effects,correct

# library(sva)
# Sample_infor <- data.frame(sample_name = colnames(exp),
#                            Sample = pdata$geo_accession,   
#                            Group = group_list )
# 
# model <- model.matrix(~factor(Sample_infor$Group))
# exp <- ComBat(dat = exp, batch = Sample_infor$Sample,mod = model)

###############################################################################################################################
###############################################################################################################################

# differentially expressed genes analysis
design=model.matrix(~group_list)
fit=lmFit(exp,design)
fit=eBayes(fit)
deg=topTable(fit,coef=2,number = Inf)
head(deg)

DEG_up <- subset(deg, P.Value<0.05&logFC>0.5&adj.P.Val < 0.2)
DEG_down <- subset(deg, P.Value<0.05&logFC<(-0.5)&adj.P.Val < 0.2)

# DEG_up ########################################################################
sorted_indices <- order(DEG_up$logFC,decreasing = T)
DEG_up_sort <- DEG_up[sorted_indices, ]

if(nrow(DEG_up_sort)>100){ 
  DEG_up_top100 <- DEG_up_sort[1:100,c("logFC","adj.P.Val")]
}else{ 
  DEG_up_top100 <- DEG_up_sort[,c("logFC","adj.P.Val")]
}

DEG_up_top100$logFC <- round(DEG_up_top100$logFC,2)
DEG_up_top100 <- cbind(rownames(DEG_up_top100),DEG_up_top100$logFC)
colnames(DEG_up_top100) <- c("Gene","logFC")

# DEG_down ####################################################################
sorted_indices <- order(DEG_down$logFC)
DEG_down_sort <- DEG_down[sorted_indices, ]

if(nrow(DEG_down_sort)>100){ 
  DEG_down_top100 <- DEG_down_sort[1:100,c("logFC","adj.P.Val")]
}else{ 
  DEG_down_top100 <- DEG_down_sort[,c("logFC","adj.P.Val")]
}

DEG_down_top100$logFC <- round(DEG_down_top100$logFC,2)
DEG_down_top100 <- cbind(rownames(DEG_down_top100),DEG_down_top100$logFC)
colnames(DEG_down_top100) <- c("Gene","logFC")


## Merge the two matrices. ##########################################################################
DEG_up_top100_merge <- apply(DEG_up_top100,2,function(x){paste(x,collapse = ";")}) %>% as.list()
DEG_up_top100_merge_mtx <- do.call(cbind,DEG_up_top100_merge)
DEG_down_top100_merge <- apply(DEG_down_top100,2,function(x){paste(x,collapse = ";")}) %>% as.list()
DEG_down_top100_merge_mtx <- do.call(cbind,DEG_down_top100_merge)
UPandDOWN <- rbind(DEG_up_top100_merge_mtx,DEG_down_top100_merge_mtx)
rownames(UPandDOWN) <- c("UP_Metastasis","UP_Primary")

write.csv(UPandDOWN, file.path("output",paste0("b_de_",GSE_Number,"_peritoneum.csv")), row.names = T)

## The second table. ##############################################################################
# X-axis: log2 fold change   Y-axis: -log10 p-value

DEG_volcano <- deg[,c(1,4)]
DEG_volcano[,2] <- round(-log10(DEG_volcano[,2]),2)
colnames(DEG_volcano) <- c("log2FC(x)","-log10Pval(y)")
write.csv(DEG_volcano, file.path("output",paste0("b_de2_",GSE_Number,"_peritoneum.csv")), row.names = T)

###############################################################################################################################
###############################################################################################################################

library(Matrix) # 1.6-4
library(clusterProfiler) # 4.8.1
library(org.Hs.eg.db) # 3.17.0

# Perform enrichment analysis using clusterProfiler.
DEG_up <- subset(deg, P.Value<0.05&logFC>0.5&adj.P.Val < 0.2)
DEG_down <- subset(deg, P.Value<0.05&logFC<(-0.5)&adj.P.Val < 0.2)
DEG_all <- subset(deg, P.Value<0.05&adj.P.Val < 0.2)

sorted_indices <- order(DEG_up$logFC,decreasing = T)
DEG_up_sort <- DEG_up[sorted_indices, ]
if(nrow(DEG_up_sort)>50){
  DEG_up <- DEG_up_sort[1:50,]}

sorted_indices <- order(DEG_down$logFC)
DEG_down_sort <- DEG_down[sorted_indices, ]
if(nrow(DEG_down_sort)>50){
  DEG_down <- DEG_down_sort[1:50,]}

sorted_indices <- order(DEG_all$P.Value,decreasing = F)
DEG_all_sort <- DEG_all[sorted_indices, ]
if(nrow(DEG_all_sort)>50){
  DEG_all <- DEG_all_sort[1:50,]}

##############################################################

DEG_up <- mapIds(org.Hs.eg.db,rownames(DEG_up),"ENTREZID","SYMBOL")
DEG_down <- mapIds(org.Hs.eg.db,rownames(DEG_down),"ENTREZID","SYMBOL")
DEG_all <- mapIds(org.Hs.eg.db,rownames(DEG_all),"ENTREZID","SYMBOL")

DEG_up <- na.omit(DEG_up)
DEG_down <- na.omit(DEG_down)
DEG_all <- na.omit(DEG_all)

###############################################################
# GO

go_up0 <- enrichGO(gene = DEG_up,
                   OrgDb = org.Hs.eg.db,
                   keyType = "ENTREZID",
                   ont = "ALL",
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                   readable = T
)
go_down0 <- enrichGO(gene = DEG_down,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                     readable = T
)
go_all0 <- enrichGO(gene = DEG_all,
                     OrgDb = org.Hs.eg.db,
                     keyType = "ENTREZID",
                     ont = "ALL",
                     pAdjustMethod = "BH",
                     pvalueCutoff = 0.05,qvalueCutoff = 0.2,
                     readable = T
)

if (!is.null(go_up0)) {
  
  go_up <- subset(go_up0@result, ONTOLOGY=="BP"&pvalue<0.05&p.adjust<0.2)
  if(nrow(go_up)<=0){
    go_up <- subset(go_up0@result,pvalue<0.05)
  }
  
  if(nrow(go_up)>0){
    sorted_indices <- order(go_up$pvalue)
    go_up_sort <- go_up[sorted_indices, ]
    
    if(nrow(go_up_sort)>10){ 
      go_up_top10 <- go_up_sort[1:10,c("Description","pvalue")]
    }else{
      go_up_top10 <- go_up_sort[,c("Description","pvalue")]
    }
    
    go_up_top10$pvalue <- round(-log10(go_up_top10$pvalue),2)
    colnames(go_up_top10) <- c("Gene Ontology","P-value(-log10)")
    
    if(all(is.infinite(go_up_top10$`P-value(-log10)`))) {
      go_up_top10$`P-value(-log10)` <- 10
    } else {
      max_value <- max(go_up_top10$`P-value(-log10)`, na.rm = TRUE)
      go_up_top10$`P-value(-log10)`[is.infinite(go_up_top10$`P-value(-log10)`)] <- max_value + 1
    }
  }else{
    go_up_top10 <- NULL
  }
}else{
  go_up_top10 <- NULL
}


##############################################################

if (!is.null(go_down0)) {
  
  go_down <- subset(go_down0@result, ONTOLOGY=="BP"&pvalue<0.05&p.adjust<0.2)
  if(nrow(go_down)<=0){
    go_down <- subset(go_down0@result,pvalue<0.05)
  }
  
  if(nrow(go_down)>0){
    sorted_indices <- order(go_down$pvalue)
    go_down_sort <- go_down[sorted_indices, ]
    
    if(nrow(go_down_sort)>10){
      go_down_top10 <- go_down_sort[1:10,c("Description","pvalue")]
    }else{
      go_down_top10 <- go_down_sort[,c("Description","pvalue")]
    }
    
    go_down_top10$pvalue <- round(-log10(go_down_top10$pvalue),2)
    colnames(go_down_top10) <- c("Gene Ontology","P-value(-log10)")
    
    if(all(is.infinite(go_down_top10$`P-value(-log10)`))) {
      go_down_top10$`P-value(-log10)` <- 10
    } else {
      max_value <- max(go_down_top10$`P-value(-log10)`, na.rm = TRUE)
      go_down_top10$`P-value(-log10)`[is.infinite(go_down_top10$`P-value(-log10)`)] <- max_value + 1
    }
  }else{
    go_down_top10 <- NULL
  }
}else{
  go_down_top10 <- NULL
}

##############################################################

if (!is.null(go_all0)) {
  
  go_all <- subset(go_all0@result, ONTOLOGY=="BP"&pvalue<0.05)
  if(nrow(go_all)<=0){
    go_all <- subset(go_all0@result,pvalue<0.05) 
  }
  
  if(nrow(go_all)>0){
    sorted_indices <- order(go_all$pvalue)
    go_all_sort <- go_all[sorted_indices, ]
    if(nrow(go_all_sort)>10){ 
      go_all_top10 <- go_all_sort[1:10,c("Description","pvalue")]
    }else{ 
      go_all_top10 <- go_all_sort[,c("Description","pvalue")]
    }
    
    go_all_top10$pvalue <- round(-log10(go_all_top10$pvalue),2)
    colnames(go_all_top10) <- c("Gene Ontology","P-value(-log10)")
    
    if(all(is.infinite(go_all_top10$`P-value(-log10)`))) {
      go_all_top10$`P-value(-log10)` <- 10
    } else {
      max_value <- max(go_all_top10$`P-value(-log10)`, na.rm = TRUE)
      go_all_top10$`P-value(-log10)`[is.infinite(go_all_top10$`P-value(-log10)`)] <- max_value + 1
    }
  }else{
    go_all_top10 <- NULL
  }
}else{
  go_all_top10 <- NULL
}

#############################################################
# KEGG

kegg_up <- enrichKEGG(DEG_up, 
                      organism = 'hsa',  
                      keyType = 'kegg', 
                      pvalueCutoff = 0.05,
                      pAdjustMethod = 'BH',
                      qvalueCutoff = 0.2)

kegg_down <- enrichKEGG(DEG_down, 
                        organism = 'hsa',
                        keyType = 'kegg', 
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        qvalueCutoff = 0.2)

kegg_all <- enrichKEGG(DEG_all, 
                        organism = 'hsa',
                        keyType = 'kegg', 
                        pvalueCutoff = 0.05,
                        pAdjustMethod = 'BH',
                        qvalueCutoff = 0.2)



if (!is.null(kegg_up)) {
  kegg_up <- subset(kegg_up@result,pvalue<0.05&p.adjust<0.2)
  
  if (nrow(kegg_up)>0){
    sorted_indices <- order(kegg_up$pvalue)
    kegg_up_sort <- kegg_up[sorted_indices, ]
    
    if(nrow(kegg_up_sort)>10){
      kegg_up_top10 <- kegg_up_sort[1:10,c("Description","pvalue")]
    }else{ 
      kegg_up_top10 <- kegg_up_sort[,c("Description","pvalue")]
    }
    
    kegg_up_top10$pvalue <- round(-log10(kegg_up_top10$pvalue),2)
    colnames(kegg_up_top10) <- c("Gene Ontology","P-value(-log10)")
    
    if(all(is.infinite(kegg_up_top10$`P-value(-log10)`))) {
      kegg_up_top10$`P-value(-log10)` <- 10
    } else {
      max_value <- max(kegg_up_top10$`P-value(-log10)`, na.rm = TRUE)
      kegg_up_top10$`P-value(-log10)`[is.infinite(kegg_up_top10$`P-value(-log10)`)] <- max_value + 1
    }
  }else{
    kegg_up_top10 <- NULL
  }
  
}else{
  kegg_up_top10 <- NULL
}

#############################################################

if (!is.null(kegg_down)) {
  kegg_down <- subset(kegg_down@result,pvalue<0.05&p.adjust<0.2)
  
  if (nrow(kegg_down)>0){
    sorted_indices <- order(kegg_down$pvalue)
    kegg_down_sort <- kegg_down[sorted_indices, ]
    
    if(nrow(kegg_down_sort)>10){ 
      kegg_down_top10 <- kegg_down_sort[1:10,c("Description","pvalue")]
    }else{ 
      kegg_down_top10 <- kegg_down_sort[,c("Description","pvalue")]
    }
    
    kegg_down_top10$pvalue <- round(-log10(kegg_down_top10$pvalue),2)
    colnames(kegg_down_top10) <- c("Gene Ontology","P-value(-log10)")
    
    if(all(is.infinite(kegg_down_top10$`P-value(-log10)`))) {
      kegg_down_top10$`P-value(-log10)` <- 10
    } else {
      max_value <- max(kegg_down_top10$`P-value(-log10)`, na.rm = TRUE)
      kegg_down_top10$`P-value(-log10)`[is.infinite(kegg_down_top10$`P-value(-log10)`)] <- max_value + 1
    }
  }else{
    kegg_down_top10 <- NULL
  }
  
}else{
  kegg_down_top10 <- NULL
}

##############################################################

if (!is.null(kegg_all)) {
  kegg_all <- subset(kegg_all@result,pvalue<0.05)
  
  if (nrow(kegg_all)>0){
    sorted_indices <- order(kegg_all$pvalue)
    kegg_all_sort <- kegg_all[sorted_indices, ]
    
    if(nrow(kegg_all_sort)>10){
      kegg_all_top10 <- kegg_all_sort[1:10,c("Description","pvalue")]
    }else{
      kegg_all_top10 <- kegg_all_sort[,c("Description","pvalue")]
    }
    
    kegg_all_top10$pvalue <- round(-log10(kegg_all_top10$pvalue),2)
    colnames(kegg_all_top10) <- c("Gene Ontology","P-value(-log10)")
    
    if(all(is.infinite(kegg_all_top10$`P-value(-log10)`))) {
      kegg_all_top10$`P-value(-log10)` <- 10
    } else {
      max_value <- max(kegg_all_top10$`P-value(-log10)`, na.rm = TRUE)
      kegg_all_top10$`P-value(-log10)`[is.infinite(kegg_all_top10$`P-value(-log10)`)] <- max_value + 1
    }
  }else{
    kegg_all_top10 <- NULL
  }
  
}else{
  kegg_all_top10 <- NULL
}

#############################################################
# Combine the six tables together.

go_up_top10_cbind <- cbind(paste(go_up_top10[,1], collapse = ";"),paste(go_up_top10[,2], collapse = ";"))
go_down_top10_cbind<- cbind(paste(go_down_top10[,1], collapse = ";"),paste(go_down_top10[,2], collapse = ";"))
go_all_top10_cbind<- cbind(paste(go_all_top10[,1], collapse = ";"),paste(go_all_top10[,2], collapse = ";"))
go_all <- do.call(rbind, list(go_up_top10_cbind, go_down_top10_cbind, go_all_top10_cbind))

kegg_up_top10_cbind<- cbind(paste(kegg_up_top10[,1], collapse = ";"),paste(kegg_up_top10[,2], collapse = ";"))
kegg_down_top10_cbind<- cbind(paste(kegg_down_top10[,1], collapse = ";"),paste(kegg_down_top10[,2], collapse = ";"))
kegg_all_top10_cbind<- cbind(paste(kegg_all_top10[,1], collapse = ";"),paste(kegg_all_top10[,2], collapse = ";"))
kegg_all <- do.call(rbind, list(kegg_up_top10_cbind, kegg_down_top10_cbind, kegg_all_top10_cbind))

GOandKEGG <- rbind(go_all,kegg_all) %>% as.data.frame()
colnames(GOandKEGG) <- c("Functions or Pathways","P-value(-log10)")
rownames(GOandKEGG) <- c("GO_Metastasis","GO_Primary","GO_All","KEGG_Metastasis","KEGG_Primary","KEGG_All")

GOandKEGG[which(GOandKEGG$`Functions or Pathways`== ""),1]<- "no significant results"
GOandKEGG[which(GOandKEGG$`Functions or Pathways`== "no significant results"),2]<- 0

write.csv(GOandKEGG, file.path("output",paste0("b_fc_",GSE_Number,"_peritoneum.csv")), row.names = T)

###############################################################################################################################
###############################################################################################################################

# https://www.jianshu.com/p/5bfb94b0e250  https://cloud.tencent.com/developer/article/1670951

### Protein-protein interaction analysis.

library(STRINGdb) # 2.14.3
library(tidyverse) # 2.0.0
load("E:/Gene_protein.Rdata")
deg <- deg[intersect(rownames(deg),Gene_protein),]

sorted_indices <- order(DEG_up$logFC,decreasing = T)
DEG_up_sort <- DEG_up[sorted_indices, ]
if(nrow(DEG_up_sort)>50){
  DEG_up <- DEG_up_sort[1:50,]}

sorted_indices <- order(DEG_down$logFC)
DEG_down_sort <- DEG_down[sorted_indices, ]
if(nrow(DEG_down_sort)>50){
  DEG_down <- DEG_down_sort[1:50,]}

DEG_all <- rbind(DEG_up,DEG_down)

#############################################################

string_db <- STRINGdb$new( version="11.0b", species=9606, 
                           score_threshold=700, input_directory="")

gene <- rownames(DEG_all)
gene %<>% bitr(fromType = "SYMBOL", 
               toType = "ENTREZID", 
               OrgDb = "org.Hs.eg.db", 
               drop = T)
			   
data_mapped <- gene %>% string_db$map(my_data_frame_id_col_names = "ENTREZID", 
                                      removeUnmappedRows = TRUE)

hit<-data_mapped$STRING_id
info <- string_db$get_interactions(hit)

links <- info %>%
  mutate(from = data_mapped[match(from, data_mapped$STRING_id), "SYMBOL"]) %>% 
  mutate(to = data_mapped[match(to, data_mapped$STRING_id), "SYMBOL"]) %>%  
  dplyr::select(from, to , last_col()) %>% 
  dplyr::rename(weight = combined_score) 

duplicated_indices <- duplicated(links)
links <- links[!duplicated_indices, ]

nodes <- links %>% { data.frame(gene = c(.$from, .$to)) } %>% distinct()

#############################################################

f_type <- c()
t_type <- c()
for(i in 1:nrow(links)){
  f_type[i] <- ifelse(links[i,1]%in%rownames(DEG_up),"metastasis","primary")
  t_type[i] <- ifelse(links[i,2]%in%rownames(DEG_up),"metastasis","primary")
}

table(c(f_type,t_type))

f_exp <- apply(exp[match(links$from,rownames(exp)),],1,mean) %>% round(.,2)
t_exp <- apply(exp[match(links$to,rownames(exp)),],1,mean) %>% round(.,2)

pi_mtx <- do.call(cbind,list(links$from,f_type,f_exp,links$to,t_type,t_exp)) %>% as.data.frame()
colnames(pi_mtx) <- c("from","f_type","f_exp","to","t_type","t_exp")


for(i in 1:nrow(pi_mtx)){
  pi_mtx$exp_cor_p[i] <- round(cor(as.numeric(as.character(exp[pi_mtx$from[i],which(group_list=="primary")])), as.numeric(as.character(exp[pi_mtx$to[i],which(group_list=="primary")])), method="pearson"),3)
  pi_mtx$exp_cor_m[i] <- round(cor(as.numeric(as.character(exp[pi_mtx$from[i],which(group_list=="metastasis")])), as.numeric(as.character(exp[pi_mtx$to[i],which(group_list=="metastasis")])), method="pearson"),3)
}

write.csv(pi_mtx,file.path("output",paste0("b_pi_",GSE_Number,"_peritoneum.csv")), row.names = F)

#############################################################

exp_met <- exp[nodes[,1],which(group_list=="metastasis")]
metgenes_cor <- cor(t(exp_met), method="pearson") %>% round(.,2)

exp_pri <- exp[nodes[,1],which(group_list=="primary")]
prigenes_cor <- cor(t(exp_pri), method="pearson") %>% round(.,2)

write.csv(metgenes_cor,file.path("output",paste0("b_pim_",GSE_Number,"_peritoneum.csv")), row.names = T)
write.csv(prigenes_cor,file.path("output",paste0("b_pip_",GSE_Number,"_peritoneum.csv")), row.names = T)

###############################################################################################################################
###############################################################################################################################
# Immune infiltration analysis

write.table(exp, file="DATA.txt", sep="\t", row.names=T, col.names=T, quote=F)
source('Cibersort.R')
cibersort_result <- CIBERSORT('LM22.txt','DATA.txt', perm = 1000, QN = T)
cib_result <- cibersort_result %<>% round(.,4)

cib_proportion <- tapply(1:nrow(cib_result), group_list, function(x){
  apply(cib_result[,1:(ncol(cib_result)-3)],2,function(y){
    round(mean(y[x]),4)
  })
})

cib_proportion_mtx <- do.call(cbind,cib_proportion) %>% as.data.frame()

write.csv(cib_proportion_mtx,file.path("output",paste0("b_ii_",GSE_Number,"_peritoneum.csv")), row.names = T)

###############################################################################################################################
###############################################################################################################################
### Survival analysis.

# https://zhuanlan.zhihu.com/p/645916220
library(survival) # 3.5-7

cox_matrix <- cbind(pdata[,c("characteristics_ch1.1","characteristics_ch1.6")],t(exp)) 
colnames(cox_matrix)[1:2] <- c("os_status","os_time")

for(i in 1:nrow(cox_matrix)){
  if(unlist(strsplit(cox_matrix[i,1],split=":"))[2]==" Dead"){
    cox_matrix[i,1] <- "1"
  }else{
    cox_matrix[i,1] <- "0"
  }
  cox_matrix[i,2] <- unlist(strsplit(cox_matrix[i,2],split=":"))[2] %>% gsub("\\s+", "", .)
}

for(i in 1:ncol(cox_matrix)){
  cox_matrix[,i] <- as.numeric(cox_matrix[,i])
}

colnames(cox_matrix) <- gsub('-', '_', colnames(cox_matrix)) 
cox_matrix_primary <- cox_matrix[which(group_list=="primary"),]

pfilter <- 0.05

cox_result_primary <- data.frame()  
for(i in colnames(cox_matrix_primary[,3:ncol(cox_matrix_primary)])){   
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ cox_matrix_primary[,i], data = cox_matrix_primary)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  if(pvalue<pfilter){ 
    cox_result_primary <- rbind(cox_result_primary,
                                cbind(gene=i,
                                      HR=round(unisum$coefficients[,2],2),
                                      L95CI=round(unisum$conf.int[,3],2),
                                      H95CI=round(unisum$conf.int[,4],2),
                                      pvalue=round(unisum$coefficients[,5],2)
                                ))
  }
}   

cox_matrix_metastasis <- cox_matrix[which(group_list=="metastasis"),]
cox_result_metastasis <- data.frame() 

for(i in colnames(cox_matrix_primary[,3:ncol(cox_matrix_primary)])){   
  unicox <- coxph(Surv(time = os_time, event = os_status) ~ cox_matrix_primary[,i], data = cox_matrix_primary)  
  unisum<- summary(unicox)   
  pvalue <- round(unisum$coefficients[,5],3) 
  if(pvalue<pfilter){ 
    cox_result_primary <- rbind(cox_result_primary,
                                cbind(gene=i,
                                      HR=round(unisum$coefficients[,2],2),
                                      L95CI=round(unisum$conf.int[,3],2),
                                      H95CI=round(unisum$conf.int[,4],2),
                                      pvalue=round(unisum$coefficients[,5],2)
                                ))
  }
}  

gene_cox <- union(cox_result_primary$gene,cox_result_metastasis$gene)
mtx_survival_primary <- cbind(cox_matrix_primary[,c(2,1)],round(cox_matrix_primary[,gene_cox],2))
colnames(mtx_survival_primary)[1:2] <- c("time","event")

mtx_survival_metastasis <- cbind(cox_matrix_metastasis[,c(2,1)],round(cox_matrix_metastasis[,gene_cox],2))
colnames(mtx_survival_metastasis)[1:2] <- c("time","event")

write.table(mtx_survival_primary,file.path("output",paste0("b_sv_",GSE_Number,"_peritoneum_p.txt")), sep = "\t",quote = FALSE,row.names = F)
write.table(mtx_survival_metastasis,file.path("output",paste0("b_sv_",GSE_Number,"_peritoneum_m.txt")),sep = "\t",quote = FALSE,row.names = F)

# example:NOS2
library(survminer) # 0.4.9

cox_matrix$group <- ifelse(cox_matrix$NOS2>=median(cox_matrix$NOS2), 'high', 'low')
cox_matrix$group <- as.factor(cox_matrix$group)

fit <- survfit(Surv(os_time,os_status) ~ group,data = cox_matrix)
P <- ggsurvplot(fit, data = cox_matrix, pval = TRUE, 
           surv.median.line = "hv",
           conf.int = TRUE,
           risk.table = TRUE,
           ncensor.plot = TRUE,
           add.all = FALSE,
           palette = "hue", 
           xlab = "Time", 
           xlim = c(0, max(cox_matrix$os_time)),
           break.time.by = 20)
		   
###############################################################################################################################
###############################################################################################################################

save(exp,pdata,group_list,file = paste0(GSE_Number,".Rdata"))
































