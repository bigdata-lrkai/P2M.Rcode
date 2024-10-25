
########## Processing of TCGA data. #########

# Cancer data with metastatic samples.
Cancer_M <- c("BLCA","BRCA","CESC","CHOL","COAD","ESCA","HNSC","KIHC","KIRC","KIRP","LIHC","LUAD","LUSC","MESO","PAAD","READ","SKCM","STAD","TGCT","THCA","UVM")
TCGA_dir <- "I:/TCGA"
TCGA_output_dir <- "I:/TCGA_output"
setwd(TCGA_dir)

sapply(file.path(TCGA_output_dir,Cancer_M,"output"), dir.create)
sapply(file.path(TCGA_output_dir,Cancer_M,"output"), function(path){
  unlink(path, recursive = TRUE)
})

sapply(file.path(TCGA_output_dir,Cancer_M), dir.create)

for(j in 1:length(Cancer_M)){
  
  ###############################################################################################################################
  ###############################################################################################################################
  
  cat("Processing: ",Cancer_M[j],"...\n")
  load(file.path(TCGA_dir,paste0(Cancer_M[j],".Rdata")))

  Gene_Transform <-TCGA_profile[,1:3]
  Gene_protein <- Gene_Transform$Symbol_Name[which(Gene_Transform$Gene_Type=="protein_coding")]
  
  TCGA_profile <- TCGA_profile[!duplicated(TCGA_profile$Symbol_Name, fromLast = FALSE),]
  rownames(TCGA_profile) <- TCGA_profile$Symbol_Name
  exp <- TCGA_profile[,-c(1:3)]
  
  exp_persent<-apply(exp,1,function(i){length(which(i==0))/length(i)})
  exp <- exp[which(exp_persent<=0.3),]
  
  ###############################################################################
  
  ajcc_pathologic_m <- c()
  for (i in 1:length(TCGA_clinical$diagnoses)) {
    if (length(TCGA_clinical$diagnoses[[i]]) == 0 || is.null(TCGA_clinical$diagnoses[[i]][[1]])) {
      ajcc_pathologic_m <- c(ajcc_pathologic_m,"NULL")
    } else {
      if(is.null(TCGA_clinical$diagnoses[[i]]$ajcc_pathologic_m))(
        ajcc_pathologic_m <- c(ajcc_pathologic_m,"NULL")
      )else{
        ajcc_pathologic_m <- c(ajcc_pathologic_m,TCGA_clinical$diagnoses[[i]]$ajcc_pathologic_m)
      }
    }
  }
  
  TCGA_clinical$M_state <- ifelse(ajcc_pathologic_m=="M0","M0",
                                  ifelse(ajcc_pathologic_m=="M1","M1","Others"))
  
  TCGA_clinical <- TCGA_clinical[which(TCGA_clinical$M_state%in% c("M0","M1")),] 
  
  exp <- exp[,which(substring(colnames(exp),1,12)%in%TCGA_clinical$submitter_id)] 
  exp <- log2(exp+1)%>%round(.,3)
  range(exp)
  
  group_list <- c()
  for(i in 1:ncol(exp)){
    if(TCGA_clinical$M_state[match(substring(colnames(exp)[i],1,12),TCGA_clinical$submitter_id)]=="M1"){
      group_list <- c(group_list,"metastasis")
    }else{
      group_list <- c(group_list,"primary")
    }
  }
  group_list = factor(group_list,levels = c("primary","metastasis")) 
  
  if(as.data.frame(table(group_list))[2,2]<2){
    cat("Cancer",Cancer_M[j]," has only one metastatic sample\n")
    next
  }
  
  ###############################################################################################################################
  ###############################################################################################################################
  
  cat("Differential gene expression analysis is in progress.\n")  
  design=model.matrix(~group_list)
  fit=lmFit(exp,design)
  fit=eBayes(fit)
  deg=topTable(fit,coef=2,number = Inf)
  head(deg)
  
  DEG_up <- subset(deg, P.Value<0.05&logFC>0.5&adj.P.Val < 0.2)
  DEG_down <- subset(deg, P.Value<0.05&logFC<(-0.5)&adj.P.Val < 0.2)
  
  # DEG_up########################################################################
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
  
  # DEG_down####################################################################
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
  
  ## Merge the two matrices ##########################################################################
  DEG_up_top100_merge <- apply(DEG_up_top100,2,function(x){paste(x,collapse = ";")}) %>% as.list()
  DEG_up_top100_merge_mtx <- do.call(cbind,DEG_up_top100_merge)
  DEG_down_top100_merge <- apply(DEG_down_top100,2,function(x){paste(x,collapse = ";")}) %>% as.list()
  DEG_down_top100_merge_mtx <- do.call(cbind,DEG_down_top100_merge)
  UPandDOWN <- rbind(DEG_up_top100_merge_mtx,DEG_down_top100_merge_mtx)
  rownames(UPandDOWN) <- c("UP_Metastasis","UP_Primary")
  
  write.csv(UPandDOWN, file.path(TCGA_output_dir,Cancer_M[j],paste0("b_de_TCGA_",Cancer_M[j],".csv")), row.names = T)
  
  ## The second table.##############################################################################
  # X-axis: log2 fold change   Y-axis: -log10 p-value
  
  DEG_volcano <- deg[,c(1,4)]
  DEG_volcano[,2] <- round(-log10(DEG_volcano[,2]),2)
  colnames(DEG_volcano) <- c("log2FC(x)","-log10Pval(y)")
  write.csv(DEG_volcano, file.path(TCGA_output_dir,Cancer_M[j],paste0("b_de2_TCGA_",Cancer_M[j],".csv")), row.names = T)
  
  ###############################################################################################################################
  ###############################################################################################################################
  
  cat("Functional enrichment analysis is in progress.\n")
  library(Matrix)
  library(clusterProfiler)
  library(org.Hs.eg.db) 
  
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
  
  ####################################################################################################################################################

  DEG_up <- mapIds(org.Hs.eg.db,rownames(DEG_up),"ENTREZID","SYMBOL")
  DEG_down <- mapIds(org.Hs.eg.db,rownames(DEG_down),"ENTREZID","SYMBOL")
  DEG_all <- mapIds(org.Hs.eg.db,rownames(DEG_all),"ENTREZID","SYMBOL")
  
  DEG_up <- na.omit(DEG_up)
  DEG_down <- na.omit(DEG_down)
  DEG_all <- na.omit(DEG_all)
  
  ###############################################################
  # GO
  # library(clusterProfiler)   
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
                      readable = T #
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
    go_all <- subset(go_all0@result, ONTOLOGY=="BP"&pvalue<0.05&p.adjust<0.2)
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
  
  write.csv(GOandKEGG, file.path(TCGA_output_dir,Cancer_M[j],paste0("b_fc_TCGA_",Cancer_M[j],".csv")), row.names = T)
  
  ###############################################################################################################################
  ###############################################################################################################################
  
  # https://www.jianshu.com/p/5bfb94b0e250  https://cloud.tencent.com/developer/article/1670951
  
  cat("Protein-protein interaction (PPI) analysis is in progress.\n")
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
  
  write.csv(pi_mtx,file.path(TCGA_output_dir,Cancer_M[j],paste0("b_pi_TCGA_",Cancer_M[j],".csv")), row.names = F)
  
  #############################################################
  
  exp_met <- exp[nodes[,1],which(group_list=="metastasis")]
  metgenes_cor <- cor(t(exp_met), method="pearson") %>% round(.,2)
  
  exp_pri <- exp[nodes[,1],which(group_list=="primary")]
  prigenes_cor <- cor(t(exp_pri), method="pearson") %>% round(.,2)
  
  write.csv(metgenes_cor,file.path(TCGA_output_dir,Cancer_M[j],paste0("b_pim_TCGA_",Cancer_M[j],".csv")), row.names = T)
  write.csv(prigenes_cor,file.path(TCGA_output_dir,Cancer_M[j],paste0("b_pip_TCGA_",Cancer_M[j],".csv")), row.names = T)
  
  ###############################################################################################################################
  ###############################################################################################################################
  
  cat("\n Immune infiltration analysis is in progress.\n")
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
  
  write.csv(cib_proportion_mtx,file.path(TCGA_output_dir,Cancer_M[j],paste0("b_ii_TCGA_",Cancer_M[j],".csv")), row.names = T)
  
  ###############################################################################################################################
  ##############################################################################################################################

  cat("Survival data is being processed\n")
  
  os_status <- TCGA_clinical[[6]]$vital_status[match(substring(colnames(exp),1,12),TCGA_clinical$submitter_id)]
  os_time <- TCGA_clinical[[6]]$days_to_death[match(substring(colnames(exp),1,12),TCGA_clinical$submitter_id)]
  
  cox_matrix <- cbind(cbind(os_status,os_time),t(exp)) 

  for(i in 1:nrow(cox_matrix)){
    if(cox_matrix[i,1]=="Dead"){
      cox_matrix[i,1] <- "1"
    }else{
      cox_matrix[i,1] <- "0"
    }
    cox_matrix[i,2] <- round((as.numeric(cox_matrix[i,2])/30),1)
  }
  colnames(cox_matrix) <- gsub('-', '_', colnames(cox_matrix))
  
  cox_matrix_primary <- cox_matrix[which(group_list=="primary"),] %>% as.data.frame()
  for(i in 1:ncol(cox_matrix_primary))(
    cox_matrix_primary[,i] <- as.numeric(cox_matrix_primary[,i])
  )
  pfilter <- 0.05
  

  cox_result_primary <- data.frame()  
  cat("Processing:cox_p\n")
  for(i in colnames(cox_matrix_primary[,3:ncol(cox_matrix_primary)])){ 
    unicox <- coxph(Surv(time = os_time, event = os_status) ~ cox_matrix_primary[,i], data = cox_matrix_primary)  
    unisum<- summary(unicox)   
    pvalue <- round(unisum$coefficients[,5],3) 
    if(is.na(pvalue)){
      next
    }
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
  
  cox_matrix_metastasis <- cox_matrix[which(group_list=="metastasis"),] %>% as.data.frame()
  for(i in 1:ncol(cox_matrix_metastasis))(
    cox_matrix_metastasis[,i] <- as.numeric(cox_matrix_metastasis[,i])
  )
  
  cox_result_metastasis <- data.frame() 
  cat("Processing:cox_m\n")

  for(i in colnames(cox_matrix_metastasis[,3:ncol(cox_matrix_metastasis)])){ 
    unicox <-coxph(Surv(time = os_time, event = os_status) ~ cox_matrix_metastasis[,i], data = cox_matrix_metastasis)  
    unisum<- summary(unicox)   
    pvalue <- round(unisum$coefficients[,5],3) 
    if(is.na(pvalue)){
      next
    }
    if(pvalue < pfilter){ 
      cox_result_metastasis <- rbind(cox_result_metastasis,
                                  cbind(gene=i,
                                        HR=round(unisum$coefficients[,2],2),
                                        L95CI=round(unisum$conf.int[,3],2),
                                        H95CI=round(unisum$conf.int[,4],2),
                                        pvalue=round(unisum$coefficients[,5],2)
                                  ))
    }
  }  
  
  gene_cox <- union(cox_result_primary$gene,cox_result_metastasis$gene)
  mtx_survival_primary <- cbind(cox_matrix_primary[,c(2,1)],cox_matrix_primary[,gene_cox])
  colnames(mtx_survival_primary)[1:2] <- c("time","event")
  
  mtx_survival_metastasis <- cbind(cox_matrix_metastasis[,c(2,1)],cox_matrix_metastasis[,gene_cox])
  colnames(mtx_survival_metastasis)[1:2] <- c("time","event")
  
  write.table(mtx_survival_primary,file.path(TCGA_output_dir,Cancer_M[j],paste0("b_sv_TCGA_",Cancer_M[j],"_p.txt")), sep = "\t",quote = FALSE,row.names = F)
  write.table(mtx_survival_metastasis,file.path(TCGA_output_dir,Cancer_M[j],paste0("b_sv_TCGA_",Cancer_M[j],"_m.txt")),sep = "\t",quote = FALSE,row.names = F)
  
}

# TGCT: There is no survival information for this cancer metastasis.
# The cancer HNSC has only one metastatic sample.


