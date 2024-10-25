####################################################################################################################################################

## Functional Enrichment Analysis
## 1) The functional enrichment results of highly expressed genes in primary and metastatic samples are presented on the sides.
## 2) In the middle, the enrichment results of all differentially expressed genes in various cell types are shown separately.

####################################################################################################################################################
####################################################################################################################################################

library(org.Hs.eg.db)
library(clusterProfiler)

dir <- file.path("E:/Data/",GSE_number)
output_dir <- file.path("E:/Data/",GSE_number,"output")

files_full<- list.files(path=dir,full.names = TRUE)
files <-basename(files_full)

## 1) The left and right sides display the functional enrichment results of highly expressed genes in primary and metastatic samples, respectively.

# Load the previously obtained differential genes between primary and metastatic samples across all cell types.

for(f in 1:1){
    
    cat("calculating: ",files[f],"\n")
    dir_annot_Method <- list.files(files_full[f],full.names = TRUE)
      
    site_all <- unique(  sapply(dir_annot_Method,function(dir){
      unlist(strsplit(basename(dir),split="_"))[2]
    }))
    
    for(s in 1:length(site_all)){
      pattern_site <- paste0(files[f],"_",site_all[s])
      matching_site_file <- grep(pattern_site, basename(dir_annot_Method), value = FALSE, ignore.case = TRUE)[1]
      dir_site_annot_Method_DEG <- file.path(dir_annot_Method[matching_site_file])
      
      pattern <- "[Aa]llcells?_diffgenes\\.csv" 
      file_allcell_DEG <- list.files(path=dir_site_annot_Method_DEG)
      matching_file <- grep(pattern, file_allcell_DEG, value = TRUE, ignore.case = TRUE)
      matching_file_path <- file.path(dir_site_annot_Method_DEG,matching_file)
      
      diffgenes <- read.csv(matching_file_path,row.names = 1,header = TRUE)
      
      DEG_up <- subset(diffgenes, p_val<0.05&avg_log2FC>0.5&p_val_adj < 0.2)
      DEG_down <- subset(diffgenes, p_val<0.05&avg_log2FC<(-0.5)&p_val_adj < 0.2)
      
      sorted_indices <- order(DEG_up$avg_log2FC,decreasing = T)
      DEG_up_sort <- DEG_up[sorted_indices, ]
      if(nrow(DEG_up_sort)>50){
        DEG_up <- DEG_up_sort[1:50,]}
      
      sorted_indices <- order(DEG_down$avg_log2FC)
      DEG_down_sort <- DEG_down[sorted_indices, ]
      if(nrow(DEG_down_sort)>50){
        DEG_down <- DEG_down_sort[1:50,]}
      
      ####################################################################################################################################################
      
      library(org.Hs.eg.db)
      DEG_up <- mapIds(org.Hs.eg.db,rownames(DEG_up),"ENTREZID","SYMBOL")
      DEG_down <- mapIds(org.Hs.eg.db,rownames(DEG_down),"ENTREZID","SYMBOL")
      
      DEG_up <- na.omit(DEG_up)
      DEG_down <- na.omit(DEG_down)
      
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
      
      if (!is.null(go_up0)) {
        
        go_up <- subset(go_up0@result, ONTOLOGY=="BP"&pvalue<0.05&p.adjust<0.2)
        if(nrow(go_up)<=0){
          go_up <- subset(go_up0@result,pvalue<0.05&p.adjust<0.2)
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
          go_down <- subset(go_down0@result,pvalue<0.05&p.adjust<0.2)
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
                              pAdjustMethod =
                              qvalueCutoff = 0.2)
      # View(kegg_down@result)
      
      if (!is.null(kegg_up)) {
        kegg_up <- subset(kegg_up@result,pvalue<0.05)
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
      
      #############################################################
      
      if (!is.null(kegg_down)) {
        kegg_down <- subset(kegg_down@result,pvalue<0.05&p.adjust<0.2)
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
      
      
      #############################################################
      # Put the four tables together
      
      go_up_top10_cbind <- cbind(paste(go_up_top10[,1], collapse = ";"),paste(go_up_top10[,2], collapse = ";"))
      go_down_top10_cbind<- cbind(paste(go_down_top10[,1], collapse = ";"),paste(go_down_top10[,2], collapse = ";"))
      go_all <- rbind(go_up_top10_cbind,go_down_top10_cbind)
      
      kegg_up_top10_cbind<- cbind(paste(kegg_up_top10[,1], collapse = ";"),paste(kegg_up_top10[,2], collapse = ";"))
      kegg_down_top10_cbind<- cbind(paste(kegg_down_top10[,1], collapse = ";"),paste(kegg_down_top10[,2], collapse = ";"))
      kegg_all <- rbind(kegg_up_top10_cbind,kegg_down_top10_cbind)
      
      GOandKEGG <- rbind(go_all,kegg_all) %>% as.data.frame()
      colnames(GOandKEGG) <- c("Functions or Pathways","P-value(-log10)")
      rownames(GOandKEGG) <- c("GO_Metastasis","GO_Primary","KEGG_Metastasis","KEGG_Primary")
      
      # GOandKEGG <- as.data.frame(apply(GOandKEGG, 2, function(row) {
      #   row[is.na(row) | row == ""] <- "no significant results"
      #   return(row)
      # }))
      
      GOandKEGG[which(GOandKEGG$`Functions or Pathways`== ""),1]<- "no significant results"
      GOandKEGG[which(GOandKEGG$`Functions or Pathways`== "no significant results"),2]<- 0
      
      name_output_GOandKEGG <- paste0(files[f],"_",site_all[s],"_GandK.csv")
      write.csv(GOandKEGG, file.path(output_dir,name_output_GOandKEGG), row.names = T)
      
      ####################################################################################################################################################
      ####################################################################################################################################################
      
      ## 2.1) The middle section presents the enrichment results of all differentially expressed genes across various cell types, categorized by the annotation results from SingleR, Sctype, and Article.
      
      # The current loop focuses on a single metastatic site. At this point, an additional annotation method corresponding to the metastatic site needs to be included.
      matching_site_files <- grep(pattern_site, basename(dir_annot_Method), value = FALSE, ignore.case = TRUE)
      for(a in 1:length(matching_site_files)){
        
        fsa <- dir_annot_Method[matching_site_files[a]]
        all_files_fsa <- list.files(path = file.path(fsa), full.names = TRUE)
        
        files_fun_go <- all_files_fsa[!grepl("^[aA]ll", basename(all_files_fsa)) & grepl("go\\.csv$", basename(all_files_fsa),ignore.case = TRUE)]
        files_fun_kegg <- all_files_fsa[!grepl("^[aA]ll", basename(all_files_fsa)) & grepl("kegg\\.csv$", basename(all_files_fsa),ignore.case = TRUE)]
        
        ## Extract the top ten significantly enriched pathway names and P-values (-log10) for each cell type in GO and KEGG.
        
        #############################################################
        # output：go_EachCellType
        go_EachCellType <- NULL
        for(i in 1:length(files_fun_go)){
          temp <- read.csv( files_fun_go[i],row.names = 1,header = TRUE)
          temp <- subset(temp, ONTOLOGY=="BP"&pvalue<0.05&p.adjust<0.2)
          
          if (nrow(temp) <= 0) {
            temp <- subset(temp, &pvalue<0.05&p.adjust<0.2)
          }
          
          if (nrow(temp) > 0) {
            sorted_indices <- order(temp$p.adjust)
            temp_sort <- temp[sorted_indices, ]
            if(nrow(temp_sort)>10){
              temp_top10 <- temp_sort[1:10,c("Description","p.adjust")]
            }else{
              temp_top10 <- temp_sort[,c("Description","p.adjust")]
            }
            temp_top10$p.adjust <- round(-log10(temp_top10$p.adjust),2)
          }else{
            temp_top10 <- NULL
          }
          
          # colnames(temp_top10) <- c("Gene Ontology","P-value(-log10)")
          go_EachCellType <- rbind(go_EachCellType,cbind(paste(temp_top10[,1], collapse = ";"),paste(temp_top10[,2], collapse = ";")))
          rownames(go_EachCellType)[nrow(go_EachCellType)] <- as.character(unlist(strsplit(basename(files_fun_go[i]),split="_GO"))[1])
        }
        go_EachCellType %<>% as.data.frame()
        colnames(go_EachCellType) <- c("Gene Ontology","P-value(-log10)")
        
        go_EachCellType[which(go_EachCellType$`Gene Ontology`== ""),1]<- "no significant results"
        go_EachCellType[which(go_EachCellType$`Gene Ontology`== "no significant results"),2]<- 0
        
        #############################################################
        # output：kegg_EachCellType
        kegg_EachCellType <- NULL
        for(i in 1:length(files_fun_kegg)){
          temp <- read.csv( files_fun_kegg[i],row.names = 1,header = TRUE)
          temp <- subset(temp, &pvalue<0.05&p.adjust<0.2)
          if (nrow(temp) != 0) {
            sorted_indices <- order(temp$p.adjust)
            temp_sort <- temp[sorted_indices, ]
            if(nrow(temp_sort)>10){
              temp_top10 <- temp_sort[1:10,c("Description","p.adjust")]
            }else{
              temp_top10 <- temp_sort[,c("Description","p.adjust")]
            }
            temp_top10$p.adjust <- round(-log10(temp_top10$p.adjust),2)
          }else{
            temp_top10 <- NULL
          }
          
          kegg_EachCellType <- rbind(kegg_EachCellType,cbind(paste(temp_top10[,1], collapse = ";"),paste(temp_top10[,2], collapse = ";")))
          rownames(kegg_EachCellType)[nrow(kegg_EachCellType)] <- as.character(unlist(strsplit(basename(files_fun_go[i]),split="_GO"))[1])
        }
        kegg_EachCellType %<>% as.data.frame()
        colnames(kegg_EachCellType) <- c("KEGG Pathways","P-value(-log10)")
        
        kegg_EachCellType[which(kegg_EachCellType$`KEGG Pathways`== ""),1]<- "no significant results"
        kegg_EachCellType[which(kegg_EachCellType$`KEGG Pathways`== "no significant results"),2]<- 0
        
        # Gets the current annotation method
        annotate_method <- substr(fsa, nchar(fsa)-1,nchar(fsa))
        
        name_go_EachCellType <- paste0(files[f],"_",site_all[s],"_",annotate_method,"_go.csv")
        name_kegg_EachCellType <- paste0(files[f],"_",site_all[s],"_",annotate_method,"_kegg.csv")
        
        write.csv(go_EachCellType, file.path(output_dir,name_go_EachCellType), row.names = T)
        write.csv(kegg_EachCellType, file.path(output_dir,name_kegg_EachCellType), row.names = T)
        
      }
    }
  }
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  