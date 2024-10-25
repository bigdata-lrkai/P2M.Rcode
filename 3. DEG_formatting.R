##################################################################################################################################################################

## Gene Differential Expression Analysis
## 1) On the left and right, display the results of highly expressed genes in primary and metastatic samples.
## 2) In the middle, showcase the differential expression results of all differentially expressed genes across various cell types.

##################################################################################################################################################################
##################################################################################################################################################################

dir <- file.path("E:/Data/",GSE_number)
output_dir <- file.path("E:/Data/",GSE_number,"output")

files_full<- list.files(path=dir,full.names = TRUE)
files <-basename(files_full)

## 1) 

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
    
    pattern <- "[Aa]llcells_diffgenes.csv" 
    file_allcell_DEG <- list.files(path=dir_site_annot_Method_DEG) 
    matching_file <- grep(pattern, file_allcell_DEG, value = TRUE, ignore.case = TRUE) 
    matching_file_path <- file.path(dir_site_annot_Method_DEG,matching_file)
    
    diffgenes <- read.csv(matching_file_path,row.names = 1,header = TRUE)
    
    DEG_up <- subset(diffgenes, p_val<0.05&avg_log2FC>0.5&p_val_adj < 0.2)
    DEG_down <- subset(diffgenes, p_val<0.05&avg_log2FC<(-0.5)&p_val_adj < 0.2)
    
    # DEG_up ########################################################################
    sorted_indices <- order(DEG_up$avg_log2FC,decreasing = T)
    DEG_up_sort <- DEG_up[sorted_indices, ]
    
    if(nrow(DEG_up_sort)>100){ 
      DEG_up_top100 <- DEG_up_sort[1:100,c("avg_log2FC","p_val_adj")]
    }else{ 
      DEG_up_top100 <- DEG_up_sort[,c("avg_log2FC","p_val_adj")]
    }
    
    DEG_up_top100$avg_log2FC <- round(DEG_up_top100$avg_log2FC,2)
    DEG_up_top100 <- cbind(rownames(DEG_up_top100),DEG_up_top100$avg_log2FC)
    colnames(DEG_up_top100) <- c("Gene","avg_log2FC")
    
    # DEG_down ####################################################################
    sorted_indices <- order(DEG_down$avg_log2FC)
    DEG_down_sort <- DEG_down[sorted_indices, ]
    
    if(nrow(DEG_down_sort)>100){
      DEG_down_top100 <- DEG_down_sort[1:100,c("avg_log2FC","p_val_adj")]
    }else{
      DEG_down_top100 <- DEG_down_sort[,c("avg_log2FC","p_val_adj")]
    }
    
    DEG_down_top100$avg_log2FC <- round(DEG_down_top100$avg_log2FC,2)
    DEG_down_top100 <- cbind(rownames(DEG_down_top100),DEG_down_top100$avg_log2FC)
    colnames(DEG_down_top100) <- c("Gene","avg_log2FC")
    
    
    ## Merge two matrices ##########################################################################
    DEG_up_top100_merge <- apply(DEG_up_top100,2,function(x){paste(x,collapse = ";")}) %>% as.list()
    DEG_up_top100_merge_mtx <- do.call(cbind,DEG_up_top100_merge)
    DEG_down_top100_merge <- apply(DEG_down_top100,2,function(x){paste(x,collapse = ";")}) %>% as.list()
    DEG_down_top100_merge_mtx <- do.call(cbind,DEG_down_top100_merge)
    UPandDOWN <- rbind(DEG_up_top100_merge_mtx,DEG_down_top100_merge_mtx)
    rownames(UPandDOWN) <- c("UP_Metastasis","UP_Primary")
    
    name_output_UPandDOWN <- paste0(files[f],"_de_",site_all[s],"_all.csv")
    write.csv(UPandDOWN, file.path(output_dir,name_output_UPandDOWN), row.names = T)
    
    ###############################################################################################################################
    
    # Each site also needs to separately output the results of various cell types annotated for the three different cell annotations.
    
    matching_site_files <- grep(pattern_site, basename(dir_annot_Method), value = FALSE, ignore.case = TRUE)
    
    for(a in 1:length(matching_site_files)){ 
      
      fsa <- dir_annot_Method[matching_site_files[a]] 
      
      all_files_fsa <- list.files(path = file.path(fsa), full.names = TRUE)
      files_de_celltype <- all_files_fsa[!grepl("^[aA]ll", basename(all_files_fsa))& grepl("diffgenes.csv$", basename(all_files_fsa))]
      
      #############################################################
      # The name for the output file is: de_EachCellType
      de_EachCellType <- NULL
      
      for(i in 1:length(files_de_celltype)){
        temp <- read.csv( files_de_celltype[i],header = TRUE)
        temp <- subset(temp, p_val<0.05&abs(avg_log2FC)>0.5&p_val_adj < 0.2)
        
        if (nrow(temp) == 0) {
          next  
        }else{
          sorted_indices <- order(temp$p_val) 
          temp_sort <- temp[sorted_indices, ]
          
          if(nrow(temp_sort)>100){
            temp_top100 <- temp_sort[1:100,c(1,3)]
          }else{
            temp_top100 <- temp_sort[,c(1,3)] }
          
          temp_top100$avg_log2FC <- round(temp_top100$avg_log2FC,2)
          
          ## The horizontal coordinates represent one unit for each cell type.
          old_min <- 0
          old_max <- 1
          
          # Target value range.
          new_min <- 0.2
          new_max <- 0.8
          
          # Original value.
          old_value <- sample(round((1:nrow(temp_top100))/nrow(temp_top100),2))
          
          # Normalize using the linear transformation formula.
          new_value <- (new_max - new_min) * (old_value - old_min) / (old_max - old_min) + new_min
          
          xlabel <- new_value+i-1
          
          de_EachCellType <- rbind(de_EachCellType,cbind(paste(temp_top100[,1], collapse = ";"),paste(temp_top100[,2], collapse = ";"),paste(xlabel, collapse = ";")))
          rownames(de_EachCellType)[nrow(de_EachCellType)] <- as.character(unlist(strsplit(basename(files_de_celltype[i]),split="_dif"))[1])
          
        }
      }
      de_EachCellType %<>% as.data.frame()
      colnames(de_EachCellType) <- c("DEG","avg_log2FC","xlabel")
      
      # Obtain the current annotation method.
      annotate_method <- substr(fsa, nchar(fsa)-1,nchar(fsa))
      name_de_EachCellType <- paste0(files[f],"_de_",site_all[s],"_",annotate_method,".csv")
      write.csv(de_EachCellType, file.path(output_dir,name_de_EachCellType), row.names = T)
      
    }
  }
}
