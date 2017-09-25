# R script for tissue-specific CpGs and genes.
# Cong Liu
# 2917/09/06

rm(list=ls())

# args = commandArgs(trailingOnly=TRUE)
# if (length(args)==0) {
#   stop("At least one arguments must be supplied.\n", call.=FALSE)
# }else{
#   sub_list_seq <- args # t test will be conducted only for these cpg subsets. 
#   cat(args)
# }

# load packages.

library(stringr)

# define parameters.
na_threshold <- 0.25 # only test CpG with not too many NAs in samples.
path <- "~/tcga_download/tissueSet"
setwd(path)


# (1) read sample.info


# added 20170919.
# add sample_type information by reading TCGA barcode.
if(!file.exists("sample_file_type.rda")){
  tissue_group_dir <- list.dirs(path = path,full.names = TRUE,recursive = F)
  for(tissue_group in tissue_group_dir){
    sample_ids <- list.dirs(path = tissue_group,full.names = FALSE,recursive = F)
    for(sample_dir in sample_ids){
      file_name <- list.files(path = paste(tissue_group,sample_dir,sep = "/"),pattern = "*gdc_hg38.txt$",full.names = F,recursive = F,all.files = F,include.dirs = F)
      sample_type <- gsub(".*TCGA-.+?-.+?-(\\d+)[A-Z]-.*","\\1",file_name,perl = T)
      tissue_type <- gsub("jhu-usc.edu_(\\w+).HumanMethylation450.*","\\1",file_name,perl = T)
      df.tmp <- data.frame(sample_id = sample_dir,file_name = file_name,tissue_type = tissue_type,sample_type = sample_type)
      if(!exists("sample_file_type")){
        sample_file_type <- df.tmp
      }else{
        sample_file_type <- rbind(sample_file_type,df.tmp)
      }
    }
  }
  tcga_project_file <- read.table(paste(path,"/tcgaProject.txt",sep=""),sep="\t",head=T)
  rownames(tcga_project_file) <- tcga_project_file[,1]
  tissue_for_each_sample <- tcga_project_file[as.character(sample_file_type$tissue_type),"Primary_Site"]  
  sample_file_type <- data.frame(sample_file_type,tissue = tissue_for_each_sample)
  #tissue_set <- levels(tcga_project_file$Primary_Site)
  rownames(sample_file_type) <- sample_file_type$sample_id
  
  # reorder sample to make it consistent with bigMatrix
  for(tissue_group in tissue_group_dir){
    sample_group_file <- paste(tissue_group,"betaMat_sample_group.txt",sep = "/")
    betaMat_sample_group <- read.table(sample_group_file,header = F,sep = "\t")
    if(!exists("sample_info")){
      sample_info <- betaMat_sample_group
    }else{
      sample_info <- rbind(sample_info,betaMat_sample_group)
    }
  }
  sample_file_type <- sample_file_type[as.character(sample_info$V1),]
  save(sample_file_type,file = "sample_file_type.rda")
  
}else{
  load("sample_file_type.rda")
}



# (2) identify CpG.
# cat("assign tissue group for each sample.\n")
# tcga_project_file <- read.table(paste(path,"/tcgaProject.txt",sep=""),sep="\t",head=T)
# rownames(tcga_project_file) <- tcga_project_file[,1]
# tissue_for_each_sample <- tcga_project_file[as.character(sample_info$tissue),]  
# tissue_for_each_sample <- tissue_for_each_sample$Primary_Site
# tissue_set <- levels(tcga_project_file$Primary_Site)
# or select your own tissues.

# perform two-group t-test
my_t.test <- function(x,tissue,tissue_info,tissue_set){
  out <- tryCatch(
    {
      # modified in 09/20/2017.
      sample_tissue_idx <- sample(which(tissue_info %in% tissue),size = 500, replace = TRUE)
      sample_remain_idx <- unlist(lapply(tissue_set[tissue_set!=tissue],function(x) sample(which(tissue_info %in% x), size = ceiling(500/length(tissue_set[tissue_set!=tissue])),replace = TRUE)))
      sample_blood_idx <- sample(which(tissue_info %in% "blood"),size = 500,replace = TRUE)
      # x1 is the tissue value.
      # x2 is the other tissue value.
      # x3 is the blood value.
      x1 <- as.numeric(x[sample_tissue_idx])
      x2 <- as.numeric(x[sample_remain_idx])
      x3 <- as.numeric(x[sample_blood_idx])
      
      # avoid the Error in if (stderr < 10 * .Machine$double.eps * max(abs(mx), abs(my))) stop("data are essentially constant") : 
      # missing value where TRUE/FALSE needed
      x1[is.infinite(x1)] <- NA 
      x2[is.infinite(x2)] <- NA
      x3[is.infinite(x3)] <- NA
      
      
      # skip too many NA CpGs.
      x1_na_percentage <- sum(is.na(x1))/length(x1)
      x2_na_percentage <- sum(is.na(x2))/length(x2)
      x3_na_percentage <- sum(is.na(x3))/length(x3)
      
      # remove tissue no CpG measured.
      if(is.na(x1_na_percentage) | is.na(x2_na_percentage) | is.na(x3_na_percentage)){
        p_value <- 9.999
        tissue_m_mean <- 9.999
        data.frame(p1_value=p_value,p2_value=p_value,tissue_m_mean=tissue_m_mean)
      }else{
        # remove Cpg less than NA threshold.
        if((x1_na_percentage > na_threshold) | (x2_na_percentage > na_threshold) | (x3_na_percentage > na_threshold)){
          p_value <- 1.001
          tissue_m_mean <- mean(x1,na.rm = TRUE)
          data.frame(p1_value=p_value,p2_value=p_value,tissue_m_mean=tissue_m_mean)
        }else{
          
          # compare tissue and other tissue.
          result_1 <- t.test(x1,x2,alternative = "greater")
          p_value_1 <- result_1$p.value
          # compare tissue and blood.
          result_2 <- t.test(x1,x3,alternative = "greater")
          p_value_2<- result_2$p.value
          
          tissue_m_mean <- as.numeric(mean(x1,na.rm = TRUE))
          
          data.frame(p1_value=p_value_1,p2_value=p_value_2,tissue_m_mean=tissue_m_mean)
        }
      }
    },
    error=function(cond) {
      #message("Here's the original error message:")
      message(paste(cond,"\n"))
      # Choose a return value in case of error
      return(data.frame(p1_value=NA,p2_value=NA,tissue_m_mean=NA))
    },
    warning=function(cond) {
      #message("Here's the original warning message:")
      message(paste(cond,"\n"))
      # Choose a return value in case of warning
      return(data.frame(p1_value=NA,p2_value=NA,tissue_m_mean=NA))
    },
    finally={
      # NOTE:
      # Here goes everything that should be executed at the end,
      # regardless of success or error.
      # If you want more than one expression to be executed, then you 
      # need to wrap them in curly brackets ({...}); otherwise you could
      # just have written 'finally=<expression>' 
      # message(paste("Processed URL:", url))
      # message("Some other message at the end")
    }
  )    
  return(out)
}

beta2M <- function(beta){
  log2(beta/(1-beta))
}

load("whole_blood_cpg_df.rda")
sub_list_seq <- c(1:98)
for(sub_list_number in sub_list_seq){
  cat("t test for cpg group: ",sub_list_number,"\n")
  tissue_group_dir <- list.dirs(path = path,full.names = TRUE,recursive = F)
  
  if(exists("beta_info")){
    rm("beta_info")
  }
  if(exists("common_cpg")){
    rm("common_cpg")
  }
  # read whole blood sample.
  for(tissue_group in tissue_group_dir){
    sub_cpg_list_file <- paste(tissue_group,"/betaMat_",sub_list_number,".txt",sep="")
    cat("reading beta matrix: ",sub_cpg_list_file,"\n")
    # betaMat_sub <- read.table(sub_cpg_list_file,header = F,row.names = 1,sep = "\t",nrows = 10) # for test only.
    betaMat_sub <- read.table(sub_cpg_list_file,header = F,row.names = 1,sep = "\t")
    if(!exists("common_cpg")){
      common_cpg <- intersect(rownames(betaMat_sub),rownames(whole_blood_cpg_df))
    }
    if(!exists("beta_info")){
      beta_info <- betaMat_sub[common_cpg,]
    }else{
      beta_info <- cbind(beta_info, betaMat_sub[common_cpg,])
    }
  }
  
  # include whole blood cpg.
  beta_info <- cbind(beta_info,whole_blood_cpg_df[common_cpg,])
  
  # convert to m_value
  m_info <- beta2M(beta_info)
  
  # tissue information.
  tissue_info <- c(as.character(sample_file_type$tissue),rep("blood",dim(whole_blood_cpg_df)[2]))
  
  # sample type information.
  sample_type_info <- c(as.character(sample_file_type$sample_type),rep("03",dim(whole_blood_cpg_df)[2]))
  
  # remove Metastatic samples. 06 and 07.
  meta_static_idx <- which(sample_type_info %in% c("06","07"))
  m_info <- m_info[,-meta_static_idx]
  tissue_info <- tissue_info[-meta_static_idx]
  sample_type_info <- sample_type_info[-meta_static_idx]
  sample_name <- as.character(sample_file_type$sample_id)[-meta_static_idx]
  
  tissue_set <- unique(tissue_info)
  
  # t-test for each tissue.
  if(exists("tTestResult")){
    rm("tTestResult")
  }
  for(tissue in tissue_set){
    cat("t-test for :",tissue,"\n")
    if(!exists("tTestResult")){
      tTestResult<- do.call(rbind,apply(m_info,1,my_t.test,tissue=tissue,tissue_info=tissue_info,tissue_set=tissue_set))
      colnames(tTestResult) <- c(paste(tissue,"pvalue_1",sep = "_"),paste(tissue,"pvalue_2",sep = "_"),paste(tissue,"meanM",sep = "_"))
    }else{
      tmp <- do.call(rbind,apply(m_info,1,my_t.test,tissue=tissue,tissue_info=tissue_info,tissue_set=tissue_set))
      colnames(tmp) <- c(paste(tissue,"pvalue_1",sep = "_"),paste(tissue,"pvalue_2",sep = "_"),paste(tissue,"meanM",sep = "_"))
      tTestResult <- cbind(tTestResult,tmp)
    }
  }
  tTestResult=(as.matrix(tTestResult))
  #tTestResult <- round(tTestResult,digits = 4)
  write.table(tTestResult,file = paste("tTestResult_",sub_list_number,".txt",sep = ""),sep="\t",row.names = T,quote = F)
}
  




