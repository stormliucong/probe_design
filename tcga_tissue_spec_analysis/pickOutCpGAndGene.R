# R script for pick out tissue-specific CpGs and genes.
# Cong Liu
# 2917/09/08
rm(list=ls())
# path = "/home/liuc/tcga_download/tTestResult_raw"
path = "~/微云同步盘/EPN/probe_design/tissue_marker/tcga_analysis"


# define functions.
tTest2binary <- function(x){
  y <- matrix(0,ncol = length(x)/2)
  pval <- x[seq(from=1,to = length(x),by = 2)]
  mean_m <- x[seq(from=2,to = length(x),by = 2)]
  sig_idx <- which(pval < 1e-5 & mean_m > 0) # pval < 1e-5 and mean beta > 0.5. 
  y[,sig_idx] <- 1
  return(y)
}

addHashByGene <- function(x){
  x = as.character(x)
  cpg_id <- x[1]
  gene_list <- unlist(strsplit(x[6],split = ";"))
  for(i in unique(gene_list)){
    if(!exists(i, envir = genehash, inherits = FALSE)){
      genehash[[i]] <- cpg_id
    }else{
      genehash[[i]] <- c(genehash[[i]],cpg_id)
    }
  }
}

scoreForGene <- function(x){
  # x is a vector of cpg ids.
  if(length(x) > 1){
    sub_y <- sig_cpg[x,]
  }else{
    sub_y <- matrix(sig_cpg[x,],nrow=1) # this is a corrected bug to avoid return 1 number.
    colnames(sub_y) <- names(sig_cpg[x,])
    rownames(sub_y) <- x
  }
  # total_sum <- max(sum(sum(sub_y)),0.1) # avoid /0.
  
  s <- apply(sub_y,2,function(x) sum(x))  
  return(s)
}

gene2tissue <- function(x){
  # 1.5 IQR + upper_q
  q <- as.numeric(summary(x))
  iqr <- q[5] - q[2]
  outlier_threshold <- 1.5*iqr + q[5]
  idx <- which(x > 5 & x == max(x) & x > outlier_threshold)
  if(length(idx) == 0){
    return(c(NA,NA,NA,NA))
  }else{
    if(length(idx) == 1){
      return(c(tissue_name[idx],x[idx],x[idx]/mean(x),T))
    }else{
      # randomly select one tissue if tie.
      idx1 <- sample(x = idx,size = 1,replace = F)
      return(c(tissue_name[idx1],x[idx1],x[idx1]/mean(x),F))
    }
  }
}

score_over_q3 <- function(x){
  y <- x
  q <- as.numeric(summary(x))
  iqr <- q[5] - q[2]
  y <- (x-q[5]+1)/(iqr+1) # add 1 psedo count.
  y[x < 5] <- NA
  return(y)
}

top20Gene <- function(x){
  idx <- which(rank(-x) <= 20)
  gene <- gene_list[idx]
  x[idx]
  return(data.frame(gene=gene,score_over_q3=x[idx]))
}


# read tTest results.
if(!file.exists("tTestResult_all.txt")){
  file_seq = c(1:99)
  file_list <- paste(path,"/tTestResult_",file_seq,".txt",sep = "")
  for(file in file_list){
    tmp <- read.table(file,header = T,sep = "\t",stringsAsFactors = F)
    if(!exists("t_test_result")){
      t_test_result <- tmp
    }else{
      t_test_result <- rbind(t_test_result,tmp)
    }
  }
  write.table(t_test_result,file = "tTestResult_all.txt",sep = "\t",col.names = T,row.names = T,quote = F)
}else{
  tTestResult_all_file <- "tTestResult_all.txt"
  t_test_result <- read.table(tTestResult_all_file,header = T,sep = "\t",stringsAsFactors = F)
}

# get significant cpg 0/1 matrix.
if(!file.exists("sig_cpg.rda")){
  tissue_wild <- colnames(t_test_result)[seq(from=1,to=dim(t_test_result)[2],by = 2)]
  tissue_name <- gsub("^(.+?)_pvalue","\\1",tissue_wild)
  sig_cpg <- apply(t_test_result,1,tTest2binary) 
  # this step is very quick( less than 1 minutes.)
  # So different criteria could be used to filter sig cpgs.
  sig_cpg <- t(sig_cpg)
  colnames(sig_cpg) <- tissue_name
  save(sig_cpg,file = "sig_cpg.rda")
}else{
  load("sig_cpg.rda")
}


# process cpg and gene relationship.
if(!file.exists("geneHash.rda")){
  cpg_anno_file <- "betaMat_cpg_annotation.txt"
  cpg_anno <- read.table(file = cpg_anno_file,header = F,sep = "\t",stringsAsFactors = F)  
  genehash <- new.env(hash=TRUE)
  apply(cpg_anno,1,addHashByGene) # takes ~ 2mins.
  # length(ls(genehash)) = 35556 
  save(genehash,file = "geneHash.rda")
  gene_list <- as.list(genehash) 
}else{
  load("geneHash.rda")
  gene_list <- as.list(genehash) 
}

# process cpg and gene relationship 2.
# only include gene_names in hg19_genes.gtf

if(!file.exists("gene_list.rda")){
  load("geneHash.rda")
  gene_list <- as.list(genehash) 
  path = "/database/hg19/Gencode"
  hg19_genes_gtf <- read.table(paste(path,"hg19_genes.gtf",sep = "/"),header = F,sep = "\t",stringsAsFactors = F)
  hg19_genes_gtf_anno <- hg19_genes_gtf$V9
  gene_id <- gsub("gene_id (.+?);.*","\\1",hg19_genes_gtf_anno,perl=T)
  gene_id <- unique(gene_id)
  # length(gene_id) = 23368
  gene_id <- intersect(gene_id,names(gene_list))
  # length(gene_id) = 19991
  gene_list <- gene_list[gene_id]
  # remove NA_list.
  save(gene_list,file = "gene_list.rda")
}else{
  load("gene_list.rda")
}

# generate gene tissue specific score.
if(!file.exists("gene_score.rda")){
  gene_score <- lapply(gene_list,scoreForGene) # takes about 15 mins.
  save(gene_score,file = "gene_score.rda")
  gene_score_df <- as.data.frame(gene_score)
}else{
  load("gene_score.rda")
  gene_score_df <- as.data.frame(gene_score)
}

# generate gene-centric list.
tissue_name <- rownames(gene_score_df)
gene2tissue_list <- t(apply(gene_score_df,2,gene2tissue))
colnames(gene2tissue_list) <- c("tissue","num_spec_cpg","over_average","unique_tissue")
gene2tissue_list <- as.data.frame(gene2tissue_list)
gene2tissue_list_spec <- gene2tissue_list[-which(is.na(gene2tissue_list$tissue)),]
table(gene2tissue_list_spec$tissue)

write.table(table(gene2tissue_list_spec$tissue),"tmp.txt",quote = F,sep = "\t")
write.table(gene2tissue_list_spec,"gene2tissue_list.txt",sep = "\t",append = F,row.names = T,col.names = T)

# generate tissue-centric list.
gene_score_over_q3 <- t(apply(gene_score_df,2,score_over_q3))
save(gene_score_over_q3,file = "gene_score_over_q3.rda")
gene_list <- rownames(gene_score_over_q3)
tissue2gene_top20 <- apply(gene_score_over_q3,2,top20Gene)
save(tissue2gene_top20,file = "tissue2gene_top20.rda")

# plot tissue-centric based heatmap
top20_gene <- (as.character(unlist(lapply(tissue2gene_top20,rownames))))
gene_score_df_top20 <- t(gene_score_df[,top20_gene])
library(plotly)
#gene_score_df_top20 <- apply(gene_score_df_top20, 1, function(x){x/mean(x)})
plot_ly(x=colnames(gene_score_df_top20), y=rownames(gene_score_df_top20), z = gene_score_df_top20, type = "heatmap",zmin=0,zmax=15)

# plot gene-centric based heatmap
tissue_sort <- names(sort(table(gene2tissue_list_spec$tissue)))
for(t in tissue_sort){
  if(!exists("gene2tissue_list_spec_sort")){
    gene2tissue_list_spec_sort <- gene2tissue_list_spec[which(gene2tissue_list_spec$tissue == t),]
  }else{
    gene2tissue_list_spec_sort <- rbind(gene2tissue_list_spec_sort,gene2tissue_list_spec[which(gene2tissue_list_spec$tissue == t),])
  }
}

gene2tissue_list_spec[match(gene2tissue_list_spec$tissue,tissue_sort),]
spec_gene <- rownames(gene2tissue_list_spec_sort)
gene_score_df_spec_gene <- t(gene_score_df[,spec_gene])
# sort tissue column by number of
gene_score_df_spec_gene <- gene_score_df_spec_gene[,tissue_sort]
#gene_score_df_spec_gene <- apply(gene_score_df_spec_gene, 1, function(x){x/mean(x)})
plot_ly(x=colnames(gene_score_df_spec_gene), y=rownames(gene_score_df_spec_gene),z = gene_score_df_spec_gene, type = "heatmap",zmin=0,zmax=15)
