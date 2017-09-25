GSE64495_series_matrix <- read.table(file = "GSE64495_series_matrix.txt",sep = "\t",header = T,fill = T,comment.char = "!",stringsAsFactors = F,)
save(GSE64495_series_matrix,file="GSE64495_series_matrix.rda")
GSE61653_series_matrix <- read.table(file = "GSE61653_series_matrix.txt",sep = "\t",header = T,fill = T,comment.char = "!",stringsAsFactors = F)
GSE61653_annotation <- read.table(file = "GSE61653_series_matrix.txt",sep = "\t",header = F,fill = T,stringsAsFactors = F,nrow=100)
GSE61653_series_matrix <- GSE61653_series_matrix[,1:65]
save(GSE61653_series_matrix,file="GSE61653_series_matrix.rda")
GSE41169_series_matrix <- read.table(file = "GSE41169_series_matrix.txt",sep = "\t",header = T,fill = T,comment.char = "!",stringsAsFactors = F)
save(GSE41169_series_matrix,file="GSE41169_series_matrix.rda")
GSE43975_series_matrix <- read.table(file = "GSE43975_series_matrix.txt",sep = "\t",header = T,fill = T,comment.char = "!",stringsAsFactors = F)
save(GSE43975_series_matrix,file="GSE43975_series_matrix.rda")
GSE44798_series_matrix <- read.table(file = "GSE44798_series_matrix.txt",sep = "\t",header = T,fill = T,comment.char = "!",stringsAsFactors = F)
save(GSE44798_series_matrix,file="GSE44798_series_matrix.rda")

c1 = GSE64495_series_matrix$ID_REF
c2 = GSE61653_series_matrix$ID_REF
c3 = GSE41169_series_matrix$ID_REF
c4 = GSE43975_series_matrix$ID_REF
c5 = GSE44798_series_matrix$ID_REF

common_cpg <- intersect(intersect(intersect(intersect(c1,c2),c3),c4),c5)
rownames(GSE64495_series_matrix) <- c1
GSE64495_series_matrix_sub <- GSE64495_series_matrix[common_cpg,-1]
rownames(GSE61653_series_matrix) <- c2
GSE61653_series_matrix_sub <- GSE61653_series_matrix[common_cpg,-1]
rownames(GSE41169_series_matrix) <- c3
GSE41169_series_matrix_sub <- GSE41169_series_matrix[common_cpg,-1]
rownames(GSE43975_series_matrix) <- c4
GSE43975_series_matrix_sub <- GSE43975_series_matrix[common_cpg,-1]
rownames(GSE44798_series_matrix) <- c5
GSE44798_series_matrix_sub <- GSE44798_series_matrix[common_cpg,-1]

GSE44798_series_matrix_sub[c(100,2000,30000),1:2]

whole_blood_cpg_matrix <- cbind(GSE64495_series_matrix_sub,GSE61653_series_matrix_sub,GSE41169_series_matrix_sub,GSE43975_series_matrix_sub,GSE44798_series_matrix_sub)
whole_blood_cpg_df <- data.frame(sapply(whole_blood_cpg_matrix,as.numeric))
rownames(whole_blood_cpg_df) = rownames(whole_blood_cpg_matrix)
save(whole_blood_cpg_matrix,file = "whole_blood_cpg_matrix.rda")
save(whole_blood_cpg_df,file = "whole_blood_cpg_df.rda")
