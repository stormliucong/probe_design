library(tidyverse)
library(readxl)
options(tibble.print_min = 5)
# read gtf.
gtf_file = read.table("/database/hg19/Gencode/hg19_genes.gtf",header = F,sep = "\t")
gene_name = gsub(pattern = "gene_id (.+?);.*",gtf_file$V9,replacement = "\\1",perl = T)
gtf_set = unique(gene_name)
length(gtf_set) # 23368

# read tcga tissue marker.
tcga_tissue_marker = read.table("../tissue_spec_gene/top20_gene_unique.txt",header = F)
# replace "." to "-"
tcga_tissue_marker = gsub("\\.","-",x = as.character(tcga_tissue_marker$V1))

# read tissue marker
cancer_maker_tissue = read_excel("cancer specific markers-refs.xlsx",sheet = 1,col_names = T)
cancer_maker_tissue = mutate(cancer_maker_tissue,tissue_count = apply(cancer_maker_tissue,1,function(x) sum(!is.na(x)))-1) %>%
  select(GENE,tissue_count)
  
# read plasma marker
cancer_maker_plasma = read_excel("cancer specific markers-refs.xlsx",sheet = 2,col_names = T)
cancer_maker_plasma = mutate(cancer_maker_plasma,plasma_count = apply(cancer_maker_plasma,1,function(x) sum(!is.na(x)))-1) %>%
  select(gene,plasma_count) %>%
  rename(GENE = gene)

# join in two table.
cancer_marker = full_join(cancer_maker_tissue,cancer_maker_plasma,by=c("GENE")) %>%
  mutate(gene_in_gtf = GENE %in% gtf_set) %>%
  mutate(gene_in_tcga = GENE %in% tcga_tissue_marker) %>%
  arrange(-gene_in_gtf,-plasma_count,-tissue_count)

tcga_tissue_marker = data.frame(GENE = as.character(tcga_tissue_marker),
                                gene_in_gtf = rep(TRUE,length(tcga_tissue_marker)),
                                gene_in_tcga = rep(TRUE,length(tcga_tissue_marker)),stringsAsFactors = F)

marker_full = full_join(cancer_marker,tcga_tissue_marker,by=c("GENE","gene_in_gtf","gene_in_tcga")) %>%
  filter(gene_in_gtf & (gene_in_tcga | (is.na(tissue_count) | tissue_count > 3))) %>%
  select(GENE) %>%
  unlist() %>%
  as.character()
# 521 gene.
write.table(marker_full,file = "candidate_gene.txt",quote = F,row.names = F,col.names = F)
marker_full = full_join(cancer_marker,tcga_tissue_marker,by=c("GENE","gene_in_gtf","gene_in_tcga"))
write_excel_csv(cancer_marker,path = "cancer_marker_tidy.csv")
