#setwd("D:/Bioinformatics_Intro/RNAseq/Arabidopsis_RNAseq/")

library(tidyverse)
list.files("03.expression",
           pattern="t_data.ctab",
           recursive = TRUE,
           full.names = TRUE) -> files
names(files)<- str_extract(files,"SRR[0-9]+")
files
#BiocManager::install("tximport")
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html
library(tximport)
tmp <- read_tsv(files[1])
tmp %>% head() %>% DT::datatable()
tx2gene <- tmp[, c("t_name", "gene_id")]
tx2gene
txi.genes <- tximport(files, type = "stringtie", tx2gene = tx2gene)

txi.genes$counts %>% as.data.frame() %>% 
  mutate(geneLength=txi.genes$length[,1]) %>% 
  round() %>% 
  rownames_to_column("gene_name") %>% 
  write_csv("gene_counts.csv")

txi.trans<-tximport(files,type = "stringtie",txOut = TRUE)
txi.trans$counts %>% 
  as.data.frame() %>% 
  mutate(geneLength=txi.trans$length[,1]) %>% 
  round() %>% 
  rownames_to_column("trans_id") %>% 
  write_csv("trans_counts.csv")

