#BiocManager::install("tximport")
#https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html

library(tximport)
library(tidyverse)

dir.create("04.expCounts")

sample_name <- list.dirs("03.stringtie",recursive=FALSE,full.names=FALSE)

files <- paste("03.stringtie",sample_name,"t_data.ctab",sep = "/")
names(files)<-sample_name

tmp <- read_tsv(files[1])

tx2gene <- tmp[, c("t_name", "gene_id")]
tx2gene
txi.genes <- tximport(files, type = "stringtie", tx2gene = tx2gene)

txi.genes$counts %>% as.data.frame() %>% 
  mutate(geneLength=txi.genes$length[,1]) %>% 
  round() %>% 
  rownames_to_column("gene_name") %>% 
  write_csv("04.expCounts/gene_counts.csv")

txi.trans<-tximport(files,type = "stringtie",txOut = TRUE)
txi.trans$counts %>% 
  as.data.frame() %>% 
  mutate(geneLength=txi.trans$length[,1]) %>% 
  round() %>% 
  rownames_to_column("trans_id") %>% 
  write_csv("04.expCounts/trans_counts.csv")

