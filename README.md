# ArabidopsisRNAseqExample

利用拟南芥的数据做转录组数据处理的一个例子，基本的流程是 
- hisat2比对 
- samtools sam文件转bam文件
- stringtie计算表达量

这三个步骤用snakemake串在了一起

- R语言tximport包获取基因和转录本水平的counts
- R语言DESeq2做差异表达分析
