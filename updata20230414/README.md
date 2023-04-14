### 第一步是安装conda

### 第二步是用conda安装需要用到的软件

```
conda install snakemake
conda install hisat2
conda install samtools
conda install stringtie
conda install r
conda install r-tidyverse
conda install bioconductor-tximport
```

### 需要准备的数据

00.raw.fq文件夹下存储的是原始的转录组测序数据 
reference/genome/ref.fa 是参考基因组，（把你的参考基因组改成ref.fa这个名字）
reference/gtf/ref.gff 是参考价基因组的注释文件 （把你的注释文件也改成ref.gff这个名字，gtf和gff都可以，统一改成ref.gff）

软件安装好，数据贮备好，运行

```
snakemake -s histat2_samtools_stringtie_rnaseq.smk --cores 12 -p
```

然后再运行

```
 Rscript getCounts.R
```

怎么把最后一步也整合到snakemake流程里，暂时遇到了问题，不知道如何解决了

最后生成4个文件夹

- 01.sam  
- 02.sorted.bam 
- 03.stringtie  
- 04.expCounts

04.expCounts 存储的是基因和转录本的count值和基因长度信息，可以用DESeq2做差异表达分析，
