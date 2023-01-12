# ArabidopsisRNAseqExample

利用拟南芥的数据做转录组数据处理的一个例子，基本的流程是 
- hisat2比对 
- samtools sam文件转bam文件
- stringtie计算表达量

这三个步骤用snakemake串在了一起

- R语言tximport包获取基因和转录本水平的counts
- R语言DESeq2做差异表达分析


## 首先是下载原始测序数据

我这里是从ENA数据库下载的，写了一个简单的python脚本

```
mkdir 00.raw.fq
cd 00.raw.fq
python download_fastq.py
```
`download_fastq.py`脚本里的内容

```
import subprocess

links = [
    "ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_1.fastq.gz",
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/003/SRR4420293/SRR4420293_2.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_1.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/004/SRR4420294/SRR4420294_2.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/005/SRR4420295/SRR4420295_1.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/005/SRR4420295/SRR4420295_2.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/006/SRR4420296/SRR4420296_1.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/006/SRR4420296/SRR4420296_2.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_1.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/007/SRR4420297/SRR4420297_2.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_1.fastq.gz',
    'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR442/008/SRR4420298/SRR4420298_2.fastq.gz']

for link in links:
    cmd = ['wget',link]
    print(' '.join(cmd))
    subprocess.run(cmd)
```

## 下载参考基因组和gff格式注释问津

```
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_chromosome_files/TAIR10_chr_all.fas.gz
dtrx TAIR10_chr_all.fas.gz
wget https://www.arabidopsis.org/download_files/Genes/TAIR10_genome_release/TAIR10_gff3/TAIR10_GFF3_genes.gff
gffread TAIR10_GFF3_genes.gff -T -o at.gtf
```

## 安装软件

我这里直接使用conda进行安装

```
conda install snakemake
conda install hisat2
conda install samtools
conda install stringtie
```
## 运行snakemake流程

```
snakemake -s 002_histat2_samtools_stringtie_rnaseq.smk --cores 8 -p
```

`002_histat2_samtools_stringtie_rnaseq.smk`文件的内容

```
SRR,FRR = glob_wildcards("00.raw.fq/{srr}_{frr}.fastq.gz")

print(SRR)

rule all:
    input:
        "reference/index/at.1.ht2",
        expand("01.sam/{srr}.sam",srr=SRR),
        expand("02.sorted.bam/{srr}.sorted.bam",srr=SRR),
        expand("02.sorted.bam/{srr}.sorted.bam.bai",srr=SRR),
        expand("03.expression/{srr}/gene_abund.tsv",srr=SRR),
        expand("03.expression/{srr}/transcripts.gtf",srr=SRR)

rule hisat2_build:
    input:
        "reference/TAIR10_chr_all.fas",
    output:
        "reference/index/at.1.ht2"
    params:
        "reference/index/at"
    threads:
        8
    resources:
        mem = 16000
    shell:
        """
        hisat2-build {input} {params}
        """
rule hista2_align:
    input:
        read01 = "00.raw.fq/{srr}_1.fastq.gz",
        read02 = "00.raw.fq/{srr}_2.fastq.gz",
        index = rules.hisat2_build.output
    output:
        sam = "01.sam/{srr}.sam"
    params:
        "reference/index/at"
    threads:
        4
    resources:
        mem = 8000
    shell:
        """
        hisat2 -p {threads} --dta -x {params} -1 {input.read01} -2 {input.read02} -S {output.sam}
        """

rule samtools_sort:
    input:
        rules.hista2_align.output.sam
    output:
        "02.sorted.bam/{srr}.sorted.bam"
    threads:
        4
    resources:
        mem = 8000
    shell:
        """
        samtools sort -@ {threads} -o {output} {input}
        """

rule samtools_index:
    input:
        rules.samtools_sort.output
    output:
        "02.sorted.bam/{srr}.sorted.bam.bai"
    threads:
        4
    resources:
        mem = 8000
    shell:
        """
        samtools index -@ {threads} {input}
        """
rule stringtie:
    input:
        gtf = "reference/at.gtf",
        bam = rules.samtools_sort.output,
        bai = rules.samtools_index.output
    output:
        trans_gtf = "03.expression/{srr}/transcripts.gtf",
        gene_abund = "03.expression/{srr}/gene_abund.tsv"
    threads:
        4
    resources:
        mem = 8000
    shell:
        """
        stringtie -p {threads} -G {input.gtf} -e -B -o {output.trans_gtf} -A {output.gene_abund} {input.bam}
        """
```

## 减小数据量
提取bam文件中比对到Chr1的双端数据

python代码
```
import pysam
import argparse

# read01 = {'read_id':[],
#           'read_seq':[],
#           'read_plus':[],
#           'read_quality':[]}
# read02 = {'read_id':[],
#           'read_seq':[],
#           'read_plus':[],
#           'read_quality':[]}




read01 = {}
read02 = {}

parser = argparse.ArgumentParser(
    formatter_class=argparse.RawDescriptionHelpFormatter,
    description='get paired-end reads from bam file',
    epilog="""
    @author: MingYan
    """
)

parser.add_argument('-b','--bam',required=True)
parser.add_argument('-c','--chr-name',required=True)
parser.add_argument('-s','--start-site',required=True,type=int)
parser.add_argument('-e','--end-site',required=True,type=int)
parser.add_argument('-r1','--read01-name',required=True)
parser.add_argument('-r2','--read02-name',required=True)

args = parser.parse_args()

bam_file = args.bam
chr_name = args.chr_name
start_site = args.start_site
end_site = args.end_site
read01_name = args.read01_name
read02_name = args.read02_name




for read in pysam.AlignmentFile(bam_file).fetch(chr_name,start_site,end_site):
    if read.is_read1 and read.is_proper_pair:
        read01[read.to_dict()['name']] = []
        read01[read.to_dict()['name']].append(read.to_dict()['seq'])
        read01[read.to_dict()['name']].append("+")
        read01[read.to_dict()['name']].append(read.to_dict()['qual'])
    if read.is_read2 and read.is_proper_pair:
        read02[read.to_dict()['name']] = []
        read02[read.to_dict()['name']].append(read.to_dict()['seq'])
        read02[read.to_dict()['name']].append("+")
        read02[read.to_dict()['name']].append(read.to_dict()['qual'])
        
fw01 = open(read01_name,'w')
fw02 = open(read02_name,'w')
        
for read_name in read01.keys():
    if read_name in read02.keys():
        fw01.write("@%s\n%s\n+\n%s\n"%(read_name,read01[read_name][0],read01[read_name][2]))
        fw02.write("@%s\n%s\n+\n%s\n"%(read_name,read02[read_name][0],read02[read_name][2]))
        
        
        
fw01.close()
fw02.close()
```

```
python get_reads_from_bam.py -b bam -c chr_name -s 1 -e 100000000 -r1 output_R1.fastq -r2 output_R2.fastq
```

得到每个fastq进行压缩后文件大小大约在20M左右

```
bgzip output_R1.fastq
bgzip output_R2.fastq
```

## seqkit 根据 id 提取fasta

```
seqkit faidx input.fasta Chr1 -r
```

## gff注释文件也只保留Chr染色体的内容

```
grep 'Chr1' input.gff > output.gff
```

## 拟南芥有参转录组 第一步：需要准备的文件

通常我们自己做转录组的实验，都是自己准备好样品，然后交给测序公司去测序，如果只是要求测序，并不要求测序公司做相关分析的话，测序公司返回给我们的就是fastq文件，现在通常都是双端测序，所以一个样本对应着是两个fastq文件，假设我们一个对照，一个处理，每个处理测3次重复，所以总共是6个样本，那我们最终会在测序公司手里拿到12个fastq文件

fastq 也是文本文件，文件里的内容按照一定的规律排列，这个规律称为fastq,fastq文件里的内容如下


```
@SRR2037320.R.2/1
TGAACTTCGAGGAATAGCAGAGACTCCGGAGCTGAAGAGACACTTTGCACGCTGAGAATCTGCAGACAGTCGTCAGGATCCAGAGATAGAATGCGCACAGT
+
@;@DBBBAHHFCDCAHBDHBEBFHIIGEHCF>CFC0?9DF3?BBDGEHGGFFBGIIHFHEHHFFFFBACCE??=B8=C@?A>ACB?CC>A>:>A@B<B>?#
@SRR2037320.R.5/1
CCGAAAGTGTACTGTATGCTGCACTGTGTGAATCACCCCAACGTCCATGTGGGAAGTCATTCCAGCCCTTTGTTCTAGATATGTAGCTCTTTGGACGGTTT
+
CCCFFFFDHHHHGJIIJJJJJJJJJJJJJIJJJJJJJJJJJJJJJJJJJHJJJJJJJIIJJIIJJJJJHHHHGFF@EDFFEEDEEDEDDDDDDDDCBDDDD
@SRR2037320.R.7/1
TTTCCATAACGCAAGAGGATAGATTTAGTATTACAAACTGAATCAATGTGCTTGACGAGGTAGACATGTTTACTAGGCAAGTACAAGAAGAAACATTTGAT
+
;??D:?DD?+CDD:1+<<+A++<CE??4C?E4*:*?CD@@DDEICDE@4BBA@?(?;B'-'.).=@=C77=CDDD);)7.;;A>;5>>,5;>=A9>A:A##
```

以四行为一个单位


fastq文件通常比较大，所以大部分情况下都是以压缩文件的形式存储

这里介绍的是有参转录组，也就是说自己的目标物种已经有了对应的参考基因组，拿到fastq文件后还需要我们自己到网上下载参考基因组和对应的注释文件，这个参考基因组具体是到哪里下载可以找到参考基因组对应的论文，通常论文里会写他的数据存储到了哪里。

大概率NCBI可以下载到参考基因组，有的物种可能研究的人比较多，自己会有专门的基因组网站，比如拟南芥 西红柿 还有蔷薇科的果实

- 1 参考基因组 fasta 格式

这是一个文本文件，文本文件的意思就是我们用电脑上的记事本软件打开这个文件可以看到里面的内容，fasta 是指文本文件的内容按照指定的格式来排列

```
>seqid_1
ATCGATCGATCG
ATCGATCG
>seqid_2
ATCGAATTTTAA
ATCGATCCC
```

- 2 参考基因组的注释文件 gff

```
##gff-version 2
##created 11/11/11 
Chr1	TAIR10	chromosome	1	30427671	.	.	.	ID=Chr1;Name=Chr1
Chr1	TAIR10	gene	3631	5899	.	+	.	ID=AT1G01010;Note=protein_coding_gene;Name=AT1G01010
Chr1	TAIR10	mRNA	3631	5899	.	+	.	ID=AT1G01010.1;Parent=AT1G01010;Name=AT1G01010.1;Index=1
Chr1	TAIR10	protein	3760	5630	.	+	.	ID=AT1G01010.1-Protein;Name=AT1G01010.1;Derives_from=AT1G01010.1
Chr1	TAIR10	exon	3631	3913	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	five_prime_UTR	3631	3759	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	3760	3913	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	3996	4276	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	3996	4276	.	+	2	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	4486	4605	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	4486	4605	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	4706	5095	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	4706	5095	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	5174	5326	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	5174	5326	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	exon	5439	5899	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	CDS	5439	5630	.	+	0	Parent=AT1G01010.1,AT1G01010.1-Protein;
Chr1	TAIR10	three_prime_UTR	5631	5899	.	+	.	Parent=AT1G01010.1
Chr1	TAIR10	gene	5928	8737	.	-	.	ID=AT1G01020;Note=protein_coding_gene;Name=AT1G01020
Chr1	TAIR10	mRNA	5928	8737	.	-	.	ID=AT1G01020.1;Parent=AT1G01020;Name=AT1G01020.1;Index=1
Chr1	TAIR10	protein	6915	8666	.	-	.	ID=AT1G01020.1-Protein;Name=AT1G01020.1;Derives_from=AT1G01020.1
Chr1	TAIR10	five_prime_UTR	8667	8737	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	8571	8666	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	8571	8737	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	8417	8464	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	8417	8464	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	8236	8325	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	8236	8325	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	7942	7987	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	7942	7987	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	7762	7835	.	-	2	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	7762	7835	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	7564	7649	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	7564	7649	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	7384	7450	.	-	1	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	7384	7450	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	7157	7232	.	-	0	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	exon	7157	7232	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	CDS	6915	7069	.	-	2	Parent=AT1G01020.1,AT1G01020.1-Protein;
Chr1	TAIR10	three_prime_UTR	6437	6914	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	exon	6437	7069	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	three_prime_UTR	5928	6263	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	exon	5928	6263	.	-	.	Parent=AT1G01020.1
Chr1	TAIR10	mRNA	6790	8737	.	-	.	ID=AT1G01020.2;Parent=AT1G01020;Name=AT1G01020.2;Index=1
Chr1	TAIR10	protein	7315	8666	.	-	.	ID=AT1G01020.2-Protein;Name=AT1G01020.2;Derives_from=AT1G01020.2
Chr1	TAIR10	five_prime_UTR	8667	8737	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	8571	8666	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	exon	8571	8737	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	8417	8464	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	exon	8417	8464	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	8236	8325	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	exon	8236	8325	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	7942	7987	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	exon	7942	7987	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	7762	7835	.	-	2	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	exon	7762	7835	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	7564	7649	.	-	0	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	exon	7564	7649	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	CDS	7315	7450	.	-	1	Parent=AT1G01020.2,AT1G01020.2-Protein;
Chr1	TAIR10	three_prime_UTR	7157	7314	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	exon	7157	7450	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	three_prime_UTR	6790	7069	.	-	.	Parent=AT1G01020.2
Chr1	TAIR10	exon	6790	7069	.	-	.	Parent=AT1G01020.2
```


数据准备好了，接下来开始分析，因为转录组的数据相对还是比较大的，

分析数据需要用到linux操作系统，因为大部分的生物信息相关软件都是linux系统下的命令行工具 ，所以linux操作系统是必须掌握的一个技能

## 简单的linux系统操作命令

windows系统的功能通常都是可以看见的，需要我们用鼠标去点对应的按钮，比如新建文件夹 在windows系统我们需要点击鼠标右键-选择新建文件夹-给文件夹命名

linux系统我们是看不见要做的事情对应的按钮的，需要我们记住一些命令 用这些命令告诉电脑我们要干什么，同样是新建一个文件夹，在linux操作系统首先要想新建文件夹的命令是什么，mkdir 然后加新建文件夹的名字 `mkdir xiaoming`
新建一个文本文件 `touch xiaoming.txt`

在文本文件里写内容 vim xiaoming.txt

查看文本文件的内容 less -S xiaoming.txt

当然这些命令很多，我们不太可能所有都记得住，我们需要做的是了解linux系统基本的使用逻辑，然后我们想做某件事不知道用什么命令来实现，我们就直接去搜索引擎里搜索。比如前面已经新建了一个文件夹，那我不想要要这个文件夹了，我要删除他，不知道用什么命令，可以用关键词 `linux delete folder`去搜索

这里只简单介绍几个基本命令，linux操作系统四一个多用，慢慢熟悉的过程

## linux安装软件

windows系统安装软件需要我们去找到对应的安装包，然后用鼠标一路点点点就可以了

linux系统安装软件同样是需要借助命令 linux有很多安装软件的方式 作为生物信息学入门安装软件 我们只需要掌握一个工具conda就可以了，掌握了这个工具基本上90%的生物信息学相关的软件就都可以安装了

下载anaconda3

## conda 添加北京外国语大学镜像

https://zhuanlan.zhihu.com/p/386843337
https://mirrors.bfsu.edu.cn/help/anaconda/

```
vi ~/.condarc
show_channel_urls: true
channels:
  - https://mirrors.bfsu.edu.cn/anaconda/cloud/bioconda/
  - https://mirrors.bfsu.edu.cn/anaconda/cloud/pytorch/
  - https://mirrors.bfsu.edu.cn/anaconda/cloud/msys2/
  - https://mirrors.bfsu.edu.cn/anaconda/cloud/conda-forge
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/main/
  - https://mirrors.bfsu.edu.cn/anaconda/pkgs/free/
  - defaults
```

新建虚拟环境，将软件安装到各自的虚拟环境里

每次需要先启动虚拟环境

安装软件

- fastqc
- fastp
- hisat2
- samtools
- stringtie

以上是在linux系统下完成，接下来用到R语言就可以在我们自己电脑上进行了

- R
- R package tximport

fastqc 是用来统计测序质量

fastp 是用来对测序数据进行过滤

hisat2是用来比对fastq到参考基因组

samtools是用来做bam sam 文件转化 这里又提到了两个文件格式 后续步骤再来介绍

stringtie 是用来计算表达量

## 下面开始实际操作

## 用命令行将服务器端的文件下载到本地

```
scp xiaoming@124.70.145.183:/tmp/dir_Gl_VKU_cp/Gl_VKU_cp.tar.gz D:\Bioinformatics_Intro\
```
