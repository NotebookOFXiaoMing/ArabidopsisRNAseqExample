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

## seqkit 根据 id 提取fasta

```
seqkit faidx input.fasta Chr1 -r
```
