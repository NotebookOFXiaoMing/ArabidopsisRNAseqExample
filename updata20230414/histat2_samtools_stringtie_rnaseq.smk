SRR,FRR = glob_wildcards("00.raw.fq/{srr}_{frr}.fastq.gz")

print(SRR)

rule all:
    input:
        "reference/index/ref.1.ht2",
        expand("01.sam/{srr}.sam",srr=SRR),
        expand("02.sorted.bam/{srr}.sorted.bam",srr=SRR),
        expand("02.sorted.bam/{srr}.sorted.bam.bai",srr=SRR),
        expand("03.stringtie/{srr}/gene_abund.tsv",srr=SRR),
        expand("03.stringtie/{srr}/transcripts.gtf",srr=SRR),
        # "04.expCounts/trans_counts.csv",
        # "04.expCounts/gene_counts.csv"

rule hisat2_build:
    input:
        "reference/genome/ref.fa",
    output:
        "reference/index/ref.1.ht2"
    params:
        "reference/index/ref"
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
        read01 = "00.raw.fq/{srr}_R1.fastq.gz",
        read02 = "00.raw.fq/{srr}_R2.fastq.gz",
        index = rules.hisat2_build.output
    output:
        sam = "01.sam/{srr}.sam"
    params:
        "reference/index/ref"
    threads:
        4
    resources:
        mem = 24000
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
        mem = 24000
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
        mem = 24000
    shell:
        """
        samtools index -@ {threads} {input}
        """
rule stringtie:
    input:
        gtf = "reference/gtf/ref.gff",
        bam = rules.samtools_sort.output,
        bai = rules.samtools_index.output
    output:
        trans_gtf = "03.stringtie/{srr}/transcripts.gtf",
        gene_abund = "03.stringtie/{srr}/gene_abund.tsv",
        t_data_ctab = "03.stringtie/{srr}/t_data.ctab",
        #log = "stringtie.log"
    threads:
        4
    resources:
        mem = 8000
    params:
        folder = "03.stringtie"
    shell:
        """
        stringtie -p {threads} -G {input.gtf} -e -B -o {output.trans_gtf} -A {output.gene_abund} {input.bam}
        """

# rule RgetCounts:
#     input:
#         gene_abund = rules.stringtie.output.gene_abund,
#         t_data_ctab = rules.stringtie.output.t_data_ctab,
#         trans_gtf = rules.stringtie.output.trans_gtf
#     output:
#         gene = "04.expCounts/gene_counts.csv",
#         trans = "04.expCounts/trans_counts.csv"
#     threads:
#         2
#     resources:
#         mem = 8000
#     params:
#         folder = rules.stringtie.params.folder
#     script:
#         "scripts/r/getCounts.R"

