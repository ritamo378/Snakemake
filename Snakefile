SAMPLES = ["A_R", "B_R", "C_R", "D_R", "E_R", "F_R", "G_R", "H_R"]
NR= ["1","2"]
GENOME = ["chr19_20Mb"]

rule fastqc:
    input:
        "reads/{sample}.fastq",
    output:
        html="qc/fastqc/{sample}.html",
        zip="qc/fastqc/{sample}_fastqc.zip"
    wrapper:
        "v1.5.0/bio/fastqc"

rule multiqc:
    input:
        expand("qc/fastqc/{sample}1_fastqc/fastqc_data.txt", sample=SAMPLES),
        expand("qc/fastqc/{sample}2_fastqc/fastqc_data.txt", sample=SAMPLES),
    output:
        "qc/multiqc.html"
    log:
        "logs/multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

rule bbduk_se:
    input:
        sample=["reads/{sample}1.fastq","reads/{sample}2.fastq"], 
        adapters = "reads/adapters.fa",
    output:
        trimmed=["trimmed/se/{sample}1.fastq.gz","trimmed/se/{sample}2.fastq.gz"],
    params:
        extra = lambda w, input: "ref={},adapters, ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10".format(input.adapters),
    threads: 7
    wrapper:
        "v1.5.0/bio/bbtools/bbduk"

rule trimmed_fastqc:
    input:
        "trimmed/se/{sample}.fastq",
    output:
        html="trimmed/qc/{sample}.html",
        zip="trimmed/qc/{sample}_fastqc.zip",
    wrapper:
        "v1.5.0/bio/fastqc"

rule trimmed_multiqc:
    input:
        expand("trimmed/qc/{sample}1_fastqc/fastqc_data.txt", sample=SAMPLES),
        expand("trimmed/qc/{sample}2_fastqc/fastqc_data.txt", sample=SAMPLES),
    output:
        "trimmed/qc/multiqc.html"
    log:
        "trimmed/logs/multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"

rule bbduk_se_polya:
    input:
        sample=["reads/{sample}1.fastq","reads/{sample}2.fastq"],
        adapters = "reads/adapters1.fa",
    output:
        trimmed=["trimmed/polya/{sample}1.fastq.gz","trimmed/polya/{sample}2.fastq.gz"],
    params:
        extra = lambda w, input: "ref={},adapters, ktrim=r k=23 mink=11 hdist=1 tpe tbo qtrim=r trimq=10".format(input.adapters),
    threads: 7
    wrapper:
        "v1.5.0/bio/bbtools/bbduk"

rule trimmed_fastqc_polya:
    input:
        "trimmed/polya/{sample}.fastq",
    output:
        html="trimmed/polya/qc/{sample}.html",
        zip="trimmed/polya/qc/{sample}_fastqc.zip",
    wrapper:
        "v1.5.0/bio/fastqc"

rule trimmed_multiqc_polya:
    input:
        expand("trimmed/polya/qc/{sample}1_fastqc/fastqc_data.txt", sample=SAMPLES),
        expand("trimmed/polya/qc/{sample}2_fastqc/fastqc_data.txt", sample=SAMPLES),
    output:
        "trimmed/polya/qc/multiqc.html"
    log:
        "trimmed/logs/multiqc.log"
    wrapper:
        "v1.5.0/bio/multiqc"


rule star_index:
    input:
        fasta=expand("annotation/{genome}.fa", genome=GENOME),
        gtf=expand("annotation/{genome}.gtf", genome=GENOME),
    output:
        directory(expand("{genome}", genome = GENOME)),
    message:
        "Testing STAR index"
    threads: 1
    params:
        extra="",
    wrapper:
        "v1.5.0/bio/star/index"
        
rule star_se:
    input:
        fq1="reads/{sample}1.fastq",
        fq2="reads/{sample}2.fastq",
        idx="chr19_20Mb",
    output:
        bam="star/se/{sample}/Aligned.sortedByCoord.out.bam",
    threads: 8
    shell:
        "STAR --runThreadN 2 --genomeDir {input.idx} --readFilesIn {input.fq1} {input.fq2} --outStd BAM_SortedByCoordinate --outSAMtype BAM SortedByCoordinate --outFileNamePrefix star/se/{wildcards.sample}/"
        
rule samtools_index:
    input:
        "{sample}_Aligned.sortedByCoord.out.bam",
    output:
        "samtools/{sample}.sorted.bam.bai",
    threads: 4 
    wrapper:
        "v1.5.0/bio/samtools/index"        
     

rule feature_counts1:
    input:
        ant="annotation/chr19_20Mb.gtf",
        bam=expand("{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
    output:
        out="counts1.txt",
    shell:
        "featureCounts -p -t exon -g gene_id -a {input.ant} -o {output.out} {input.bam} -s 1"

rule feature_counts2:
    input:
        ant="annotation/chr19_20Mb.gtf",
        bam=expand("{sample}_Aligned.sortedByCoord.out.bam", sample=SAMPLES),
    output:
        out="counts2.txt"
    shell:
        "featureCounts -p -t exon -g gene_id -a {input.ant} -o {output.out} {input.bam} -s 2"

