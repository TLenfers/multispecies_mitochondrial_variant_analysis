# mapping of the reads to the reference genome using Burrows-Wheeler Aligner
rule bwa_mem:
    input:
        amb="data/reference/{reference}.fa.amb",
        ref="data/reference/{reference}.fa",
        fastq_r1="data/samples/{sample}_R1.fastq.gz",
        fastq_r2="data/samples/{sample}_R2.fastq.gz",
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    output:
        "results/mapped/{reference}/{sample}.bam",
    conda:
        "workflow/envs/bwa.yaml"
    log:
        "logs/{reference}/{sample}.bwa.log",
    shell:
        "bwa mem {input.ref} {input.fastq_r1} {input.fastq_r2} | samtools sort| samtools view -b > {output}"


# creates a bam index
rule idx_bam:
    input:
        "results/mapped/{reference}/{sample}.bam",
    output:
        "results/mapped/{reference}/{sample}.bam.bai",
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    conda:
        "workflow/envs/samtools.yaml"
    log:
        "logs/{reference}/{sample}.bam.bai.log",
    shell:
        "samtools index {input}"


# creates index for the reference sequence in fasta format
rule idx_fasta:
    input:
        "data/reference/{reference}.fa",
    output:
        "data/reference/{reference}.fa.fai",
    conda:
        "workflow/envs/samtools.yaml"
    log:
        "logs/{reference}/idx_fasta.log",
    shell:
        "samtools faidx {input}"
