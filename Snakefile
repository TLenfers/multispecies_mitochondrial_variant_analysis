rule create_idx:
    input: "data/reference/{reference}.fasta"
    output: "data/reference/{reference}.fasta.amb", "data/{reference}.fasta.ann", "data/{reference}.fasta.bwt", "data/{reference}.fasta.pac", "data/{reference}.fasta.sa"
    conda: "workflow/envs/bwa.yaml"
    shell: "bwa index {input}"

rule map_reads:
    input:
        "data/reference/{reference}.fasta",
        "data/samples/{sample}.fastq.gz"
    output:
        "results/mapped/{sample}.bam"
    conda:
        "workflow/envs/mapping.yaml"
    shell:
        "bwa mem {input} | samtools view -Sb -> {output}"