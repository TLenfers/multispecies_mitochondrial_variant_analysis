# TODO: wildcards: [name]R1_001.fastq.gz
# TODO: doku: name for samples: sampleid_L001_R1 / R2
# TODO: doku: Referenz Datei immer .fa
# TODO: sampleid kommen in config
# TODO: geneierte datein: [sampleid].vcf

# TODO: Ziel Pipeline: fastqc, vcf (für ein sample), dann generisch mit wildcards) 

configfile: "config/config.yaml"

rule all:
    input:
        expand("results/calls_bcftools/{reference}/{sample}.vcf",sample=config["samples"],reference=config["reference"])   

rule get_ref:
    input:
        "data/reference/{reference}.fa" 
    output:
        "data/reference/{reference}.fa"
    run:
        if wildcards.{reference} == "mouse":
            shell("curl http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.MT.fa.gz | zcat > {output}")

rule bwa_idx:
    input:
        "data/reference/{reference}.fa" 
    output:
        "data/reference/{reference}.fa.amb", 
        "data/reference/{reference}.fa.ann", 
        "data/reference/{reference}.fa.bwt", 
        "data/reference/{reference}.fa.pac", 
        "data/reference/{reference}.fa.sa" 
    conda:
        "workflow/envs/bwa.yaml"
    shell: 
        "bwa index {input}"

rule bwa_mem:
    input:  
        amb="data/reference/{reference}.fa.amb",
        ref="data/reference/{reference}.fa",
        fastq_r1 = "data/samples/{sample}_L001_R1_001.fastq.gz",
        fastq_r2 = "data/samples/{sample}_L001_R2_001.fastq.gz"
    output: 
        "results/mapped/{reference}/{sample}.bam"
    conda: 
        "workflow/envs/bwa.yaml"
    shell: "bwa mem {input.ref} {input.fastq_r1} {input.fastq_r2} | samtools sort| samtools view -b > {output}"


rule idx_bam:
    input:
        "results/mapped/{reference}/{sample}.bam"
    output:
        "results/mapped/{reference}/{sample}.bam.bai"
    conda:
        "workflow/envs/bwa.yaml"
    shell:
        "samtools index {input}"

rule idx_fasta:
    input:
        "data/reference/{reference}.fa"
    output:
        "data/reference/{reference}.fa.fai"
    conda:
        "workflow/envs/samtools.yaml"
    shell:
        "samtools faidx {input}"

# variant calling (bcftools)

#rule call_variants:
#    input:
#        fa=expand("data/reference/{reference}.fa", reference=config["reference"]),
#        bam=expand("results/mapped/{sample}.sorted.bam", sample=config["samples"])
#    output:
#        expand("results/calls/{sample}.vcf", sample=config["samples"])
#    conda:
#        "workflow/envs/bcftools.yaml"
#    shell:
#        "bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output}"

rule call_variants:
    input:
        bam="results/mapped/{reference}/{sample}.bam",
        bamidx = "results/mapped/{reference}/{sample}.bam.bai",
        ref="data/reference/{reference}.fa",
        index="data/reference/{reference}.fa.fai"
    output:
        "results/calls_bcftools/{reference}/{sample}.vcf"
    wildcard_constraints: 
        reference="[a-z]+"
    shell:
        "bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv --ploidy 1 -> {output}"

# variant calling (mutserve)

rule mutserve_all:
    input:
        expand("results/calls/{sample}.vcf",sample=config["samples"])

rule mutserve:
    input:
        bam="results/mapped/{sample}.bam",
        bamidx = "results/mapped/{sample}.bam.bai",
        reference="data/reference/{reference}.fa",
        index="data/reference/{reference}.fa.fai"
    output:
        "results/calls_mutserve/{reference}/{sample}.vcf"
    shell:
        "./mutserve/mutserve call --reference {input.reference} --output {output} {input.bam}"
        #"java -jar ./mutserve/mutserve.jar call --reference {input.reference} --output {output} {input.bam}"


#rule fastqc_german:
#    input: "data/raw_german/{sample}_L001_{mate}_001.fastq.gz"
#    output: "fastqc_german/{sample}_L001_{mate}_001_fastqc.html",
#            "fastqc_german/{sample}_L001_{mate}_001_fastqc/summary.txt"
#    conda: "envs/fastqc.yaml"
#    shell: "fastqc --extract --outdir=fastqc_german/ {input}"
#
#rule fastqc_german_all:
#    input: expand("fastqc_german/{sample}_L001_{mate}_001_fastqc.html",sample=INDIVIDUALS_GERMAN2021,mate=["R1","R2"])
#
#rule fastqc_german_summary:
#    input: expand("fastqc_german/{sample}_L001_{mate}_001_fastqc/summary.txt",sample=INDIVIDUALS_GERMAN2021,mate=["R1","R2"])
#    output: "fastqc_german/fastqc_summary.txt"
#    run:
#        with open(output[0],"w") as f_out:
#            header = ["","Basic Statistics","Per base sequence quality",\
#            "Per tile sequence quality","Per sequence quality scores", \
#            "Per base sequence content","Per sequence GC content", \
#            "Per base N content","Sequence Length Distribution", \
#            "Sequence Duplication Levels","Overrepresented sequences", \
#            "Adapter Content"]
#            f_out.write("\t".join(header)+"\n")
#            for filename in input:
#                f_out.write(filename.split("/")[1]+"\t")
#                with open(filename,"r") as f_in:
#                    i = 0
#                    for line in f_in:
#                        if i<10:
#                            f_out.write(line.split("\t")[0]+"\t")
#                            i += 1
#                        else:
#                            f_out.write(line.split("\t")[0]+"\n")