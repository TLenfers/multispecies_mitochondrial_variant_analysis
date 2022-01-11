# TODO: wildcards: [name]R1_001.fastq.gz
# TODO: doku: name for samples: sampleid_L001_R1 / R2
# TODO: doku: Referenz Datei immer .fa
# TODO: sampleid kommen in config
# TODO: geneierte datein: [sampleid].vcf

# TODO: Ziel Pipeline: fastqc, vcf (für ein sample), dann generisch mit wildcards) 

configfile: "config/config.yaml"

rule all:
    input:
        expand("calls/{sample}.vcf",sample=config["samples"])   

rule create_idx:
    input:
        ref=expand("data/reference/{reference}.fa", reference=config["reference"])
    output:
         "{input.ref}.amb", "{input.ref}.ann", "{input.ref}.bwt", "{input.ref}.pac", "{input.ref}.sa"
    conda:
        "workflow/envs/bwa.yaml"
    shell: 
        "bwa index {input.ref}"

rule bwa_mem:
    input:  
        ref=expand("data/reference/{reference}.fa", reference=config["reference"]),
        fastq_r1 = expand("data/samples/{sample}_L001_R1_001.fastq.gz",sample=config["samples"]),
        fastq_r2 = expand("data/samples/{sample}_L001_R2_001.fastq.gz",sample=config["samples"])
    output: 
        expand("results/mapped/{sample}.bam", sample=config["samples"])
    conda: 
        "workflow/envs/bwa.yaml"
    shell: "bwa mem {input.ref} {input.fastq_r1} {input.fastq_r2} | samtools view -b > {output}"


rule sort_alignments:
    input:
        #"results/mapped/{sample}.bam"
        expand("results/mapped/{sample}.bam", sample=config["samples"])
    output:
        #"results/mapped/{sample}.sorted.bam"
        expand("results/mapped/{sample}.sorted.bam", sample=config["samples"])
    conda:
        "envs/mapping.yaml"
    shell:
        "samtools sort -o {output} {input}"

# variant calling (bcftools)

rule call_variants:
    input:
        fa=expand("data/reference/{reference}.fa", reference=config["reference"]),
        bam=expand("results/mapped/{sample}.sorted.bam", sample=config["samples"])
    output:
        expand("results/calls/{sample}.vcf", sample=config["samples"])
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.fa} {input.bam} | bcftools call -mv - > {output}"


# variant calling (mutserve)

rule mutserve:
    input:
        bam=expand("results/mapped/{sample}.sorted.bam", sample=config["samples"]),
        #reference=expand("data/reference/{reference}.fa", reference=config["reference"]),
        reference="mutserve/rCRS.fasta.fai",
        index=expand("data/reference/{reference}.fa.fai", reference=config["reference"])
    output:
        expand("results/calls/{sample}.vcf",sample=config["samples"]),
        #expand("results/calls/{sample}.txt",sample=config["samples"]),
        #expand("results/calls/{sample}_raw.txt",sample=config["samples"])
    shell:
        "./mutserve/mutserve call --reference {input.reference} --output {output} {input.bam}"


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