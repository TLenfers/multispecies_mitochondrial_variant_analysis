configfile: "config/config.yaml"

rule all:
    input:
        expand("results/calls_bcftools/{reference}/norm_{sample}.vcf.gz",sample=config["samples"],reference=config["reference"]),
        expand("results/calls_bcftools/{reference}/merged.vcf",reference=config["reference"]),
        expand("results/calls_bcftools/{reference}/commonVariants.tsv",reference=config["reference"]),
        expand("results/plots/{reference}/ref_heatmap.pdf",reference=config["reference"]),
        expand("results/plots/{reference}/ref_heatmap_clusterrow.pdf",reference=config["reference"]),
        expand("results/plots/{reference}/alt_heatmap.pdf",reference=config["reference"]),
        expand("results/plots/{reference}/alt_heatmap_clusterrow.pdf",reference=config["reference"]),
        expand("results/sequences/{reference}/{sample}.fa",reference=config["reference"],sample=config["samples"])
        #expand("results/calls_mutserve/{reference}/{sample}.vcf.gz",sample=config["samples"],reference=config["reference"])#,
        #expand("results/calls_mutserve/{reference}/merged.vcf",reference=config["reference"])

# TODO : add dog reference, fix download
#rule get_ref:
#    input:
#        "data/reference/{reference}.fa"
#    output:
#        "data/reference/{reference}.fa"
#    run:
#        if wildcards.{reference} == "mouse":
#            shell("curl http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.MT.fa.gz | zcat > {output}")


rule bwa_idx:
    input:
        "data/reference/{reference}.fa" 
    output:
        "data/reference/{reference}.fa.amb", 
        "data/reference/{reference}.fa.ann", 
        "data/reference/{reference}.fa.bwt", 
        "data/reference/{reference}.fa.pac", 
        "data/reference/{reference}.fa.sa" 
    wildcard_constraints:
        reference="[A-Za-z0-9]+"
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
    wildcard_constraints:
        reference="[A-Za-z0-9]+"
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
    wildcard_constraints:
        reference="[A-Za-z0-9]+"
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
rule call_variants:
    input:
        bam="results/mapped/{reference}/{sample}.bam",
        bamidx = "results/mapped/{reference}/{sample}.bam.bai",
        ref="data/reference/{reference}.fa",
        index="data/reference/{reference}.fa.fai"
    output:
        "results/calls_bcftools/{reference}/{sample}.vcf"
    wildcard_constraints:
        reference="[A-Za-z0-9]+"
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv --ploidy 1 -> {output}"

rule normalize_variants:
    input:
        vcf="results/calls_bcftools/{reference}/{sample}.vcf.gz",
        ref="data/reference/{reference}.fa"
    output:
        "results/calls_bcftools/{reference}/norm_{sample}.vcf"
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools norm  --fasta-ref {input.ref} --check-ref -m {input.vcf} | bcftools view -Ov -o {output}"


# variant calling (mutserve)

#rule mutserve_all:
#    input:
#        expand("results/calls/{sample}.vcf",sample=config["samples"])

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

##################################

# merge vcf
rule zip_vcf:
    input:
        "results/calls_bcftools/{reference}/{sample}.vcf"
    output:
        "results/calls_bcftools/{reference}/{sample}.vcf.gz"
    conda:
        "workflow/envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+"
    shell:
        "bgzip -f {input}; tabix -f -p vcf {output}"



rule merge_vcf:
    input:
        expand("results/calls_bcftools/{{reference}}/{sample}.vcf.gz",sample=config["samples"])
    output:
        "results/calls_bcftools/{reference}/merged.vcf"
    conda: "workflow/envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+"
    shell: "bcftools merge -m none -O v {input} > {output}"


rule common_variants:
    input:
        in_file= "results/calls_bcftools/{reference}/merged.vcf"
    output:
        "results/calls_bcftools/{reference}/commonVariants.tsv"
    conda:
        "workflow/envs/r-heatmap.yaml"
    script:
        "workflow/scripts/common_variants.R"


rule plot_variant_heatmap:
    input:
        "results/calls_bcftools/{reference}/commonVariants.tsv"
    output:
        "results/plots/{reference}/ref_heatmap.pdf",
        "results/plots/{reference}/ref_heatmap_clusterrow.pdf",
        "results/plots/{reference}/alt_heatmap.pdf",
        "results/plots/{reference}/alt_heatmap_clusterrow.pdf"
    conda:
        "workflow/envs/r-heatmap.yaml"
    script:
        "workflow/scripts/plot_variant_heatmap.R"

rule prepare_for_seq:
    input:
            vcf="results/calls_bcftools/{reference}/{sample}.vcf.gz"
    output:
        "results/sequences/{reference}/{sample}.vcf.gz",
        "results/sequences/{reference}/{sample}.vcf.gz.tbi"
    conda: "envs/vcftools.yaml"
    shell: "vcftools --gzvcf {input.vcf} --remove-indels --minQ 200.0 --recode --recode-INFO-all --stdout |bgzip > {output[0]}; tabix -p vcf {output[0]}"

# Replace '*', which is the gap symbol in mutserv VCF, with '-', the appropriate 
# gap symbol in fasta format
# Further, replace newlines not in the header line
# Further, replace the header line to update the 'sequence name'
rule vcf_to_fasta:
    input: 
           ref="data/reference/{reference}.fa" ,
            vcf="results/sequences/{reference}/{sample}.vcf.gz",
            index="results/sequences/{reference}/{sample}.vcf.gz.tbi"
    output: "results/sequences/{reference}/{sample}.fa"
    conda: "envs/bcftools.yaml"
    shell: "bcftools consensus {input.vcf} < {input.ref} > {output} "