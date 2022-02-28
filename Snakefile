# TODO: doku: name for samples: sampleid_L001_R1 / R2
# TODO: doku: Referenz Datei immer .fa

configfile: "config/config.yaml"

rule all:
    input:
        expand("results/calls_bcftools/{reference}/{sample}.vcf",sample=config["samples"],reference=config["reference"]),
        expand("results/calls_bcftools/{reference}/merged.vcf",reference=config["reference"]),
        expand("results/calls_bcftools/{reference}/normalized/{sample}.vcf.gz",sample=config["samples"],reference=config["reference"]),
expand("results/calls_bcftools/{reference}/normalized/merged_normalized.vcf",reference=config["reference"])#,
    #todo uncomment
        #expand("results/plots/{reference}/heatmap.pdf",reference=config["reference"])

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
        reference="[A-Za-z0-9]+",
        sample="[A-Za-z0-9]+"
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
        reference="[A-Za-z0-9]+",
        sample="[A-Za-z0-9]+"
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
        "results/calls_bcftools/{reference}/normalized/{sample}.vcf"
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools norm  --fasta-ref {input.ref} --check-ref -m {input.vcf} | bcftools view -Ov -o {output}"


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

rule idx_normalized:
    input:
        "results/calls_bcftools/{reference}/normalized/{sample}.vcf"
    output:
        "results/calls_bcftools/{reference}/normalized/{sample}.vcf.gz"
    conda:
        "workflow/envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
        sample="[A-Za-z0-9]+"
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

rule merge_vcf_normalized:
    input:
        expand("results/calls_bcftools/{{reference}}/normalized/{sample}.vcf.gz",sample=config["samples"])
    output:
        "results/calls_bcftools/{reference}/normalized/merged_normalized.vcf"
    conda: "workflow/envs/bcftools.yaml"
    shell: "bcftools merge -m none -O v {input} > {output}"




# FIXME conda
#rule plot_variant_heatmap:
#    input:
#        "results/calls_bcftools/{reference}/merged.vcf"
#    output:
#        "results/plots/{reference}/heatmap.pdf",
#        "results/plots/{reference}/heatmap_clusterrow.pdf",
#        "results/plots/{reference}/heatmap_clusterrow_removedCommonVariants.pdf",
#        "results/plots/{reference}/heatmap_removedCommonVariants.pdf"
#    conda:
#        "workflow/envs/r-heatmap.yaml"
#    script:
#        "workflow/scripts/plot_variant_heatmap.R"