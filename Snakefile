configfile: "config/config.yaml"

# Executing all rules for complete analysis
# Works for Human, Mouse and dog mitochondrial data
rule all:
    input:
        expand(
            "results/calls_bcftools/{reference}/norm_{sample}.vcf.gz",
            caller="bcftools",
            sample=config["samples"],
            reference=config["reference"],
        ),
        expand(
            "results/calls_{caller}/{reference}/merged.vcf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/calls_{caller}/{reference}/commonVariants.tsv",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/ref_heatmap.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/ref_heatmap_clusterrow.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/alt_heatmap.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/alt_heatmap_clusterrow.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/sequences/{caller}/{reference}/{sample}.fa",
            caller="bcftools",
            reference=config["reference"],
            sample=config["samples"],
        ),


# Executing all rules for complete analysis for Human mitochondrial data
# adds extra the step of calling the mitochondrial variants with mutserve
# mutserve executable hast to be installed and in the path mutserve/
rule all_human:
    input:
        expand(
            "results/calls_{caller}/{reference}/norm_{sample}.vcf.gz",
            caller="bcftools",
            sample=config["samples"],
            reference=config["reference"],
        ),
        expand(
            "results/calls_{caller}/{reference}/merged.vcf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/calls_{caller}/{reference}/commonVariants.tsv",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/ref_heatmap.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/ref_heatmap_clusterrow.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/alt_heatmap.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/plots/{caller}/{reference}/alt_heatmap_clusterrow.pdf",
            caller="bcftools",
            reference=config["reference"],
        ),
        expand(
            "results/sequences/{caller}/{reference}/{sample}.fa",
            caller="bcftools",
            reference=config["reference"],
            sample=config["samples"],
        ),
        # mutserve:
        expand(
            "results/calls_{caller}/{reference}/{sample}.vcf.gz",
            caller="mutserve",
            sample=config["samples"],
            reference=config["reference"],
        ),
        expand(
            "results/calls_{caller}/{reference}/merged.vcf",
            caller="mutserve",
            reference=config["reference"],
        ),
        #TODO add plots & common variants for mutserve | r-scrpips for bcftools cannot be used


#        expand("results/calls_{caller}/{reference}/commonVariants.tsv",caller="mutserve",reference=config["reference"])#,
# expand("results/plots/{caller}/{reference}/ref_heatmap.pdf",caller="mutserve",reference=config["reference"]),
# expand("results/plots/{caller}/{reference}/ref_heatmap_clusterrow.pdf",caller="mutserve",reference=config["reference"]),
# expand("results/plots/{caller}/{reference}/alt_heatmap.pdf",caller="mutserve",reference=config["reference"]),
# expand("results/plots/{caller}/{reference}/alt_heatmap_clusterrow.pdf",caller="mutserve",reference=config["reference"]),

# TODO add consensus sequence for mutserve
# expand("results/sequences/{caller}/{reference}/{sample}.fa",caller="mutserve",reference=config["reference"],sample=config["samples"])

# downloading the mitochondrial reference genome for mouse, dog and human data
rule get_ref:
    output:
        "data/reference/{reference}.fa",
    run:
        if config["reference"] == "mouse":
            shell(
                "wget -O- http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.MT.fa.gz | zcat > {output}"
            )
        if config["reference"] == "dog":
            shell(
                "wget -O- https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=U96639.2&db=nuccore&report=fasta| zcat > {output}"
            )
        if config["reference"] == "human":
            shell(
                "wget -O- http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz |zcat > {output} && sed -i '1c\>chrM' {output}"
            )


# creates the index for the reference genome
rule bwa_idx:
    input:
        "data/reference/{reference}.fa",
    output:
        "data/reference/{reference}.fa.amb",
        "data/reference/{reference}.fa.ann",
        "data/reference/{reference}.fa.bwt",
        "data/reference/{reference}.fa.pac",
        "data/reference/{reference}.fa.sa",
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    conda:
        "workflow/envs/bwa.yaml"
    shell:
        "bwa index {input}"


# mapping of the reads to the reference genome using Burrows-Wheeler Aligner
rule bwa_mem:
    input:
        amb="data/reference/{reference}.fa.amb",
        ref="data/reference/{reference}.fa",
        fastq_r1="data/samples/{sample}_L001_R1_001.fastq.gz",
        fastq_r2="data/samples/{sample}_L001_R2_001.fastq.gz",
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    output:
        "results/mapped/{reference}/{sample}.bam",
    conda:
        "workflow/envs/bwa.yaml"
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
    shell:
        "samtools faidx {input}"


###############################################################################
# variant calling (bcftools)
###############################################################################
# calls the variants using bcftools from the bam file
rule call_variants:
    input:
        bam="results/mapped/{reference}/{sample}.bam",
        bamidx="results/mapped/{reference}/{sample}.bam.bai",
        ref="data/reference/{reference}.fa",
        index="data/reference/{reference}.fa.fai",
    output:
        "results/calls_bcftools/{reference}/{sample}.vcf",
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools mpileup -f {input.ref} {input.bam} | bcftools call -mv --ploidy 1 -> {output}"


# left-align and normalization of the variants
rule normalize_variants:
    input:
        vcf="results/calls_bcftools/{reference}/{sample}.vcf.gz",
        ref="data/reference/{reference}.fa",
    output:
        "results/calls_bcftools/{reference}/norm_{sample}.vcf",
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools norm  --fasta-ref {input.ref} --check-ref -m {input.vcf} | bcftools view -Ov -o {output}"


###############################################################################
# variant calling (mutserve)
###############################################################################


# calls the variants using bcftools from the bam file
# works only for human mitochondrial data
rule mutserve:
    input:
        bam="results/mapped/{reference}/{sample}.bam",
        bamidx="results/mapped/{reference}/{sample}.bam.bai",
        ref="data/reference/{reference}.fa",
        index="data/reference/{reference}.fa.fai",
    output:
        "results/calls_mutserve/{reference}/{sample}.vcf",
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    shell:
        "./mutserve/mutserve call --reference {input.ref} --output {output} {input.bam}"
        #"java -jar ./mutserve/mutserve.jar call --reference {input.ref} --output {output} {input.bam}"


###############################################################################
# analysis

# compress the vcf file
rule zip_vcf:
    input:
        "results/calls_{caller}/{reference}/{sample}.vcf",
    output:
        "results/calls_{caller}/{reference}/{sample}.vcf.gz",
    conda:
        "workflow/envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    shell:
        "bgzip -f {input}; tabix -f -p vcf {output}"


# merges the variants
rule merge_vcf:
    input:
        expand(
            "results/calls_{{caller}}/{{reference}}/{sample}.vcf.gz",
            sample=config["samples"],
        ),
    output:
        "results/calls_{caller}/{reference}/merged.vcf",
    conda:
        "workflow/envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    shell:
        "bcftools merge -m none -O v {input} > {output}"


# extracts the common variants out of the merged vcf file
rule common_variants:
    input:
        in_file="results/calls_{caller}/{reference}/merged.vcf",
    output:
        "results/calls_{caller}/{reference}/commonVariants.tsv",
    conda:
        "workflow/envs/r-heatmap.yaml"
    script:
        "workflow/scripts/common_variants.R"


# creates the heatmap plots for the common variants of dog and mouse data
rule plot_variant_heatmap:
    input:
        "results/calls_{caller}/{reference}/commonVariants.tsv",
    output:
        "results/plots/{caller}/{reference}/ref_heatmap.pdf",
        "results/plots/{caller}/{reference}/ref_heatmap_clusterrow.pdf",
        "results/plots/{caller}/{reference}/alt_heatmap.pdf",
        "results/plots/{caller}/{reference}/alt_heatmap_clusterrow.pdf",
    conda:
        "workflow/envs/r-heatmap.yaml"
    script:
        "workflow/scripts/plot_variant_heatmap.R"

# prepares variants data for creating the consensus sequence
# removes indels and SNPs with low quality
# quality score is set to min of 200 given out by bcftools
rule prepare_for_seq:
    input:
        vcf="results/calls_{caller}/{reference}/{sample}.vcf.gz",
    output:
        "results/sequences/{caller}/{reference}/{sample}.vcf.gz",
        "results/sequences/{caller}/{reference}/{sample}.vcf.gz.tbi",
    conda:
        "workflow/envs/vcftools.yaml"
    shell:
        "vcftools --gzvcf {input.vcf} --remove-indels --minQ 200.0 --recode --recode-INFO-all --stdout |bgzip > {output[0]}; tabix -p vcf {output[0]}"


# creates the consensus sequence
rule vcf_to_fasta:
    input:
        ref="data/reference/{reference}.fa",
        vcf="results/sequences/{caller}/{reference}/{sample}.vcf.gz",
        index="results/sequences/{caller}/{reference}/{sample}.vcf.gz.tbi",
    output:
        "results/sequences/{caller}/{reference}/{sample}.fa",
    conda:
        "workflow/envs/bcftools.yaml"
    shell:
        "bcftools consensus {input.vcf} < {input.ref} > {output} "
