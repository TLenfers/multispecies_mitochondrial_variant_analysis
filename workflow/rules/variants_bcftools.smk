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
        "envs/bcftools.yaml"
    log:
        "logs/{reference}/{sample}.bcftools.log",
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
        "../envs/bcftools.yaml"
    log:
        "logs/{reference}/{sample}.log",
    shell:
        "bcftools norm  --fasta-ref {input.ref} --check-ref -m {input.vcf} | bcftools view -Ov -o {output}"
