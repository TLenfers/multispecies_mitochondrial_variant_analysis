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
    container:
        "docker://stephenturner/mutserve"
    log:
        "logs/calls_mutserve/{reference}/{sample}.log",
    shell:
        "mutserve call --reference {input.ref} --output {output} {input.bam}"
