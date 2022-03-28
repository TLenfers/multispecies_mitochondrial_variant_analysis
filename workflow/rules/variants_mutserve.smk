###############################################################################
# variant calling (mutserve)
###############################################################################


# calls the variants using bcftools from the bam file
# works only for human mitochondrial data
# mutserve is not available in conda, therefore a container will be used
rule mutserve:
    input:
        bam="results/mapped/{reference}/{sample}.bam",
        bamidx="results/mapped/{reference}/{sample}.bam.bai",
        ref="{{reference_path}}/{reference}.fa",
        index="{{reference_path}}/{reference}.fa.fai",
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
