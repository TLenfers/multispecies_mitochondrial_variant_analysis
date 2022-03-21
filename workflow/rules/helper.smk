# compress the vcf file
rule zip_vcf:
    input:
        "results/calls_{caller}/{reference}/{sample}.vcf",
    output:
        "results/calls_{caller}/{reference}/{sample}.vcf.gz",
    conda:
        "../envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    log:
        "logs/{caller}/{reference}/zip_vcf/{sample}.log",
    shell:
        "bgzip -f {input}; tabix -f -p vcf {output}"
