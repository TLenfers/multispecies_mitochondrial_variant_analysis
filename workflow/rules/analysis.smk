# merges the variants
rule merge_vcf:
    input:
        expand(
            "results/calls_{{caller}}/{{reference}}/{sample}.vcf.gz",
            sample=samples,
        ),
    output:
        "results/calls_{caller}/{reference}/merged.vcf",
    conda:
        "../envs/bcftools.yaml"
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    log:
        "logs/{caller}/{reference}/merge_vcf/merged.log",
    shell:
        "bcftools merge -m none -O v {input} > {output}"


# extracts the common variants out of the merged vcf file
rule common_variants:
    input:
        in_file="results/calls_{caller}/{reference}/merged.vcf",
    output:
        out="results/calls_{caller}/{reference}/commonVariants.tsv",
    conda:
        "../envs/r-heatmap.yaml"
    log:
        "logs/calls_{caller}/{reference}/commonVariants.log",
    script:
        "../scripts/common_variants.R"


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
        "../envs/r-heatmap.yaml"
    log:
        "logs/{caller}/{reference}/plot_variant_heatmap.log",
    script:
        "../scripts/plot_variant_heatmap.R"
