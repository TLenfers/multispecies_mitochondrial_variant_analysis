# prepares variants data for creating the consensus sequence
# removes indels and SNPs with low quality
# quality score is set to min of 200 given out by bcftools
rule prepare_for_seq:
    input:
        vcf="results/calls_{caller}/{reference}/{sample}.vcf.gz",
    output:
        vcf="results/sequences/{caller}/{reference}/{sample}.vcf.gz",
        tbi="results/sequences/{caller}/{reference}/{sample}.vcf.gz.tbi",
    conda:
        "../envs/vcftools.yaml"
    log:
        "logs/{caller}/{reference}/{sample}.prepare_for_seq.log",
    shell:
        "vcftools --gzvcf {input.vcf} --remove-indels --minQ 200.0 --recode --recode-INFO-all --stdout |bgzip > {output.vcf}; tabix -p vcf {output.vcf}"


# creates the consensus sequence
rule vcf_to_fasta:
    input:
        ref=Path.joinpath(Path(config["reference_path"]), "{reference}.fa"),
        vcf="results/sequences/{caller}/{reference}/{sample}.vcf.gz",
        index="results/sequences/{caller}/{reference}/{sample}.vcf.gz.tbi",
    output:
        fa="results/sequences/{caller}/{reference}/{sample}.fa",
    conda:
        "../envs/bcftools.yaml"
    log:
        "logs/vcf_to_fasta/{caller}/{reference}/{sample}.log",
    shell:
        "bcftools consensus {input.vcf} < {input.ref} > {output.fa} "
