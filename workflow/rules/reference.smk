# downloading the mitochondrial reference genome for mouse, dog and human data
rule get_ref:
    output:
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa"),
    log:
        "logs/{reference}/get_ref.log",
    run:
        if config["reference"] == "mouse":
            shell(
                "wget -O- http://ftp.ensembl.org/pub/current_fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.chromosome.MT.fa.gz | zcat > {output}"
            )
        if config["reference"] == "dog":
            shell(
                "wget -O- 'https://www.ncbi.nlm.nih.gov/sviewer/viewer.fcgi?id=U96639.2&db=nuccore&report=fasta' > {output}"
            )
        if config["reference"] == "human":
            shell(
                "wget -O- http://ftp.ensembl.org/pub/current_fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.MT.fa.gz |zcat > {output} && sed -i '1c\>chrM' {output}"
            )


# creates the index for the reference genome
rule bwa_idx:
    input:
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa"),
    output:
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa.amb"),
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa.ann"),
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa.bwt"),
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa.pac"),
        Path.joinpath(Path(config["reference_path"]), "{reference}.fa.sa"),
    wildcard_constraints:
        reference="[A-Za-z0-9]+",
    conda:
        "../envs/bwa.yaml"
    log:
        "logs/{reference}/bwa_idx.log",
    shell:
        "bwa index {input}"
