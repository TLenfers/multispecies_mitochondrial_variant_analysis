# multispecies_mitochondrial_variant_calling
- all to be analyzed samples have to be put in /data/samples
- all sample files have to be of type .fastaq.gz
- the reference file has to be put in /dat/reference as {reference}.fa
- in the config in /config/config.yaml the names of all to be analyzed files must be entered under *samples:[]*, as well as the name of the reference file in *reference:*
- to execute the pipeline snakemake has to be installed inside a conda environment