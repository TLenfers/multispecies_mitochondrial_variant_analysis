# Multispecies Mitochondrial Variant Analysis 

This workflow performs a variant analysis on mitochondrial genomes using the bcftools variant caller.
For human samples, the workflow also performs a variant analysis using mutserve.

## Installation
1. Clone this repo
```bash
git clone https://github.com/tlenfers/multispecies_mitochondrial_variant_analysis.git
cd multispecies_mitochondrial_variant_analysis
```

2. Install dependencies
```bash
# download Miniconda3 installer
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O miniconda.sh
# install Conda (respond by 'yes')
bash miniconda.sh
# update Conda
conda update -y conda
# setup channels 
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
# create & activate new env with installed deps
conda env create -n wf -f environment.yaml
conda activate wf
```
## Configuration
**Config files**:

  - [`config.yaml`](config/config.yaml) - analysis-specific settings 
  - [`environment.yaml`](environment.yaml) - software dependencies and versions
  - [`samples.tsv`](config/samples.tsv) - list of (paired) samples

**Samples:**

  - Put all sample names in a single column in [`samples.tsv`](/config/samples.tsv).
  - Put all sample files in `workflow/data/samples/`
  - assumed naming convention:
    - sampleName_R1.fastq.gz
    - sampleName_R2.fastq.gz

**Reference:**

To analyse dog, mouse or human samples the corresponding reference will be downloaded.
Define the to be analysed species in [`config.yaml`](/config/config.yaml) under reference.

If you want to analyse a different species or use your own reference, enter the name of the file in [`config.yaml`](config/config.yaml).
The reference file (name.fa) should be put in `workflow/data/reference`.

## Execute the workflow
```bash
cd workflow
# 'dry' run only checks I/O files
snakemake -np

# To run mutlipecies variant analysis
snakemake -j n all --use-conda --use-singularity
# where n is the numer of cores to use

# To run human variant analysis with mutserve
snakemake -j n all_human --use-conda --use-singularity
# where n is the numer of cores to use
```


## Output
All output files are put in `/results` and in their own subfolder regarding the used reference and caller.


If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.
