# Multispecies Mitochondrial Variant Analysis 
 [![Actions Status](https://github.com/TLenfers/multispecies_mitochondrial_variant_analysis/workflows/Linting/badge.svg)](https://github.com/TLenfers/multispecies_mitochondrial_variant_analysis/actions)
 
This workflow performs a variant analysis on mitochondrial genomes using the bcftools variant caller.
For human samples, the workflow also performs a variant analysis using mutserve.
 
![mt_analyses](https://user-images.githubusercontent.com/14835042/160881351-98898821-e286-49ea-9f10-6547a1798ecd.png) 

## Local installation
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
source activate wf
```
## Configuration
**Config files**:

  - [`config.yaml`](config/config.yaml) - analysis-specific settings 
  - [`environment.yaml`](environment.yaml) - software dependencies and versions
  - [`samples.tsv`](config/samples.tsv) - list of (paired) samples

**Samples:**

  - Put all sample names in a single column in [`samples.tsv`](/config/samples.tsv).
  - Add the data folder to [`config.yaml`](config/config.yaml) where all files to be analysed are located
    - standard path is set to `data/` 
  - assumed naming convention:
    - sampleName_R1.fastq.gz
    - sampleName_R2.fastq.gz

**Reference:**

To analyse dog, mouse or human samples the corresponding reference will be downloaded.
Define the to be analysed species in [`config.yaml`](/config/config.yaml) under reference.

If you want to analyse a different species or use your own reference, enter the name of the file and it's path in [`config.yaml`](config/config.yaml).

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
The results are in a sub-folder corresponding to the name of the reference file used.  
- `/results/calls_bcftools` contains all called variants using bcftools. The variants are normalized and saved as `sample_name.vcf.gz`. In addition, the file `mergerd.vcf` is created in which all variants are merged together.
  - Example file without header:
```bash
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	results/mapped/human/sample_name.bam
    chrM	73	.	A	G	225.417	.	DP=253;VDB=3.59147e-17;SGB=-0.693147;MQSBZ=0;FS=0;MQ0F=0;AC=1;AN=1;DP4=0,0,240,5;MQ=60	GT:PL	1:255,0
    chrM	146	.	T	C	225.422	.	DP=242;VDB=0.795672;SGB=-0.693147;MQSBZ=0;FS=0;MQ0F=0;AC=1;AN=1;DP4=0,0,165,52;MQ=60	GT:PL	1:255,0
```    
- `/results/calls_mutserve` contains all called variants using mutserve.
- `/results/mapped` contains all aligned reads as `sample_name.bam` ad their index file `sample_name.bam.bai`. 
- `/results/plots` contains the created heatmap plots for the bctools caller. Example plots:
  - [ref_heatmap.pdf](https://github.com/TLenfers/multispecies_mitochondrial_variant_analysis/files/8379593/ref_heatmap.pdf)
  - [ref_heatmap_clusterrow.pdf](https://github.com/TLenfers/multispecies_mitochondrial_variant_analysis/files/8379594/ref_heatmap_clusterrow.pdf)
  - The name of the samples is on the X-axis, the variants on the Y-axis
  - The values of the heatmap refer to the Phred-scaled likelihood for homomorphic reference allele (scale 0-255; 255: reference is very unlikely -> alternative more likely).
  - The plots of `alt_heatmap` are containing the Phred-scaled likelihood for homomorphic alternative allele, i.e. that the variant is present at this position (scale 0-255; 0: variant is present).
- `/results/sequences` contains the created consensus sequences for each sample in regard to the used reference in fasta format as `sample_name.fa`.

## Snakedeploy usage
The usage of this workflow is described in the [Snakemake Workflow Catalog](https://snakemake.github.io/snakemake-workflow-catalog?usage=TLenfers/multispecies_mitochondrial_variant_analysis).
```bash
# To run human variant analysis with mutserve
snakemake --cores all all_human --use-conda --use-singularity 
```

If you use this workflow in a paper, don't forget to give credits to the authors by citing the URL of this (original) repository.

