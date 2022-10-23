# Project #1 “What causes antibiotic resistance?” Alignment to reference, variant calling.

## Prepare the workspace

##### Create a directory for files and results

```ruby
mkdir WorkshpoBI/Project1/raw_data
cd WorkshpoBI/Project1/raw_data
```

##### Download raw data

1. Reference sequence of the parental _E. coli_ strain (sequence in `.fna` and annotation in `.gff`):
```ruby
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
```
2. Raw Illumina sequencing reads of resistant to the ampicillin _E. coli_ strain:

Just go to the [link](https://doi.org/10.6084/m9.figshare.10006541.v3) and click "download".

##### Checking the data

```ruby
ls
```
- `GCF_000005845.2_ASM584v2_genomic.fna.gz`
- `GCF_000005845.2_ASM584v2_genomic.gff.gz`
- `amp_res_1.fastq.gz`
- `amp_res_2.fastq.gz`
