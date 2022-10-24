# Project #1 “What causes antibiotic resistance?” Alignment to reference, variant calling.

## Prepare the workspace

### Create a directory for files and results

```ruby
mkdir WorkshpoBI/Project1/raw_data
cd WorkshpoBI/Project1/raw_data
```

### Download raw data

1. Reference sequence of the parental _E. coli_ strain:
```ruby
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gff.gz
```
2. Raw Illumina sequencing reads of resistant to the ampicillin _E. coli_ strain:

Just go to the [link](https://doi.org/10.6084/m9.figshare.10006541.v3) and click "download".

### Checking the data

```ruby
ls
```
- `GCF_000005845.2_ASM584v2_genomic.fna.gz` — reference sequence
- `GCF_000005845.2_ASM584v2_genomic.gff.gz` — annotation
- `amp_res_1.fastq.gz` — forward rerads
- `amp_res_2.fastq.gz` — reverse rerads

## Inspect raw sequencing data manually

### Checking the files structure

```ruby
zcat amp_res_1.fastq.gz | head -20
zcat amp_res_2.fastq.gz | head -20
```
Each read has 4 lines of information, and then the next read starts on the following line. The first line starts with the `@` symbol, and contains identifiers and information about this read. The next line contains the actual sequence, then on line three there is a `+`, which may sometimes have the identifier and info repeated. Line 4 contains the quality string, where ASCII characters encode the quality score
for each base. 

| Meaning  | Line |
| ------------- |-------------|
| Identifiers       | `@SRR1363257.37 GWZHISEQ01:153:C1W31ACXX:5:1101:14027:2198 length=101`    |
| Sequence     | `GGTTGCAGATTCGCAGTGTCGCTGTTCCAGCGCATCACATCTTTGATGTTCACGCCGTGGCGTACACG`     |
|  Separator  | `+`    |
| Quality     | `@?:=:;DBFADH;CAECEE@@:FFHGAE4?C?DE<BFGEC>?>FHE4BFFIIFHIBABEECA83;>>@`     |


### Counting reads in files

```ruby
zcat amp_res_1.fastq.gz | wc -l
zcat amp_res_2.fastq.gz | wc -l
```
Each file contains `1 823 504` lines, which means `455 876` reads.

## Inspect raw sequencing data with fastqc. Filtering the reads.

### Installing fastqc via mamba

```ruby
wget -O Mambaforge.sh "https://github.com/conda-forge/miniforge/releases/latest/download/Mambaforge-$(uname)-$(uname -m).sh"
bash Mambaforge.sh -b -p "${HOME}/conda"
source "${HOME}/conda/etc/profile.d/conda.sh”
conda activate
mamba install -c bioconda fastqc
```
### Checking fastqc

```ruby
fastqc -h
```
If we see the [manual page](https://home.cc.umanitoba.ca/~psgendb/doc/fastqc.help), everything is fine.

### Run the program fastqc
```ruby
fastqc -o ./QCbefore ./raw_data/amp_res_1.fastq ./raw_data/amp_res_2.fastq
```
### Basic statistics of raw reads
| Measure  | Value (Forward)|Value (Reverse)|
| ------------- |-------------|-------------|
| Total Sequences | 455876    | 455876 |
|  Sequence length  | 101   | 101  |
| %GC    | 50  | 50 |
| Low quality    | Per base sequence quality,  Per tile sequence quality |Per base sequence quality |

