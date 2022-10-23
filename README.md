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

## Inspect raw sequencing data manually

##### Checking the files structure

```ruby
zcat amp_res_1.fastq.gz | head -20
zcat amp_res_2.fastq.gz | head -20
```
Each read has 4 lines of information, and then the next read starts on the following line. The first line starts with the `@` symbol, and contains identifiers and information about this read. The next line contains the actual sequence, then on line three there is a `+`, which may sometimes have the identifier and info repeated. Line 4 contains the quality string, where ASCII characters encode the quality score
for each base. 

| Meaning  | Line |
| ------------- |:-------------:|
| Identifiers       | `@SRR1363257.37 GWZHISEQ01:153:C1W31ACXX:5:1101:14027:2198 length=101`    |
| Sequence     | `GGTTGCAGATTCGCAGTGTCGCTGTTCCAGCGCATCACATCTTTGATGTTCACGCCGTGGCGTACACG`     |
|  Separator  | `+`    |
| Quality     | `@?:=:;DBFADH;CAECEE@@:FFHGAE4?C?DE<BFGEC>?>FHE4BFFIIFHIBABEECA83;>>@`     |

