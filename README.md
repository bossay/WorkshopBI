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
> - `GCF_000005845.2_ASM584v2_genomic.fna.gz` — reference sequence
> - `GCF_000005845.2_ASM584v2_genomic.gff.gz` — annotation
> - `amp_res_1.fastq.gz` — forward rerads
> - `amp_res_2.fastq.gz` — reverse rerads

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

### Run fastqc

```ruby
fastqc -o ./QCbefore amp_res_1.fastq amp_res_2.fastq
```
### Basic statistics of raw reads

| Measure  | Value (Forward)|Value (Reverse)|
| ------------- |-------------|-------------|
| Total Sequences | 455 876    | 455 876 |
|  Sequence length  | 101   | 101  |
| %GC    | 50  | 50 |
| Low quality    | Per base sequence quality,  Per tile sequence quality |Per base sequence quality |

## Filtering the reads
### Installing trimmomatic via conda

Install a program for trimming adapters and low-quality reads - [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic)

```ruby
conda install -c bioconda trimmomatic
```
### Checking trimmomatic

```ruby
trimmomatic
```
If we see the manual page, everything is fine.

### Run trimmomatic

Run Trimmomatic in paired end mode, with following parameters:
+ Cut bases off the start of a read if quality below 20 (`LEADING`)
+ Cut bases off the end of a read if quality below 20 (`TRAILING`)
+ Trim reads using a sliding window approach, with window size 10 and average quality within the window 20 (`SLIDINGWINDOW:10:20`)
+ Drop the read if it is below length 20 (`MINLEN`)

```ruby
trimmomatic PE -threads 4 amp_res_1.fastq.gz amp_res_2.fastq.gz amp_1.trimmed.fastq amp_1un.trimmed.fastq amp_2.trimmed.fastq amp_2un.trimmed.fastq LEADING:20 TRAILING:20 SLIDINGWINDOW:10:20 MINLEN:20
```
The command is executed with a comment:
> Quality encoding detected as phred33
>
> Input Read Pairs: 455876 Both Surviving: 446259 (97.89%) Forward Only Surviving: 9216 (2.02%) Reverse Only Surviving: 273 (0.06%) Dropped: 128 (0.03%)

### Counting reads in files

```ruby
wc -l amp_1.trimmed.fastq
```
Each file contains `1 785 036` lines, which means `446 259` reads.

### Repeat the fastqc analysis 

```ruby
fastqc -o ./QCafter amp_1.trimmed.fastq amp_2.trimmed.fastq
```
### Basic statistics of raw reads

| Measure  | Value (Forward)|Value (Reverse)|
| ------------- |-------------|-------------|
| Total Sequences | 446 259    | 446 259 |
|  Sequence length  | 20-101   | 20-101  |
| %GC    | 50  | 50 |
| Low quality    | Per tile sequence quality |- |

## Aligning sequences to reference
### Checking bwa

```ruby
bwa
```
If we see the manual page, everything is fine.

## Run bwa index

```ruby
bwa index -p index GCF_000005845.2_ASM584v2_genomic.fna
```
There are 4 new files in the directory:

> - `GCF_000005845.2_ASM584v2_genomic.fna.amb`
> - `GCF_000005845.2_ASM584v2_genomic.fna.ann`
> - `GCF_000005845.2_ASM584v2_genomic.fna.bwt`
> - `GCF_000005845.2_ASM584v2_genomic.fna.pac`
> - `GCF_000005845.2_ASM584v2_genomic.fna.sa`

## Align your reads
### Run bwa mem

```ruby
bwa mem GCF_000005845.2_ASM584v2_genomic.fna amp_1.trimmed.fastq amp_2.trimmed.fastq > GCF_000005845.2_ASM584v2_genomic.fna.sam
```
The names of the files `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` must match the name of the file with the reference genome.

## Compress SAM file
### Run samtools

```ruby
samtools view -S -b GCF_000005845.2_ASM584v2_genomic.fna.sam > GCF_000005845.2_ASM584v2_genomic.fna.bam
```
### Checking basic statistics
```ruby
samtools flagstat GCF_000005845.2_ASM584v2_genomic.fna.bam
```
>- 892776 + 0 in total (QC-passed reads + QC-failed reads)
> - 0 + 0 secondary
> - 258 + 0 supplementary
> - 0 + 0 duplicates
> - 891649 + 0 mapped (99.87% : N/A)
> - 892518 + 0 paired in sequencing
> - 446259 + 0 read1
> - 446259 + 0 read2
> - 888554 + 0 properly paired (99.56% : N/A)
> - 890412 + 0 with itself and mate mapped
> - 979 + 0 singletons (0.11% : N/A)
> - 0 + 0 with mate mapped to a different chr
> - 0 + 0 with mate mapped to a different chr (mapQ>=5)

## Sort and index BAM file

### Sort bam file by sequence coordinate on reference
```ruby
samtools sort GCF_000005845.2_ASM584v2_genomic.fna.bam -o GCF_000005845.2_ASM584v2_genomic.fna.bam_sorted.bam
```

### Index bam file
```ruby
samtools index GCF_000005845.2_ASM584v2_genomic.fna.bam_sorted.bam
```

### Visualize the data

Visualize what the data actually looks like in the [IGV browser](https://software.broadinstitute.org/software/igv/)

Select `Genomes`, `Load Genome from File` and select reference genome. Then select `File`, `Open from file` and select BAM file.

## Variant calling

### Create an intermediate mpileup file

```ruby
samtools mpileup -f GCF_000005845.2_ASM584v2_genomic.fna GCF_000005845.2_ASM584v2_genomic.fna.bam_sorted.bam > GCF_000005845.2_ASM584v2_genomic.mpileup
```
### VarScan installation
Download `VarScan.v2.4.0.jar` from [GitHub repository](https://github.com/dkoboldt/varscan) to your working directory.

### Checking VarScan

```ruby
java -jar VarScan.v2.4.0.jar mpileup2snp -h
```
If we see the manual page, everything is fine.

### Run VarScan

Run the program with the threshold 0.5 (minimum 50% of non-reference bases at a position required to
call it a mutation in the sample):

```ruby
java -jar VarScan.v2.4.0.jar mpileup2snp GCF_000005845.2_ASM584v2_genomic.mpileup --min-var-freq 0.5 --variants --output-vcf 1 > VarScan_results.vcf
```
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.5
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from GCF_000005845.2_ASM584v2_genomic.mpileup
> - 4641343 bases in pileup file
> - 9 variant positions (6 SNP, 3 indel)
> - 1 were failed by the strand-filter
> - 5 variant positions reported (5 SNP, 0 indel)

More about [Variant Call Format](https://samtools.github.io/hts-specs/VCFv4.2.pdf) format

## Automatic SNP annotation

### SnpEff installation

Use automatic tools for snp annotation `SnpEff`

```ruby
conda install -c bioconda snpEff
```
### Download sequence and annotation of reference

```ruby
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.gbff.gz
```
### Create database

Create empty text file `snpeff.config`, and add there just one string: 
>k12.genome : ecoli_K12

Create folder for the database

```ruby
mkdir -p data/k12
```
Put there .gbk file (unzip and rename to genes.gbk):

```ruby
gunzip GCF_000005845.2_ASM584v2_genomic.gbff.gz
cp GCF_000005845.2_ASM584v2_genomic.gbff data/k12/genes.gbk
```

Create database:
```ruby
snpEff build -genbank -v k12
```
### Annotation
```ruby
snpEff ann k12 VarScan_results.vcf > annotation_k12.vcf
```

### Visualize .vcf file

Go back to IGV browser and add two more tracks - `annotation_k12.vcf` and annotation `GCF_000005845.2_ASM584v2_genomic.gff.gz`.


## Ready! You`re incredible!
