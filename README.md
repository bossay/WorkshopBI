# Project #2 “Why did I get the flu?”. Deep sequencing, error control, p-value, viral evolution.

## 1. Prepare the workspace

### 1.1. Create a directory for files and results

```ruby
mkdir WorkshpoBI/Project2/raw_data
cd WorkshpoBI/Project2/raw_data
```

### 1.2. Download raw data

Just go to the [link](http://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/001/SRR1705851/) and click to `SRR1705851.fastq.gz` to get the sequence being studied.

Go to the [link](https://www.ncbi.nlm.nih.gov/nuccore/KF848938.1?report=fasta), click to `Send to`, select `File` and format `FASTA` to get the reference sequence. I called it `reference.fasta`.


## 2. Inspect raw sequencing data

### 2.1. Checking the files structure

```ruby
zcat SRR1705851.fastq.gz | head -20
```
The data has a standard .fasta file structure.

### 2.2. Counting reads in files

```ruby
zcat SRR1705851.fastq.gz | wc -l
```
Each file contains `1 433 060` lines, which means `358 265` reads.

### 2.3. Run fastqc

```ruby
mkdir ../QC
gunzip SRR1705851.fastq.gz
fastqc -o ../QC SRR1705851.fastq
```
### 2.4. Basic statistics of raw reads

| Measure  | Value (Forward)|
| ------------- |-------------|
| Total Sequences | 358 265    |
|  Sequence length  | 35-151  |
| %GC    | 42  |
| Low quality    | Per base sequence content, Sequence Duplication Levels |

## 3. Aligning sequences to reference

### 3.1. Run bwa index

```ruby
bwa index -p reference.fasta reference.fasta
```
There are 5 new files in the directory:

> - `reference.fasta.amb`
> - `reference.fasta.ann`
> - `reference.fasta.bwt`
> - `reference.fasta.pac`
> - `reference.fasta.sa`

### 3.2. Run pipe with bwa mem and samtools

```ruby
bwa mem reference.fasta SRR1705851.fastq | samtools view -S -b - | samtools sort -o SRR1705851_aligned_sorted.bam
```
The names of the files `.amb`, `.ann`, `.bwt`, `.pac`, `.sa` must match the name of the file with the reference genome.

### 3.3. Index bam file
```ruby
samtools index SRR1705851_aligned_sorted.bam
```
### 3.4. Create an intermediate mpileup file

Set depth limit with the -d flag.

```ruby
samtools mpileup -d 0 -f reference.fasta SRR1705851_aligned_sorted.bam > SRR1705851.mpileup
```
## 4. Look for common variants with VarScan

### 4.1. VarScan installation

Download `VarScan.v2.4.0.jar` from [GitHub repository](https://github.com/dkoboldt/varscan) to your working directory.

### 4.2. Run VarScan (0.95)

Run the program with the threshold 0.95:

```ruby
java -jar VarScan.v2.4.0.jar mpileup2snp SRR1705851.mpileup --min-var-freq 0.95 --variants --output-vcf 1 > VarScan_results_0_95.vcf
```
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.95
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from SRR1705851.mpileup
> - 1665 bases in pileup file
> - 5 variant positions (5 SNP, 0 indel)
> - 0 were failed by the strand-filter
> - 5 variant positions reported (5 SNP, 0 indel)

### 4.3. Extract the basic information (0.95)

Extract the basic variant information from the `.vcf` file:

```ruby
cat VarScan_results_0_95.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > Variants_0_95.txt
```

### 4.4. Visualize .vcf file
Go to `IGV` browser and examine the mutations you find.

##  5. Look for rare variants with VarScan

### 5.1. Run VarScan (0.001)

Run the program with the threshold 0.001:

```ruby
java -jar VarScan.v2.4.0.jar mpileup2snp SRR1705851.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results_0_001.vcf
```
> - Only SNPs will be reported
> - Warning: No p-value threshold provided, so p-values will not be calculated
> - Min coverage:   8
> - Min reads2:     2
> - Min var freq:   0.001
> - Min avg qual:   15
> - P-value thresh: 0.01
> - Reading input from SRR1705851.mpileup
> - 1665 bases in pileup file
> - 23 variant positions (21 SNP, 2 indel)
> - 0 were failed by the strand-filter
> - 21 variant positions reported (21 SNP, 0 indel)


### 5.2. Extract the basic information (0.001)

Extract the basic variant information from the `.vcf` file:

```ruby
cat VarScan_results_0_001.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > Variants_0_001.txt
```

## 6. Inspect and align the control sample sequencing data

### 6.1. Download raw data

```ruby
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/008/SRR1705858/SRR1705858.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/009/SRR1705859/SRR1705859.fastq.gz
wget ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR170/000/SRR1705860/SRR1705860.fastq.gz

gunzip SRR1705858.fastq.gz
gunzip SRR1705859.fastq.gz
gunzip SRR1705860.fastq.gz
```
### 6.2. Repeat all the steps for new files

```ruby
cat SRR1705858.fastq | wc -l
cat SRR1705859.fastq | wc -l
cat SRR1705860.fastq | wc -l
```
The files contain `256 586`, `233 327` and `249 964` reads, respectively.

```ruby
fastqc -o ../QC SRR1705858.fastq
fastqc -o ../QC SRR1705859.fastq
fastqc -o ../QC SRR1705860.fastq
```
The per base sequence quality was sufficient for all samples.

```ruby
bwa mem reference.fasta SRR1705858.fastq | samtools view -S -b - | samtools sort -o SRR1705858_aligned_sorted.bam

bwa mem reference.fasta SRR1705859.fastq | samtools view -S -b - | samtools sort -o SRR1705859_aligned_sorted.bam

bwa mem reference.fasta SRR1705860.fastq | samtools view -S -b - | samtools sort -o SRR1705860_aligned_sorted.bam
```

```ruby
samtools index SRR1705858_aligned_sorted.bam
samtools index SRR1705859_aligned_sorted.bam
samtools index SRR1705860_aligned_sorted.bam
```

```ruby
samtools mpileup -d 0 -f reference.fasta SRR1705858_aligned_sorted.bam > SRR1705858.mpileup

samtools mpileup -d 0 -f reference.fasta SRR1705859_aligned_sorted.bam > SRR1705859.mpileup

samtools mpileup -d 0 -f reference.fasta SRR1705860_aligned_sorted.bam > SRR1705860.mpileup
```

```ruby
java -jar VarScan.v2.4.0.jar mpileup2snp SRR1705858.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results58_0_001.vcf

java -jar VarScan.v2.4.0.jar mpileup2snp SRR1705859.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results59_0_001.vcf

java -jar VarScan.v2.4.0.jar mpileup2snp SRR1705860.mpileup --min-var-freq 0.001 --variants --output-vcf 1 > VarScan_results60_0_001.vcf
```

```ruby
cat VarScan_results58_0_001.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > Variants58_0_001.txt

cat VarScan_results59_0_001.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > Variants59_0_001.txt

cat VarScan_results60_0_001.vcf | awk 'NR>24 {print $1, $2, $4, $5, $10}' > Variants60_0_001.txt
```
### 6.3. Calculate the average and standard deviation

 Calculate the average and standard deviation of the frequencies reported within each list.

## Ready! You`re incredible!
