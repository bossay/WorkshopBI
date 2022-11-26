# Project #3 E.coli outbreak investigation 

## 1. Prepare the workspace

### 1.1. Create a directory for files and results

```ruby
mkdir WorkshpoBI/Project3/raw_data
cd WorkshpoBI/Project3/raw_data
```

### 1.2. Download raw data

We provide three libraries from the TY2482:

> `SRR292678` - paired end, insert size 470 bp
>
> `SRR292862` - mate pair, insert size 2 kb
>
> `SRR292770` - mate pair, insert size 6 kb

Download the sequences:

```ruby
wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R1_001.fastq.gz

wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292678sub_S1_L001_R2_001.fastq.gz

wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R1_001.fastq.gz

wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292862_S2_L001_R2_001.fastq.gz

wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R1_001.fastq.gz

wget https://d28rh4a8wq0iu5.cloudfront.net/bioinfo/SRR292770_S1_L001_R2_001.fastq.gz
```
## 2. Inspect raw sequencing data

### 2.1. Create a directory and unpack the files

```ruby
mkdir ../QC
gunzip *.fastq.gz
```
### 2.2. Run fastqc

```ruby
fastqc -o ../QC *.fastq
```

### 2.4. Basic statistics of raw reads

| **Measure**         | **SRR292678 (F)** | **SRR292678 (R)** | **SRR292770 (F)** | **SRR292770 (R)** | **SRR292862 (F)** | **SRR292862 (R)** |
|---------------------|-------------------|-------------------|-------------------|-------------------|-------------------|-------------------|
| **Total Sequences** | 5 499 346           | 5 499 346           | 5 102 041           | 5 102 041           | 5 102 041           | 5 102 041           |
| **Sequence length** | 90                | 90                | 49                | 49                | 49                | 49                |
| **%GC**             | 49                | 49                | 50                | 49                | 50                | 49                |
| **Poor quality**   | 0                 | 0                 | 0                 | 0                 | 0                 | 0                 |


### 2.5. K-mer profile and genome size estimation

Install Jellyfish: a fast k-mer counter
```ruby
sudo apt-get install jellyfish
```

Run Jellyfish (k-mer sizes of 31, length ~ 5.5M):
```ruby
jellyfish count -m 31 -s 6M -C SRR292678sub_S1_L001_R1_001.fastq
jellyfish histo mer_counts.jf > hist_j.txt
```
### 2.6. Visualize k-mer distribution
Genome size can be calculated by counting k-mer frequency of the read data. Using the Jellyfish output we plot the graph with R:

```r
setwd('C:\\Users\\pmbog\\source\\WorkshpoBI\\Project3\\raw_data')
hist_k_mer <- read.table("hist_j.txt")

plot(hist_k_mer[4:150,],type="l")
points(hist_k_mer[4:150,])
```
Calculate the total number of k-mer in the distribution:

```r
sum(as.numeric(hist_k_mer[1:817,1]*hist_k_mer[1:817,2]))
```
In this case there are about `329 960 760` k-mers in the histogram.

Next, we want to know the peak position. From the graph, we can see its close to 60. Thus we examine the number close to 60 and find the maximum value
```
hist_k_mer[55:65,]
```
|    |       |
|----|-------|
| 55 | 86843 |
| 56 | 87293 |
| 57 | 87712 |
| 58 | 87824 |
| 59 | 88556 |
| 60 | 88875 |
| 61 | 89428 |
| 62 | 89747 |
| 63 | 88592 |
| 64 | 88192 |
| 65 | 87854 |

In this case, the peak is at 62. Then, the genome size can be estimated as:

```
sum(as.numeric(hist_k_mer[1:817,1]*hist_k_mer[1:817,2]))/62
```
This reads as `5 321 948` - `5.32 Gb`.

### 2.7. Estimate the genome size 

Let's calculate the size of the genome using the formula: 
`N = (M*L)/(L-K+1)`
`Genome_size = T/N`
> M = 62: k-mer peak
>
> K = 31: k-mer-size
>
> L - 90: avg read length
>
> T = 5499346*90 = 494941140: total bases
>
> N = (62*90)/(90-31+1) = 93: depth of coverage
>
> **G = 494941140/93 = 5321948**

## Assembling E. coli X genome from paired reads

### Installing SPAdes via conda

```ruby
conda install spades -c bioconda
```

To verify that the software was installed correctly, run it in the test mode:
```ruby
spades.py --test
```
If you do everything right, the last line of the output will be `Thank you for using SPAdes!`. More about [SPAdes](https://cab.spbu.ru/software/spades/)

### Run SPAdes
Run SPAdes in the paired-end mode, providing paired reads of E. coli X. from the library SRR292678 (forward and reverse):
 ```ruby
spades.py -1 SRR292678sub_S1_L001_R1_001.fastq -2 SRR292678sub_S1_L001_R2_001.fastq -o spades
```
### Installing QUAST via conda

```ruby
conda install quast -c bioconda
```
### Run QUAST
```ruby
cd spades
quast.py contigs.fasta scaffolds.fasta
```

## Effect of read correction
Repeat step 2.5 for corrected reads.

```ruby
cd spades/corrected
jellyfish count -m 31 -s 6M -C SRR292678sub_S1_L001_R1_001.00.0_0.cor.fastq
jellyfish histo mer_counts.jf > hist_j_cor.txt
```
Repeat the visualization and calculation of the genome size in R:
```r
hist_k_mer_cor <- read.table("spades_my\\corrected\\hist_j_cor.txt")
plot(hist_k_mer_cor[4:150,],type="l")
points(hist_k_mer_cor[4:150,])

sum(as.numeric(hist_k_mer_cor[1:839,1]*hist_k_mer_cor[1:839,2]))
hist_k_mer_cor[50:70,]
sum(as.numeric(hist_k_mer_cor[1:839,1]*hist_k_mer_cor[1:839,2]))/64
```
In this case total number of k-mer is `329 960 760`, peak position - `64`, and genome size `5 155 220` - `5.16 Gb`.

Compare                   | **Uncorrected reads** | **Corrected reads** 
---------------------------|:---------------------:|:-------------------:
**Total number of k-mers** | 329 960 760           | 329 934 081         
**Peak position**           | 62                    | 64                  
**Peak value**              | 89 747                | 88 759              
**Genome size**             | 5 321 948             | 5 155 220           

## Impact of reads with large insert size

Run SPAdes providing all three libraries: SRR292678 as a paired ends,  SRR292862 and SRR292770 as a mate pairs:
```ruby
spades.py --pe1-1 SRR292678sub_S1_L001_R1_001.fastq --pe1-2 SRR292678sub_S1_L001_R2_001.fastq --mp1-1 SRR292770_S1_L001_R1_001.fastq --mp1-2 SRR292770_S1_L001_R2_001.fastq --mp2-1 SRR292862_S2_L001_R1_001.fastq --mp2-2 SRR292862_S2_L001_R2_001.fastq -o spades_three
```
Run QUAST:
```ruby
cd spades_three
quast.py contigs.fasta scaffolds.fasta
```
## Genome Annotation
