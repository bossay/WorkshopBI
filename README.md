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
gunzip SRR292678sub_S1_L001_R1_001.fastq.gz SRR292678sub_S1_L001_R2_001.fastq.gz SRR292770_S1_L001_R1_001.fastq.gz SRR292770_S1_L001_R2_001.fastq.gz SRR292862_S2_L001_R1_001.fastq.gz SRR292862_S2_L001_R2_001.fastq.gz
```
### 2.2. Run fastqc

```ruby
fastqc -o ../QC SRR292678sub_S1_L001_R1_001.fastq SRR292678sub_S1_L001_R2_001.fastq SRR292770_S1_L001_R1_001.fastq SRR292770_S1_L001_R2_001.fastq SRR292862_S2_L001_R1_001.fastq SRR292862_S2_L001_R2_001.fastq
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
```
sudo apt-get install jellyfish
```

Run Jellyfish (k-mer sizes of 31, length ~ 5.5M):
```
jellyfish count -m 31 -s 6M -C SRR292678sub_S1_L001_R1_001.fastq
jellyfish histo mer_counts.jf > hist_j.txt
```
### Visualize k-mer distribution
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

### Estimate the genome size 

Let's calculate the size of the genome using the formula: 
`N = (M*L)/(L-K+1)`
`Genome_size = T/N`
> M = 62: k-mer peak
>
> K = 31: k-mer-size
>
> L - 90: avg read length
>
> T = 5499346: total bases
>
> N = (62*90)/(90-31+1) = 93: depth of coverage
>
> **G = 5499346/93 = 59132.75**



## Assembling E. coli X genome from paired reads
