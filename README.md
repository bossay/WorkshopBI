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
| **Poor quality**    | 0                 | 0                 | 0                 | 0                 | 0                 | 0                 |


### 2.5. K-mer profile and genome size estimation

Install Jellyfish: a fast k-mer counter
```ruby
sudo apt-get install jellyfish
```

Run Jellyfish (k-mer sizes of 31, length ~ 5.5M):
```ruby
jellyfish count -m 31 -s 6M -C SRR292678sub_S1_L001_R1_001.fastq SRR292678sub_S1_L001_R2_001.fastq -o SRR292678_v2.jf
jellyfish histo SRR292678_v2.jf > hist_j_v2.txt
```
### 2.6. Visualize k-mer distribution
Genome size can be calculated by counting k-mer frequency of the read data. Using the Jellyfish output we plot the graph with R:

```r
setwd('C:\\Users\\pmbog\\source\\WorkshpoBI\\Project3\\raw_data')
hist_k_mer <- read.table("hist_j_v2.txt")

plot(hist_k_mer[4:400,],type="l", xlab = "Unique k-mers", ylab = "Count")
points(hist_k_mer[4:400,])
```
Calculate the total number of k-mer in the distribution:

```r
sum(as.numeric(hist_k_mer[1:1569,1]*hist_k_mer[1:1569,2]))
```
In this case there are about `659 921 520` k-mers in the histogram.

Next, we want to know the peak position. From the graph, we can see its close to 120-130. Thus we examine the number close to 120-130 and find the maximum value
```
hist_k_mer[115:135,]
```
| **Count** | **Unique k-mers** |
|-----------|-------------------|
| ...       | ...               |
| 123       | 46072             |
| 124       | 45509             |
| **125**   | **46189**         |
| 126       | 45664             |
| 127       | 45458             |
| ...       | ...               |

In this case, the peak is at 125. Then, the genome size can be estimated as:

```
sum(as.numeric(hist_k_mer[1:1569,1]*hist_k_mer[1:1569,2]))/125
```
This reads as `5 279 372` - `5.28 Gb`.

### 2.7. Estimate the genome size 

Let's calculate the size of the genome using the formula: 
`N = (M×L)/(L-K+1)`
`Genome_size = T/N`
> M = 125: k-mer peak
>
> K = 31: k-mer-size
>
> L - 90: avg read length
>
> T = 5499346×2×90 = 989882280: total bases
>
> N = (12590)/(90-31+1) = 187.5: depth of coverage
>
> **G = 989882280/187.5 = 5279372**

## 3. Assembling E. coli X genome from paired reads

### 3.1. Installing SPAdes via conda

```ruby
conda install spades -c bioconda
```

To verify that the software was installed correctly, run it in the test mode:
```ruby
spades.py --test
```
If you do everything right, the last line of the output will be `Thank you for using SPAdes!`. More about [SPAdes](https://cab.spbu.ru/software/spades/)

### 3.2. Run SPAdes
Run SPAdes in the paired-end mode, providing paired reads of E. coli X. from the library SRR292678 (forward and reverse):
 ```ruby
spades.py -1 SRR292678sub_S1_L001_R1_001.fastq -2 SRR292678sub_S1_L001_R2_001.fastq -o spades
```
### 3.3. Installing QUAST via conda

```ruby
conda install quast -c bioconda
```
### 3.4. Run QUAST
```ruby
cd spades
quast.py contigs.fasta scaffolds.fasta
```

## 4. Effect of read correction
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

Compare                    | **Uncorrected reads** | **Corrected reads** 
---------------------------|:---------------------:|:-------------------:
**Total number of k-mers** | 659 921 520           | 329 934 081         
**Peak position**          | 125                   | 64                  
**Peak value**             | 46 189                | 88 759              
**Genome size**            | 5 279 372             | 5 155 220           

## 5. Impact of reads with large insert size

Run SPAdes providing all three libraries: SRR292678 as a paired ends,  SRR292862 and SRR292770 as a mate pairs:
```ruby
spades.py --pe1-1 SRR292678sub_S1_L001_R1_001.fastq --pe1-2 SRR292678sub_S1_L001_R2_001.fastq --mp1-1 SRR292770_S1_L001_R1_001.fastq --mp1-2 SRR292770_S1_L001_R2_001.fastq --mp2-1 SRR292862_S2_L001_R1_001.fastq --mp2-2 SRR292862_S2_L001_R2_001.fastq -o spades_three
```
Run QUAST:
```ruby
cd spades_three
quast.py contigs.fasta scaffolds.fasta
```
## 6. Genome Annotation

### 6.1. Installing Prokka
More about [Prokka](https://github.com/tseemann/prokka)
```ruby
 conda create -n prokka_env -c conda-forge -c bioconda prokka
 conda activate prokka_env
```
To check that the installation was successful
```ruby
prokka --version
```
That's right, if you see `prokka 1.14.6`

### 6.2. Run Prokka
Run Prokka on the `scaffolds.fasta` file from the SPAdes output with default parameters. 
```ruby
prokka --outdir ./prokka --centre X --compliant spades_three/scaffolds.fasta
```
## 7. Finding the closest relative of E. coli X
Our goal is to find the known genome most similar to the pathogenic strain (and to infer the properties of E. coli X from it). An efficient approach is to select one important and evolutionarily conserved gene to compare with all other sequenced genomes. The gene we will use is 16S ribosomal RNA.

To find the 16S rRNA in the collected E. coli X genome we will use the `Barrnap` rRNA gene prediction tool. 

### 7.1. Installing Barrnap via conda
More about [Barrnap](https://github.com/tseemann/barrnap)
```ruby
conda install -c bioconda -c conda-forge barrnap
```
### 7.2. Run Barrnap
 
```ruby
barrnap -o rrna.fa < ./spades_three/contigs.fasta > rrna.gff
head -n 3 rrna.fa
```
>16S_rRNA::NODE_14_length_114134_cov_78.973238:45-1583(-)
TTGAAGAGTTTGATCATGGCTCAGATTGAACGCTGGCGGCAGGCCTAACACATGCAAGTCGAACGGTAACAGGAAACAGCTTGCTGTTTCGCTGACGAGTGGCGGACGGGTGAGTAATGTCTGGGAAACTGCCTGATGGAGGGGGATAACTACTGGAAACGGTAGCTAATACCGCATAACGTCGCAAGACCAAAGAGGGGGACCTTCGGGCCTCTTGCCATCGGATGTGCCCAGATGGGATTAGCTTGTTGGTGGGGTAACGGCTCACCAAGGCGACGATCCCTAGCTGGTCTGAGAGGATGACCAGCCACACTGGAACTGAGACACGGTCCAGACTCCTACGGGAGGCAGCAGTGGGGAATATTGCACAATGGGCGCAAGCCTGATGCAGCCATGCCGCGTGTATGAAGAAGGCCTTCGGGTTGTAAAGTACTTTCAGCGGGGAGGAAGGGAGTAAAGTTAATACCTTTGCTCATTGACGTTACCCGCAGAAGAAGCACCGGCTAACTCCGTGCCAGCAGCCGCGGTAATACGGAGGGTGCAAGCGTTAATCGGAATTACTGGGCGTAAAGCGCACGCAGGCGGTTTGTTAAGTCAGATGTGAAATCCCCGGGCTCAACCTGGGAACTGCATCTGATACTGGCAAGCTTGAGTCTCGTAGAGGGGGGTAGAATTCCAGGTGTAGCGGTGAAATGCGTAGAGATCTGGAGGAATACCGGTGGCGAAGGCGGCCCCCTGGACGAAGACTGACGCTCAGGTGCGAAAGCGTGGGGAGCAAACAGGATTAGATACCCTGGTAGTCCACGCCGTAAACGATGTCGACTTGGAGGTTGTGCCCTTGAGGCGTGGCTTCCGGAGCTAACGCGTTAAGTCGACCGCCTGGGGAGTACGGCCGCAAGGTTAAAACTCAAATGAATTGACGGGGGCCCGCACAAGCGGTGGAGCATGTGGTTTAATTCGATGCAACGCGAAGAACCTTACCTGGTCTTGACATCCACGGAAGTTTTCAGAGATGAGAATGTGCCTTCGGGAACCGTGAGACAGGTGCTGCATGGCTGTCGTCAGCTCGTGTTGTGAAATGTTGGGTTAAGTCCCGCAACGAGCGCAACCCTTATCCTTTGTTGCCAGCGGTCCGGCCGGGAACTCAAAGGAGACTGCCAGTGATAAACTGGAGGAAGGTGGGGATGACGTCAAGTCATCATGGCCCTTACGACCAGGGCTACACACGTGCTACAATGGCGCATACAAAGAGAAGCGACCTCGCGAGAGCAAGCGGACCTCATAAAGTGCGTCGTAGTCCGGATTGGAGTCTGCAACTCGACTCCATGAAGTCGGAATCGCTAGTAATCGTGGATCAGAATGCCACGGTGAATACGTTCCCGGGCCTTGTACACACCGCCCGTCACACCATGGGAGTGGGTTGCAAAAGAAGTAGGTAGCTTAACCTTCGGGAGGGCGCTTACCACTTTGTGATTCATGACTGGGGTGAAGTCGTAACAAGGTAACCGTAGGGGAACCTGCGGTTGGATCACCTCCTT
>16S_rRNA::NODE_21_length_75978_cov_83.376460:45-1583(-)

### 7.3. Search for the relative genome in the RefSeq
We will now use `BLAST` to search for the genome in the `RefSeq` database with 16S rRNA that is most similar to the 16S rRNA that we just found.

+ Open the  [NCBI BLAST homepage](http://blast.ncbi.nlm.nih.gov)
+ Select `Nucleotide blast`
+ Select `Reference Genome Database (refseq_genomes)` in the `Database` field
+ Select `Escherichia coli` in the `Organism` field
+ Set the time range using parameter `PDAT` in the `Entrez Query` field `1900/01/01:2011/01/01[PDAT]`
+ Other parameters should be specified as default

In the results, we found probably the closest relative of E. coli X, which can be used as a reference genome.

>[Escherichia coli 55989, Sequence ID: NC_011748.1](https://www.ncbi.nlm.nih.gov/nucleotide/NC_011748.1?report=genbank&log$=nucltop&blast_rank=1&RID=S6AWDH15013)

Save the found genome in `.fasta` format.

## 8. What is the genetic cause of HUS?

### 8.1. Installing Mauve

Download the [Mauve](https://darlinglab.org/mauve/download.html) version for your system. You may also need to install [Java x64](https://www.java.com/en/download/manual.jsp).

### 8.2. Compare the E. coli X with the reference genome

+ Open `Mauve`
+ Select `File` → `Align with progressiveMauve...`
+ Press `Add sequences`
+ Select the reference genome `NC_011748_1.fasta`, then the annotated E. coli X genome (`scaffolds.gbk` or `scaffolds.gbf`, depending on version)
+ Start the alignment

Look for specific genes using `Sequence Navigator`. Select `Product` in the left window and enter the name of the desired gene in the right window.

We found 2 shiga toxin-related genes stxA and stxB.

## 9. Tracing the source of toxin genes in E. coli X
We have discovered that the E. coli X genome encodes Shiga-like toxin genes (stxA, stxB). Now let's figure out how this strain has acquired these weapons.

Consider the genes that are located next to the Shiga-type toxin genes. Let's use [BLAST](http://blast.ncbi.nlm.nih.gov) to find which organism they belong to. Most likely our bacterium acquired the toxin genes from this organism.

>Escherichia coli O103:H2 str. 12009
>
>Escherichia coli O157:H7 str. EC508 ???
>
>Escherichia coli S88 ???

## 10. Antibiotic resistance detection
To search for genes responsible for antibiotic resistance, we will use [ResFinder](https://cge.food.dtu.dk/services/ResFinder/), which specifically searches a database of genes implicated in antibiotic resistance, identifying similarities between the sequenced genome and this database using local alignment.

+ Visit the [ResFinder](https://cge.food.dtu.dk/services/ResFinder/) homepage
+ Upload the `scaffolds.fasta` file from the SPAdes output
+ In the field `Select Antimicrobial configuration` select `All`

For comparison, do the same for the reference strain.

| **Antimicrobial**            | **Class**                   | **E. coli X** | **E. coli 55989** |
|------------------------------|------------------------------|---------------|--------------------|
| **Amoxicillin**              | Beta-lactam                  | Resistant     | No resistance      |
| **Ampicillin**               | Beta-lactam                  | Resistant     | No resistance      |
| **Aztreonam**                | Beta-lactam                  | Resistant     | No resistance      |
| **Benzylkonium chloride**    | Quaternary ammonium compound | Resistant     | No resistance      |
| **Cefepime**                 | Beta-lactam                  | Resistant     | No resistance      |
| **Cefotaxime**               | Beta-lactam                  | Resistant     | No resistance      |
| **Ceftazidime**              | Beta-lactam                  | Resistant     | No resistance      |
| **Ceftriaxone**              | Beta-lactam                  | Resistant     | No resistance      |
| **Cephalothin**              | Beta-lactam                  | Resistant     | No resistance      |
| **Cetylpyridinium chloride** | Quaternary ammonium compound | Resistant     | No resistance      |
| **Chlorhexidine**            | Quaternary ammonium compound | Resistant     | No resistance      |
| **Doxycycline**              | Tetracycline                 | Resistant     | Resistant          |
| **Ethidium bromide**         | Quaternary ammonium compound | Resistant     | No resistance      |
| **Minocycline**              | Tetracycline                 | No resistance | Resistant          |
| **Piperacillin**             | Beta-lactam                  | Resistant     | No resistance      |
| **Streptomycin**             | Aminoglycoside               | Resistant     | No resistance      |
| **Sulfamethoxazole**         | Folate pathway antagonist    | Resistant     | No resistance      |
| **Tetracycline**             | Tetracycline                 | Resistant     | Resistant          |
| **Ticarcillin**              | Beta-lactam                  | Resistant     | No resistance      |
| **Trimethoprim**             | Folate pathway antagonist    | Resistant     | No resistance      |

## 11. Antibiotic resistance mechanism

We have discovered that the E. coli X genome encodes Shiga-like toxin genes (stxA, stxB). Now let's figure out how this strain has acquired these weapons.

We found that E. coli X acquired additional resistance to antibiotics. Some bacteria acquire specific enzymes that can degrade the antibiotic molecule itself. In our case, one of these enzymes is β-lactamase.
Let's look at the genes that are located next to the β-lactamase genes. Use [BLAST](http://blast.ncbi.nlm.nih.gov) to find which organism they belong to.

>Escherichia coli MS 182-1
>
>Escherichia coli E22

## Ready! You`re incredible!
