!!! If you are trying to find my lab logs, it has probably been a very long time since the deadline. I have moved all the lab logs to another repository - [BI_2022](https://github.com/bossay/BI_2022)

# Project #4 Tardigrades: from genestealers to space marines

## 1. Prepare the workspace

### 1.1. Create a directory for files and results

```ruby
mkdir WorkshpoBI/Project4/raw_data
cd WorkshpoBI/Project4/raw_data
```

### 1.2. Download raw data

For this project we will be using a sequence of the *Ramazzottius varieornatus*, the YOKOZUNA-1 strain

Download the sequences:

```ruby
wget ftp.ncbi.nlm.nih.gov/genomes/all/GCA/001/949/185/GCA_001949185.1_Rvar_4.0/GCA_001949185.1_Rvar_4.0_genomic.fna.gz
```
## 2. Structural annotation

### 2.1. Installing AUGUSTUS

```ruby
sudo apt install augustus augustus-data augustus-doc
```
More about [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus) in the project's repository.

### 2.2. Runing AUGUSTUS

```ruby
augustus --species=nasonia GCA_001949185.1_Rvar_4.0_genomic.fna --progress=true > augustus.whole.gff
```

### 2.3. Extraction protein sequences

We extract protein sequences (fasta) from prediction results (usually GFF/GTF) using the getAnnoFasta.pl script.

Download this script and make it executable.

```ruby
wget http://augustus.gobics.de/binaries/scripts/getAnnoFasta.pl
chmod +x getAnnoFasta.pl
```

Count the number of obtained proteins:

```ruby
grep '>' augustus.whole.aa | wc -l
```
> 16435

## 3. Physical localization

Tandem mass spectrometry can be used to analyze the chromatin fraction and provide a list of peptides that have been linked to DNA:

|                       |                        |                            |                   |
|-----------------------|------------------------|----------------------------|-------------------|
| >1 VFSALDVLR          | >11 LTDDELNQAMK        | >21 TNAAGQLIGPGGIAINDAGLTR | >41 GPEDEAGGPPK   |
| >2 SPLQTDEIR          | >12 AATGAVQSSASK       | >22 VFSALDVLR              | >42 SSGGNSPDVVVR  |
| >3 AATGAVLSSADSR      | >13 VLLDNQDDYELK       | >23 SPLQTDEIR              | >43 STVGQTPPQNLQR |
| >4 LLDDENDYELK        | >14 VFSALDVLR          | >24 QFQDYSNSR              |                   |
| >5 VTGSSQGAINQQQAK    | >15 EIYANAQPGK         | >25 AATGAVLSSADSR          |                   |
| >6 VTGSSAQQIDINQAR    | >16 QIMDNEGLR          | >26 VFSALDVLR              |                   |
| >7 LTSSGTGAGSAPAAAK   | >17 SAPLSASPISAR       | >27 EIYANAQPGK             |                   |
| >8 NIPVGGVNTEATGDNYIR | >18 SSSSSSQESQSSSSQVR  | >28 QFQDYSNSR              |                   |
| >9 QGGMGMSGGMGGADR    | >19 STVGQTPPQNLQR      | >29 NIPVGGVNTEATGDNYIR     |                   |
| >10 GPEDEAGGPPK       | >20 APVNGIVTDANGNQIQVR | >30 AATGAVQSSASK           |                   |

You can download this list [here](https://disk.yandex.ru/d/xJqQMGX77Xueqg).

### 3.1. Installing BLAST+

```ruby
conda activate
mamba install -c bioconda blast
```

### 3.2. Runing BLAST+

Creating database:

```ruby
makeblastdb -in augustus.whole.aa -dbtype prot -out tardigrade_db
```
Blastp with [outfmt6 output](https://www.metagenomics.wiki/tools/blast/blastn-output-format-6). Added additional column with coverage:

```ruby
blastp -db tardigrade_db -query peptides.fa -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qcovs" > blastp_of_tardi_proteins.txt
```

Parsing for whole protein sequences:

```ruby
xargs samtools faidx augustus.whole.aa < names.txt > proteins.faa
```
> `names.txt` - file with the names of unique proteins

## 4. Localization prediction

### 4.1. WoLF PSORT

Submit our set of proteins to the [WoLF PSORT](https://wolfpsort.hgc.jp/) server:

> g10513.t1 details nucl: 20, cyto_nucl: 14.5, cyto: 7, extr: 3, E.R.: 1, golg: 1
>
>g10514.t1 details nucl: 19, cyto_nucl: 15, cyto: 9, extr: 3, mito: 1
>
>g11806.t1 details nucl: 18, cyto_nucl: 11.8333, mito: 5, extr: 4, cyto: 3.5, cyto_pero: 2.66667, cysk_plas: 1
>
>g11960.t1 details nucl: 32
>
>g14472.t1 details nucl: 28, plas: 2, cyto: 1, cysk: 1
>
>g15484.t1 details nucl: 17.5, cyto_nucl: 15.3333, cyto: 12, cyto_mito: 6.83333, plas: 1, golg: 1
>
>g16318.t1 details nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1
>
>g16368.t1 details nucl: 20.5, cyto_nucl: 13, extr: 5, cyto: 4.5, E.R.: 1, golg: 1
>
>g5927.t1 details nucl: 30.5, cyto_nucl: 16.5, cyto: 1.5
>
>g7861.t1 details nucl: 16, cyto_nucl: 14, cyto: 8, plas: 5, pero: 1, cysk: 1, golg: 1
>
>g8100.t1 details nucl: 16.5, cyto_nucl: 12.5, cyto: 7.5, plas: 5, extr: 2, E.R.: 1
>
>g8312.t1 details nucl: 15.5, cyto_nucl: 15.5, cyto: 12.5, mito: 2, plas: 1, golg: 1

Parsing protein sequences localized in the nucleus

```ruby
xargs samtools faidx augustus.whole.aa < nucl_prot.txt > proteins_nucleo.faa
```

### 4.2. TargetP 1.1 Server

Submit our set to the [TargetP server](https://services.healthtech.dtu.dk/service.php?TargetP-2.0). For several proteins were found signal peptides in their sequences. We are not interested in secretory proteins, we will remove them from further analysis.

> Select Organism group: Non-plant

| **# ID**  | **Prediction** | **OTHER** | **SP**   | **mTP**  |
|-----------|----------------|-----------|----------|----------|
| g10513.t1 | OTHER          | 0.999999  | 0.000001 | 0.000000 |
| g10514.t1 | OTHER          | 0.999543  | 0.000349 | 0.000107 |
| g11513.t1 | OTHER          | 0.999434  | 0.000401 | 0.000164 |
| g11806.t1 | OTHER          | 0.998977  | 0.000887 | 0.000136 |
| g11960.t1 | OTHER          | 0.999996  | 0.000002 | 0.000002 |
| g12510.t1 | OTHER          | 0.999738  | 0.000099 | 0.000163 |
| g14472.t1 | OTHER          | 0.999999  | 0.000001 | 0.000000 |
| g15484.t1 | OTHER          | 0.999980  | 0.000010 | 0.000010 |
| g16318.t1 | OTHER          | 0.997047  | 0.002953 | 0.000000 |
| g16368.t1 | OTHER          | 0.996693  | 0.003307 | 0.000000 |
| g2203.t1  | OTHER          | 0.999869  | 0.000031 | 0.000100 |
| g3428.t1  | OTHER          | 0.999903  | 0.000033 | 0.000064 |
| g4106.t1  | OTHER          | 0.729658  | 0.266917 | 0.003425 |
| g4970.t1  | OTHER          | 0.999996  | 0.000003 | 0.000001 |
| g5237.t1  | OTHER          | 0.999545  | 0.000345 | 0.000111 |
| g5443.t1  | OTHER          | 0.952853  | 0.043784 | 0.003363 |
| g5510.t1  | OTHER          | 0.999108  | 0.000016 | 0.000876 |
| g5927.t1  | OTHER          | 0.999995  | 0.000001 | 0.000004 |
| g7861.t1  | OTHER          | 0.999975  | 0.000004 | 0.000022 |
| g8100.t1  | OTHER          | 0.999955  | 0.000024 | 0.000021 |
| g8312.t1  | OTHER          | 0.999930  | 0.000065 | 0.000004 |


## 5. BLAST search

| **Protein**   | **Description**                                        | **Query coverage** | **E-value** | **Ident** | **Accession** |
|---------------|--------------------------------------------------------|--------------------|-------------|-----------|----------------|
| **g10513.t1** |                                                        |                    |             |           |                |
| **g10514.t1** |                                                        |                    |             |           |                |
| **g11806.t1** |                                                        |                    |             |           |                |
| **g11960.t1** | E3 ubiquitin-protein ligase BRE1B                      | 0,96               | 0           | 26.96%    | Q8CJB9.1       |
| **g14472.t1** |                                                        |                    |             |           |                |
| **g15484.t1** | Vacuolar protein sorting-associated protein 51 homolog | 78%                | 0.0         | 45.03%    | Q155U0.1       |
| **g16318.t1** | Eukaryotic translation initiation factor 3 subunit A   | 40%                | 4,00E-08    | 36.11%    | A2VD00.1       |
| **g16368.t1** | Eukaryotic translation initiation factor 3 subunit A   | 40%                | 4,00E-08    | 36.11%    | A2VD00.1       |
| **g5927.t1**  | Glucosamine 6-phosphate N-acetyltransferase            | 0,14               | 0           | 38.64%    | Q17427.1       |
| **g7861.t1**  | Inositol monophosphatase 2                             | 0,99               | 0           | 37.21%    | B4F769.1       |
| **g8100.t1**  | Inositol monophosphatase 3                             | 0,22               | 0           | 36.04%    | Q2YDR3.1       |
| **g8312.t1**  | Vacuolar protein sorting-associated protein 41         | 0,84               | 0.0         | 40.84%    | Q5KU39.1       |

## 6. Pfam prediction

We can predict the function of proteins even if we cannot find orthologous sequences in the databases using Blast. 
Let's use [HMMER](https://www.ebi.ac.uk/Tools/hmmer/).

+ select `Search`
+ select `hmmscan` tool
+ select the `Pfam` database
+ insert protein

| **Protein**   | **Id**        | **Accession** | **Clan** | **Description**                                           | **Ind.** | **Cond.** |
|---------------|---------------|---------------|----------|-----------------------------------------------------------|----------|-----------|
| **g10513.t1** | None          |               |          |                                                           |          |           |
| **g10514.t1** | None          |               |          |                                                           |          |           |
| **g11806.t1** | None          |               |          |                                                           |          |           |
| **g11960.t1** | zf-C3HC4      | PF00097.28    | CL0229   | Zinc finger, C3HC4 type (RING finger)                     | 4.2e-05  | 2.1e-09   |
| **g14472.t1** | None          |               |          |                                                           |          |           |
| **g15484.t1** | Vps51         | PF08700.14    | CL0295   | Vps51/Vps67                                               | 1.3e-23  | 3.2e-27   |
| **g16318.t1** | None          |               |          |                                                           |          |           |
| **g16368.t1** | None          |               |          |                                                           |          |           |
| **g5927.t1**  | None          |               |          |                                                           |          |           |
| **g7861.t1**  | SNF2-rel_dom  | PF00176.26    | CL0023   | SNF2-related domain                                       | 1.2e-28  | 1.8e-32   |
| **g7861.t1**  | HARP          | PF07443.16    | n/a      | HepA-related protein (HARP)                               | 2.6e-10  | 4.0e-14   |
| **g8100.t1**  | Inositol_P    | PF00459.28    | CL0171   | Inositol monophosphatase family                           | 1.9e-37  | 2.0e-41   |
| **g8100.t1**  | MKLP1_Arf_bdg | PF16540.8     | n/a      | Arf6-interacting domain of mitotic kinesin-like protein 1 | 5.1e-27  | 5.2e-31   |
| **g8312.t1**  | Clathrin      | PF00637.23    | CL0020   | Region in Clathrin and VPS                                | 5.4e-23  | 2.7e-27   |

## 7. Integrate pieces of evidence

| **Protein**   | **Best BLAST hit**                                     | **E value** | **Relative organism**  | **Pham domains**           | **Probable localization (WoLF PSORT)** | **Probable localization (TargetP)** |
|---------------|--------------------------------------------------------|-------------|------------------------|----------------------------|----------------------------------------|-------------------------------------|
| **g10513.t1** | None                                                   |             |                        |                            | Nuclear                                | Other                               |
| **g10514.t1** | None                                                   |             |                        |                            | Nuclear                                | Other                               |
| **g11806.t1** | None                                                   |             |                        |                            | Nuclear                                | Other                               |
| **g11960.t1** | E3 ubiquitin-protein ligase BRE1B                      | 0           | Rattus norvegicus      | zf-C3HC4                   | Nuclear                                | Other                               |
| **g14472.t1** | None                                                   |             |                        |                            | Nuclear                                | Other                               |
| **g15484.t1** | Vacuolar protein sorting-associated protein 51 homolog | 0.0         | Danio rerio            | Vps51                      | Nuclear                                | Other                               |
| **g16318.t1** | Eukaryotic translation initiation factor 3 subunit A   | 4,00E-08    | Xenopus laevis         |                            | Nuclear                                | Other                               |
| **g16368.t1** | Eukaryotic translation initiation factor 3 subunit A   | 4,00E-08    | Xenopus laevis         |                            | Nuclear                                | Other                               |
| **g5927.t1**  | Glucosamine 6-phosphate N-acetyltransferase            | 0           | Caenorhabditis elegans |                            | Nuclear                                | Other                               |
| **g7861.t1**  | Inositol monophosphatase 2                             | 0           | Rattus norvegicus      | SNF2-rel_dom, HARP         | Nuclear                                | Other                               |
| **g8100.t1**  | Nuclear                                                | Other       |                        | Inositol_P, MKLP1_Arf_bdg  | Nuclear                                | Other                               |
| **g8312.t1**  | Inositol monophosphatase 3                             | 0           | Danio rerio            | Clathrin                   | Nuclear                                | Other                               |
