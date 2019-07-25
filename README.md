# Ripper: A tool to process Illumina RIP-Seq/HITS-CLIP data

## Purpose and General Information

This pipeline will perform every task to analyze **prokaryotic** RIP-Seq/HITS-CLIP Illumina data. It was developed to find regions transcribing RNAs that bind to a coimmunoprecipitated protein and check if these regions are within Insertion Sequence (IS) regions.  

A few downstream analyses are also performed, for example: finding the position of interaction regions relative to the annotated IS and genes; computing the correlation between AT-content and interaction region density; creating circular representations of AT-content and interaction regions density for each replicon.  

The pipeline steps can be briefly described as:  

1. Pre-trimming quality control  
2. Trimming  
3. Post-trimming quality control  
4. Downloading reference genome from NCBI RefSeq; building HISAT2 index and aligning to ref. genome  
5. Converting SAM files to BAM and sorting them by read name (performed by SAMtools)  
6. Adjusting position of multi-mappers using MMR  
7. Creating coverage files  
8. Creating annotation of interaction regions  
9. Finding the position of interaction regions relative to annotated insertion sequences and genes  
10. Computing correlation between AT-content and density of interaction regions  
11. Plotting circular representations of AT-content and density of interaction regions for each replicon  

**Warning**  

* This pipeline was designed to work only with Illumina RIP-Seq/HITS-CLIP single-end data.
* It was designed as an in-house tool of analysis, so there are many requirements to be satisfied that need to be installed manually. The program will check all the requirements before start processing.

## Usage

```
bash ripper.sh <spp> <control_lib> <positionAnalysis> <threads> <genomeURL> <readsize> <strandspecific> <invertstrand> <additionalPlots> <windowsize> <stepsize> <log2fcthreshold> <minlength>

spp [VARCHAR]:               prefix used by the script to read and write target files.
                             I advise the use of species name. e.g.: Hsalinarum
                             the hisat2 indexes will be set as misc/Hsalinarum.1.ht2 .
                             by the same means, the IS annotation file should be
                             Hsalinarum-ISSaga-checked.gff3

control_lib [VARCHAR]:       DEPRECATED. (you must assign any character, but it will take no effect).
                             this was an option before, but now is hard coded.
                             e.g.: control

positionAnalysis [VARCHAR]:  flag indicating whether position analysis will be performed
                             for IS and genes; do not use it if you do not have IS annotation files
                             e.g.: y

threads [INT]:               number of threads passed to R in quality control,
                             trimmomatic, hisat2, samtools and mmr
                             e.g.: 8

genomeURL [VARCHAR]:         URL to the NCBI RefSeq FTP.
                             it is used to get the reference genome and annotation files
                             e.g.: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz

readsize [INT]:              mean read size passed to mmr.
                             e.g.: 50

strandspecific [VARCHAR]:    flag indicating whether RNA-Seq data is stranded or not
                             e.g.: y

invertstrand [VARCHAR]:      flag indicating whether RNA-Seq data is from dUTP libraries
                             and should be inverted when computing genome-wide read depth
                             e.g.: y

additionalPlots [VARCHAR]:   DEPRECATED (you must assign any character but it will take no effect).
                             flag indicating whether additional read depth plots should be created;
                             if additionalPlots = y, the program will create additional coverage files
                             with non-normalized counts and normalized log2 counts
                             e.g.: n

windowsize [INT]:            size of windows while computing genome AT content
                             and RNA Binding Protein (RBP) density
                             e.g.: 250

stepsize [INT]:              size of step while computing genome AT content and RBP density
                             e.g.: 250

log2fcthreshold [INT]:       threshold used to identify which regions interact with given RBP;
                             for example using log2fcthreshold = 2, means that the program will
                             consider a base bound to RBP only if a position satisfies
                             rip count/control count >= 4.
                             e.g.: 2

minlength [INT]:             minimum length of transcript fragments interacting with RBP;
                             for example, if minlength = 10, the program will create fragments
                             only if there are at least 10 consecutive bases satisfying log2fcthreshold
                             e.g.: 10
```

e.g.:  

```
bash ripper.sh Hsalinarum control y 8 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz 50 y y n 250 250 2 10
```

if something goes wrong, one may remove all the directories and files created by this script by doing  

```
rm -r 1stQC 2ndQC trimmed sam bam cov positionAnalysis positionAnalysisGenes gccontent correlationAnalysis circlize'
```

## Required programs and files

### Programs:  

* bedtools v2.27.1 @ PATH  
* curl v7.47.0 @ PATH  
* emboss v6.6.0.0  
* hisat2 v2.1.0 @ PATH  
* MMR default version @ PATH (https://github.com/ratschlab/mmr)  
* samtools v1.9 @ PATH  
* trimmomatic v0.36 @ /opt/Trimmomatic-0.36/trimmomatic-0.36.jar  
* R @ PATH (please, check the required packages below)  

### R packages:  

* Rqc (BioConductor)  
* rtracklayer (BioConductor)  
* gplots (CRAN)  
* MSG (CRAN)  
* VennDiagram (CRAN)  
* ggplot2 (CRAN)  
* circlize (CRAN)  

### Directory structure and files:  

```
yourDirectory
├── misc
│   ├── adap.fa
│   └── <spp>-ISSaga-checked.gff3
├── raw
│   ├── control.fastq.gz
│   └── rip.fastq.gz
├── ripper.sh
└── scripts
    ├── 1stQC.R
    ├── 2ndQC.R
    ├── check-genes-with-lsm.sh
    ├── circlize.R
    ├── computeGC.sh
    ├── correlationAnalysis.R
    ├── lsm-positionGenes.R
    └── lsm-position.R
```

**adap.fa** must be a fasta file containing what adapters would look like if they are sequenced. For more information, check https://support.illumina.com/bulletins/2016/12/what-sequences-do-i-use-for-adapter-trimming.html  

The file must look like:  

```
>TruSeq_indexed_adapter
AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC
```

**\<spp\>-ISSaga-checked.gff3** must be a gff3 annotation file containing only insertion sequences. In order to work properly, the attributes field must be exactly the same order as the example below (ID= ; Name= ; rpt_family= ). You can also disable this step by setting positionAnalysis=n.  

```
NC_001869.1 artemis mobile_genetic_element  1718    5019    .   -   .   ID=GenBank:repeat_region:NC_001869:1718:5019.repeat_region;Name=ISH7A;rpt_family=ISNCY
```

several directories and files will be created during the execution and they will be placed majorly in "yourDirectory" and "misc".  


```
                                     __..--"".          .""--..__             
                               _..-``        /\        /\        ``-.._       
                           _.-`           __/\c\      / /\__           `-._   
                        .-`     __..---```    \o\    / /    ```---..__     `-.
                      .`  _.--``               \b\  / /               ``--._  `.
                     / .-`                      \h\/ /                      `-. \
                    /.`                          \c\/                          `.\
                    ´                            /\ \                            `

```

[ASCII Art Credits: David Palmer](http://www.ascii-art.de/ascii/s/scythe.txt)

