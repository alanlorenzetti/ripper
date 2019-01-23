#!/bin/bash

# alorenzetti 

####################################
# Ripper's brief description
####################################
# this script is the scaffold for running RIP-Seq analysis
# and getting results about RNAs derived from Insertion Sequences and Genes
# that interact with this protein in vivo
#
# it relies on many scripts, files and programs, and therefore the setting up
# is quite complicated. Nevertheless, the tool is aware of required files and
# should not work before everything is set up
####################################

version=0.3.0
lastupdate=20190123

# please, check the README.md file before using this script
# there is also a version of the manual on the end of this file
# which can be accessed using ./frtc.sh --help

# starting an if statement to show the help
# which is presented on the end of the file
# this if statement only ends on the end of
# this script
if [ "$1" != "--help" ] ; then

# showing usage hints if no arguments are supplied
if [ $# -ne 13 ] ; then echo "

Ripper:
A tool to process Illumina RIP-Seq/HITS-CLIP data.
Version: $version
Last update: $lastupdate

Usage:

bash ripper.sh <spp> <control_lib> <positionAnalysis> <threads> <genomeURL> <readsize> <strandspecific> <invertstrand> <additionalPlots> <windowsize> <stepsize> <log2fcthreshold> <minlength>

spp [VARCHAR]:               prefix used by the script to read and write target files.
                             I advise the use of species name. e.g.: Hsalinarum
                             the hisat2 indexes will be set as misc/Hsalinarum.1.ebwt .
                             by the same means, the IS annotation file should be
                             Hsalinarum-ISSaga-checked.gff3

control_lib [VARCHAR]:       prefix used to read the control fastq file.
                             e.g.: control

positionAnalysis [VARCHAR]:  flag indicating whether position analysis will be performed
                             for IS and genes; do not use it if you does not have IS annotation files
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

additionalPlots [VARCHAR]:   DEPRECATED (you can assign any character but it will take no effect)
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

e.g.:
bash ripper.sh Hsalinarum control y 8 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz 50 y y n 250 250 2 10

if something goes wrong, one may remove all the directories and files created by this script by doing
rm -r 1stQC 2ndQC trimmed sam bam cov positionAnalysis positionAnalysisGenes gccontent correlationAnalysis circlize

" && exit 1

fi

####################################
# ARGUMENT VARIABLES
####################################
spp=${1}
control=${2}
positionAnalysis=${3} # choose if IS position analysis will be performed
threads=${4} # choose how many threads to be used

genomeURL=${5} # genome url from ncbi ftp
annotURL=${genomeURL/.fna.gz/.gff.gz}

readsize=${6} # mean readsize to input in mmr

if [ ${7} == "y" ] ; then
    strandspecific="-S" ; else
    strandspecific=" "
    hisatstrand="FR"
fi

invertstrand=${8}

if [ $invertstrand == "y" ] ; then
    hisatstrand="RF" ; else
    hisatstrand="FR"
fi

# additionalPlots=${9} # DEPRECATED
additionalPlots="n" # not read from user input anymore
windowsize=${10} # these two parameters for GC content computation
stepsize=${11}

## parameters for finding interaction regions
# threshold=1.58496 is log2(3)
threshold=${12}
# minimum size of interaction region
minsize=${13} 

####################################
# HARD CODED VARIABLES
####################################
rawdir="raw"
miscdir="misc"
scriptsdir="scripts"
firstqcdir="1stQC"
trimmeddir="trimmed"
secondqcdir="2ndQC"
samdir="sam"
bamdir="bam"
covdir="cov"
positionanalysisdir="positionAnalysis"
positionanalysisgenesdir="positionAnalysisGenes"
gccontentdir="gccontent"
correlationanalysisdir="correlationAnalysis"
circlizedir="circlize"

####################################
# PROGRAM STAMP
####################################
echo "Ripper:
A tool to process Illumina RIP-Seq/HITS-CLIP data.
Version: $version
Last update: $lastupdate"

####################################
# CALL AND DATE
####################################
dateAndTime=`date`
echo "$dateAndTime"
echo "Call: $0 $@"

####################################
# CHECKING DEPENDENCIES
####################################

# list of program dependencies
dependencies="curl bedtools R infoseq hisat2 samtools mmr infoseq convert"

echo "Checking dependencies..."

# checking programs
for i in $dependencies ; do
    command -v $i > /dev/null >&1 || { echo >&2 "$i isn't installed. Aborting" ; exit 1; }
done 

if [ ! -e /opt/Trimmomatic-0.36/trimmomatic-0.36.jar ] ; then echo >&2 "trimmomatic isn't installed. Aborting" ; exit 1; fi

# checking R packages
rpackages="Rqc gplots MSG VennDiagram ggplot2 circlize rtracklayer"

for i in $rpackages ; do
    R --slave -e 'args=commandArgs(trailingOnly=T); if(!require(args[1], quietly=T, character.only=T)){quit(save="no", status=1)}else{quit(save="no", status=0)}' \
      --args $i > /dev/null 2>&1
    if [ $? == 1 ] ; then echo >&2 "$i R package is not installed. Aborting" ; exit 1 ; fi
done

# checking files
if [ ! -e $miscdir/adap.fa ] ; then echo >&2 "Adapter sequences not found. Aborting" ; exit 1 ; fi
if [ ! -e $rawdir/$control.fastq.gz ] ; then echo >&2 "Control file not found. Aborting" ; exit 1 ; fi

if [ "$positionAnalysis" == "y" ] ; then
    if [ ! -e $miscdir/$spp-ISSaga-checked.gff3 ] ; then echo >&2 "IS annotation not found. Aborting" ; exit 1 ; fi
fi

# checking scripts
scripts="1stQC.R 2ndQC.R correlationAnalysis.R lsm-position.R lsm-positionGenes.R computeGC.sh circlize.R"

for i in $scripts ; do
    if [ ! -e $scriptsdir/$i ] ; then
        echo >&2 "Missing $i script. Aborting"
        exit 1
    fi
done

echo "Done!"

####################################
# PROCESSING STARTS HERE
####################################

####################################
# checking quality
####################################
if [ ! -d $firstqcdir ] ; then
    echo "First round of quality control"
    mkdir $firstqcdir
    R -q -f $scriptsdir/1stQC.R --args $threads $firstqcdir $rawdir > /dev/null 2>&1

    echo "Done!"
fi

####################################
# trimming with trimmomatic
####################################
# adapter sequence example AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG
if [ ! -d $trimmeddir ] ; then
    echo "Trimming fastq files"
    mkdir $trimmeddir
    for i in $rawdir/*.fastq.gz ; do
        out=$(echo $i | sed "s/^$rawdir/$trimmeddir/g")
        log=$(echo $i | sed "s/^$rawdir/$trimmeddir/g;s/.fastq.gz$/.log/")
        java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads $threads $i $out ILLUMINACLIP:$miscdir/adap.fa:2:30:10 SLIDINGWINDOW:5:28 MINLEN:15 > $log 2>&1
    done

    echo "Done!"
fi

####################################
# checking quality after trimming
####################################

if [ ! -d $secondqcdir ] ; then
    echo "Second round of quality control"
    mkdir $secondqcdir
    R -q -f $scriptsdir/2ndQC.R --args $threads $secondqcdir $trimmeddir > /dev/null 2>&1

    echo "Done!"
fi

####################################
# downloading genome and creating idx
####################################

if [ ! -f $miscdir/$spp.fa ] ; then
    echo >&2 "Downloading genome and creating HISAT2 index"
    curl -u anonymous: $genomeURL 2> /dev/null | zcat > $miscdir/$spp".fa" || { echo >&2 "Genome file download failed. Aborting" ; exit 1; }

    # building index
    hisat2-build $miscdir/$spp".fa" $miscdir/$spp > /dev/null 2>&1
    echo "Done!"
fi

####################################
# aligning to reference genome
####################################
if [ ! -d $samdir ] ; then
    echo "Aligning reads to reference genome"
    mkdir $samdir

    ## aligning to genome
    for i in $trimmeddir/*.fastq.gz ; do
        prefix=$(echo $i | sed "s/^$trimmeddir/$samdir/;s/.fastq.gz$//")
        out=$(echo $i | sed "s/^$trimmeddir/$samdir/;s/fastq.gz$/sam/")
        log=$(echo $i | sed "s/^$trimmeddir/$samdir/;s/fastq.gz$/log/")
        hisat2 --no-spliced-alignment \
                       --no-softclip \
               --rdg 1000,1 \
               --rfg 1000,1 \
                       --rna-strandness $hisatstrand \
                       -k 1000 \
                       -p $threads \
                       -x $miscdir/$spp \
                       -U $i \
                       --summary-file $log > $out 
    done

    echo "Done!"
fi

####################################
# processing bam files
####################################
if [ ! -d $bamdir ] ; then
    echo "Processing BAM files"
    mkdir $bamdir

    # sam2bam ; sort bam
    for i in $samdir/*.sam; do
        out=$(echo $i | sed "s/^$samdir/$bamdir/;s/sam$/bam/")
        log=$(echo $i | sed "s/^$samdir/$bamdir/;s/sam$/log/")
        samtools sort -@ $threads -n -o $out $i
    done

    # to do correction of multi-mapping reads
    # we are going to use MMR
    for i in $bamdir/*.bam; do
        out=$(echo $i | sed 's/.bam$/-mmr.bam/')
        log=$(echo $i | sed 's/.bam$/-mmr.log/')
        mmr -o $out -t $threads $strandspecific -b -R $readsize $i
    done

    # we need to resort bam using default mode in order to proceed
    # we are also indexing them for visualization
    for i in $bamdir/*-mmr.bam; do
        samtools sort -@ $threads -o $i $i
        samtools index -b $i
    done

    # I will also sort and index the input bam in order to compare mmr and non-mmr versions
    for i in `ls $bamdir/*.bam | grep -v mmr`; do
        samtools sort -@ $threads -o $i $i
        samtools index -b $i
    done

    echo "Done!"
fi

####################################
# Creating coverage files
####################################
# creating coverage file for rev and fwd
# ATTENTION HERE:
# for Truseq Stranded mRNA library preparation kit (Hsalinarum library)
# single-end reads are actually the reverse complement of original input RNA
# therefore the strands must be inverted to compute the genome coverage
if [ ! -d $covdir ] ; then
    echo "Creating read depth files"
    mkdir $covdir

    if [ "$invertstrand" == "y" ]
    then
        strandfwd="-"
        strandrev="+"
    else
        strandfwd="+"
        strandrev="-"
    fi

    # note the inversion of strand on the output file name if  invertstrand == y
    for i in $bamdir/*-mmr.bam ; do
      prefix=$(echo $i | sed "s/^$bamdir/$covdir/;s/-mmr.bam$//")
      bedtools genomecov -ibam $i -d -strand $strandrev > $prefix"-counts-rev.txt" 
      bedtools genomecov -ibam $i -d -strand $strandfwd > $prefix"-counts-fwd.txt" 
    done

    # creating correction factor for normalization
    # creating file to store number of alignments for each seq lib
    if [ ! -f $miscdir/readCounts.txt ] ; then
        touch $miscdir/readCounts.txt
    else
        rm $miscdir/readCounts.txt
        touch $miscdir/readCounts.txt
    fi

    for i in $bamdir/*-mmr.bam ; do
        prefix=$(echo $i | sed "s/^$bamdir.//;s/-mmr.bam$//")
        alnReads=`samtools view -@ $threads $i | wc -l`
        echo -e "$prefix\t$alnReads" >> $miscdir/readCounts.txt
    done

    # storing number of maximum aligned reads of all libs
    maxReads=`awk 'OFS="\t" {if($2 >= n){n = $2}} END {print n}' $miscdir/readCounts.txt`
    # creating correction factor file
    awk -v maxReads=$maxReads 'OFS="\t" {print $1, $2, maxReads/$2}' $miscdir/readCounts.txt > awk.tmp
    # renaming it
    mv awk.tmp $miscdir/readCounts.txt

    ## CONTROL LIB

    # creating fwd count file for GGB
    for i in $covdir/$control*counts-fwd.txt ; do
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-fwd.txt//")
        corFactor=`grep -w $idxLib $miscdir/readCounts.txt | awk '{print $3}'`
    
        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) >= 1)\
                            {print $1, "+", $2, $3*corFactor}\
                            else\
                            {print $1, "+", $2, "1"}}' $i > $name"-normalized-ggb.txt"

        # normalized count (bedgraph format)
        bedtools genomecov -ibam $bamdir/$idxLib-mmr.bam -bga -strand $strandfwd -scale $corFactor > $covdir/$idxLib"-counts-fwd-normalized.bedgraph"
    done

    # creating rev count file for GGB
    for i in $covdir/$control*counts-rev.txt ; do
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-rev.txt//")
        corFactor=`grep -w $idxLib $miscdir/readCounts.txt | awk '{print $3}'`

        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) >= 1)\
                            {print $1, "-", $2, $3*corFactor}\
                            else\
                            {print $1, "-", $2, "1"}}' $i > $name"-normalized-ggb.txt"

        # normalized count (bedgraph format)
        bedtools genomecov -ibam $bamdir/$idxLib-mmr.bam -bga -strand $strandrev -scale $corFactor > $covdir/$idxLib"-counts-rev-normalized.bedgraph"

    done

    ## NON CONTROL LIBS

    # creating fwd count file for GGB
    nonControl=`ls $covdir/*counts-fwd.txt | grep -v $control`
    for i in $nonControl ; do 
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-fwd.txt//")
        corFactor=`grep -w $idxLib $miscdir/readCounts.txt | awk '{print $3}'`

        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {print $1, "+", $2, $3*corFactor}' $i > $name"-normalized-ggb.txt" 

        # normalized count (bedgraph format)
        bedtools genomecov -ibam $bamdir/$idxLib-mmr.bam -bga -strand $strandfwd -scale $corFactor > $covdir/$idxLib"-counts-fwd-normalized.bedgraph"
    done

    # creating rev count file for GGB
    nonControl=`ls $covdir/*counts-rev.txt | grep -v $control`
    for i in $nonControl ; do
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-rev.txt//")
        corFactor=`grep -w $idxLib $miscdir/readCounts.txt | awk '{print $3}'`
   
        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {print $1, "-", $2, $3*corFactor}' $i > $name"-normalized-ggb.txt"

        # normalized count (bedgraph format)
        bedtools genomecov -ibam $bamdir/$idxLib-mmr.bam -bga -strand $strandrev -scale $corFactor > $covdir/$idxLib"-counts-rev-normalized.bedgraph"
    done

    ####################################
    # Computing fold changes
    ####################################

    # fwd normalized count
    nonControl=`ls $covdir/*-fwd-normalized-ggb.txt | grep -v $control`
    for i in $nonControl ; do
        name=$(echo $i | sed 's/.txt//')
        paste $i $covdir/$control-counts-fwd-normalized-ggb.txt |\
        awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
                      {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
                       else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
                       else {print $1, $2, $3, "0"}}\
                       else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"

        # converting to bedgraph
        outname=$(echo $name | sed 's/-ggb//')
        awk -v FS="\t" -v OFS="\t" '{print $1,$3-1,$3,$4}' $name"-log2FC.txt" > $outname"-log2FC.bedgraph"

    done

    # rev normalized count
    nonControl=`ls $covdir/*-rev-normalized-ggb.txt | grep -v $control`
    for i in $nonControl ; do
        name=$(echo $i | sed 's/.txt//')
        paste $i $covdir/$control-counts-rev-normalized-ggb.txt |\
        awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
                      {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
                       else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
                       else {print $1, $2, $3, "0"}}\
                       else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"

        # converting to bedgraph
        outname=$(echo $name | sed 's/-ggb//')
        awk -v FS="\t" -v OFS="\t" '{print $1,$3-1,$3,$4}' $name"-log2FC.txt" > $outname"-log2FC.bedgraph"

    done

    echo "Done!"
fi

####################################
# ENTERING NEW PROGRAM: createInTracks.sh 
####################################
echo "Creating regions of interaction"

# this program creates regions of interaction (GFF3)
# based on the log2FC files

## functional variables
flag=0
n=0

# fwd processing
strandname=fwd
for i in $covdir/*fwd-normalized-ggb-log2FC.txt ; do
    sort -Vk1,1 -Vk3,3 $i |\
    awk -v flag=$flag -v minsize=$minsize -v strandname=$strandname -v threshold=$threshold 'OFS="\t" {\
            if($4 >= threshold){\
                if(flag == 0){\
                    oldacc=$1;\
                    oldstrand=$2;\
                    oldpos=$3;\
                    oldscore=$4;\
                    startpos=$3;\
                    flag+=1\
                }else{\
                    if(oldacc == $1 && $3 == oldpos+1){\
                        oldacc=$1;\
                        oldstrand=$2;\
                        oldpos=$3;\
                        oldscore+=$4\
                    }else{\
                        if(oldpos-startpos+1 >= minsize){\
                            meanscore=oldscore/(oldpos-startpos+1);\
                            print oldacc,"in_house","misc_feature",startpos,oldpos,meanscore,oldstrand,".","ID=""Lsm_interaction_"strandname"_"flag;\
                            oldacc=$1;\
                            oldstrand=$2;\
                            oldpos=$3;\
                            oldscore=$4;\
                            startpos=$3;\
                            flag += 1\
                        }else{\
                            oldacc=$1;\
                            oldstrand=$2;\
                            oldpos=$3;\
                            oldscore=$4;\
                            startpos=$3\
                        }
                    }
                }
            }
        }'
done > "interaction-regions-entire-genome-"$strandname".gff3"

# rev processing
strandname=rev
for i in $covdir/*rev-normalized-ggb-log2FC.txt ; do
    sort -Vk1,1 -Vk3,3 $i |\
    awk -v flag=$flag -v minsize=$minsize -v strandname=$strandname -v threshold=$threshold 'OFS="\t" {\
        if($4 >= threshold){\
                if(flag == 0){\
                    oldacc=$1;\
                    oldstrand=$2;\
                    oldpos=$3;\
                    oldscore=$4;\
                    startpos=$3;\
                    flag+=1\
                }else{\
                    if(oldacc == $1 && $3 == oldpos+1){\
                        oldacc=$1;\
                        oldstrand=$2;\
                        oldpos=$3;\
                        oldscore+=$4\
                    }else{\
                        if(oldpos-startpos+1 >= minsize){\
                            meanscore=oldscore/(oldpos-startpos+1);\
                            print oldacc,"in_house","misc_feature",startpos,oldpos,meanscore,oldstrand,".","ID=""Lsm_interaction_"strandname"_"flag;\
                            oldacc=$1;\
                            oldstrand=$2;\
                            oldpos=$3;\
                            oldscore=$4;\
                            startpos=$3;\
                            flag += 1\
                        }else{\
                            oldacc=$1;\
                            oldstrand=$2;\
                            oldpos=$3;\
                            oldscore=$4;\
                            startpos=$3\
                        }
                    }
                }
            }
        }'
done > "interaction-regions-entire-genome-"$strandname".gff3"

echo "Done!"

####################################
# performing position analysis
####################################
# for insertion sequences and genes
if [ "$positionAnalysis" == "y" ] ; then
    echo "Performing position analysis"

    if [ ! -d $positionanalysisdir ] ; then
        mkdir $positionanalysisdir

        # getting intersections between IS file and interaction regions
        # this will reduce time for the next step
        for i in interaction-regions-entire-genome-*.gff3 ; do
            name=$(echo $i | sed s/-entire-genome//)
            bedtools intersect -wa -a $i -b $miscdir/$spp-ISSaga-checked.gff3 | uniq > $name
        done

        # generating heatmaps for lsm position
        R --slave -q -f $scriptsdir/lsm-position.R --args $spp $miscdir $positionanalysisdir > /dev/null 2>&1
    fi

    if [ ! -d $positionanalysisgenesdir ] ; then
        mkdir $positionanalysisgenesdir
        bash $scriptsdir/check-genes-with-lsm.sh $spp $annotURL $miscdir $scriptsdir $positionanalysisgenesdir
    fi

    echo "Done!"
fi

####################################
# computing genome-wide GC content
####################################
if [ ! -d $gccontentdir ] ; then
    echo "Computing GC content"

    mkdir $gccontentdir
    bash $scriptsdir/computeGC.sh $spp $miscdir $windowsize $stepsize $gccontentdir

    echo "Done!"
fi

####################################
# correlation analysis
####################################
if [ ! -d $correlationanalysisdir ] ; then
    echo "Performing correlation analysis"

    mkdir $correlationanalysisdir
    R --slave -q -f $scriptsdir/correlationAnalysis.R --args $windowsize $stepsize $threshold $gccontentdir $covdir $correlationanalysisdir > /dev/null 2>&1

    echo "Done!"
fi

####################################
# creating circlize ideograms
####################################
if [ ! -d $circlizedir ] ; then
    echo "Generating circlize ideograms"

    mkdir $circlizedir

    R --slave -q -f $scriptsdir/circlize.R --args $spp $miscdir $gccontentdir $correlationanalysisdir $circlizedir > /dev/null 2>&1
fi

echo "Done."

# continuing the help if-else statement
else
echo '
####################################
# USAGE MANUAL
####################################

bash ripper.sh <spp> <control_lib> <positionAnalysis> <threads> <genomeURL> <readsize> <strandspecific> <invertstrand> <additionalPlots> <windowsize> <stepsize> <log2fcthreshold> <minlength>

spp [VARCHAR]:               prefix used by the script to read and write target files.
                             I advise the use of species name. e.g.: Hsalinarum
                             the hisat2 indexes will be set as misc/Hsalinarum.1.ebwt .
                             by the same means, the IS annotation file should be
                             Hsalinarum-ISSaga-checked.gff3

control_lib [VARCHAR]:       prefix used to read the control fastq file.
                             e.g.: control

positionAnalysis [VARCHAR]:  flag indicating whether position analysis will be performed
                             for IS and genes; do not use it if you does not have IS annotation files
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

additionalPlots [VARCHAR]:   DEPRECATED (you can assign any character but it will take no effect)
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

e.g.:
bash ripper.sh Hsalinarum control y 8 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz 50 y y n 250 250 2 10

if something goes wrong, one may remove all the directories and files created by this script by doing
rm -r 1stQC 2ndQC trimmed sam bam cov positionAnalysis positionAnalysisGenes gccontent correlationAnalysis circlize

####################################
REQUIRED FILES AND PROGRAMS
####################################

programs:

trimmomatic v0.36 @ /opt/Trimmomatic-0.36/trimmomatic-0.36.jar
hisat2 v1.1.2 @ PATH
samtools v1.3.1 @ PATH
bedtools v2.21.0 @ PATH
MMR default version @ PATH
emboss v6.6.0.0
R @ PATH (please, check the required packages below)

R packages: 

Rqc (BioConductor)
rtracklayer (BioConductor)
gplots (CRAN)
MSG (CRAN)
VennDiagram (CRAN)
ggplot2 (CRAN)
circlize (CRAN)

files:

IS annotation contained in misc directory
scripts contained in misc directory
adap fasta adap.fa contained in misc directory
raw data (below) contained in the raw directory
rip fastq rip.fastq.gz (single-end Illumina sequencing)
control fastq control.fastq.gz (single-end Illumina sequencing)

required directory structure and files:

yourDirectory
├── misc
│   ├── adap.fa
│   └── Hsalinarum-ISSaga-checked.gff3
├── raw
│   ├── <control_lib>.fastq.gz
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

several directories and files will be created during the execution and they will be placed
majorly in "yourDirectory" and "misc"

'

fi
