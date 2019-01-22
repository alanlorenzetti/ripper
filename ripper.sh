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
lastupdate=20190122

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

additionalPlots [VARCHAR]:   flag indicating whether additional read depth plots should be created;
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
rm -r 1stQC 2ndQC trimmed sam bam cov positionAnalysis positionAnalysisGenes gccontent correlationAnalysis circos

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

additionalPlots=${9}
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
circosdir="circos"

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
dependencies="curl bedtools R infoseq hisat2 samtools mmr circos"

echo "Checking dependencies..."

# checking programs
for i in $dependencies ; do
    command -v $i > /dev/null >&1 || { echo >&2 "$i isn't installed. Aborting" ; exit 1; }
done 

if [ ! -e /opt/Trimmomatic-0.36/trimmomatic-0.36.jar ] ; then echo >&2 "trimmomatic isn't installed. Aborting" ; exit 1; fi

# checking files
if [ ! -e $miscdir/adap.fa ] ; then echo >&2 "Adapter sequences not found. Aborting" ; exit 1 ; fi
if [ ! -e $rawdir/$control.fastq ] ; then echo >&2 "Control file not found. Aborting" ; exit 1 ; fi

if [ "$positionAnalysis" == "y" ] ; then
    if [ ! -e $miscdir/$spp-ISSaga-checked.gff3 ] ; then echo >&2 "IS annotation not found. Aborting" ; exit 1 ; fi
fi

# checking scripts
scripts="1stQC.R 2ndQC.R correlationAnalysis.R lsm-position.R lsm-positionGenes.R computeGC.sh createCircosFiles.sh circos.conf"

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
    R -q -f $scriptsdir/1stQC.R --args $threads > /dev/null 2>&1

    echo "Done!"
fi

####################################
# trimming with trimmomatic
####################################
# adapter sequence example AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG
if [ ! -d $trimmeddir ] ; then
    echo "Trimming fastq files"
    mkdir $trimmeddir
    for i in $rawdir/*.fastq ; do
        out=$(echo $i | sed "s/^$rawdir/$trimmeddir/g")
        log=$(echo $i | sed "s/^$rawdir/$trimmeddir/g;s/.fastq$/.log/")
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
    R -q -f $scriptsdir/2ndQC.R --args $threads > /dev/null 2>&1

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

    # creating file to store number of alignments
    # for each seq lib
    touch $samdir/readsAligned.txt

    ## aligning to genome
    for i in $trimmeddir/*.fastq ; do
        prefix=$(echo $i | sed "s/^$trimmeddir/$samdir/;s/.fastq$//")
        out=$(echo $i | sed "s/^$trimmeddir/$samdir/;s/fastq$/sam/")
        log=$(echo $i | sed "s/^$trimmeddir/$samdir/;s/fastq$/log/")
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
        uniqReads=`grep "aligned exactly 1 time" $log | sed 's/^    \(.*\) (.*$/\1/'` # extracting number of uniq aligned reads from log
        multiReads=`grep "aligned >1 times" $log | sed 's/^    \(.*\) (.*$/\1/'` # extracting number of multi aligned reads from log
        alnReads=`echo "$uniqReads+$multiReads" | bc` # computing how many reads have been aligned
        printf "$prefix-multi\t$alnReads\n" >> $samdir/readsAligned.txt # storing it to file
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
    for i in $bamdir/*-mmr.bam; do
        samtools sort -@ $threads -o $i $i
    done

    # I will also sort the input bam in order to compare mmr and non-mmr versions
    for i in `ls $bamdir/*.bam | grep -v mmr`; do
        samtools sort -@ $threads -o $i $i
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
        for i in $bamdir/*-mmr.bam ; do
            prefix=$(echo $i | sed "s/^$bamdir/$covdir/;s/-mmr.bam$//")
            bedtools genomecov -ibam $i -d -strand + > $prefix"-counts-rev.txt" # note the inversion of strand on the output file name
            bedtools genomecov -ibam $i -d -strand - > $prefix"-counts-fwd.txt" # note the inversion of strand on the output file name
        done
    else
        for i in $bamdir/*-mmr.bam ; do
            prefix=$(echo $i | sed "s/^$bamdir/$covdir/;s/-mmr.bam$//")
            bedtools genomecov -ibam $i -d -strand - > $prefix"-counts-rev.txt"
            bedtools genomecov -ibam $i -d -strand + > $prefix"-counts-fwd.txt"
        done
    fi

    # creating correction factor for normalization
    maxReads=`awk 'OFS="\t" {if($2 >= n){n = $2}} END {print n}' $samdir/readsAligned.txt` # storing number of maximum aligned reads of a lib
    awk -v maxReads=$maxReads 'OFS="\t" {print $1, $2, maxReads/$2}' $samdir/readsAligned.txt > awk.tmp # creating correction factor file
    mv awk.tmp $samdir/readsAligned.txt # renaming it

    ## CONTROL LIB

    # creating fwd count file for GGB
    for i in $covdir/$control*counts-fwd.txt ; do
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-fwd.txt//")
        corFactor=`grep -w $idxLib $samdir/readsAligned.txt | awk '{print $3}'`
    
        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) >= 1)\
                            {print $1, "+", $2, $3*corFactor}\
                            else\
                            {print $1, "+", $2, "1"}}' $i > $name"-normalized-ggb.txt" 
        if [ "$additionalPlots" == "y" ] ; then
            # absolut count (prepared to compute fold change)
            awk -v corFactor=$corFactor 'OFS="\t" {if(($3) >= 1)\
                                {print $1, "+", $2, $3}\
                                else\
                                {print $1, "+", $2, "1"}}' $i > $name"-absolut-ggb.txt"

            # normalized log2 count
            awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) != 0)\
                                  {print $1, "+", $2, log($3*corFactor)/log(2)}\
                                   else\
                                  {print $1, "+", $2, "0"}}' $i > $name"-normalized-log2-ggb.txt"
        fi
    done

    # creating rev count file for GGB
    for i in $covdir/$control*counts-rev.txt ; do
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-rev.txt//")
        corFactor=`grep -w $idxLib $samdir/readsAligned.txt | awk '{print $3}'`

        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) >= 1)\
                            {print $1, "-", $2, $3*corFactor}\
                            else\
                            {print $1, "-", $2, "1"}}' $i > $name"-normalized-ggb.txt"
        if [ "$additionalPlots" == "y" ] ; then
            # absolut count (prepared to compute fold change)
            awk -v corFactor=$corFactor 'OFS="\t" {if(($3) >= 1)\
                                {print $1, "-", $2, $3}\
                                else\
                                {print $1, "-", $2, "1"}}' $i > $name"-absolut-ggb.txt" 

            # normalized log2 count
            awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) != 0)\
                      {print $1, "-", $2, log($3*corFactor)/log(2)}\
                       else\
                      {print $1, "-", $2, "0"}}' $i > $name"-normalized-log2-ggb.txt"
        fi
    done

    ## NON CONTROL LIBS

    # creating fwd count file for GGB
    nonControl=`ls $covdir/*counts-fwd.txt | grep -v $control`
    for i in $nonControl ; do 
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-fwd.txt//")
        corFactor=`grep -w $idxLib $samdir/readsAligned.txt | awk '{print $3}'`

        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {print $1, "+", $2, $3*corFactor}' $i > $name"-normalized-ggb.txt" 

        if [ "$additionalPlots" == "y" ] ; then
            # absolut count (prepared to compute fold change)
            awk -v corFactor=$corFactor 'OFS="\t" {print $1, "+", $2, $3}' $i > $name"-absolut-ggb.txt" 
    
            # normalized log2 count
            awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) != 0)\
                      {print $1, "+", $2, log($3*corFactor)/log(2)}\
                       else\
                      {print $1, "+", $2, "0"}}' $i > $name"-normalized-log2-ggb.txt"
        fi
    done

    # creating rev count file for GGB
    nonControl=`ls $covdir/*counts-rev.txt | grep -v $control`
    for i in $nonControl ; do
        name=$(echo $i | sed 's/.txt//')
        idxLib=$(echo $i | sed "s/^$covdir\///;s/-counts-rev.txt//")
        corFactor=`grep -w $idxLib $samdir/readsAligned.txt | awk '{print $3}'`
   
        # normalized count (prepared to compute fold change)
        awk -v corFactor=$corFactor 'OFS="\t" {print $1, "-", $2, $3*corFactor}' $i > $name"-normalized-ggb.txt"
    
        if [ "$additionalPlots" == "y" ] ; then
            # absolut count (prepared to compute fold change)
            awk -v corFactor=$corFactor 'OFS="\t" {print $1, "-", $2, $3}' $i > $name"-absolut-ggb.txt"
    
            # normalized log2 count
            awk -v corFactor=$corFactor 'OFS="\t" {if(($3*corFactor) != 0)\
                      {print $1, "-", $2, log($3*corFactor)/log(2)}\
                       else\
                      {print $1, "-", $2, "0"}}' $i > $name"-normalized-log2-ggb.txt"
        fi
    done

    ####################################
    # Computing fold changes
    ####################################
    
    if [ "$additionalPlots" == "y" ] ; then
        # fwd absolut count
        nonControl=`ls $covdir/*-fwd-absolut-ggb.txt | grep -v $control`
        for i in $nonControl ; do
            name=$(echo $i | sed 's/.txt//')
            paste $i $covdir/$control-multi-counts-fwd-absolut-ggb.txt |\
            awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
                          {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
                           else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
                           else {print $1, $2, $3, "0"}}\
                           else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"
        done
    
        # rev absolut count
        nonControl=`ls $covdir/*-rev-absolut-ggb.txt | grep -v $control`
        for i in $nonControl ; do
            name=$(echo $i | sed 's/.txt//')
            paste $i $covdir/$control-multi-counts-rev-absolut-ggb.txt |\
            awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
                          {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
                           else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
                           else {print $1, $2, $3, "0"}}\
                           else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"
        done
    fi

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
        R --slave -q -f $scriptsdir/lsm-position.R --args $spp > /dev/null 2>&1
    fi

    if [ ! -d $positionanalysisgenesdir ] ; then
        mkdir $positionanalysisgenesdir
        bash $scriptsdir/check-genes-with-lsm.sh $spp $annotURL
    fi
fi

####################################
# computing genome-wide GC content
####################################
if [ ! -d $gccontentdir ] ; then
    echo "Computing GC content"
    mkdir $gccontentdir
    bash $scriptsdir/computeGC.sh $miscdir/$spp.fa $windowsize $stepsize
fi

####################################
# correlation analysis
####################################
if [ ! -d $correlationanalysisdir ] ; then
    echo "Performing correlation analysis"
    mkdir $correlationanalysisdir
    R --slave -q -f $scriptsdir/correlationAnalysis.R --args $windowsize $stepsize $threshold > /dev/null 2>&1
fi

####################################
# creating circos ideograms
####################################
if [ ! -d $circosdir ] ; then
    echo "Generating Circos Ideograms"
    mkdir $circosdir
    bash $scriptsdir/createCircosFiles.sh $spp $positionAnalysis
    
    infoseq -only -name -length $miscdir/$spp.fa 2> /dev/null | tail -n +2 |\
    while read replicon length ; do
        name=$(echo $replicon | sed 's/\.[0-9]//')
        unit=${#length}
        if [ $unit -ge 7 ]
        then
            circos -conf $scriptsdir/circos.conf -param karyotype=$circosdir/$name"-karyotype.txt" -outputfile $circosdir/$name -silent
        else
            circos -conf $scriptsdir/circos.conf \
                   -param karyotype=$circosdir/$name"-karyotype.txt" \
                   -param ticks/multiplier="1e-3" \
                   -param ticks/format="%.2f Kb" \
                   -param ticks/tick/spacing="25000u" \
                   -outputfile $circosdir/$name \
                   -silent
        fi
    done
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

additionalPlots [VARCHAR]:   flag indicating whether additional read depth plots should be created;
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
rm -r 1stQC 2ndQC trimmed sam bam cov positionAnalysis positionAnalysisGenes gccontent correlationAnalysis circos

####################################
REQUIRED FILES AND PROGRAMS
####################################

programs:

trimmomatic v0.36 @ /opt/Trimmomatic-0.36/trimmomatic-0.36.jar
hisat2 v1.1.2 @ PATH
samtools v1.3.1 @ PATH
bedtools v2.21.0 @ PATH
R @ PATH ; also it has many library dependencies that should be installed automatically
MMR default version @ PATH
circos v0.69-6 @ PATH

files:

IS annotation contained in misc directory
scripts contained in misc directory
adap fasta adap.fa contained in misc directory
raw data (below) contained in the raw directory
rip fastq rip.fastq (single-end Illumina sequencing)
control fastq control.fastq (single-end Illumina sequencing)

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
    ├── circos.conf
    ├── computeGC.sh
    ├── correlationAnalysis.R
    ├── createCircosFiles.sh
    ├── lsm-positionGenes.R
    └── lsm-position.R

several directories and files will be created during the execution and they will be placed
majorly in "yourDirectory" and "misc"

'
fi
