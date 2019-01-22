#!/bin/bash

####################
# alorenzetti 201708
####################

####################
# Brief description
####################
# this script is the scaffold for running RIP-Seq analysis
# and getting results about RNAs derived from Insertion Sequences and Genes
# that interact with this protein in vivo
#
# it relies on many scripts, files and programs, and therefore the setting up
# is quite complicated. Nevertheless, the tool is aware of required files and
# should not work before everything is set up
#
# alpha version 0.2
####################

####################
# Usage
####################
# ./workflow.sh <spp> <control_lib> <positionAnalysis> <threads> <genomeURL> <readsize> <strandspecific> <invertstrand> <additionalPlots> <windowsize> <stepsize> <log2fcthreshold> <minlength>
# e.g.:
# ./workflow.sh Hsalinarum control y 20 ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/006/805/GCF_000006805.1_ASM680v1/GCF_000006805.1_ASM680v1_genomic.fna.gz 50 y y n 250 250 2 10
#
# if something goes wrong, one may remove all the directories and files created by this script by doing
# rm -r 1stQC 2ndQC trimmed sam bam cov positionAnalysis positionAnalysisGenes gccontent correlationAnalysis circos
#
# brief description of every argument:
#
# <spp> [CHAR] prefix used by the script to read and write target files; I advise the use of species name
# 	       for example: if spp=Hsalinarum, the bowtie indexes will be set as misc/Hsalinarum.1.ebwt;
#	       by the same means, the IS annotation file should be misc/Hsalinarum-ISSaga-checked.gff3
#
# <control_lib> [CHAR] prefix used to read the control fastq file. for example raw/control.fastq
#
# <positionAnalysis> [y/n] flag indicating whether position analysis will be performed for IS and genes
#		           don't use it if you doesn't have IS annotation files
#
# <threads> [INT] number of threads passed to R in quality control, trimmomatic, bowtie and mmr
#
# <genomeURL> [CHAR] URL to the NCBI RefSeq FTP; it's used to get the reference genome file and annotation file
#
# <readsize> [INT] mean read size passed to mmr
#
# <strandspecific> [y/n] flag indicating whether RNA-Seq data is stranded or not; required by mmr
#
# <invertstrand> [y/n] flag indicating whether RNA-Seq data is from dUTP libraries and should be inverted when
#                      computing genome-wide read depth
#
# <additionalPlots> [y/n] flag indicating whether additional read depth plots should be created;
#			  if additionalPlots = y, the program will create additional coverage files with
#			  non-normalized counts and normalized log2 counts
#
# <windowsize> [INT] size of windows while computing genome AT content and RNA Binding Protein (RBP) density
# 
# <stepsize> [INT] size of step while computing genome AT content and RBP density
#
# <log2fcthreshold> [INT] threshold used to identify which regions interact with given RBP;
#			  for example using log2fcthreshold = 2, means that the program will
# 			  consider a base bound to RBP only if a position satisfies rip count/control count >= 4
#
# <minlength> [INT] minimum length of transcript fragments interacting with RBP;
#		    for example, if minlength = 10, the program will create fragments only if there
#		    are at least 10 consecutive bases satisfying log2fcthreshold
####################

####################
# Required files and programs
####################
#
# required programs:
# trimmomatic v0.36 @ /opt/Trimmomatic-0.36/trimmomatic-0.36.jar
# bowtie v1.1.2 @ PATH
# samtools v1.3.1 @ PATH
# bedtools v2.21.0 @ PATH
# R @ PATH ; also it has many library dependencies that should be installed automatically
# MMR @ PATH
# circos v0.69-6 @ PATH
#
# required files:
# IS annotation contained in misc directory
# scripts contained in misc directory
# adap fasta adap.fa contained in misc directory
# raw data (below) contained in the raw directory
# rip fastq rip.fastq (single-end Illumina sequencing)
# control fastq control.fastq (single-end Illumina sequencing)
#
# required directory structure and files
#
# _WORKINGDIRECTORY
#	|
#	|-- workflow.sh
#	|-- _raw
#	     |-- <control>.fastq
#	     |-- rip.fastq
#	|-- _misc
#	     |-- 1stQC.R
#	     |-- 2ndQC.R
#	     |-- correlationAnalysis.R
#	     |-- lsm-positionGenes.R
#	     |-- lsm-position.R
#	     |-- check-genes-with-lsm.sh
#	     |-- computeGC.sh
#	     |-- createCircosFiles.sh
#	     |-- circos.conf
#	     |-- <spp>-ISSaga-checked.gff3
#	     |-- adap.fa
#
# several directories and files will be created during the execution and they will be placed
# majorly in _WORKINGDIRECTORY and _misc
####################

######################
# checking dependencies
######################
echo >&2 "Checking dependencies"

command -v bedtools > /dev/null >&1 || { echo >&2 "bedtools isn't installed. Aborting" ; exit 1; }
command -v R > /dev/null >&1 || { echo >&2 "R isn't installed. Aborting" ; exit 1; }
command -v infoseq > /dev/null >&1 || { echo >&2 "emboss isn't installed. Aborting" ; exit 1; }
command -v bowtie > /dev/null >&1 || { echo >&2 "bowtie isn't installed. Aborting" ; exit 1; }
command -v samtools > /dev/null >&1 || { echo >&2 "samtools isn't installed. Aborting" ; exit 1; }
command -v mmr > /dev/null >&1 || { echo >&2 "mmr isn't installed. Aborting" ; exit 1; }
command -v circos > /dev/null >&1 || { echo >&2 "circos isn't installed. Aborting" ; exit 1; }
if [ ! -e /opt/Trimmomatic-0.36/trimmomatic-0.36.jar ] ; then echo >&2 "trimmomatic isn't installed. Aborting" ; exit 1; fi

# checking files
if [ ! -e misc/adap.fa ] ; then echo >&2 "Adapter sequences not found. Aborting" ; exit 1 ; fi
if [ ! -e raw/$2.fastq ] ; then echo >&2 "Control file not found. Aborting" ; exit 1 ; fi

if [ "$3" == "y" ] ; then
	if [ ! -e misc/$1-ISSaga-checked.gff3 ] ; then echo >&2 "IS annotation not found. Aborting" ; exit 1 ; fi
fi

# checking scripts
if [ ! -e misc/1stQC.R ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/2ndQC.R ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/correlationAnalysis.R ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/lsm-position.R ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/lsm-positionGenes.R ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/computeGC.sh ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/createCircosFiles.sh ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi
if [ ! -e misc/circos.conf ] ; then echo >&2 "Missing script. Aborting" ; exit 1 ; fi

######################
# getting parameters
######################
spp=${1}
control=${2}
positionAnalysis=${3} # choose if IS position analysis will be performed
threads=${4} # choose how many threads to be used
genomeURL=${5} # genome url from ncbi ftp
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

## parameters for interaction regions finding
# threshold=1.58496 is log2(3)
threshold=${12}
minsize=${13} # minimum size of interaction region

######################
# checking quality
######################
if [ ! -d 1stQC ] ; then
	echo >&2 "First round of quality control"
	mkdir 1stQC
	R -q -f misc/1stQC.R --args $threads > /dev/null 2>&1
fi

######################
# trimming with trimmomatic
######################
# adapter sequence example AGATCGGAAGAGCACACGTCTGAACTCCAGTCACAGTCAACTCTCGTATGCCGTCTTCTGCTTG
if [ ! -d trimmed ] ; then
	echo >&2 "Trimming fastq files"
	mkdir trimmed
	for i in raw/*.fastq ; do
		out=$(echo $i | sed 's/^raw/trimmed/g')
		log=$(echo $i | sed 's/^raw/trimmed/g;s/.fastq$/.log/')
		java -jar /opt/Trimmomatic-0.36/trimmomatic-0.36.jar SE -threads $threads $i $out ILLUMINACLIP:misc/adap.fa:2:30:10 SLIDINGWINDOW:5:28 MINLEN:15 > $log 2>&1
	done
fi

######################
# checking quality after trimming
######################

if [ ! -d 2ndQC ] ; then
	echo >&2 "Second round of quality control"
	mkdir 2ndQC
	R -q -f misc/2ndQC.R --args $threads > /dev/null 2>&1
fi

######################
# downloading genome and creating idx
######################

if [ ! -f misc/$spp.fa ] ; then
	echo >&2 "Downloading genome and creating Bowtie index"
	wget -O misc/$spp.fa.gz $genomeURL
	gunzip -f -k misc/$spp.fa.gz
fi	

# building index
hisat2-build misc/$spp".fa" misc/$spp > /dev/null 2>&1

######################
# aligning to reference genome
######################
if [ ! -d sam ] ; then
	echo >&2 "Aligning reads to reference genome"
	mkdir sam

	# creating file to store number of alignments
	# for each seq lib
	touch sam/readsAligned.txt

	## aligning to genome
	for i in trimmed/*.fastq ; do
		prefix=$(echo $i | sed 's/^trimmed/sam/;s/.fastq$//')
		out=$(echo $i | sed 's/^trimmed/sam/;s/fastq$/sam/')
		log=$(echo $i | sed 's/^trimmed/sam/;s/fastq$/log/')
		hisat2 --no-spliced-alignment \
                       --no-softclip \
		       --rdg 1000,1 \
		       --rfg 1000,1 \
                       --rna-strandness $hisatstrand \
                       -k 1000 \
                       -p $threads \
                       -x misc/$spp \
                       -U $i \
                       --summary-file $log > $out 
		uniqReads=`grep "aligned exactly 1 time" $log | sed 's/^    \(.*\) (.*$/\1/'` # extracting number of uniq aligned reads from log
		multiReads=`grep "aligned >1 times" $log | sed 's/^    \(.*\) (.*$/\1/'` # extracting number of multi aligned reads from log
		alnReads=`echo "$uniqReads+$multiReads" | bc` # computing how many reads have been aligned
		printf "$prefix-multi\t$alnReads\n" >> sam/readsAligned.txt # storing it to file
	done
fi

######################
# processing bam files
######################
if [ ! -d bam ] ; then
	echo >&2 "Processing BAM files"
	mkdir bam
	# sam2bam ; sort bam
	for i in sam/*.sam; do
		out=$(echo $i | sed 's/^sam/bam/;s/sam$/bam/')
		log=$(echo $i | sed 's/^sam/bam/;s/sam$/log/')
		samtools sort -@ $threads -n -o $out $i
	done

	# to do correction of multi-mapping reads
	# we are going to use MMR
	for i in bam/*.bam; do
		out=$(echo $i | sed 's/.bam$/-mmr.bam/')
		log=$(echo $i | sed 's/.bam$/-mmr.log/')
		mmr -o $out -t $threads $strandspecific -b -R $readsize $i
	done

	# we need to resort bam using default mode in order to proceed
	for i in bam/*-mmr.bam; do
		samtools sort -@ $threads -o $i $i
	done

	# I will also sort the input bam in order to compare mmr and non-mmr versions
	for i in `ls bam/*.bam | grep -v mmr`; do
		samtools sort -@ $threads -o $i $i
	done
fi

######################
# Creating coverage files
######################
# creating coverage file for rev and fwd
# ATTENTION HERE:
# for Truseq Stranded mRNA library preparation kit (Hsalinarum library)
# single-end reads are actually the reverse complement of original input RNA
# therefore the strands must be inverted to compute the genome coverage
if [ ! -d cov ] ; then
	echo >&2 "Creating read depth files"
	mkdir cov
	if [ "$invertstrand" == "y" ]
	then
		for i in bam/*-mmr.bam ; do
			prefix=$(echo $i | sed 's/^bam/cov/;s/-mmr.bam$//')
			bedtools genomecov -ibam $i -d -strand + > $prefix"-counts-rev.txt" # note the inversion of strand on the output file name
			bedtools genomecov -ibam $i -d -strand - > $prefix"-counts-fwd.txt" # note the inversion of strand on the output file name
		done
	else
		for i in bam/*-mmr.bam ; do
			prefix=$(echo $i | sed 's/^bam/cov/;s/-mmr.bam$//')
			bedtools genomecov -ibam $i -d -strand - > $prefix"-counts-rev.txt"
			bedtools genomecov -ibam $i -d -strand + > $prefix"-counts-fwd.txt"
		done
	fi

	# creating correction factor for normalization
	maxReads=`awk 'OFS="\t" {if($2 >= n){n = $2}} END {print n}' sam/readsAligned.txt` # storing number of maximum aligned reads of a lib
	awk -v maxReads=$maxReads 'OFS="\t" {print $1, $2, maxReads/$2}' sam/readsAligned.txt > awk.tmp # creating correction factor file
	mv awk.tmp sam/readsAligned.txt # renaming it

	## CONTROL LIB

	# creating fwd count file for GGB
	for i in cov/$control*counts-fwd.txt; do
		name=$(echo $i | sed 's/.txt//')
		idxLib=$(echo $i | sed 's/^cov\///;s/-counts-fwd.txt//')
		corFactor=`grep -w $idxLib sam/readsAligned.txt | awk '{print $3}'`
	
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
	for i in cov/$control*counts-rev.txt; do
		name=$(echo $i | sed 's/.txt//')
		idxLib=$(echo $i | sed 's/^cov\///;s/-counts-rev.txt//')
		corFactor=`grep -w $idxLib sam/readsAligned.txt | awk '{print $3}'`

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
	nonControl=`ls cov/*counts-fwd.txt | grep -v $control`
	for i in $nonControl; do 
		name=$(echo $i | sed 's/.txt//')
		idxLib=$(echo $i | sed 's/^cov\///;s/-counts-fwd.txt//')
		corFactor=`grep -w $idxLib sam/readsAligned.txt | awk '{print $3}'`

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
	nonControl=`ls cov/*counts-rev.txt | grep -v $control`
	for i in $nonControl; do
		name=$(echo $i | sed 's/.txt//')
		idxLib=$(echo $i | sed 's/^cov\///;s/-counts-rev.txt//')
		corFactor=`grep -w $idxLib sam/readsAligned.txt | awk '{print $3}'`
	
		
	
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

	######################
	# Computing fold changes
	######################
	
	if [ "$additionalPlots" == "y" ] ; then
		# fwd absolut count
		nonControl=`ls cov/*-fwd-absolut-ggb.txt | grep -v $control`
		for i in $nonControl; do
			name=$(echo $i | sed 's/.txt//')
			paste $i cov/$control-multi-counts-fwd-absolut-ggb.txt |\
			awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
	        		      {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
	        		       else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
	        		       else {print $1, $2, $3, "0"}}\
	        		       else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"
		done
	
		# rev absolut count
		nonControl=`ls cov/*-rev-absolut-ggb.txt | grep -v $control`
		for i in $nonControl; do
			name=$(echo $i | sed 's/.txt//')
			paste $i cov/$control-multi-counts-rev-absolut-ggb.txt |\
			awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
	        		      {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
	        		       else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
	        		       else {print $1, $2, $3, "0"}}\
	        		       else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"
		done
	fi

	# fwd normalized count
	nonControl=`ls cov/*-fwd-normalized-ggb.txt | grep -v $control`
	for i in $nonControl; do
		name=$(echo $i | sed 's/.txt//')
		paste $i cov/$control-counts-fwd-normalized-ggb.txt |\
		awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
	        	      {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
	        	       else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
        		       else {print $1, $2, $3, "0"}}\
        		       else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"
	done

	# rev normalized count
	nonControl=`ls cov/*-rev-normalized-ggb.txt | grep -v $control`
	for i in $nonControl; do
		name=$(echo $i | sed 's/.txt//')
		paste $i cov/$control-counts-rev-normalized-ggb.txt |\
		awk 'OFS="\t" {if($4 == 0 || $8 == 0)\
        		      {if($4 > $8){print $1, $2, $3, log($4)/log(2)}\
        		       else if($8 > $4){print $1, $2, $3, -1 * (log($8)/log(2))}\
        		       else {print $1, $2, $3, "0"}}\
        		       else {print $1, $2, $3, (log($4/$8))/log(2)}}' > $name"-log2FC.txt"
	done
fi

###########################################
# ENTERING NEW PROGRAM: createInTracks.sh #
###########################################
echo >&2 "Creating regions of interaction"

# this program creates regions of interaction (GFF3)
# based on the log2FC files

## functional variables
flag=0
n=0

# fwd processing
strandname=fwd
for i in cov/*fwd-normalized-ggb-log2FC.txt
do
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
for i in cov/*rev-normalized-ggb-log2FC.txt
do
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

######################
# performing position analysis
######################
# for insertion sequences and genes
if [ "$positionAnalysis" == "y" ] ; then
	echo >&2 "Performing position analysis"

	if [ ! -d positionAnalysis ] ; then
		mkdir positionAnalysis

		# getting intersections between IS file and interaction regions
		# this will reduce time for the next step
		for i in interaction-regions-entire-genome-*.gff3 ; do
			name=$(echo $i | sed s/-entire-genome//)
			bedtools intersect -wa -a $i -b misc/$spp-ISSaga-checked.gff3 | uniq > $name
		done

		# generating heatmaps for lsm position
		R --slave -q -f misc/lsm-position.R --args $spp > /dev/null 2>&1
	fi

	if [ ! -d positionAnalysisGenes ] ; then
		gfffile=$(echo $genomeURL | sed 's/.fna.gz$/.gff.gz/')
		mkdir positionAnalysisGenes
		bash misc/check-genes-with-lsm.sh $spp $gfffile
	fi
fi

######################
# computing genome-wide GC content
######################
if [ ! -d gccontent ] ; then
	echo >&2 "Computing GC content"
	mkdir gccontent
	bash misc/computeGC.sh misc/$spp.fa $windowsize $stepsize
fi

######################
# correlation analysis
######################
if [ ! -d correlationAnalysis ] ; then
	echo >&2 "Performing correlation analysis"
	mkdir correlationAnalysis
	R --slave -q -f misc/correlationAnalysis.R --args $windowsize $stepsize $threshold > /dev/null 2>&1
fi

######################
# creating circos ideograms
######################
if [ ! -d circos ] ; then
	echo >&2 "Generating Circos Ideograms"
	mkdir circos
	bash misc/createCircosFiles.sh $spp $positionAnalysis
	
	infoseq -only -name -length misc/$spp.fa 2> /dev/null | tail -n +2 |\
	while read replicon length ; do
		name=$(echo $replicon | sed 's/\.[0-9]//')
		unit=${#length}
		if [ $unit -ge 7 ]
		then
			circos -conf misc/circos.conf -param karyotype=circos/$name"-karyotype.txt" -outputfile circos/$name -silent
		else
			circos -conf misc/circos.conf \
			       -param karyotype=circos/$name"-karyotype.txt" \
			       -param ticks/multiplier="1e-3" \
			       -param ticks/format="%.2f Kb" \
			       -param ticks/tick/spacing="25000u" \
			       -outputfile circos/$name \
			       -silent
		fi
	done
fi

echo >&2 "Done."

