#!/bin/bash

# alorenzetti 201705
# this script will take lsm-interaction-positions
# in gff3 format and then process it to find genes
# interacting with this protein in H. salinarum NRC-1

# usage
# ./check-genes-with-lsm.sh

# program requirements
# bedtools v2.25.0 or superior
# R v3.4.0 or superior
# emboss v6.6.0.0

# file requirements
# interaction-regions-entire-genome-fwd.gff3 - interaction regions for LSm protein fwd
# interaction-regions-entire-genome-rev.gff3 - interaction regions for LSm protein fwd
# Hsalinarum.fa

spp=$1
gfffile=$2

##########
# 1st part - downloading and parsing annotation file from ncbi
##########

if [ ! -f "misc/$spp.gff3" ] ; then
	wget -O misc/$spp.gff3.gz $gfffile > /dev/null 2>&1
	gunzip -f -k misc/$spp.gff3.gz
	
	sed -i '/^#/d' misc/$spp.gff3.gz
	awk -v OFS="\t" -v FS="\t" '{if($3 == "gene"){print}}' misc/$spp.gff3 > tmp1
	mv tmp1 misc/$spp.gff3
fi

##########
# 2nd part - taking regions interacting with genes
##########

bedtools intersect -wa -a interaction-regions-entire-genome-fwd.gff3 -b misc/$spp.gff3 |\
uniq > positionAnalysisGenes/interaction-regions-genes.gff3

bedtools intersect -wa -a interaction-regions-entire-genome-rev.gff3 -b misc/$spp.gff3 |\
uniq >> positionAnalysisGenes/interaction-regions-genes.gff3

bedtools intersect -wa -a misc/$spp.gff3 -b interaction-regions-entire-genome-fwd.gff3 |\
uniq > tmp1

bedtools intersect -wa -a misc/$spp.gff3 -b interaction-regions-entire-genome-rev.gff3 |\
uniq > tmp2

cat tmp1 tmp2 | sort -Vk1,1 -k4,4 | uniq > positionAnalysisGenes/$spp-annot-genesWithInteraction.gff3

rm tmp1 tmp2

##########
# 3rd part - downstream analysis using R
##########

R --slave -q -f misc/lsm-positionGenes.R --args positionAnalysisGenes/interaction-regions-genes.gff3 positionAnalysisGenes/$spp-annot-genesWithInteraction.gff3 > /dev/null 2>&1


