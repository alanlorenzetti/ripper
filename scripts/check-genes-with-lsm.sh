#!/bin/bash

# alorenzetti jan 2019
# this script will take lsm-interaction-positions
# in gff3 format and then process it to find genes
# interacting with this protein in H. salinarum NRC-1

# usage
# ./check-genes-with-lsm.sh

# program requirements
# bedtools v2.25.0 or superior
# R v3.4.0 or superior

# file requirements
# interaction-regions-entire-genome-fwd.gff3 - interaction regions for LSm protein fwd
# interaction-regions-entire-genome-rev.gff3 - interaction regions for LSm protein fwd
# Hsalinarum.fa

spp=$1
annotURL=$2
miscdir=$3
scriptsdir=$4
positionanalysisgenesdir=$5

##########
# 1st part - downloading and parsing annotation file from ncbi
##########

if [ ! -f $miscdir/$spp".gff" ] ; then
    curl -u anonymous: $annotURL 2> /dev/null |\
    zcat | awk -v FS="\t" -v OFS="\t" '{if($3 != "region"){print}}' \
    > $miscdir/$spp"-full-annotation.gff" || { echo >&2 "Annotation file download failed. Aborting" ; exit 1; }

    cp $miscdir/$spp"-full-annotation.gff" $miscdir/$spp".gff"
	sed -i '/^#/d' $miscdir/$spp".gff"
	awk -v OFS="\t" -v FS="\t" '{if($3 == "gene"){print}}' $miscdir/$spp".gff" > tmp1
	mv tmp1 $miscdir/$spp".gff"
fi

##########
# 2nd part - taking regions interacting with genes
##########

bedtools intersect -wa -a interaction-regions-entire-genome-fwd.gff3 -b $miscdir/$spp.gff |\
uniq > $positionanalysisgenesdir/interaction-regions-genes.gff3

bedtools intersect -wa -a interaction-regions-entire-genome-rev.gff3 -b $miscdir/$spp.gff |\
uniq >> $positionanalysisgenesdir/interaction-regions-genes.gff3

bedtools intersect -wa -a $miscdir/$spp.gff -b interaction-regions-entire-genome-fwd.gff3 |\
uniq > tmp1

bedtools intersect -wa -a $miscdir/$spp.gff -b interaction-regions-entire-genome-rev.gff3 |\
uniq > tmp2

cat tmp1 tmp2 | sort -Vk1,1 -k4,4 | uniq > $positionanalysisgenesdir/$spp-annot-genesWithInteraction.gff3

rm tmp1 tmp2

##########
# 3rd part - downstream analysis using R
##########

R --slave -q -f $scriptsdir/lsm-positionGenes.R \
  --args $positionanalysisgenesdir/interaction-regions-genes.gff3 \
  $positionanalysisgenesdir/$spp-annot-genesWithInteraction.gff3 \
  $positionanalysisgenesdir > /dev/null 2>&1
