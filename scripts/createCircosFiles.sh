#!/bin/bash

# this script will generate circos files
# it takes as arg the spp

spp=$1
positionAnalysis=$2

# karyotypes
infoseq misc/$spp.fa 2> /dev/null | tail -n +2 | awk '{print $3,$6}' |\
while read replicon length ; do
	prefix=$(echo $replicon | sed 's/\.[0-9]//')
	echo "chr - $replicon $replicon 1 $length grey" > circos/$prefix"-karyotype.txt"
done

# IShighlights
if [ "$positionAnalysis" == "y" ]
then
	awk '{print $1,$4,$5,fill_color=black}' misc/$spp-ISSaga-checked.gff3 > circos/IShighlights.txt
else
	touch circos/IShighlights.txt
fi
