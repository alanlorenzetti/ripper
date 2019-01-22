#!/bin/bash

# alorenzetti 201708

# require emboss
# require bedtools

# usage ./computeGC.sh <spp> <miscdir> <window_size> <step_size> <gccontent_dir>
# e.g ./computeGC.sh Hsalinarum.fa 50 50

spp=$1
miscdir=$2
windowsize=$3
stepsize=$4
gccontentdir=$5

# genome fasta file = $1
infoseq -only -length -name $miscdir/$spp.fa | tail -n +2 | sed "s/ \+/	/" > gtmp 2> /dev/null

# creating windows
bedtools makewindows -g gtmp -w $windowsize -s $stepsize > wtmp

# computing GC content
bedtools nuc -fi $miscdir/$spp.fa -bed wtmp |\
tail -n +2 |\
awk -v OFS="\t" -v FS="\t" '{print $1, $2+1, $3, $5}' > $gccontentdir/gccontent.txt

rm gtmp wtmp

# computing average GC content for LSm-bound transcripts and replicons
# for LSm-bound transcripts
echo replicon_name Avg_GC_content_genome AVG_GC_content_transcript > tmp2

for i in `grep "^>" $1 | sed 's/^>//;s/ .*$//'` ; do

	cat interaction-regions-entire-genome-fwd.gff3 interaction-regions-entire-genome-rev.gff3 |\
	awk -v replicon=$i -v FS="\t" -v OFS="\t" '{if($1 == replicon){print}}' > tmp1
	GCtranscript=`bedtools nuc -fi $1 -bed tmp1 | awk -v sum=0 '{sum=sum+$11}END{print sum/NR*100}'`
	GCgenome=`infoseq $1 2> /dev/null | grep $i | awk '{print $7}'`

	echo $i $GCgenome $GCtranscript
done >> tmp2

rm tmp1

tail -n +2 tmp2 |\
while read replicon avggenomegc avgtxgc ; do
	awk -v OFS="\t" -v replicon=$replicon -v avggenomegc=$avggenomegc '{if($1 == replicon){print $1,$2,$3,avggenomegc/100}}' $gccontentdir/gccontent.txt
done > $gccontentdir/avggccontent.txt

rm tmp2
