# alorenzetti 201705

# this script was adapted from a previous one
# which analyzed lsm interaction positions in IS
# in this way some variables might have nonsense names
# and further unnecessary processing

# this one is intended to find interaction regions in genes

# usage R -f lsm-position.R --args <lsm-interaction-file> <geneAnnotation-file>
# e.g. R -f lsm-position.R --args <interaction-regions-genes.gff3> <Hsalinarum-annot-genesWithInteraction.gff3

# it finds in which position there is a putative lsm interaction
# five sense and
# five antisense regions are possible

# sense
# |  1  |  2  |  3  |  4  |  5  |
# -------------------------------
# ·            Gene             ·
# -------------------------------
# |  1  |  2  |  3  |  4  |  5  |
# antisense

# takes as input the gff gene annotation
# and also the gff containing lsm positions

# getting args
args = commandArgs(trailingOnly=TRUE)
positionanalysisgenesdir=args[3]

# libs
# gplots
if(!require("gplots")){install.packages("gplots", repos="https://vps.fmvz.usp.br/CRAN/")}; library("gplots")
# MSG
if(!require("MSG")){install.packages("MSG", repos="https://vps.fmvz.usp.br/CRAN/")}; library("MSG")
# VennDiagram
if(!require("VennDiagram")){install.packages("VennDiagram", repos="https://vps.fmvz.usp.br/CRAN/")}; library("VennDiagram")

## misc settings
# gffheader = c("chr", "source", "key", "start", "end", "score", "strand", "phase", "att")
# dividing the gene in five parts; each block is 20 per cent of seq
pct = 0.2

## loading lsm gff3s
lsm = read.delim(args[1], header = FALSE, row.names = NULL)
#lsm = read.delim("interaction-regions-genes.gff3", header = FALSE, row.names = NULL)

## loading NCBI gff3
IS = read.table(args[2], header=FALSE, row.names = NULL, sep="\t")
#IS = read.table("Hsalinarum-annot.gff3", header=FALSE, row.names = NULL, sep="\t")

# subsetting genes
IS = IS[IS[,3] == "gene",]; rownames(IS) = NULL
IS[,4] = as.numeric(IS[,4])
IS[,5] = as.numeric(IS[,5])

# adjusting att col
IS$V9 = gsub("\\+", " ", IS$V9)

# creating table to store results
results = NULL

# creating table to get the name of a gene
# on which each lsm interaction relies on
genes = NULL

# if results.csv exists
# it will just skip the following analysis
if(!file.exists(paste0(positionanalysisgenesdir,"/results.csv")){
  
  # iterating over IS to find lsm interactions lying on it
  for(i in 1:dim(IS)[1]){
    # parsing gene table
    ISline = IS[i,]
    ISchr = as.character(ISline$V1)
    ISstart = as.numeric(ISline$V4)
    ISend = as.numeric(ISline$V5)
    ISstrand = as.character(ISline$V7)
    ISname = sub(".*;Name=(.*?)\\;.*$", "\\1", as.character(ISline$V9), perl = TRUE)
  
    # setting pct
    ISpct = (ISend - ISstart + 1) * pct ; ISpct = round(ISpct)
    
    # setting regions for an IS
    one = seq(ISstart, ISstart + ISpct)
    two = seq(ISstart + ISpct + 1, ISstart + ISpct*2)
    three = seq(ISstart + 2*ISpct + 1, ISstart + ISpct*3)
    four = seq(ISstart + 3*ISpct + 1, ISstart + ISpct*4)
    five = seq(ISstart + 4*ISpct + 1, ISend)
    
    # setting counters for actual IS
    sone = 0 ; stwo = 0 ; sthree = 0 ; sfour = 0 ; sfive = 0 # sense
    asone = 0 ; astwo = 0; asthree = 0 ; asfour = 0 ; asfive = 0 # antisense
    
    # reading lsm positions
    for(j in seq(1, dim(lsm)[1])){
      lsmline = lsm[j,]
      lsmstart = as.numeric(lsmline$V4)
      lsmend = as.numeric(lsmline$V5)
      lsmrange = as.numeric(seq(lsmstart, lsmend))
      if(lsmline[1] == ISchr){ # only compare if chr matches
        if(lsmline[7] == ISstrand){ # makes sense here
          # position conditionals here
          if(sum(lsmrange %in% one) > 0){
            sone = sone + 1
          }
          if(sum(lsmrange %in% two) > 0){
            stwo = stwo + 1
          }
          if(sum(lsmrange %in% three) > 0){
            sthree = sthree + 1
          }
          if(sum(lsmrange %in% four) > 0){
            sfour = sfour + 1
          }
          if(sum(lsmrange %in% five) > 0){
            sfive = sfive + 1
          }
        }else{ # antisense here
          # position conditionals here
          if(sum(lsmrange %in% one) > 0){
            asone = asone + 1
          }
          if(sum(lsmrange %in% two) > 0){
            astwo = astwo + 1
          }
          if(sum(lsmrange %in% three) > 0){
            asthree = asthree + 1
          }
          if(sum(lsmrange %in% four) > 0){
            asfour = asfour + 1
          }
          if(sum(lsmrange %in% five) > 0){
            asfive = asfive + 1
          }
        }
      }else{next}
    }
    
    # for gene on minus strand I must invert antisense counts
    # gamb
    if(ISline$V7 == "+"){
      results = rbind(results, c(ISchr, ISstart, ISend, ISstrand, ISname,
                                 sone, stwo, sthree, sfour, sfive,
                                 asone, astwo, asthree, asfour, asfive))
    }else{
      results = rbind(results, c(ISchr, ISstart, ISend, ISstrand, ISname,
                                 sfive, sfour, sthree, stwo, sone,
                                 asfive, asfour, asthree, astwo, asone))
    }
  }
  
# closing the conditional checking results.csv existence
# it will load results.csv in case it exists
}else{results = read.csv(paste0(positionanalysisgenesdir, "/results.csv", header=TRUE)}

# transforming results table type; adding header to results
results = as.data.frame(results)
results[,c(2,3,6:15)] = apply(results[,c(2,3,6:15)], 2, function(x) as.numeric(x))
colnames(results) = c("chr", "start", "end", "strand", "name", 
                      "one_sense", "two_sense", "three_sense", "four_sense", "five_sense",
                      "one_as", "two_as", "three_as", "four_as", "five_as")

# writing results to table
write.csv(results, file = paste0(positionanalysisgenesdir, "/results.csv", quote = FALSE, row.names = FALSE)

# creating a matrix to further analysis
M = results[,6:15]
rownames(M) = paste(results$name, rownames(M), sep="_")
rownames(M) = gsub(" ", "_", rownames(M))

# removing genes that doesn't have any lsm
M = M[apply(M, 1, sum) != 0,]

# number of lsm regions per gene doesn't make sense
# to compute binary distance
# switching to presence/absence mode
MM = (results[,6:15] >= 1) * 1
rownames(MM) = paste(results$name, 1:dim(MM)[1], sep="_")
rownames(MM) = gsub(" ", "_", rownames(MM))

# removing genes that doesn't have any lsm interaction
MM = MM[apply(MM, 1, sum) != 0,]

# making clusters
dist = dist(MM, method = "binary", upper = TRUE, diag = TRUE)
clust = hclust(dist, method = "ward.D")

# plotting heatmap
MM = as.matrix(MM)

svgpaste0(positionanalysisgenesdir, "/heatmap.svg", width = 10, height = 60)
heatmap.2(MM,
          scale="none",
          trace = "none",
          key = FALSE,
          Rowv = as.dendrogram(clust),
          Colv = "none",
          srtCol=45,
          margins=c(10,25),
          lhei = c(0.2, 4),
          dendrogram = "row",
          cexRow = 0.8,
          cexCol = 1,
          col = c("white", "red"),
          sepwidth=c(0.001, 0.001),
          sepcolor="black",
          colsep=1:ncol(MM),
          rowsep=1:nrow(MM))
dev.off()

##############
# generating report for lsm-interacting genes
##############

# genes with lsm interaction on sense strand
filter = apply(results[,6:10], 1, sum) != 0
genesWithLsmSense = as.character(results[filter,5])
write.table(genesWithLsmSense, paste0(positionanalysisgenesdir, "/genesWithLsmSense.txt", sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# genes with lsm interaction on antisense strand
filter = apply(results[,11:15], 1, sum) != 0
genesWithLsmAntiSense = as.character(results[filter,5])
write.table(genesWithLsmAntiSense, paste0(positionanalysisgenesdir, "/genesWithLsmAntiSense", sep = " ",
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# drawing venn diagram
vennInput=list(Sense=genesWithLsmSense, Antisense=genesWithLsmAntiSense)
vennDiag = venn.diagram(vennInput, filename = paste0(positionanalysisgenesdir, "/vennDiagram.png", imagetype = "png",
                        fill=c(2,5), alpha=0.5, cex = 1.5)

# genes with lsm interaction only on sense strand
#filter = apply(results[,6:10], 1, sum) != 0 & apply(results[,11:15], 1, sum) == 0
#genesOnlyWithLsmSense = as.character(results[filter,5])

# genes with lsm interaction only on antisense strand
#filter = apply(results[,6:10], 1, sum) == 0 & apply(results[,11:15], 1, sum) != 0
#genesOnlyWithLsmAntiSense = as.character(results[filter,5])

# genes with lsm interaction on sense and antisense
#filter = apply(results[,6:10], 1, sum) != 0 & apply(results[,11:15], 1, sum) != 0
#genesOnlyWithLsmSenseAndAntiSense = as.character(results[filter,5])
