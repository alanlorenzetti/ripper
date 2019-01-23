# alorenzetti jan 2019
# v0.3

# it finds where is a lsm interaction within one IS
# five sense/antisense regions are possible:

# sense
# |  1  |  2  |  3  |  4  |  5  |
# -------------------------------
# ·             IS              ·
# -------------------------------
# |  1  |  2  |  3  |  4  |  5  |
# antisense

# takes as input the gff IS annotation with family
# and also the gff containing lsm interaction regions

# getting args
args = commandArgs(trailingOnly=TRUE)

spp=args[1]
miscdir=args[2]
positionanalysisdir=args[3]

# libs
if(!require("gplots")){install.packages("gplots", repos="https://vps.fmvz.usp.br/CRAN/"); library("gplots")}
if(!require("MSG")){install.packages("MSG", repos="https://vps.fmvz.usp.br/CRAN/"); library("MSG")}

## misc settings
# gffheader = c("chr", "source", "key", "start", "end", "score", "strand", "phase", "att")
# dividing the IS in five parts; each block is 20 per cent of seq
pct = 0.2

## loading lsm gff3s
if(file.size("interaction-regions-fwd.gff3") != 0 & file.size("interaction-regions-rev.gff3") != 0){
	lsmfwd = read.delim("interaction-regions-fwd.gff3", header=FALSE, row.names = NULL)
	lsmrev = read.delim("interaction-regions-rev.gff3", header=FALSE, row.names = NULL)
	# just one object for all interaction regions
	lsm = rbind(lsmfwd, lsmrev)
}

if(file.size("interaction-regions-fwd.gff3") != 0 & file.size("interaction-regions-rev.gff3") == 0){
	lsm = read.delim("interaction-regions-fwd.gff3", header=FALSE, row.names = NULL)
}

if(file.size("interaction-regions-fwd.gff3") == 0 & file.size("interaction-regions-rev.gff3") != 0){
	lsm = read.delim("interaction-regions-fwd.gff3", header=FALSE, row.names = NULL)
}

if(file.size("interaction-regions-fwd.gff3") == 0 & file.size("interaction-regions-rev.gff3") == 0){
	quit(save = "no", status = 1, runLast = FALSE)
}

## loading ISSaga gff3
IS = read.delim(paste0(miscdir, "/", spp,"-ISSaga-checked.gff3"), header=FALSE, row.names = NULL)

# subsetting mobile_elements
IS = IS[IS[,3] == "mobile_element",]; rownames(IS) = NULL
IS[,4] = as.numeric(IS[,4])
IS[,5] = as.numeric(IS[,5])

# adjusting att col
IS$V9 = gsub("\\+", " ", IS$V9)

# creating table to store results
results = NULL

# iterating over IS to find interactions lying on it
for(i in 1:dim(IS)[1]){
  # parsing IS table
  ISline = IS[i,]
  ISchr = as.character(ISline$V1)
  ISstart = as.numeric(ISline$V4)
  ISend = as.numeric(ISline$V5)
  ISstrand = as.character(ISline$V7)
  ISname = sub(".*;Name=(.*?)\\;.*$", "\\1", as.character(ISline$V9), perl = TRUE)
  ISfamily = sub(".*rpt_family=", "", as.character(ISline$V9), perl = TRUE)
  
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
  
  # reading interaction regions
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
  
  # for IS on minus strand I must invert antisense counts
  # gamb
  if(ISline$V7 == "+"){
    results = rbind(results, c(ISchr, ISstart, ISend, ISstrand, ISname, ISfamily,
                               sone, stwo, sthree, sfour, sfive,
                               asone, astwo, asthree, asfour, asfive))
  }else{
    results = rbind(results, c(ISchr, ISstart, ISend, ISstrand, ISname, ISfamily,
                               sfive, sfour, sthree, stwo, sone,
                               asfive, asfour, asthree, astwo, asone))
  }
}

# transforming results table type; adding header to results
results = as.data.frame(results)
results[,c(2,3,7:16)] = apply(results[,c(2,3,7:16)], 2, function(x) as.numeric(x))
colnames(results) = c("chr", "start", "end", "strand", "name", "family", 
                      "one_sense", "two_sense", "three_sense", "four_sense", "five_sense",
                      "one_as", "two_as", "three_as", "four_as", "five_as")

## getting only ISs with at least one interaction region
#ISwithLsm = apply(results[,7:18], 1, sum) > 0
#results = results[ISwithLsm,]

# number of interaction regions per IS doesn't matter
# switching to presence/absence mode
results[,7:16] = (results[,7:16] >= 1) * 1

# writing results to table
write.csv(results, file = paste0(positionanalysisdir, "/results.txt"), quote = FALSE, row.names = FALSE)

## This block for redundant names

# making clusters
M = results[,7:16] # all redundant IS
rownames(M) = paste(results$name, results$family, rownames(M), sep = "_")
rownames(M) = gsub(" ", "_", rownames(M))
dist = dist(M, method = "binary", upper = TRUE, diag = TRUE)
clust = hclust(dist, method = "ward.D")

# creating vector of colors for each family
x = as.factor(sub(" .*$", "", results$family))
col = vec2col(x, name = "Spectral")

# plotting heatmap
M = as.matrix(M)

svg(paste0(positionanalysisdir, "/redundantIS.svg"), width = 8.267, height = 11.692)
heatmap.2(M,
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
          colsep=1:ncol(M),
          rowsep=1:nrow(M),
          RowSideColors = col)
dev.off()

## This block for non-redundant names

# setting non redundant results by name
NRresults = results[!duplicated(results[,5]),]
MM = NRresults[,7:16]
rownames(MM) = paste(NRresults$name, NRresults$family, rownames(MM), sep = "_")
rownames(MM) = gsub(" ", "_", rownames(MM))
dist = dist(MM, method = "binary", upper = TRUE, diag = TRUE)
clust = hclust(dist, method = "ward.D")

# creating vector of colors for each family
x = as.factor(sub(" .*$", "", NRresults$family))
col = vec2col(x, name = "Spectral")

# plotting heatmap
MM = as.matrix(MM)
svg(paste0(positionanalysisdir, "/nonredundantIS.svg"), width = 8.267, height = 11.692)
heatmap.2(MM,
          scale="none",
          trace = "none",
          key = FALSE,
          Rowv = as.dendrogram(clust),
          Colv = "none",
          srtCol=45,
          margins=c(10,25),
          lhei = c(0.2, 4),
          lwid = c(0.4, 2, 8),
          dendrogram = "row",
          cexRow = 1,
          cexCol = 1,
          col = c("white", "red"),
          sepwidth=c(0.001, 0.001),
          sepcolor="black",
          colsep=1:ncol(MM),
          rowsep=1:nrow(MM),
          RowSideColors = col)
dev.off()

