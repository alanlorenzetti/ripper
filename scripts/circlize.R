# alorenzetti jan 2019

library(circlize)
library(rtracklayer)

args = commandArgs(trailingOnly = T)

spp = args[1]
miscdir = args[2]
gccontentdir = args[3]
correlationanalysisdir = args[4]
circlizedir = args[5]

# svg parameters?
ht=10
wh=10

# min, mean and max gc should be computed for entire genome?
# if "n" the program will compute min, mean and max separately for chr
# and the combination of plasmids
genomewidegc="y"

# default pal to match ggplot2 colors
# defaultPal=hue_pal()(10)

# reading genome info
genomeInfo=read.delim(paste0(miscdir, "/", spp, ".fa.fai"), header=F)
genomeInfo=genomeInfo[,c(1,2)]

# creating a df to represent the layout of genome
# two last columns are required but not meaningful
df=data.frame(acc=genomeInfo[,1],
              start=0,
              acc=genomeInfo[,2],
              foo=genomeInfo[,1],
              bar=genomeInfo[,1],
              stringsAsFactors = F)
df$acc = as.character(df$acc)

# reading and preparing gccontent
gc = read.delim(file = paste0(gccontentdir, "/", "gccontent.txt"), header=F, stringsAsFactors = F)
gc$V2 = gc$V2 - 1

# reading and preparing lsm density
lsm = read.delim(file = paste0(correlationanalysisdir, "/", "lsmDensity.txt"), header=F, stringsAsFactors = F)
lsm$V2 = lsm$V2 - 1

for(repliconidx in 1:length(df$acc)){
  # setting parameters
  circos.clear()
  circos.par("start.degree" = 90,
             "track.height" = 0.075)
  
  # saving
  svg(paste0(circlizedir, "/", df[repliconidx,1], ".svg"), height = ht, width = wh)
  par(cex=1.5)
  
  #### init plot with the genome layout ####
  circos.initializeWithIdeogram(df[repliconidx,], plotType = c("labels","axis")[2], chromosome.index = df[repliconidx,1])
  
  ####regions IS plasmids#####
  is = as.data.frame(rtracklayer::import(paste0(miscdir, "/", spp, "-ISSaga-checked.gff3")))
  is$rpt_family = as.character(sub("\\+.*$", "", is$rpt_family))
  isfamilies=unique(is$rpt_family)
  bed_list = vector("list", length = length(isfamilies))
  n=1
  for(i in isfamilies){
    bed_list[[n]] = is[is$rpt_family == i,c("seqnames", "start", "end")]
    n=n+1
  }
  
  # plotting annotated IS regions
  circos.genomicTrackPlotRegion(bed_list, ylim = c(0,1), track.height=0.075, bg.col = adjustcolor("white", alpha.f = 0.9),
                                panel.fun = function(region, value, ...) {
                                  i=getI(...)
                                  circos.genomicRect(region, value, col = "gray")
                                }
  )
  
  #####plotting gc content heat map#####
  if(genomewidegc == "y"){
    mingc=min(gc$V4)
    maxgc=max(gc$V4)
    meangc=mean(gc$V4)
  }else{
    mingc=min(gc[gc$V1 == df[repliconidx,1], 4])
    maxgc=max(gc[gc$V1 == df[repliconidx,1], 4])
    meangc=mean(gc[gc$V1 == df[repliconidx,1], 4])
  }
  
  col=colorRamp2(c(mingc,meangc,maxgc), c("green", "black", "red"))
  circos.genomicHeatmap(gc, col = col, numeric.column = 4, connection_height = 0.0001, heatmap_height = 0.1)
  
  #####lsm density heat map#####
  col=colorRamp2(c(0,1), c("white", "red"))
  circos.genomicHeatmap(lsm, col = col, numeric.column = 4, connection_height = 0.0001, heatmap_height = 0.1)
  # need a border down here
  
  # closing dev to save plot
  dev.off()
}
