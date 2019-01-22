# alorenzetti 201708
# this script will compute LSm density and compute the correlation with AT content

# getting args
args = commandArgs(trailingOnly=TRUE)
windowsize = args[1]
stepsize = args[2]
lfc = args[3]

# loading libs
# ggplot2
if(!require("ggplot2")){install.packages("ggplot2", repos="https://vps.fmvz.usp.br/CRAN/"); library("ggplot2")}

# getting AT content
gc = read.delim("gccontent/gccontent.txt", header=FALSE)
at = gc
at[,4] = 1 - at[,4]
at[,1] = factor(at[,1], levels=unique(at[,1]))

meandf=aggregate(at[,4]~at[,1], FUN=mean)
colnames(meandf) = c("V1", "meanAT")

# getting LSm density
lsmfwd = read.delim("cov/rip-counts-fwd-normalized-ggb-log2FC.txt", header=FALSE)
lsmrev = read.delim("cov/rip-counts-rev-normalized-ggb-log2FC.txt", header=FALSE)
lsm = lsmfwd
lsm[,4] = ((lsmfwd[,4] >= lfc) * 1) + ((lsmrev[,4] >= lfc) * 1)
lsm[,4] = (lsm[,4] >= 1) * 1
lsm[,1] = as.character(lsm[,1])
lsmMean = at
indexes = table(lsmMean[,1])

for(i in 1:length(indexes)){
	if(i == 1){
		first=1
		last=indexes[1]

		for(j in first:last){
			lsmMean[j,4] = mean(lsm[lsmMean[j,2]:lsmMean[j,3],4])
		}

	} else {
		first=sum(indexes[1:i-1])+1
		last=sum(indexes[1:i])

		for(j in first:last){
			lsmsub = subset(lsm, V1 == names(indexes[i]))
			lsmMean[j,4] = mean(lsmsub[lsmMean[j,2]:lsmMean[j,3],4])
		}
	}
}

# writing lsm density file
write.table(lsmMean, "correlationAnalysis/lsmDensity.txt", row.names=F, col.names=F, quote=F, sep="\t")

# joining AT content with lsm density
at[,5] = lsmMean[,4]
at = subset(at, V5 != 0)

# scatter plot function
scatterPlot = function(df){
  ggplot(df, aes(x=V4, y=V5)) +
    scale_x_continuous(limits=c(0.2, 0.8), breaks=seq(0.2, 0.8, length.out = 5)) +
    scale_y_continuous(limits=c(0, 1), breaks=seq(0, 1, length.out = 5)) +
    xlab(paste("AT content per", windowsize, "bp window")) +
    ylab(paste("RBP density per", windowsize, "bp window")) +
#   geom_point(alpha = 0.5, shape=1) +
    geom_density2d(col = "black") +
    geom_smooth(method = "lm", formula = y ~ x, se = T, col = "black") +
    facet_grid(. ~ V1) +
    geom_vline(aes(xintercept = meanAT), col="red", data=meandf) +
    geom_text(aes(x=meanAT-.01, y=.75),
              data=meandf,
              label="Mean AT content",
              angle=90,
              size=3.5) +
    theme_bw()
}

png("correlationAnalysis/scatterPlot.png", width=600, height=300)
scatterPlot(at)
dev.off()

svg("correlationAnalysis/scatterPlot.svg", width=10, height=4)
scatterPlot(at)
dev.off()

for(i in 1:length(indexes)){
	atsub = subset(at, V1 == names(indexes[i]))
	corTest = cor.test(atsub[,4], atsub[,5], method="spearman", alternative="greater")
	filename = paste0("correlationAnalysis/", sub(".[0-9]$", "", names(indexes[i])), "-correlationTest.txt")
	write(capture.output(corTest), filename)
}


