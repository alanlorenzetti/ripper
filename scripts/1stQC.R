# retrieving variables from command line
args = commandArgs(trailingOnly=TRUE)
threads = args[1]

# loading Rqc
if(!require('Rqc')){
	source('https://bioconductor.org/biocLite.R')
	biocLite('Rqc')
	library('Rqc')	
}

# workaround
.readFrequency <- function (chunk)
{
	tbl <- table(as.character(sread(chunk)))
	count <- as.integer(tbl)
	hash <- names(tbl)
	data.frame(hash, count, stringsAsFactors = FALSE)
}
assignInNamespace(".readFrequency", .readFrequency, "Rqc")

# listing files
files = list.files(path='raw', pattern='.fastq', full.names=T)

# generating report
qual = rqcQA(files, workers=threads)
rqcReport(qual, outdir = "1stQC", file = "raw-quality-Rqc")
