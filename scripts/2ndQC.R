# retrieving variables from command line
args = commandArgs(trailingOnly=TRUE)
threads = args[1]
secondqcdir = args[2]
trimmeddir = args[3]

# loading Rqc
library("Rqc")

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
files = list.files(path=trimmeddir, pattern='.fastq', full.names=T)

# generating report
qual = rqcQA(files, workers=threads)
rqcReport(qual, outdir = secondqcdir, file = "trimmed-quality-Rqc")
