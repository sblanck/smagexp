#!/usr/bin/env Rscript
# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("optparse")

##### Read options
option_list=list(
		make_option("--id",type="character",default=NULL,help="GSE ID from GEO databse (required)"),
		make_option("--report",type="character",default=NULL,help="Text file summarizing conditions of the experiment")

);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$id)){
	print_help(opt_parser)
	stop("Recount id required.", call.=FALSE)
}

#loading libraries
suppressPackageStartupMessages(require(GEOquery))

studyID=opt$id
reportFile=opt$report

dir.create("./split", showWarnings = TRUE, recursive = FALSE)

url <- download_study(studyID)
load(file.path(studyID, 'rse_gene.Rdata'))
rse <- scale_counts(rse_gene)
counts=assay(rse)
conditions=rse$title

for (i in 1:ncol(counts))
{
	currentCount=as.data.frame(counts[,i])
	sampleID=colnames(counts)[i]
	colnames(currentCount)=sampleID
	write.table(x=currentCount,file=paste0("./split/",sampleID,"_",conditions[i]),sep="\t",row.names = T, col.names = T)
}

write.table(as.data.frame(cbind(sampleID=colnames(counts),title=conditions)),quote = FALSE,col.names =TRUE, row.names=FALSE,file=reportFile,sep="\t")