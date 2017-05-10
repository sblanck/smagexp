#!/usr/bin/env Rscript
# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("optparse")

##### Read options
option_list=list(
		make_option("--id",type="character",default=NULL,help="GSE ID from GEO databse (required)"),
		make_option("--transformation",type="character",default=NULL,help="log2 transformation (required)"),
		make_option("--data",type="character",default=NULL,help="A table containing the expression data"),
		make_option("--rdata",type="character",default="NULL",help="rdata object containing eset object"),
		make_option("--conditions",type="character",default=NULL,help="Text file summarizing conditions of the experiment")
		
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$id)){
	print_help(opt_parser)
	stop("GEOdata id required.", call.=FALSE)
}

#loading libraries
suppressPackageStartupMessages(require(GEOquery))

GEOQueryID=opt$id
GEOQueryData=opt$data
GEOQueryRData=opt$rdata
conditionFile=opt$conditions
transformation=opt$transformation

data1=getGEO(GEOQueryID)
eset=data1[[1]]

#check if datas are in log2 space
normalization<-function(data){
	ex <- exprs(data)
	qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
	LogC <- (qx[5] > 100) ||
			(qx[6]-qx[1] > 50 && qx[2] > 0) ||
			(qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
	if (LogC) { ex[which(ex <= 0)] <- NaN
		return (log2(ex)) } else {
		return (ex)
	}
}

if (transformation=="auto"){
	exprs(eset)=normalization(eset)
} else if (transformation=="yes"){
	exprs(eset)=log2(exprs(eset))			
}

matrixData=exprs(eset)
write.table(matrixData,col.names=NA,row.names=TRUE,sep="\t",file=GEOQueryData)

#Construcion of condition file
#if there is data in "source_name_ch1" field, we keep this data as a condition
#else we keep the "description" field data.
if (length(unique(tolower(pData(data1[[1]])["source_name_ch1"][,1])))>1)
{
	conditions=pData(data1[[1]])["source_name_ch1"]
	description=paste0(as.vector(pData(data1[[1]])["geo_accession"][,1]), " ",as.vector(pData(data1[[1]])["title"][,1]), " ", as.vector(conditions[,1]))
}	else
{
	conditions=pData(data1[[1]])["description"]
	description=paste0(as.vector(pData(data1[[1]])["geo_accession"][,1]), " ",as.vector(pData(data1[[1]])["title"][,1]), " ", as.vector(conditions[,1]))
}

conditions[,1]=tolower(conditions[,1])
pData(eset)["source_name_ch1"]=conditions

write.table(cbind(conditions,description),quote = FALSE,col.names = FALSE, row.names=TRUE,file=conditionFile,sep="\t")
save(eset,conditions,file=GEOQueryRData)