#!/usr/bin/env Rscript
# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("optparse")

##### Read options
option_list=list(
		make_option("--input",type="character",default="NULL",help="rdata object containing eset object"),
		make_option("--conditions",type="character",default=NULL,help="Text file summarizing conditions of the experiment (required)"),
		make_option("--normalization",type="character",default=NULL,help="log2 transformation"),
		make_option("--annotations",type="character",default="NULL",help="rdata object containing eset object"),
		make_option("--rdataoutput",type="character",default="NULL",help="output rdata object containing eset object"),
		make_option("--htmloutput",type="character",default=NULL,help="Output html report"),
		make_option("--htmloutputpath",type="character",default="NULL",help="Path of output html report"),
		make_option("--htmltemplate",type="character",default="NULL",help="html template)")


);


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input)){
	print_help(opt_parser)
	stop("input required.", call.=FALSE)
}

if(is.null(opt$conditions)){
	print_help(opt_parser)
	stop("conditions input required.", call.=FALSE)
}


#loading libraries
suppressPackageStartupMessages(require(GEOquery))

suppressPackageStartupMessages(require(Biobase))
suppressPackageStartupMessages(require(GEOquery))
suppressPackageStartupMessages(require(GEOmetadb))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(jsonlite))
suppressPackageStartupMessages(require(affy))
suppressPackageStartupMessages(require(dplyr))
suppressPackageStartupMessages(require(affyPLM))

dataFile=opt$input
normalization=opt$normalization
conditionsFile=opt$conditions
annotation=opt$annotations
result_export_eset=opt$rdataoutput
result=opt$htmloutput
result.path=opt$htmloutputpath
result.template=opt$htmltemplate

dir.create(result.path, showWarnings = TRUE, recursive = FALSE)

data=as.matrix(read.table(file = dataFile,row.names=1,header=TRUE))
#conditions=read.table(file=conditionsFile,sep = "\t",row.names=1)
htmlfile=readChar(result.template, file.info(result.template)$size)

colnames(conditions)=c("source_name_ch1","description")
#phenodata<-new("AnnotatedDataFrame",data=conditions)

eset=ExpressionSet(assayData=data,annotation=annotation)

if (normalization == "quantile") {
	eset <- normalize.ExpressionSet.quantiles(eset, transfn="log")
} else if (normalization == "log2") {
	exprs(eset) = log2(exprs(eset)) 
} 

boxplotnorm="boxplotnorm.png"
png(boxplotnorm,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(data.frame(exprs(eset)),las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOTNORM###",replacement = boxplotnorm, fixed = TRUE)
file.copy(boxplotnorm,result.path)

plotMAnorm="plotMAnorm.png"
nblines=length(colnames(data))%/%3 + as.numeric((length(colnames(data))%%3)!=0) 
png(plotMAnorm,width=800,height =300*nblines )
par(mfrow=c(nblines,3))
##for (i in 1:length(colnames(data))){
	MAplot(eset)
#}

dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###PLOTMANORM###",replacement = plotMAnorm, fixed = TRUE)
file.copy(plotMAnorm,result.path)
#write.table(tolower(c(condition1Name,condition2Name)),quote = FALSE,col.names = FALSE, row.names=FALSE,file=result_export_conditions)
#saveConditions=c(condition1Name,condition2Name)
save(eset,file=result_export_eset)
write(htmlfile,result)

#l=list()
#for(i in 1:length(esets))
#{
#	l[[paste("study",i,sep="")]]<-res[[i]]
#}
#l[["Meta"]]=res[[length(res)-1]]
#showVenn(res,file.path(temp.files.path,"venn.png"))
#writeLines(c("<h2>Venn diagram</h2>"),file.conn)
#writeLines(c("<img src='venn.png'><br/><br/>"),file.conn)
#writeLines(c("</body></html>"),file.conn)
#close(file.conn)