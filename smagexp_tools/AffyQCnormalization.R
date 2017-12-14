#!/usr/bin/env Rscript
# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("optparse")

##### Read options
option_list=list(
		make_option("--input",type="character",default="NULL",help="rdata object containing eset object"),
		make_option("--normalization",type="character",default=NULL,help="normalization method"),
		make_option("--nbresult",type="character",default=NULL,help="number of result displayed results"),
		make_option("--rdataoutput",type="character",default="NULL",help="output rdata object containing eset object"),
		make_option("--htmloutput",type="character",default=NULL,help="Output html report"),
		make_option("--htmloutputpath",type="character",default="NULL",help="Path of output html report"),
		make_option("--htmltemplate",type="character",default=NULL,help="html template)")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if(is.null(opt$input)){
	print_help(opt_parser)
	stop("input required.", call.=FALSE)
}

#loading libraries

suppressPackageStartupMessages(require(Biobase))
suppressPackageStartupMessages(require(GEOquery))
suppressPackageStartupMessages(require(GEOmetadb))
suppressPackageStartupMessages(require(limma))
suppressPackageStartupMessages(require(jsonlite))
suppressPackageStartupMessages(require(affy))
suppressPackageStartupMessages(require(affyPLM))
suppressPackageStartupMessages(require(dplyr))

listInput <- trimws( unlist( strsplit(trimws(opt$input), ",") ) )

celList=vector()
celFileNameList=vector()

for (i in 1:length(listInput))
{
	inputFileInfo <- unlist( strsplit( listInput[i], ';' ) )
	celList=c(celList,inputFileInfo[1])
	celFileNameList=c(celFileNameList,inputFileInfo[2])
}


normalization=opt$normalization
result_export_eset=opt$rdataoutput
result=opt$htmloutput
result.path=opt$htmloutputpath
result.template=opt$htmltemplate

dir.create(result.path, showWarnings = TRUE, recursive = TRUE)
for(i in 1:length(celList))
{
	file.copy(celList[i],paste0("./",celFileNameList[i]))
}

data <- ReadAffy(filenames=celFileNameList, celfile.path=".")
htmlfile=readChar(result.template, file.info(result.template)$size)

boxplot="boxplot.png"
png(boxplot,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(data,las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOT###",replacement = boxplot, fixed = TRUE)
file.copy(boxplot,result.path)

images="images.png"
nblines=length(celList)%/%4 + as.numeric((length(celList)%%4)!=0) 
png(images,width=800,height = 200*nblines)
par(mfrow=c(nblines,4))
image(data)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###IMAGES###",replacement = images, fixed = TRUE)
file.copy(images,result.path)


plotMA="plotMA.png"
nblines=length(celList)%/%3 + as.numeric((length(celList)%%3)!=0) 
png(plotMA,width=800,height =300*nblines )
par(mfrow=c(nblines,3))
MAplot(data)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###PLOTMA###",replacement = plotMA, fixed = TRUE)
file.copy(plotMA,result.path)


if (normalization == "rma") {
	eset <- rma(data)
} else if (normalization == "quantile") {
	eset = rma(data,background = FALSE,normalize = TRUE)
} else if (normalization == "background"){
	eset = rma(data,background = TRUE ,normalize = FALSE)
} else if (normalization == "log2") {
	eset = rma(data,background = FALSE ,normalize = FALSE)
}
	

boxplotnorm="boxplotnorm.png"
png(boxplotnorm,width=800,height = 400)
par(mar=c(7,5,1,1))
boxplot(data.frame(exprs(eset)),las=2,outline=FALSE)
dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###BOXPLOTNORM###",replacement = boxplotnorm, fixed = TRUE)
file.copy(boxplotnorm,result.path)

plotMAnorm="plotMAnorm.png"
nblines=length(celList)%/%3 + as.numeric((length(celList)%%3)!=0) 
png(plotMAnorm,width=800,height =300*nblines )
par(mfrow=c(nblines,3))
#for (i in 1:length(celList)){
MAplot(eset)
#}

dev.off()
htmlfile=gsub(x=htmlfile,pattern = "###PLOTMANORM###",replacement = plotMAnorm, fixed = TRUE)
file.copy(plotMAnorm,result.path)
save(eset,file=result_export_eset)
write(htmlfile,result)

