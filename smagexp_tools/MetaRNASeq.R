#!/usr/bin/env Rscript
# setup R error handling to go to stderr
options( show.error.messages=F, error = function () { cat( geterrmessage(), file=stderr() ); q( "no", 1, F ) } )

# we need that to not crash galaxy with an UTF8 error on German LC settings.
loc <- Sys.setlocale("LC_MESSAGES", "en_US.UTF-8")

library("optparse")

##### Read options
option_list=list(
		make_option("--input",type="character",default="NULL",help="list of rdata objects containing eset objects"),
		make_option("--fdr",type="character",default=NULL,help="Adjusted p-value threshold to be declared differentially expressed"),
		make_option("--result",type="character",default=NULL,help="text file containing result of the meta-analysis"),
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

suppressPackageStartupMessages(require(affy))
suppressPackageStartupMessages(require(annaffy))
suppressPackageStartupMessages(require(VennDiagram))
suppressPackageStartupMessages(require(GEOquery))

listInput <- trimws( unlist( strsplit(trimws(opt$input), ",") ) )

listfiles=vector()
listfilenames=vector()
nbreplicates=vector()

for (i in 1:length(listInput))
{
	inputFileInfo <- unlist( strsplit( listInput[i], ';' ) )
	listfiles=c(listfiles,inputFileInfo[1])
	listfilenames=c(listfilenames,inputFileInfo[2])
	nbreplicates[i]=as.numeric(inputFileInfo[3])
}


outputfile <- opt$result
result.html = opt$htmloutput
html.files.path=opt$htmloutputpath
result.template=opt$htmltemplate

alpha=as.numeric(opt$fdr)

listData=lapply(listfiles,read.table)
orderData=lapply(listData, function(x) x[order(x[1]), ])
rawpval=lapply(orderData,function(x) x[6])
rawpval=lapply(rawpval, function(x) as.numeric(unlist(x)))

DE=list()
DE=lapply(orderData, function(x) ifelse(x[7]<=0.05,1,0))

FC=list()
FC=lapply(orderData, function(x) x[3])

DE=as.data.frame(DE)
colnames(DE)=listfilenames
FC=as.data.frame(FC)
colnames(FC)=listfilenames

# combinations
library(metaRNASeq)
library(UpSetR)
fishcomb<-fishercomb(rawpval, BHth=alpha)
warning(length(rawpval))
invnormcomb<-invnorm(rawpval, nrep=nbreplicates, BHth=alpha)
#DE[["fishercomb"]]<-ifelse(fishcomb$adjpval<=alpha,1,0)
#DE[["invnormcomb"]]<-ifelse(invnormcomb$adjpval<=alpha,1,0)

signsFC<-mapply(FC,FUN=function(x) sign(x))
sumsigns<-apply(signsFC,1,sum)
commonsgnFC<-ifelse(abs(sumsigns)==dim(signsFC)[2],sign(sumsigns),0)

DEresults <- data.frame(DE=DE,"DE.fishercomb"=ifelse(fishcomb$adjpval<=alpha,1,0),"DE.invnorm"=ifelse(invnormcomb$adjpval<=alpha,1,0))

unionDE <- unique(c(fishcomb$DEindices,invnormcomb$DEindices))
FC.selecDE <- data.frame(DEresults[unionDE,],FC[unionDE,],signFC=commonsgnFC[unionDE])
keepDE <- FC.selecDE[which(abs(FC.selecDE$signFC)==1),]

fishcomb_de <- rownames(keepDE)[which(keepDE[,"DE.fishercomb"]==1)]
invnorm_de <- rownames(keepDE)[which(keepDE[,"DE.invnorm"]==1)]
indstudy_de = list()
for (i in 1:length(listfiles)) {
	currentIndstudy_de = rownames(keepDE)[which(keepDE[,i]==1)]
	indstudy_de[[listfilenames[i]]]=currentIndstudy_de
}

IDDIRRfishcomb=IDD.IRR(fishcomb_de,indstudy_de)
IDDIRRinvnorm=IDD.IRR(invnorm_de,indstudy_de)

conflits<-data.frame(ID=listData[[1]][rownames(DEresults),1],DE=DEresults,FC=FC,signFC=commonsgnFC)
#write DE outputfile
write.table(conflits, outputfile,sep="\t",,row.names=FALSE)
library(VennDiagram)
DE_num=apply(keepDE[,1:(length(listfiles)+2)], 2, FUN=function(x) which(x==1))
#DE_num=apply(DEresults, 2, FUN=function(x) which(x==1))
temp.venn.plot = file.path( html.files.path, paste("venn.png"))
if (length(listfiles)<=2) {
	title="VENN DIAGRAM"
	width=500
	venn.plot<-venn.diagram(x=as.list(DE_num),filename=NULL, col="black", fill=1:length(DE_num)+1,alpha=0.6)
	png(temp.venn.plot,width=width,height=500)
	grid.draw(venn.plot)
	dev.off()
} else {
	title="UPSETR DIAGRAM"
	width=1000
	png(temp.venn.plot,width=width,height=500)
	upset(fromList(as.list(DE_num)), order.by = "freq")
	dev.off()
	
}


library(jsonlite)
matrixConflits=as.matrix(conflits)
datajson=toJSON(matrixConflits,pretty = TRUE)
summaryFishcombjson=toJSON(as.matrix(t(IDDIRRfishcomb)),pretty = TRUE)
summaryinvnormjson=toJSON(as.matrix(t(IDDIRRinvnorm)),pretty = TRUE)

vennFilename="venn.png"
vennFile=file.path(html.files.path,vennFilename)
htmlfile=readChar(result.template, file.info(result.template)$size)
htmlfile=gsub(x=htmlfile,pattern = "###DATAJSON###",replacement = datajson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###FISHSUMMARYJSON###",replacement = summaryFishcombjson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###INVSUMMARYJSON###",replacement = summaryinvnormjson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###VENN###",replacement = vennFilename, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###WIDTH###",replacement = as.character(width), fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###TITLE###",replacement = title, fixed = TRUE)
write(htmlfile,result.html)

