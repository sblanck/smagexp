cargs <- commandArgs()
cargs <- cargs[(which(cargs == "--args")+1):length(cargs)]
nbargs=length(cargs)
listfiles=vector()
listfilenames=vector()
for (i in seq(1,nbargs-6,2)) {
	listfiles=c(listfiles,cargs[[i]])
	listfilenames=c(listfilenames,cargs[[i+1]])
}
#mod<-cargs[[length(cargs) - 6]]
outputfile <- cargs[[length(cargs) - 5]]
result.html = cargs[[length(cargs) - 4]]
html.files.path=cargs[[length(cargs) - 3]]
result.template=cargs[[length(cargs) - 2]]

alpha=0.05

#print(comparison)

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
# the comparison must only have two values and the conds must
# be a vector from those values, at least one of each.

#if (length(comparison) != 2) {
#	stop("Comparison type must be a tuple: ", cargs[length(cargs) - 8])
#}

sink("/dev/null")
dir.create(html.files.path, recursive=TRUE)
#library(DESeq)
#library(HTSFilter)

#DE=list()
#FC=list()
#i=1

# Open the html output file
#file.conn = file(diag.html, open="w")

#writeLines( c("<html><body>"), file.conn)

# Perform deseq analysis on each study
#for(i in 1:length(listfiles))
#{
#  f=listfiles[i]
#  fname=listfilenames[i]
#  study_name=unlist(strsplit(fname,"[.]"))[1]
#  print(paste0("study.name ",study_name))
#  d <- read.table(f, sep=" ", header=TRUE, row.names=1)
#  conds<-sapply(strsplit(colnames(d),"[.]"),FUN=function(x) x[1])
#  if (length(unique(conds)) != 2) {
#  	warning(as.data.frame(strsplit(colnames(d),"[.]")))
#  	stop("You can only have two columns types: ", paste(conds,collapse=" "))
#  }
#  if (!identical(sort(comparison), sort(unique(conds)))) {
#  	stop("Column types must use the two names from Comparison type, and vice versa.  Must have at least one of each in the Column types.\nColumn types: ", cargs[2], "\n", "Comparison type: ", cargs[3])
#  }
#  if (length(d) != length(conds)) {
#  	stop("Number of total sample columns in counts file must correspond to the columns types field.  E.g. if column types is 'kidney,kidney,liver,liver' then number of sample columns in counts file must be 4 as well.")
#  }
#  
#  cds <- newCountDataSet(d, conds)
#  cds <- estimateSizeFactors(cds)
#  
#  cdsBlind <- estimateDispersions( cds, method="blind" )
#  
#  if (length(conds) != 2) {
#  	cds <- estimateDispersions( cds )
#  	norep = FALSE
#  }
#  
#  if (length(conds) == 2) {
#  	cds <- estimateDispersions( cds, method=method, sharingMode=mod, fitType="parametric" )
#  	norep = TRUE
#  }
#  
#  filter<-HTSFilter(cds, plot=FALSE)
#  cds.filter<-filter$filteredData
#  on.index<-which(filter$on==1)
#
#  res<-as.data.frame(matrix(NA,nrow=nrow(cds),ncol=ncol(cds)))
#  nbT <- nbinomTest(cds.filter, comparison[1], comparison[2])
#  colnames(res)<-colnames(nbT)
#  res[on.index,]<-nbT
#  #write.table(res[order(res$padj), ], file=outputfile, quote=FALSE, row.names=FALSE, sep="\t")
#  
#
#  temp.pval.plot = file.path( html.files.path, paste("PvalHist",i,".png",sep=""))
#  png( temp.pval.plot, width=500, height=500 )
#  hist(res$pval, breaks=100, col="skyblue", border="slateblue", main="")
#  dev.off()
#  
#  writeLines( c("<h2>P-value histogram for ",study_name,"</h2>"), file.conn)
#  writeLines( c("<img src='PvalHist",i,".png'><br/><br/>"), file.conn)
#  
#  #on enregistre la p-value
#  rawpval[[study_name]]<-res$pval
#  DE[[study_name]]<-ifelse(res$padj<=alpha,1,0)
#  FC[[study_name]]<-res$log2FoldChange 
#  
#  i=i+1
#}


# combinations
library(metaRNASeq)
fishcomb<-fishercomb(rawpval, BHth=alpha)
warning(length(rawpval))
invnormcomb<-invnorm(rawpval, nrep=c(8,8), BHth=alpha)
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

#conflits<-data.frame(ID=listData[[1]][rownames(DEresults),1],Fishercomb=DEresults[["DE.fishercomb"]],Invnormcomb=DEresults[["DE.invnorm"]],sign=commonsgnFC)
conflits<-data.frame(ID=listData[[1]][rownames(DEresults),1],DE=DEresults,FC=FC,signFC=commonsgnFC)
#write DE outputfile
write.table(conflits, outputfile,sep="\t",,row.names=FALSE)
library(VennDiagram)
DE_num=apply(DEresults, 2, FUN=function(x) which(x==1))
venn.plot<-venn.diagram(x=as.list(DE_num),filename=NULL, col="black", fill=1:length(DE_num)+1,alpha=0.6)
temp.venn.plot = file.path( html.files.path, paste("venn.png"))
png(temp.venn.plot,width=500,height=500)
grid.draw(venn.plot)
dev.off()

library(jsonlite)
matrixConflits=as.matrix(conflits)
datajson=toJSON(matrixConflits,pretty = TRUE)
summaryFishcombjson=toJSON(as.matrix(t(IDDIRRfishcomb)),pretty = TRUE)
summaryinvnormjson=toJSON(as.matrix(t(IDDIRRinvnorm)),pretty = TRUE)


#vennsplit=strsplit(result.venn,split="/")[[1]]
#venn=paste0("./",vennsplit[length(vennsplit)])


vennFilename="venn.png"
vennFile=file.path(html.files.path,vennFilename)
htmlfile=readChar(result.template, file.info(result.template)$size)
htmlfile=gsub(x=htmlfile,pattern = "###DATAJSON###",replacement = datajson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###FISHSUMMARYJSON###",replacement = summaryFishcombjson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###INVSUMMARYJSON###",replacement = summaryinvnormjson, fixed = TRUE)
htmlfile=gsub(x=htmlfile,pattern = "###VENN###",replacement = vennFilename, fixed = TRUE)
write(htmlfile,result.html)

#library(VennDiagram)
#flog.threshold(ERROR)
#
##venn.plot<-venn.diagram(x = c(res[c(1:(length(res)-3))],meta=list(res$Meta)),filename = v, col = "black", fill = c(1:(length(res)-2)), margin=0.05, alpha = 0.6,imagetype = "png")
#dir.create(result.path, showWarnings = TRUE, recursive = FALSE)
#
#showVenn<-function(liste,file)
#{
#	venn.plot<-venn.diagram(x = liste,
#			filename = vennFilename, col = "black",
#			fill = 1:length(liste)+1,
#			margin=0.05, alpha = 0.6,imagetype = "png")
##	png(file);
##	grid.draw(venn.plot);
##	dev.off();
#	
#}
#
#l=list()
#for(i in 1:length(esets))
#{
#	l[[paste("study",i,sep="")]]<-res[[i]]
#}
#l[["Meta"]]=res[[length(res)-1]]
#showVenn(l,vennFile)
#file.copy(vennFilename,result.path)


#writeLines( c("<h2>Venn Plot</h2>"), file.conn)
#writeLines( c("<img src='venn.png'><br/><br/>"), file.conn)
#writeLines( c("</body></html>"), file.conn)
#close(file.conn)
#print("passe6")
#sink(NULL)
