#remove everything after the '.'
#looks like N####.*****
#return N###
truncateID <- function(x){
	pat <- "(N[0-9]+)(.*)"	
	id<-sub(pat, '\\1', x)
	print(id)
	junk<-sub(pat,'\\2',x)	
	print(junk)
	return(id)
}


library(Biobase)
library(impute)
library(hgu133a2.db) # i need to find the correct library to process these samples (should be affy)
setwd("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/")
data.dir <- "./RawData/"

#dataset should already include gene name and entrez id, entrez ids with -1 mean missing, we want to filter those first
## Load Paul's dataset
paul.raw <- read.delim(file="./RawData/NormalOrganoidAnnotatedDat.txt", header=TRUE, sep="\t",check.names=FALSE, stringsAsFactors=FALSE)

dim(paul.raw)
## 25927    44
colnames(paul.raw)

#columns 2-31 = sample data
#column 34 = GeneName
#column 35 = GeneSymbol
#column 37 = EntrezID (what we will match on)

paul.exprs<-as.matrix(paul.raw[,(2:31)])


na.frac.genes <- apply(paul.exprs,1,
                       function(x){sum(is.na(x))/length(x)})
                       
hist(na.frac.genes,500)
dim(paul.exprs)
## 25927    30

#apparently no missing data...
paul.exprs <- paul.exprs[-which(na.frac.genes > .2),]
dim(paul.exprs)
#25927    30

#no missing data, so no need to impute

#set the row names to the entrez ids
entrezID <- paul.raw[(1:length(paul.exprs)),37]

#rownames(paul.exprs) <- entrezID

#read in demographic data
paul.clin <- read.xls("./RawData/ethnicity.xlsx",header=TRUE,
                     check.names=FALSE,
                     stringsAsFactors=FALSE)
 
paul.pData <- paul.clin[,1:ncol(paul.clin)]    
rownames(paul.pData) <- paul.clin$ID  

#make sure the names match
colnames(paul.exprs) <- truncateID(colnames(paul.exprs))               


#make the Annotated data frame with entrez id (37)
paul.fData <- paul.raw[,34:37]
paul.fData<- cbind(paul.raw[,1],paul.fData)
rownames(paul.fData)<-paul.raw[,1] #first row (not sure what this identifier is)

paul.featureData <- new("AnnotatedDataFrame", data = paul.fData)

paul.phenoData <- new("AnnotatedDataFrame", data = paul.pData)
paul.eset <- new("ExpressionSet",exprs = paul.exprs, phenoData= paul.phenoData,
                   featureData = paul.featureData)
                   
dim(paul.eset)
#Features  Samples 
#   25927       30          

#save it!
save(paul.eset,file="paul_eset.RData")