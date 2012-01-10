library(Biobase)
library(impute)

## Load NKI295 dataset
nki295.raw <- read.delim(file="e:/UNC09/NKI295_Data/NKI295_LogRatio_data.TXT", header=TRUE,
                  sep="\t",check.names=FALSE, stringsAsFactors=FALSE)
dim(nki295.raw)
## [1] 24481   305

## First 10 columns are:
## [1] "ID"            "Systemic_name" "probe_seq"     "geneid"
## [5] "symbol"        "accession"     "refseq"        "chromosome"
## [9] "map_location"  "description"
## Systematic_name, and ID are unique indicators - use Systematic_name as rownames for expression matrix
## Samples begin at column 11, all names are unique
length(colnames(nki295.raw[11:ncol(nki295.raw)])) == length(unique(colnames(nki295.raw[11:ncol(nki295.raw)])))

nki295.exprs <-  as.matrix(nki295.raw[1:nrow(nki295.raw),11:ncol(nki295.raw)])
rownames(nki295.exprs) <- nki295.raw[,"Systemic_name"]

## Remove flaged spots
load("nki295_gene_flags.RData")
dim(flag.data)
colnames(flag.data) <- gsub("Sample ","",colnames(flag.data))
all(rownames(flag.data) == rownames(nki295.exprs))
all(colnames(flag.data) == colnames(nki295.exprs))

nki295.exprs <- nki295.exprs+flag.data

na.frac.genes <- apply(nki295.exprs,1,
                       function(x){sum(is.na(x))/length(x)})
hist(na.frac.genes,500)
dim(nki295.exprs)
## [1] 24481   295
nki295.exprs <- nki295.exprs[-which(na.frac.genes > .2),]
dim(nki295.exprs)

na.frac.arrays <- apply(nki295.exprs,2,
                       function(x){sum(is.na(x))/length(x)})
hist(na.frac.arrays)

## impute the NA's introduced from the flagged spots
sum(is.na(nki295.exprs))/length(nki295.exprs) ## about .5%, impute
nki295.exprs.imputed <- impute.knn(nki295.exprs,k=10, rng.seed=98765)$data



nki295.clin <- read.xls("e:/UNC09/NKI295_Data/Clinical_Data_Supplement295-chang.xls",
                     check.names=FALSE,
                     stringsAsFactors=FALSE)

nki295.pData <- nki295.clin[,2:ncol(nki295.clin)]
rownames(nki295.pData) <- nki295.clin$ID
nki295.pData <- nki295.pData[colnames(nki295.exprs.imputed),]
all(rownames(nki295.pData) == colnames(nki295.exprs.imputed))

nki295.fData <- nki295.raw[,1:9]
nki295.fData <- sapply(nki295.fData,function(x){x[grep("^---$",x)] <- NA;x})
nki295.fData <- data.frame(nki295.fData,stringsAsFactors=FALSE)
rownames(nki295.fData) <- nki295.raw[,"Systemic_name"]
nki295.fData <- nki295.fData[rownames(nki295.exprs.imputed),]

nki295.featureData <- new("AnnotatedDataFrame", data = nki295.fData)

nki295.phenoData <- new("AnnotatedDataFrame", data = nki295.pData)


nki295.eset <- new("ExpressionSet",exprs = nki295.exprs.imputed, phenoData= nki295.phenoData,
                   featureData = nki295.featureData)
save(nki295.eset,file="nki295_eset.RData")
