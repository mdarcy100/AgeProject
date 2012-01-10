library(Biobase)

## Caldas Data
caldas.raw <- read.delim(file="E-UCON-1_data_avg_135_raw_bySymbol.TXT",sep="\t")
caldas.exprs <- caldas.raw[1:nrow(caldas.raw),2:ncol(caldas.raw)]
colnames(caldas.exprs) <- gsub("X","caldas_",colnames(caldas.exprs))
rownames(caldas.exprs) <- caldas.raw[1:nrow(caldas.raw),1]

caldas.clin.data <- read.xls("caldas_info_135_with_predictions.XLS",check.names=FALSE)
rownames(caldas.clin.data) <- paste("caldas_",caldas.clin.data$"patient No_",sep="")

caldas.clin.data <- caldas.clin.data[colnames(caldas.exprs),]
all(colnames(caldas.exprs) == rownames(caldas.clin.data))

caldas.phenoData <- new("AnnotatedDataFrame", data = caldas.clin.data)

## Remove the "(includes EG:" annotation from the rownames and add as featureData
geneid <- gsub(" \\(includes EG:([A-Za-z0-9]+)\\)","",rownames(caldas.exprs))
names(geneid) <- rownames(caldas.exprs)
geneid <- data.frame(geneid,stringsAsFactors=FALSE)
geneid <- new("AnnotatedDataFrame", data = geneid)


caldas.eset <- new("ExpressionSet", exprs = as.matrix(caldas.exprs),
                  phenoData = caldas.phenoData,featureData=geneid)

save(caldas.eset,file="caldas_eset.RData")
