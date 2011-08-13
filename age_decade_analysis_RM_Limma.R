# Originally done by Jason Pirone in Fall 2010
#################################################################################
##
## Set working and data directories
##
#################################################################################

setwd("/Users/mdarcy100/Desktop/MTroester/AgePaper/AgeManuscript_062011/")
data.dir <- "./RawData/"

#################################################################################
##
## Load required libraries
##
#################################################################################
#MD: needed to install biobase prior to using the library Biobase
#source("http://bioconductor.org/biocLite.R")
#biocLite("Biobase")
#########

library(Biobase)
## library(genefilter)
library(hgug4112a.db)
library(gdata)
## library(pamr)
library(qvalue)
## library(cluster)
library(survival)
library(gplots)
library(limma)
## library(samr)
library(RColorBrewer)
## library(GOstats)
## library(clValid)

##MD: can't seem to install made4 package - get error with the first command.  the second works
#install.packages("made4", repos = "http://www.bioconductor.org")
#need to install made4 prior to using it
#source("http://bioconductor.org/biocLite.R")
#biocLite("made4")
library(made4)
###MD: needed to also install mclust, flexmix, modeltools,multcomp,mvtnorm for fpc to install
library(fpc)

source("age_paper_scripts2.R")

#################################################################################
##
## Load processed microarray data sets
##
#################################################################################

load("all_rm_eset.RData") ## process_All_RMsamples.R)
load("caldas_eset.RData") ## process_CaldasSamples.RRe
load("nki295_eset.RData") ## process_NKI295samples.R

#################################################################################
##
## Analyze the RM dataset to find genes that show changes with age
##
#################################################################################

#################################################################################
#   MD NOTE FOR HERSELF:  ALL THIS CODE USES A DATASTRUCTURE called EXPRESSIONSET
####  (ExpressionSet {Biobase})
#################################################################################
## Restrict analysis to subjects with age >= 20 who are pre- or peri-
## menopausal (i.e., menopause = 1 or 2)

preperi.age.eset <- all.rm.eset[,all.rm.eset$age >= 20 &
                                all.rm.eset$menopause %in% c(1,2)]
dim(preperi.age.eset)
## Features  Samples
##    31456       76

## Basic prefiltering
## Reduce probes to those having an ENTREZ Id
entrezIds <- mget(featureNames(preperi.age.eset), envir=hgug4112aENTREZID)
haveEntrezIds <- names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(preperi.age.eset)) - length(haveEntrezIds)

age.eset <- preperi.age.eset[haveEntrezIds,]

dim(age.eset)
##Below were Jason's numbers.
## Features  Samples
##    23729       76

## Monica had more since there are more ids in July 2011
#Features  Samples 
#   23941       76 



## Create a new variable 'age.decade' that discretizes age by decade
## (age.decade = 2 if 20 <= age < 30, etc.)

age.decade<- cut(age.eset$age,
                breaks = seq(10,(max(age.eset$age)+10)-max(age.eset$age)%%10,by=10),
                include.lowest=TRUE,labels=FALSE,right=FALSE)

pData(age.eset)$age.decade <- age.decade

## Determine the interquartile range for each gene - will be used to filter results
## after applying the empirical Bayes step in Limma
## With limma

iqr.age.eset <- esApply(age.eset,1,IQR)



## Regression analysis using Limma

## Regression on age.decade
design.age.decade <- model.matrix(~ age.eset$age.decade)
colnames(design.age.decade) <- c("Intercept","age.decade.limma")

age.decade.limma.fit1 <- lmFit(age.eset,design.age.decade)
age.decade.limma.ebayes <- eBayes(age.decade.limma.fit1)

## Apply IQR cutoff of median(iqr.age.eset)
age.decade.limma.filtered <- age.decade.limma.ebayes[iqr.age.eset>median(iqr.age.eset),]
##MD checking dimension
dim(age.decade.limma.filtered)
#[1] 11970     2

#this is a summary of the results of the regression - much nicer than what Monica did
age.decade.limma.summary <- topTable(age.decade.limma.filtered,coef="age.decade.limma",adjust.method="fdr",num=Inf)

##MD did this to figure out what the data structure looked like:
#colnames(age.decade.limma.summary)
#[1] "ID"        "logFC"     "AveExpr"   "t"         "P.Value"   "adj.P.Val" "B"        
 

## Distribution of p-values
pval.summary(age.decade.limma.summary$P.Value)
## JASON's values below:
## <1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1
##     13    120   1075   2154   3293   4748  11864

##Monica's values slightly different:
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1 
#    14    119   1086   2179   3323   4797  11970 

pval.summary(age.decade.limma.summary$adj.P.Val)
## Jason's values:
## <1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1
##      0      0      0      2      2    770  11864

## Monica's values slightly different
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1 
#     0      0      0      2      2    802  11970 

## Compute qvalue if desired, however the correctness of using this method to
## Validity of computing q-value based on filtered p-values?
q.value.age.decade <- qvalue(p=age.decade.limma.summary$P.Value)
summary(q.value.age.decade)
## Jason's values below
##         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value     13    120  1075   2154  3293 4748 11864
## q-value      0      0     2      2  1726 4797 11864
##Monica's values slightly different
#        <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#	p-value     14    119  1086   2179  3323 4797 11970
#	q-value      0      0     2      9  1744 4855 11970

age.decade.limma.summary$qvalue <- q.value.age.decade$qvalue

plot(age.decade.limma.summary$logFC,
     -log10(age.decade.limma.summary$P.Value),pch=19,col='blue',
     xlab='Regression Coefficient',ylab='log10(raw P value)')
abline(h=-log10(.001),col='green')
abline(h=-log10(.05),col='purple')
abline(h=-log10(.1),col='orange')

## Regression on age
design.age <- model.matrix(~ age.eset$age)
colnames(design.age) <- c("Intercept","age.limma")

age.limma.fit1 <- lmFit(age.eset,design.age)
age.limma.ebayes <- eBayes(age.limma.fit1)

age.limma.filtered <- age.limma.ebayes[iqr.age.eset>median(iqr.age.eset),]
age.limma.summary <- topTable(age.limma.filtered,coef="age.limma",adjust.method="fdr",num=Inf)


## Distribution of p-values
pval.summary(age.limma.summary$P.Value)
## <1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1
##      5     41    574   1391   2465   3925  11864
##Monica's values:
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1 
#     5     42    582   1406   2496   3964  11970 

pval.summary(age.limma.summary$adj.P.Val)
## <1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1
##      0      0      0      0      0      0  11864
##Monica's values:
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1 
#     0      0      0      0      0      0  11970 

q.value.age <- qvalue(p=age.limma.summary$P.Value)
summary(q.value.age)

##         <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
## p-value      5     41   574   1391  2465 3925 11864
## q-value      0      0     0      0     0 1509 11864

##Monica's values:
#        <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#p-value      5     42   582   1406  2496 3964 11970
#q-value      0      0     0      0     0 1682 11970

age.limma.summary$qvalue <- q.value.age$qvalue

plot(age.limma.summary$logFC,-log10(age.limma.summary$P.Value))

#################################################################################
##
## Make a gene list from the resuls of the regression on age.decade
##
#################################################################################

## Select genes that show a greater than 1 fold change and have an
## adjusted p-value less than 0.1, there should be fewer results than in the non-restricted

sig.age.decade.fold <- age.decade.limma.summary[abs(age.decade.limma.summary$logFC*diff(range(age.eset$age.decade)))>=1.0
                                    & age.decade.limma.summary$adj.P.Val < 0.1,]
                                    
sig.age.decade.fold.probe <- sig.age.decade.fold$ID
sig.age.decade.fold.entrez <- unlist(mget(sig.age.decade.fold$ID,env=hgug4112aENTREZID))
sig.age.decade.fold.symbol <- unlist(mget(sig.age.decade.fold$ID,env=hgug4112aSYMBOL))



###MD added - write out a file of the decade gene list
toWrite_SigAge<-cbind(sig.age.decade.fold.probe,sig.age.decade.fold.entrez,sig.age.decade.fold.symbol)
write.table(toWrite_SigAge,file="MD_SigGenes_080311.txt",sep="\t",col.names=NA)

#MD - try again without the fold change parameter because I had gar fewer results than the 
# original list
sig.age.decade <- age.decade.limma.summary[age.decade.limma.summary$adj.P.Val < 0.1,]
sig.age.decade.probe <- sig.age.decade$ID
sig.age.decade.entrez <- unlist(mget(sig.age.decade$ID,env=hgug4112aENTREZID))
sig.age.decade.symbol <- unlist(mget(sig.age.decade$ID,env=hgug4112aSYMBOL))



toWrite_SigAgeAllFold<-cbind(sig.age.decade.probe,sig.age.decade.entrez,sig.age.decade.symbol)
write.table(toWrite_SigAgeAllFold,file="MD_SigGenesAllFold_080311.txt",sep="\t",col.names=NA)


## Plot a few genes
## Green circles are log2(expr), orange line is the fitted regression model, and the blue line
## is a loess fit. For the plots, age, not age.decade, was used.

par(mfrow=c(4,4))
for(i in 33:(33+15)){
    y <- exprs(age.eset)[sig.age.decade$ID[i],]-median(exprs(age.eset)[sig.age.decade$ID[i],])
    plot(age.eset$age,y,pch=19,col='green',
         xlab="Age (years)",ylab="log2(exprs)",main=unlist(mget(sig.age.decade$ID[i],env=hgug4112aSYMBOL)))
    abline(lm(y~age.eset$age),col='orange')
    lines(loess.smooth(age.eset$age,y),col="blue")
}


#################################################################################
##
## Visualize Cluster/partition the samples
##
#################################################################################

## Get expressions corresponding to genes determined above
## Get expression with genes that have fold change
# sig.age.decade.entrez
exprs.age.clust <- exprs(age.eset)[sig.age.decade.fold$ID,]

exprs.age.clust.nki <- exprs(age.eset)[sig.age.decade.fold$ID,]

#############################################################################################
### MD: SUMMER 2011 
###   Save files for future use and manipulation when doing Creighton method
###   Write clusterable files out, ordered by age
############################################################################################
sig.fold.rm.exprs <- age.eset[sig.age.decade.fold$ID,] #only includes those with at least a certain fold change (~165 in Aug 2011)
sig.full.rm.exprs <-age.eset[sig.age.decade$ID,] #includes all 10% FDR genes (~802 in Aug 2011)

save(sig.fold.rm.exprs,file="sig_fold_rm_eset.RData")
save(sig.full.rm.exprs,file="sig_full_rm_eset.RData")##

#order age before writing out to the file
sig.fold.rm.exprs.order<-exprs(sig.fold.rm.exprs[,order(age.eset$age)])
sig.full.rm.exprs.order<-exprs(sig.full.rm.exprs[,order(age.eset$age)])

colname_age_order<-paste("Age:",sig.fold.rm.exprs[,order(age.eset$age)]$age)
colnames(sig.fold.rm.exprs.order)<-colname_age_order
colnames(sig.full.rm.exprs.order)<-colname_age_order

sig_fold_identifier<-paste(sig.age.decade.fold.probe,sig.age.decade.fold.entrez,sig.age.decade.fold.symbol,sep="|")
sig_full_identifier<-paste(sig.age.decade.probe,sig.age.decade.entrez,sig.age.decade.symbol,sep="|")

clusterFoldFile<-cbind(sig_fold_identifier,as.matrix(sig.fold.rm.exprs.order))
clusterFullFile<-cbind(sig_full_identifier,as.matrix(sig.full.rm.exprs.order))
write.table(clusterFoldFile,file=paste("ClusterableAgeSigFold_", Sys.Date(),".txt",sep=""),sep="\t",col.names=NA)
write.table(clusterFullFile,file=paste("ClusterableAgeSigFull_", Sys.Date(),".txt",sep=""),,sep="\t",col.names=NA)

############################################################################################
## Age Cluster
medians.exprs.age.clust <- apply(exprs.age.clust,1,median)
exprs.age.clust.cent <- t(scale(t(exprs.age.clust),medians.exprs.age.clust))

age.clust.1 <- heatplot(exprs.age.clust.cent,classvec=age.eset$age.decade,
                        scale='none',method='centroid',returnSampleTree=TRUE)

## Make a matrix with symbols as rownames to use for matching
exprs.age.clust.cent.symb <- exprs.age.clust.cent
rownames(exprs.age.clust.cent.symb) <- sig.age.decade.symbol
col.cols <- ifelse(age.eset$age<25,'yellow','blue')


age.order.ind <- sort.int(age.eset$age,index.return=TRUE)$ix
exprs.age.order <- exprs.age.clust[,age.order.ind]
age.decade.order <- age.eset$age.decade[age.order.ind]
age.cols <- brewer.pal(9,"BuPu")[6:9]
names(age.cols) <- as.character(unique(age.decade.order))
col.cols <- age.cols[as.character(age.decade.order)]

heatmap.2(exprs.age.order,
          Colv="none",
          dendrogram="row",
          scale='row',
          col=greenred(32),
          trace="none",
          ColSideColors=col.cols)

################################################################################
##
## Use the new gene list to predict on the NKI295 data
##
################################################################################

## Prepare NKI295 data
## The nki295 expression ratios are most likely (need to verify) log10

## Determine which probes do not have an Entrez ID either NAs or ''
## more of a problem with 'symbol'

dim(nki295.eset)
nki295.eset <- nki295.eset[!is.na(featureData(nki295.eset)$geneid),]
dim(nki295.eset)



if(any(featureData(nki295.eset)$symbol == "")){
    nki295.eset <- nki295.eset[-which(featureData(nki295.eset)$symbol==""),]
}
dim(nki295.eset)

## Median center the NKI295 data
medians.nki295 <- esApply(nki295.eset,1,median)
nki295.exprs.scale <- t(scale(t(exprs(nki295.eset)),center=medians.nki295,scale=FALSE))

## Create a new eset
nki295.scale.eset <- nki295.eset
exprs(nki295.scale.eset) <- nki295.exprs.scale
dim(nki295.scale.eset)

## Remove low variability probes - this may not be necessary
iqr.val.nki295.scale <- esApply(nki295.scale.eset,1,function(x){IQR(x,na.rm = TRUE)})

hist(iqr.val.nki295.scale,sqrt(length(iqr.val.nki295.scale)))


high.iqr.inds.nki295 <- which(iqr.val.nki295.scale > median(iqr.val.nki295.scale))

nki295.scale.eset <- nki295.scale.eset[high.iqr.inds.nki295,]
dim(nki295.scale.eset)
# MD did this
#Features  Samples 
#    8800      295 


high.iqr.probes.nki295 <- mod.findLargest(featureNames(nki295.scale.eset),
                                          iqr.val.nki295.scale[high.iqr.inds.nki295],
                                         featureData(nki295.scale.eset)$symbol)

#MD did this:
#length(high.iqr.probes.nki295)
# 8034

nki295.scale.eset <- nki295.scale.eset[high.iqr.probes.nki295,]
dim(nki295.scale.eset)
#MD did this:
#Features  Samples 
#    8034      295 

## Remove genes with more than one value in the geneid -
## This is very confusing. In some cases a probe matches to two different symbols/ids,
## both of which are distinct and valid. In other
## cases one or more of the ids/symbols are invalid.
## To avoid confusion, conduct analysis without the mult. ids (approx 317)
all(grep("\\|",featureData(nki295.scale.eset)$symbol) == grep("\\|",featureData(nki295.scale.eset)$geneid))

nki295.mult.inds <- grep("\\|",featureData(nki295.scale.eset)$symbol)
nki295.mult.ids <- cbind(featureData(nki295.scale.eset)$geneid[nki295.mult.inds],featureData(nki295.scale.eset)$symbol[nki295.mult.inds])

nki295.scale.eset <- nki295.scale.eset[-nki295.mult.inds,]
dim(nki295.scale.eset)
# MD results:
#Features  Samples 
#    7717      295 

## Map probes found from Limma analysis of RM data to entrez ids for the NKI295 data
## I have not fully QC'd the matching - could be a better way of doing this.
lm.age.list.entrez <- na.omit(sig.age.decade.fold.entrez)
lm.age.list.symbol <- na.omit(sig.age.decade.fold.symbol)


nki.match.vec <- vector("numeric",length=dim(nki295.scale.eset)["Features"])
names(nki.match.vec) <- featureNames(nki295.scale.eset)

for(e.id in unique(lm.age.list.entrez)){
    e.id.pattern <- paste("(^",e.id,"$|\\|",e.id,"$|^",e.id,"\\|$|\\|",e.id,"\\|)",sep="")
    nki.match.ind <- grepl(e.id.pattern,featureData(nki295.scale.eset)$geneid)
    if(any(nki.match.ind)){
        nki.match.vec[nki.match.ind] <- nki.match.vec[nki.match.ind]+1
    }
}

nki.gene.match.full <- names(nki.match.vec[nki.match.vec !=0])
nkiSub.eset <- nki295.scale.eset[nki.gene.match.full,]

#### MD testing something  exprs(age.eset)[sig.age.decade.fold$ID,]
dim(nkiSub.eset)
###MD results:
#Features  Samples 
#      67      295 


###### AFTER we have compared the NKI data with signature, then we can find subset to genes also in the caldas dataset
	## Which genes are also in the caldas eset?
	## Caldas does not have Entrez Ids. Match by symbol
	common.genes.nki.caldas <- intersect(featureNames(caldas.eset),featureData(nkiSub.eset)$symbol)
	##MD:
	length(common.genes.nki.caldas)
	different.genes.nki.caldas <- setdiff(featureNames(caldas.eset),featureData(nkiSub.eset)$symbol)
	##MD:
	length(different.genes.nki.caldas)
	nkiSub.eset <- nkiSub.eset[featureData(nkiSub.eset)$symbol %in% common.genes.nki.caldas,]
	dim(nkiSub.eset)
	##MD: 8/11
	#Features  Samples 
	#      52      295


 

## MD did this
### featureData(nkiSub.eset)$symbol
# [1] "ABCA1"   "ABCA8"   "ABLIM2"  "AKAP12"  "ALDH1A1" "APBB1IP" "BCAT1"   "BST1"    "CD28"    "CDKN2A"  "CDO1"    "CISH"    "CRTAP"  
#[14] "CTSG"    "CXCL14"  "CYP4B1"  "DDR1"    "DIO3"    "DLC1"    "DOCK11"  "EPB41L3" "FGF18"   "FZD4"    "GHR"     "IER5"    "JAK2"   
#[27] "LAMA2"   "LUM"     "MS4A2"   "MS4A4A"  "NOX4"    "NPR3"    "NPY1R"   "NPY5R"   "NT5E"    "NXN"     "ODZ3"    "PALMD"   "PDE4DIP"
#[40] "PGM5"    "PLOD2"   "RBP4"    "RHOQ"    "RSAD2"   "SDCBP"   "SIX1"    "SNX10"   "SP5"     "SVEP1"   "TNFSF11" "TRAF4"   "WISP2"  


################################################################################
##
## Cluster the NKI295 data using the RM derived gene list
##
################################################################################

nki.pam.2.daisy <- pam(daisy(t(exprs(nkiSub.eset))),2)
## nki.pam.2.dist <- pam(as.dist(1-cor((exprs(nkiSub.eset)))),2)
## nki.pam.eucl <- pam(t(nkiSub.exprs.scale),2)
clusplot(nki.pam.2.daisy, labels = 2, col.p = nki.pam.2.daisy$clustering)

all(colnames(nkiSub.eset) ==  names(nki.pam.2.daisy$clustering))

colnames(exprs(nkiSub.eset))

nki.clust.ind <- sort.int(nki.pam.2.daisy$clustering, index.return=TRUE)$ix

pData(nkiSub.eset)$clustering <- nki.pam.2.daisy$clustering

###MD is getting an error below;  will try with nki.clust.ind instead of clust.ind
#nki.pam.clust.sort <- names(nki.pam.2.daisy$clustering[clust.ind])
nki.pam.clust.sort <- names(nki.pam.2.daisy$clustering[nki.clust.ind])
nkiSub.clust.eset <- nkiSub.eset[,nki.clust.ind]

heatmap.2(exprs(nkiSub.clust.eset),
          dendrogram="row",
          trace = "none",
          Rowv = TRUE,
          Colv = FALSE,
          scale="none",
          col=greenred(11),
          ColSideColors = ifelse(nkiSub.clust.eset$clustering == 1, "blue","red"))

## Overall survival for the NKI patitents
nki.surv.clust.death <- coxph(Surv(nkiSub.clust.eset$"survival(death)",nkiSub.clust.eset$"event_death" ) ~
                    +factor(nkiSub.clust.eset$clustering))
summary(nki.surv.clust.death)
## MD got the same results as below
## Call:
## coxph(formula = Surv(nkiSub.clust.eset$"survival(death)", nkiSub.clust.eset$event_death) ~
##     +factor(nkiSub.clust.eset$clustering))

##   n= 295, number of events= 79

##                                         coef exp(coef) se(coef)     z Pr(>|z|)
## factor(nkiSub.clust.eset$clustering)2 1.2847    3.6136   0.2431 5.284 1.26e-07

## factor(nkiSub.clust.eset$clustering)2 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##                                       exp(coef) exp(-coef) lower .95 upper .95
## factor(nkiSub.clust.eset$clustering)2     3.614     0.2767     2.244      5.82

## Rsquare= 0.099   (max possible= 0.941 )
## Likelihood ratio test= 30.67  on 1 df,   p=3.061e-08
## Wald test            = 27.92  on 1 df,   p=1.262e-07
## Score (logrank) test = 31.86  on 1 df,   p=1.658e-08


#below is the pvalue
summary(nki.surv.clust.death)$coef[5]
plot(survfit(Surv(nkiSub.clust.eset$"survival(death)",nkiSub.clust.eset$"event_death" ) ~
                    +factor(nkiSub.clust.eset$clustering)),
     lty = 1:2,col=c("blue","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time",main="NKI survival")
legend(15, 1, legend = c(paste("cluster 1, N = ",length(nkiSub.clust.eset[,nkiSub.clust.eset$clustering==1]),sep=""), paste("cluster 2, N = ",length(nkiSub.clust.eset[,nkiSub.clust.eset$clustering==2]),sep="")), fill = c("blue", "red"),lty = c(1,2)
       , bty = "n", cex = .8)
       text(15,.2,paste("p = ", round(summary(nki.surv.clust.death)$coef[5],8),sep=''), cex=.7)
      

################################################################################
##
## Cluster the Caldas data using the RM derived gene list
##
################################################################################

## Predict on Caldas - already processed by Fathi, not sure what was done. Can we
## dig this up?

## medians.caldas <- esApply(caldasSub.eset,1,median)
## caldas.exprs.cent <-  t(scale(t(exprs(caldas.eset)),center=medians.caldas,scale=FALSE))
## exprs(caldasSub.eset) <- caldas.exprs.cent

caldasSub.eset <- caldas.eset[which(featureNames(caldas.eset) %in% common.genes.nki.caldas),]
dim(caldasSub.eset)
##MD below
#Features  Samples 
#      52      135

caldas.pam.2.daisy <- pam(daisy(t(exprs(caldasSub.eset))),2)
## caldas.pam.2.dist <- pam(as.dist(1-cor((exprs(caldasSub.eset)))),2)
## caldas.pam.eucl <- pam(t(caldasSub.exprs.scale),2)
clusplot(caldas.pam.2.daisy, labels = 2, col.p = caldas.pam.2.daisy$clustering)

all(colnames(caldasSub.eset) ==  names(caldas.pam.2.daisy$clustering))

caldas.clust.ind <- sort.int(caldas.pam.2.daisy$clustering, index.return=TRUE)$ix

pData(caldasSub.eset)$clustering <- caldas.pam.2.daisy$clustering
caldas.pam.clust.sort <- names(caldas.pam.2.daisy$clustering[caldas.clust.ind])

caldasSub.clust.eset <- caldasSub.eset[,caldas.clust.ind]

heatmap.2(exprs(caldasSub.clust.eset),
          dendrogram="row",
          trace = "none",
          Rowv = TRUE,
          Colv = FALSE,
          scale="none",
          col=greenred(11),
          ColSideColors = ifelse(caldasSub.clust.eset$clustering == 1, "blue","red"))

## Overall survival for the NKI patitents
## Restrict to samples where 'Death' == 0 or 1

caldas.test.ind <- which(caldasSub.clust.eset$Dead == 0 | caldasSub.clust.eset$Dead == 1)
caldas.surv.clust.death <- coxph(Surv(caldasSub.clust.eset[,caldas.test.ind]$Survival/12,caldasSub.clust.eset[,caldas.test.ind]$Dead ) ~
                    +factor(caldasSub.clust.eset[,caldas.test.ind]$clustering))
summary(caldas.surv.clust.death)

###NOTE: MD got the same results as below

## Call:
## coxph(formula = Surv(caldasSub.clust.eset[, caldas.test.ind]$Survival/12,
##     caldasSub.clust.eset[, caldas.test.ind]$Dead) ~ +factor(caldasSub.clust.eset[,
##     caldas.test.ind]$clustering))

##   n= 122, number of events= 35

##                                                                coef exp(coef)
## factor(caldasSub.clust.eset[, caldas.test.ind]$clustering)2 -0.9549    0.3848
##                                                             se(coef)      z
## factor(caldasSub.clust.eset[, caldas.test.ind]$clustering)2   0.3455 -2.764
##                                                             Pr(>|z|)
## factor(caldasSub.clust.eset[, caldas.test.ind]$clustering)2  0.00571 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

##                                                             exp(coef)
## factor(caldasSub.clust.eset[, caldas.test.ind]$clustering)2    0.3848
##                                                             exp(-coef)
## factor(caldasSub.clust.eset[, caldas.test.ind]$clustering)2      2.598
##                                                             lower .95 upper .95
## factor(caldasSub.clust.eset[, caldas.test.ind]$clustering)2    0.1955    0.7575

## Rsquare= 0.062   (max possible= 0.929 )
## Likelihood ratio test= 7.84  on 1 df,   p=0.005101
## Wald test            = 7.64  on 1 df,   p=0.00571
## Score (logrank) test = 8.23  on 1 df,   p=0.004114


plot(survfit(Surv(caldasSub.clust.eset[,caldas.test.ind]$Survival/12,caldasSub.clust.eset[,caldas.test.ind]$Dead) ~
                    +factor(caldasSub.clust.eset[,caldas.test.ind]$clustering)),
     lty = 1:2,col=c("blue","orange"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time")
legend(15, 1, legend = c("cluster 1", "cluster 2"), lty = c(1,2)
       , bty = "n")
