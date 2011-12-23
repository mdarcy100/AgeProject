#M D'Arcy
# final analysis for the age paper
# # Slightly modified by M D'Arcy in Summer 2011 (Creighton method added to look at correlation with age profile)
# 	1) retrieve clusterable gene list without including fold change (~800)
#   2) use fold change restriction to get a smaller set of genes (~165) for use:
#			a)  in survival analysis
#			b)  to generate age profile
#				* calculate difference between median of old vs young expression in 165 genes, 
#				* set to -1 or 1 for positive or negative (up or down regulated relative)
#	3) Creighton method to test NKI and Caldas dataset again age-set
#			a) calculate the correlation for both datasets against the 'profile' generated above.
# 			b) see if it corresponds to good (hopefully old) vs bad (hopefully young) signatures.
#   4) Test against "Paul's" data in the same way
# will need to process his data a bit first however
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

load("sig_fold_rm_eset.RData") ## all significant probes with > certain fold change (164 probes)
load("sig_full_rm_eset.RData") ## all significant probes, will use this one (802 probes) for crieghton analysis since too few overlap with the two other datasets
load("caldas_eset.RData") ## process_CaldasSamples.R
load("nki295_eset.RData") ## process_NKI295samples.R
load("paul_eset.RData")

dim(sig.fold.rm.eset)
dim(sig.full.rm.eset)
dim(caldas.eset)
dim(nki295.eset)
dim(paul.eset)

#####################################################################################
#Prior to any analysis, calculate how many genes generally went up and down with age
exprs.age.full.young<-sig.full.rm.eset[,sig.full.rm.eset$age<30]
dim(exprs.age.full.young) #20

exprs.age.full.old<-sig.full.rm.eset[,sig.full.rm.eset$age>39]
dim(exprs.age.full.old) #23

medians.exprs.full.young <- apply(exprs.age.full.young,1,median)
medians.exprs.full.old <- apply(exprs.age.full.old,1,median)

#will need to set to -1 if < 0 or 1 if > 0 #if young - old < 0 (-1), then expression goes up with age, otherwise down
age_signature_full<-as.matrix(medians.exprs.full.young-medians.exprs.full.old) 
age_signature_cast<-as.matrix(apply(age_signature_full,1,castDirection))

##### the file of up/down probes to write to file
entrez.full<-as.vector(unlist(mget(rownames(exprs(sig.full.rm.eset)),env=hgug4112aENTREZID)))
up_down<-cbind(entrez.full,age_signature_cast,medians.exprs.full.young,medians.exprs.full.old)
colnames(up_down)<-c("entrezid","direction","median_young","median_old")
write.table(up_down,file=paste("UpDown_", Sys.Date(),".txt",sep=""),sep="\t",col.names=NA)
#write out this matrix of 1/-1 values to get relative increases/decreases with age
#age_names<-rownames(age_signature_full) #vector of probe'ids
#####################################################################################


#rownames(exprs(sig.fold.rm.eset)) #this should be the probe id
#rownames(exprs(sig.full.rm.eset)) #this should be the probe id, will retreive the entrez id from this value
rownames(exprs(caldas.eset[1:2,])) #these are symbols (will need to get entrez id id)
rownames(exprs(nki295.eset[1:2,])) #ACCESSION Number (need entrez id)
rownames(exprs(paul.eset[1:2,]))
#RefSeq accession number of sequence used to design oligonucleotide probe
#################################################################################
###		MD - Aug 2011 - for Creighton method and making an age profile need those 
### 	
###		< 30 and > 40
### from data generated after main analysis (see age_decade_analysis_RM.R)

#featureData(sig.fold.rm.eset)

###first collapse the entrez ids to single genes
	
##### First collapseID based on entrez ids...
#     there are 145 of these (from 164) from the fold change ones
#     there are 719/802 unique probes
   
#xd<-apply(xd,2,as.numeric)
   
# the entrez ids associated with the probes of our data
entrez.fold<-as.vector(unlist(mget(rownames(exprs(sig.fold.rm.eset)),env=hgug4112aENTREZID)))
entrez.full<-as.vector(unlist(mget(rownames(exprs(sig.full.rm.eset)),env=hgug4112aENTREZID)))

   	
collapsed.sig.fold.rm.eset<-collapseIDs(exprs(sig.fold.rm.eset),entrez.fold,"mean")
collapsed.sig.full.rm.eset<-collapseIDs(exprs(sig.full.rm.eset),entrez.full,"mean")

#prior to calculating the age signature, determine which NKI genes are in these collapsed lists, then do survival analysis on these genes

## Prepare NKI295 data (stolen almost exactly from age_decade_analysis_RM.R)
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
#list.entrez <- na.omit(unlist(mget(rownames(exprs(sig.fold.rm.eset)),env=hgug4112aENTREZID)))
#list.symbol <- na.omit(unlist(mget(rownames(exprs(sig.fold.rm.eset)),env=hgug4112aSYMBOL)))

nki.match.vec <- vector("numeric",length=dim(nki295.scale.eset)["Features"])
names(nki.match.vec) <- featureNames(nki295.scale.eset)

for(e.id in unique(entrez.fold)){ #should be 719 of these
    e.id.pattern <- paste("(^",e.id,"$|\\|",e.id,"$|^",e.id,"\\|$|\\|",e.id,"\\|)",sep="")
    nki.match.ind <- grepl(e.id.pattern,featureData(nki295.scale.eset)$geneid)
    if(any(nki.match.ind)){
        nki.match.vec[nki.match.ind] <- nki.match.vec[nki.match.ind]+1
    }
}

nki.gene.match.full <- names(nki.match.vec[nki.match.vec !=0])
nkiSub.eset <- nki295.scale.eset[nki.gene.match.full,]

#### MD testing something 
dim(nkiSub.eset)
###MD results:
#Features  Samples 
#      380      295

#################################################################################
#
# FINDING AGE SIGNATURE FOR COMMON GENES IN NKI DATASET  
# to get age signature for NKI data (use data where there is overlap)
#
#################################################################################
RM_NKI_OVERLAP.eset <- collapsed.sig.fold.rm.eset[rownames(collapsed.sig.fold.rm.eset)%in% featureData(nkiSub.eset)$geneid,]
 dim(RM_NKI_OVERLAP.eset )
 #this value is only true for the full significant set - far fewer in the fold list
#[1] 380  76
	
#NOTE: will need to sort the NKI data prior to correlation analysis

exprs.age.fold.young<-RM_NKI_OVERLAP.eset[,sig.fold.rm.eset$age<30]
dim(exprs.age.fold.young) #20

exprs.age.fold.old<-RM_NKI_OVERLAP.eset[,sig.fold.rm.eset$age>39]
dim(exprs.age.fold.old) #23

medians.exprs.young <- apply(exprs.age.fold.young,1,median)
medians.exprs.old <- apply(exprs.age.fold.old,1,median)


young_signature<-as.matrix(medians.exprs.young)
old_signature<-as.matrix(medians.exprs.old)

#will need to set to -1 if < 0 or 1 if > 0
age_signature<-as.matrix(medians.exprs.young-medians.exprs.old)
age_signature_cast<-as.matrix(apply(age_signature,1,castDirection))

age_names<-rownames(age_signature) #vector of probe'ids

#################################################################################
#
### NKI SURVIVAL ANALYSIS to get good/bad prognosis people, then do correlation
#
#################################################################################

# SORT PRIOR TO CLUSTERING/SURVIVAL ANALYSIS !!!!sort data by entrez id
nkiSub.eset.age.clust<-exprs(nkiSub.eset) #should be 380 elements in Aug 2011

#set the rownames of expression data to entrez id so the same as our dataset
rownames(nkiSub.eset.age.clust) <- featureData(nkiSub.eset)$geneid

#sort by entrez ids since our dataset is sorted by that rowname - need to be in order to do correlation
nkiSub.eset.age.clust<-nkiSub.eset.age.clust[order(rownames(nkiSub.eset.age.clust)),]

#reset column name of expression array
exprs(nkiSub.eset)<-nkiSub.eset.age.clust

nki.pam.2.daisy <- pam(daisy(t(nkiSub.eset.age.clust)),2)

clusplot(nki.pam.2.daisy, labels = 2, col.p = nki.pam.2.daisy$clustering)

all(colnames(nkiSub.eset) ==  names(nki.pam.2.daisy$clustering))

colnames(exprs(nkiSub.eset))

nki.clust.ind <- sort.int(nki.pam.2.daisy$clustering, index.return=TRUE)$ix

pData(nkiSub.eset)$clustering <- nki.pam.2.daisy$clustering


nki.pam.clust.sort <- names(nki.pam.2.daisy$clustering[nki.clust.ind])
nkiSub.clust.eset <- nkiSub.eset[,nki.clust.ind]


## Overall survival for the NKI patitents (380 people who overlap genes with full dataset)
nki.surv.clust.death <- coxph(Surv(nkiSub.clust.eset$"survival(death)",nkiSub.clust.eset$"event_death" ) ~
                    +factor(nkiSub.clust.eset$clustering))
summary(nki.surv.clust.death)

plot(survfit(Surv(nkiSub.clust.eset$"survival(death)",nkiSub.clust.eset$"event_death" ) ~
                    +factor(nkiSub.clust.eset$clustering)),
     lty = 1:2,col=c("blue","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time")


###   SUBSET data for good and bad survival
nkiSub.clust.eset.bad<-nkiSub.clust.eset[,nkiSub.clust.eset$clustering==2]
dim(nkiSub.clust.eset.bad)
#Features  Samples 
#     380      182 	

nkiSub.clust.eset.good<-nkiSub.clust.eset[,nkiSub.clust.eset$clustering==1]
dim(nkiSub.clust.eset.good)
#Features  Samples 
#     380      113 	


###########################################################################################
#     CREIGHTON METHOD for NKI data
####  MD: to get correlation matrix for nkiSub.est from AgeSignature
#### need to get those entrez ids that appear in age signature
###########################################################################################

#look at characteristics of good prognosis correlated with the signature

	nki_castgood_cor<-cor(exprs(nkiSub.clust.eset.good),age_signature_cast)
	nki_castgood_cor<-sort(nki_castgood_cor)
	hist(nki_castgood_cor)
	length(nki_castgood_cor)
	median(nki_castgood_cor)
	mean(nki_castgood_cor)
	quantile(nki_castgood_cor)

	
	nki_good_cor<-cor(exprs(nkiSub.clust.eset.good),age_signature)
	nki_good_cor<-sort(nki_good_cor)
	hist(nki_good_cor)
	length(nki_good_cor)
	median(nki_good_cor)
	mean(nki_good_cor)
	quantile(nki_good_cor)
	

	nki_goodyoung_cor<-cor(exprs(nkiSub.clust.eset.good),young_signature)
	nki_goodyoung_cor<-sort(nki_goodyoung_cor)
	hist(nki_goodyoung_cor)
	length(nki_goodyoung_cor)
	median(nki_goodyoung_cor)
	mean(nki_goodyoung_cor)
	quantile(nki_goodyoung_cor)


	nki_goodold_cor<-cor(exprs(nkiSub.clust.eset.good),old_signature)
	nki_goodold_cor<-sort(nki_goodold_cor)
	hist(nki_goodold_cor)
	length(nki_goodold_cor)
	median(nki_goodold_cor)
	mean(nki_goodold_cor)
	quantile(nki_goodold_cor)



#look at characteristics of bad prognosis correlated with the signature

	nki_castbad_cor<-cor(exprs(nkiSub.clust.eset.bad),age_signature_cast)
	nki_castbad_cor<-sort(nki_castbad_cor)
	hist(nki_castbad_cor)
	length(nki_castbad_cor)
	median(nki_castbad_cor)
	mean(nki_castbad_cor)
	quantile(nki_castbad_cor)

	nki_bad_cor<-cor(exprs(nkiSub.clust.eset.bad),age_signature)
	 nki_bad_cor<-sort(nki_bad_cor)
	length(nki_bad_cor)
	median(nki_bad_cor)
	mean(nki_bad_cor)
	quantile(nki_bad_cor)

	nki_badyoung_cor<-cor(exprs(nkiSub.clust.eset.bad),young_signature)
	nki_badyoung_cor<-sort(nki_badyoung_cor)
	length(nki_badyoung_cor)
	median(nki_badyoung_cor)
	mean(nki_badyoung_cor)
	quantile(nki_badyoung_cor)
	
	nki_badold_cor<-cor(exprs(nkiSub.clust.eset.bad),old_signature)
	nki_badold_cor<-sort(nki_badold_cor)
	length(nki_badold_cor)
	median(nki_badold_cor)
	mean(nki_badold_cor)
	quantile(nki_badold_cor)
	
	
#image of correlations
#Image 1	
	plot(nki_bad_cor,pch='.',col='red', main="NKI: correlation with age signatures", ylab="correlation",xlab="sample")
	points(nki_castbad_cor,pch=20,col='orange')
	points(nki_badyoung_cor, pch='*', col='red')
	points(nki_badold_cor, pch='+', col='red')
	
	points(nki_castgood_cor,pch=20,col='purple')
	points(nki_good_cor, pch='.', col='blue')
	points(nki_goodyoung_cor, pch='*', col='blue')
	points(nki_goodold_cor, pch='+', col='blue')
	
	text(50,.35,"orange=creighton correlation with bad dx", cex=.9)
	text(50,.30,"purple=creighton correlation with good dx", cex=.9)
	text(25,.2, "blue=good", cex=.9)
	text(25,.15, "red=bad", cex=.9)
	text(25,.25,"*=young, +=old")


	#image 2 - only include the orange and purple lines
	#write these images out to pdf/jpeg files
	pdf('Figures/NKICreighton.pdf')
	plot(nki_bad_cor,pch=20,col='orange', main="NKI: correlation with age signatures", ylab="correlation",xlab="samples")
	points(nki_castgood_cor,pch=20,col='purple')
	legend(5,.3, c(paste("Good Prognosis,","N = ",length(nki_castgood_cor),sep=""), paste("Poor Prognosis,","N = ",length(nki_castbad_cor),sep="")),fill = c("purple", "orange"),title=paste(" "))
	text(120,-.25,paste("p = ",round(t.test(nki_castgood_cor,nki_castbad_cor)$p.value,6),sep=""))
	dev.off()
	
	jpeg('Figures/NKICreighton.jpeg')
	plot(nki_bad_cor,pch=20,col='orange', main="NKI: correlation with age signatures", ylab="correlation",xlab="samples")
	points(nki_castgood_cor,pch=20,col='purple')
	legend(5,.3, c(paste("Good Prognosis,","N = ",length(nki_castgood_cor),sep=""), paste("Poor Prognosis,","N = ",length(nki_castbad_cor),sep="")),fill = c("purple", "orange"),title=paste(" "))
	text(120,-.25,paste("p = ",round(t.test(nki_castgood_cor,nki_castbad_cor)$p.value,6),sep=""))
	dev.off()

	
		#image 3
	pdf('Figures/NKICreightonBP.pdf')
	bp_names<-c("Good Prognosis","Poor Prognosis")
	boxplot(nki_castgood_cor, nki_castbad_cor,names=bp_names,fill=c("purple","orange"))
	title("NKI: correlation with age signature")
	#legend(1,2, c("Good Prognosis", "Bad Prognosis"),fill = c("purple", "orange"))
	text(1,median(nki_castgood_cor)+.05,paste("median=",round(median(nki_castgood_cor),2)),cex=.7)
	text(2,median(nki_castbad_cor)+.05,paste("median=",round(median(nki_castbad_cor),2)),cex=.7)
	text(1.5,-.2,paste("p = ",round(t.test(nki_castgood_cor,nki_castbad_cor)$p.value,9),sep=""))
	dev.off()
	
	jpeg('Figures/NKICreightonBP.jpeg')
	bp_names<-c("Good Prognosis","Poor Prognosis")
	boxplot(nki_castgood_cor, nki_castbad_cor,names=bp_names,fill=c("purple","orange"))
	title("NKI: correlation with age signature")
	#legend(1,2, c("Good Prognosis", "Bad Prognosis"),fill = c("purple", "orange"))
	text(1,median(nki_castgood_cor)+.05,paste("median=",round(median(nki_castgood_cor),2)),cex=.7)
	text(2,median(nki_castbad_cor)+.05,paste("median=",round(median(nki_castbad_cor),2)),cex=.7)
	text(1.5,-.2,paste("p = ",round(t.test(nki_castgood_cor,nki_castbad_cor)$p.value,9),sep=""))
	dev.off()

	
	
###########################################################################################
###########################################################################################
####
####   CALDAS DATASET ANALYSIS
####
###########################################################################################
###########################################################################################


#
# FINDING AGE SIGNATURE FOR COMMON GENES IN CALDAS DATASET  
# to get age signature for caldas data (use data where there is overlap)
# No entrez ids for this dataset - need to use the the symbol and then collapse on that
#  We will be looking for the signature in those genes with a large fold change, not just all genes with a good FDR
#################################################################################
#featureNames(caldas.eset)

#get the symbols associated with the probe names
symbol.fold <- as.vector(unlist(mget(rownames(exprs(sig.fold.rm.eset)),env=hgug4112aSYMBOL)))
symbol.full<-as.vector(unlist(mget(rownames(exprs(sig.full.rm.eset)),env=hgug4112aSYMBOL)))

#collapse on the symbol so that they are unique	(145/164 and 719/802)
collapsed.symbol.fold.rm.eset<-collapseIDs(exprs(sig.fold.rm.eset),symbol.fold,"mean")
collapsed.symbol.full.rm.eset<-collapseIDs(exprs(sig.full.rm.eset),symbol.full,"mean")

##find overlap with the CALDAS dataset
intersect(rownames(collapsed.symbol.fold.rm.eset),featureNames(caldas.eset))
intersect(featureNames(caldas.eset),rownames(collapsed.symbol.fold.rm.eset))
## 

RM_CALDAS_OVERLAP.eset <- collapsed.symbol.fold.rm.eset[rownames(collapsed.symbol.fold.rm.eset)%in% featureNames(caldas.eset),]
dim(RM_CALDAS_OVERLAP.eset )
#[1] 111  76
	
#NOTE: will not need to sort prior to correlation analysis (the two datasets appear to be sorted the same way already (by symbol name))

exprs.age.fold.young<-RM_CALDAS_OVERLAP.eset[,sig.fold.rm.eset$age<30]
dim(exprs.age.fold.young) #20

exprs.age.fold.old<-RM_CALDAS_OVERLAP.eset[,sig.fold.rm.eset$age>39]
dim(exprs.age.fold.old) #23

medians.exprs.young <- apply(exprs.age.fold.young,1,median)
medians.exprs.old <- apply(exprs.age.fold.old,1,median)

young_signature<-as.matrix(medians.exprs.young)
old_signature<-as.matrix(medians.exprs.old)

#will need to set to -1 if < 0 or 1 if > 0
age_signature<-as.matrix(medians.exprs.young-medians.exprs.old)
age_signature_cast<-as.matrix(apply(age_signature,1,castDirection))
age_names<-rownames(age_signature) #vector of probe'ids


####CLUSTER THE GENES (first subset the original caldas dataset to those genes that overlap with original list)
caldasSub.eset <- caldas.eset[featureNames(caldas.eset) %in% rownames(RM_CALDAS_OVERLAP.eset),]
dim(caldasSub.eset)
#Features  Samples 
#     111      135 

caldasSub.eset.clust<-exprs(caldasSub.eset)

caldas.pam.2.daisy <- pam(daisy(t(caldasSub.eset.clust)),2)

clusplot(caldas.pam.2.daisy, labels = 2, col.p = caldas.pam.2.daisy$clustering)

all(colnames(caldasSub.eset) ==  names(caldas.pam.2.daisy$clustering))

caldas.clust.ind <- sort.int(caldas.pam.2.daisy$clustering, index.return=TRUE)$ix
colnames(exprs(caldasSub.eset))

pData(caldasSub.eset)$clustering <- caldas.pam.2.daisy$clustering
caldas.pam.clust.sort <- names(caldas.pam.2.daisy$clustering[caldas.clust.ind])

caldasSub.clust.eset <- caldasSub.eset[,caldas.clust.ind]

caldas.test.ind <- which(caldasSub.clust.eset$Dead == 0 | caldasSub.clust.eset$Dead == 1)

####SURVIVAL for 135 people, not as significant as other group
caldas.surv.clust.death <- coxph(Surv(caldasSub.clust.eset[,caldas.test.ind]$Survival/12,caldasSub.clust.eset[,caldas.test.ind]$Dead ) ~
                    +factor(caldasSub.clust.eset[,caldas.test.ind]$clustering))
summary(caldas.surv.clust.death)


plot(survfit(Surv(caldasSub.clust.eset[,caldas.test.ind]$Survival/12,caldasSub.clust.eset[,caldas.test.ind]$Dead) ~
                    +factor(caldasSub.clust.eset[,caldas.test.ind]$clustering)),
     lty = 1:2,col=c("blue","red"), mark.time = TRUE, ylab = "Probability",
     xlab = "Survival Time")
legend(15, 1, legend = c("cluster 1", "cluster 2"), lty = c(1,2)
       , bty = "n")


###   SUBSET data for good and bad survival
caldasSub.clust.eset.bad<-caldasSub.clust.eset[,caldasSub.clust.eset$clustering==2]
dim(caldasSub.clust.eset.bad)
#Features  Samples 
#     111      66 	

caldasSub.clust.eset.good<-caldasSub.clust.eset[,caldasSub.clust.eset$clustering==1]
dim(caldasSub.clust.eset.good)
#Features  Samples 
#     111     69	


###########################################################################################
#     CREIGHTON METHOD for CALDAS data
####  MD: to get correlation matrix for caldasSub.eset from AgeSignature
#### these are matched and sorted by symbol
###########################################################################################

#look at characteristics of good prognosis correlated with the signature

	caldas_castgood_cor<-cor(exprs(caldasSub.clust.eset.good),age_signature_cast)
	caldas_castgood_cor<-sort(caldas_castgood_cor)
	hist(caldas_castgood_cor)
	length(caldas_castgood_cor)
	median(caldas_castgood_cor)
	mean(caldas_castgood_cor)
	quantile(caldas_castgood_cor)

	caldas_good_cor<-cor(exprs(caldasSub.clust.eset.good),age_signature)
	caldas_good_cor<-sort(caldas_good_cor)
	hist(caldas_good_cor)
	length(caldas_good_cor)
	median(caldas_good_cor)
	mean(caldas_good_cor)
	quantile(caldas_good_cor)


	caldas_goodyoung_cor<-cor(exprs(caldasSub.clust.eset.good),young_signature)
	caldas_goodyoung_cor<-sort(caldas_goodyoung_cor)
	hist(caldas_goodyoung_cor)
	length(caldas_goodyoung_cor)
	median(caldas_goodyoung_cor)
	mean(caldas_goodyoung_cor)
	quantile(caldas_goodyoung_cor)


	caldas_goodold_cor<-cor(exprs(caldasSub.clust.eset.good),old_signature)
	caldas_goodold_cor<-sort(caldas_goodold_cor)
	hist(caldas_goodold_cor)
	length(caldas_goodold_cor)
	median(caldas_goodold_cor)
	mean(caldas_goodold_cor)
	quantile(caldas_goodold_cor)


#look at characteristics of bad prognosis correlated with the signature

	caldas_castbad_cor<-cor(exprs(caldasSub.clust.eset.bad),age_signature_cast)
	caldas_castbad_cor<-sort(caldas_castbad_cor)
	hist(caldas_castbad_cor)
	length(caldas_castbad_cor)
	median(caldas_castbad_cor)
	mean(caldas_castbad_cor)
	quantile(caldas_castbad_cor)

	caldas_bad_cor<-cor(exprs(caldasSub.clust.eset.bad),age_signature)
	 caldas_bad_cor<-sort(caldas_bad_cor)
	length(caldas_bad_cor)
	median(caldas_bad_cor)
	mean(caldas_bad_cor)
	quantile(caldas_bad_cor)
	
	caldas_badyoung_cor<-cor(exprs(caldasSub.clust.eset.bad),young_signature)
	caldas_badyoung_cor<-sort(caldas_badyoung_cor)
	length(caldas_badyoung_cor)
	median(caldas_badyoung_cor)
	mean(caldas_badyoung_cor)
	quantile(caldas_badyoung_cor)

	caldas_badold_cor<-cor(exprs(caldasSub.clust.eset.bad),old_signature)
	caldas_badold_cor<-sort(caldas_badold_cor)
	length(caldas_badold_cor)
	median(caldas_badold_cor)
	mean(caldas_badold_cor)
	quantile(caldas_badold_cor)
	
	
#image of correlations	
	plot(caldas_bad_cor,pch='.',col='red', main="CALDAS: correlation with age signatures", ylab="correlation",xlab="sample")
	points(caldas_castbad_cor,pch=20,col='orange')
	points(caldas_badyoung_cor, pch='*', col='red')
	points(caldas_badold_cor, pch='+', col='red')
	
	points(caldas_castgood_cor,pch=20,col='purple')
	points(caldas_good_cor, pch='.', col='blue')
	points(caldas_goodyoung_cor, pch='*', col='blue')
	points(caldas_goodold_cor, pch='+', col='blue')
	
	
	text(50,.25,"orange=creighton correlation with bad dx", cex=.8)
	text(50,.30,"purple=creighton correlation with good dx", cex=.8)
	text(25,.20, "blue=good", cex=.7)
	text(25,.15, "red=bad", cex=.7)
	text(25,.25,"*=young, +=old")
	
	#image 2 - only include the orange and purple lines
	# write out to pdf and jpeg files
	pdf('Figures/CaldasCreighton.pdf')
	plot(caldas_bad_cor,pch=20,col='orange', main="CALDAS: correlation with age signatures", ylab="correlation",xlab="samples")
	points(caldas_castgood_cor,pch=20,col='purple')
	legend(5,.25, c(paste("Good Prognosis,","N=",length(caldas_castgood_cor),sep=""), paste("Poor Prognosis,","N=",length(caldas_castbad_cor),sep="")),fill = c("purple", "orange"),title=paste(" "))
	text(50,-.2,paste("p = ",round(t.test(caldas_castgood_cor,caldas_castbad_cor)$p.value,9),sep=""))
	dev.off()
	
	jpeg('Figures/CaldasCreighton.jpeg')
	plot(caldas_bad_cor,pch=20,col='orange', main="CALDAS: correlation with age signatures", ylab="correlation",xlab="samples")
	points(caldas_castgood_cor,pch=20,col='purple')
	legend(5,.25, c(paste("Good Prognosis,","N=",length(caldas_castgood_cor),sep=""), paste("Poor Prognosis,","N=",length(caldas_castbad_cor),sep="")),fill = c("purple", "orange"),title=paste(" "))
	text(50,-.2,paste("p = ",round(t.test(caldas_castgood_cor,caldas_castbad_cor)$p.value,9),sep=""))
	dev.off()

	
#BOXPLOTS
	#image 3
	pdf('Figures/CaldasCreightonBP.pdf')
	bp_names<-c("Good Prognosis","Poor Prognosis")
	boxplot(caldas_castgood_cor, caldas_castbad_cor,names=bp_names,fill=c("purple","orange"))
	title("CALDAS: correlation with age signature")
	#legend(1,2, c("Good Prognosis", "Bad Prognosis"),fill = c("purple", "orange"))
	text(1,median(caldas_castgood_cor)+.02,paste("median=",round(median(caldas_castgood_cor),2)),cex=.7)
	text(2,median(caldas_castbad_cor)+.02,paste("median=",round(median(caldas_castbad_cor),2)),cex=.7)
	text(1.5,-.2,paste("p = ",round(t.test(caldas_castgood_cor,caldas_castbad_cor)$p.value,9),sep=""))
	dev.off()
	
	jpeg('Figures/CaldasCreightonBP.jpeg')
	bp_names<-c("Good Prognosis","Poor Prognosis")
	boxplot(caldas_castgood_cor, caldas_castbad_cor,names=bp_names,fill=c("purple","orange"))
	title("CALDAS: correlation with age signature")
	#legend(1,2, c("Good Prognosis", "Bad Prognosis"),fill = c("purple", "orange"))
	text(1,median(caldas_castgood_cor)+.02,paste("median=",round(median(caldas_castgood_cor),2)),cex=.7)
	text(2,median(caldas_castbad_cor)+.02,paste("median=",round(median(caldas_castbad_cor),2)),cex=.7)
	text(1.5,-.2,paste("p = ",round(t.test(caldas_castgood_cor,caldas_castbad_cor)$p.value,9),sep=""))
	dev.off()


#################################################################################
#
# FINDING AGE SIGNATURE FOR COMMON GENES IN PAUL's DATASET  
# to get age signature for PAUL data (use data where there is overlap)
#
#################################################################################
RM_PAUL_OVERLAP.eset <- collapsed.sig.fold.rm.eset[rownames(collapsed.sig.fold.rm.eset)%in% featureData(paul.scale.eset)$EntrezID,]
 dim(RM_PAUL_OVERLAP.eset )
#[1] 93  76
	
#NOTE: will need to sort the NKI data prior to correlation analysis

exprs.age.fold.young<-RM_PAUL_OVERLAP.eset[,sig.fold.rm.eset$age<30]
dim(exprs.age.fold.young) #20

exprs.age.fold.old<-RM_PAUL_OVERLAP.eset[,sig.fold.rm.eset$age>39]
dim(exprs.age.fold.old) #23

medians.exprs.young <- apply(exprs.age.fold.young,1,median)
medians.exprs.old <- apply(exprs.age.fold.old,1,median)


young_signature<-as.matrix(medians.exprs.young)

old_signature<-as.matrix(medians.exprs.old)

#will need to set to -1 if < 0 or 1 if > 0
age_signature<-as.matrix(medians.exprs.young-medians.exprs.old)
age_signature_cast<-as.matrix(apply(age_signature,1,castDirection))


#need to subset the paul data so that it only includes the genes that are common in both
paulSub.eset <- paul.eset[featureData(paul.eset)$EntrezID %in% rownames(RM_PAUL_OVERLAP.eset),]
paulSub.eset.order<-exprs(paulSub.eset[,order(paulSub.eset$Age)])


length(featureData(paulSub.eset)$EntrezID)


length(unique(featureData(paulSub.eset)$EntrezID))
#93

# we need to collapse this data, because there are two non-unique entrez id's
collapsed.paul.eset<-collapseIDs(exprs(paulSub.eset),featureData(paulSub.eset)$EntrezID,"mean")


#sort by entrez ids since our dataset is sorted by that rowname - need to be in order to do correlation
paul.collapsed.entrezsorted<-collapsed.paul.eset[order(rownames(collapsed.paul.eset)),]

#this data needs to be sorted prior to normalization

#normalize affymatrix data with RMA


paul_cor<-cor(paul.collapsed.entrezsorted,age_signature_cast)

#now sort to age prior to correlation


dim(paulSub.eset)

	

###########################################################################################################
###########################################################################################################
###########################################################################################################
###  Look at demographic data of the different samples, need to set those with 999 bmi to NA

# there were 30 people missing
#what is the max age:
max(sig.fold.rm.eset$age)
#[1] 55
#what is the min age:
min(sig.fold.rm.eset$age)
#[1] 20

sig.fold.rm.eset$bmi<-ifelse(sig.fold.rm.eset$bmi > 900, NA, sig.fold.rm.eset$bmi)
mean(sig.fold.rm.eset$bmi,na.rm = TRUE) #29.99565
bmisummary<-aov(sig.fold.rm.eset$bmi~factor(floor(sig.fold.rm.eset$age/10)))
summary(bmisummary)
print(model.tables(bmisummary,"means"),digits=4) 

#factor(floor(sig.fold.rm.eset$age/10)) 
#        2     3     4    5
#    30.54 29.36 30.08 30.5
#rep 15.00 16.00 13.00  2.0

#1 = pre-menopause, 2=peri-menopause
table(sig.fold.rm.eset$menopause)
#1  2 
#62 14 

table(factor(floor(sig.fold.rm.eset$age/10)))
# 2  3  4  5 
#20 33 17  6
 
table(sig.fold.rm.eset$race)
# 1  2  3  4 
#48  6 19  3 

###########################################################################################################
###########################################################################################################
#For two groups survival (NKI/CALDAS) and tells the pathological characteristics (ER status, size,grade, node status) and 
#subtype (n, percentage, or average value forquantitative variables).  Also give the age of the patients in the two groups?
#  Don't know anything about the phenotypic data available so need to find out
# colnames(pData(caldas.eset))
#AGE, Meno, ER of list 135 details
#Nodes, "Meno", Grade, Size, Call = subtype
quantile(caldas.eset$AGE)
#     0%      25%      50%      75%     100% 
#32.00000 50.08966 58.00000 63.73854 70.00000

#ER status
table(caldas.eset$`ER of list 135 details`)
# 0  1 
#40 93 

#Size
table(factor(caldas.eset$Size))
mean(mean(caldas.eset$Size, na.rm=TRUE)) # there is missing data
#[1] 1.864179

#Grade
table(factor(caldas.eset$Grade))
# 1  2  3 
# 35 50 49 

#Nodes
table(factor(caldas.eset$`Total nodes`))
caldas.eset$`Total nodes`  #missing and . are the same?

#subtype
table(caldas.eset$Call)
#Basal   Her2   LumA   LumB Normal 
#    17     21     59     25     13 

mean(caldas.eset$AGE)
#56.85066
table(factor(floor(caldas.eset$AGE/10)))
 # 3  4  5  6  7 
 # 4 28 43 57  3
 hist(caldas.eset$AGE)
 
table(caldas.eset$Meno)
# 1  2 
# 44 90 
