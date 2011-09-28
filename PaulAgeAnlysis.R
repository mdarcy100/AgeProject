#M D'Arcy
# final analysis for the age paper
# # Slightly modified by M D'Arcy in Summer 2011 (Creighton method added to look at correlation with age profile)
#    1) Use smaller gene set that includes fold change in addition to 10% FDR (N=164)
#			b)  to generate age profile
#				* calculate difference between median of old vs young expression in 164 genes, 
#				* set to -1 or 1 for positive or negative (up or down regulated relative)
#	
#   2) Test against Paul's data using signature associated with age
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

load("sig_fold_rm_eset.RData") ## all significant probes with > certain fold change (164 probes)
load("sig_full_rm_eset.RData") ## all significant probes, will use this one (802 probes) for crieghton analysis since too few overlap with the 
load("paul_eset.RData")
rownames(exprs(paul.eset[1:2,]))
dim(paul.eset)
#Features  Samples 
#   25927       30 
# the entrez ids associated with the probes of our data
entrez.fold<-as.vector(unlist(mget(rownames(exprs(sig.fold.rm.eset)),env=hgug4112aENTREZID)))
entrez.full<-as.vector(unlist(mget(rownames(exprs(sig.full.rm.eset)),env=hgug4112aENTREZID)))

#need to collapse probe id based on entrez id (there are duplicate gene probes)
   	
collapsed.sig.fold.rm.eset<-collapseIDs(exprs(sig.fold.rm.eset),entrez.fold,"mean")
collapsed.sig.full.rm.eset<-collapseIDs(exprs(sig.full.rm.eset),entrez.full,"mean")


###########################################################################################################
###########################################################################################################
#							DEALING WITH PAUL'S DATA.
###########################################################################################################
###########################################################################################################
#  we will want to:
#		1)  test his data against the age signature
#		2)   do a 1-D sort of genes that are in the list and do a cluster
# NOTE: we don't need to do a survival analysis on these people (because they are alive...)

##Median center Paul's affy data
medians.paul <-esApply(paul.eset,1,median)
paul.exprs.scale <- t(scale(t(exprs(paul.eset)),center=medians.paul,scale=FALSE))


#create a new dataset for paul's data
paul.scale.eset <- paul.eset
exprs(paul.scale.eset) <- paul.exprs.scale
dim(paul.scale.eset)
#Features  Samples 
#   25927       30 

#make sure centered around zero
hist(apply(exprs(paul.scale.eset),1,median))


## Remove low variability probes from paul's data
## DONT DO THIS - JUST FIND THE OVERLAP
#iqr.val.paul.scale <- esApply(paul.scale.eset,1,function(x){IQR(x,na.rm = TRUE)})

# they seem to be centered around .5
#hist(iqr.val.paul.scale,sqrt(length(iqr.val.paul.scale)))

#high.iqr.inds.paul<- which(iqr.val.paul.scale > median(iqr.val.paul.scale))

#paul.scale.eset <- paul.scale.eset[high.iqr.inds.paul,]
#dim(paul.scale.eset)
###Features  Samples 
###   12963       30 

hist(exprs(paul.scale.eset))

#first narrow down 164 genes to those that exist in Paul's data
RM_PAUL_OVERLAP.eset <- collapsed.sig.fold.rm.eset[rownames(collapsed.sig.fold.rm.eset)%in% featureData(paul.scale.eset)$EntrezID,]
 dim(RM_PAUL_OVERLAP.eset )
#138  76

hist(apply(exprs(paul.scale.eset),1,sd))
median(apply(exprs(paul.scale.eset),1,sd))
#0.2443691

######CREIGHTON METHOD STARTS HERE: 
###############################################
#find the expression of young in the overlapped genes
exprs.age.fold.young<-RM_PAUL_OVERLAP.eset[,sig.fold.rm.eset$age<30] #young defined as < 30 years of age
dim(exprs.age.fold.young) #20

#find the expression of the old
exprs.age.fold.old<-RM_PAUL_OVERLAP.eset[,sig.fold.rm.eset$age>39] #old defined as > 39 years of age
dim(exprs.age.fold.old) #23
medians.exprs.young <- apply(exprs.age.fold.young,1,median)
medians.exprs.old <- apply(exprs.age.fold.old,1,median)

young_signature<-as.matrix(medians.exprs.young)
old_signature<-as.matrix(medians.exprs.old)

#will need to set to -1 if < 0 or 1 if > 0
age_signature<-as.matrix(medians.exprs.young-medians.exprs.old)
age_signature_cast<-as.matrix(apply(age_signature,1,castDirection)) #cast direction is a function MD wrote (example: -.52 -> -1, .52 ->1)

#####CREIGHTON METHODS STOP - BELOW WE DO THE CORRELATION WITH TEST DATASET
paulSub.eset <- paul.scale.eset[featureData(paul.scale.eset)$EntrezID %in% rownames(RM_PAUL_OVERLAP.eset),]

paulSub.eset.order<-exprs(paulSub.eset[,order(paulSub.eset$Age)])
colname_age_order<-paste("Age:",paulSub.eset[,order(paulSub.eset$Age)]$Age)


colnames(paulSub.eset.order)<-colname_age_order

exprs(paulSub.eset)<-paulSub.eset.order

length(featureData(paulSub.eset)$EntrezID)
#140 - there must be duplicates, will need to collapse

length(unique(featureData(paulSub.eset)$EntrezID))
#138

# we need to collapse this data, because there are two non-unique entrez id's
collapsed.paul.eset<-collapseIDs(exprs(paulSub.eset),featureData(paulSub.eset)$EntrezID,"mean")

#reset expression data
exprs(paulSub.eset)<-collapsed.paul.eset

#sort by entrez ids since our dataset is sorted by that rowname - need to be in order to do correlation
paul.collapsed.entrezsorted<-collapsed.paul.eset[order(rownames(collapsed.paul.eset)),]
colnames(paul.collapsed.entrezsorted)<-colname_age_order
exprs(paulSub.eset)<-paul.collapsed.entrezsorted

#write out a clusterable file - this is sorted by age
write.table(paul.collapsed.entrezsorted,file=paste("PaulClusterableFile_", Sys.Date(),".txt",sep=""),sep="\t",col.names=NA)


######FIND CORRELATION BETWEEN PAUL'S DATA AND OUR SIGNATURE.
paul_cor<-cor(paul.collapsed.entrezsorted,age_signature_cast)
median(paul_cor[1:9]) #young people
median(paul_cor[24:30]) #old people


barplot(paul_cor[1:30], main ="Correlation  - sorted by age", names.arg=names.arg)
 median_high=median(paul_cor[20:30])
 median_low=median(paul_cor[1:9])
 median_middle=median(paul_cor[10:19])
 text(7,.08, "median: age < 30 = .07",col='red',cex=.9,font=2)
 text(17,-.04, "median: 30 <= age < 40 = -.04",col='blue',cex=.9, font=2)
 text(27,-.09, "median: age >= 40 = -.09",col='purple',cex=.9,font=2)
 
 ages=sort(paulSub.eset$Age)
 group=c("L","L","L","L","L","L","L","L","L","M","M","M","M","M","M","M","M","M","M","H","H","H","H","H","H","H","H","H","H","H")
 cor_data<-cbind(ages, group,paul_cor)
 reg_line<-lm(paul_cor~ages)
 anova_summary<-aov(paul_cor~factor(group))
 
 kruskal.test(paul_cor~factor(group)) #p=.285
  ltrend.test(paul_cor, factor(group), score=NULL, option = "median")
  #entrez id is not sorted anymore
#featureData(paulSub.eset)$EntrezID<-rownames(paul.collapsed.entrezsorted)

