####MD Sept 2011 - this analysis may have been done before, but i don't know where it is
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

## DO analysis on everyone (all menopausal statuses (1-4))

#1/2 are pre-peri menopausal status, 3 is post menopausal status. 4= unknown status
menopause.eset <- all.rm.eset[,all.rm.eset$age >= 20 &
                                all.rm.eset$menopause %in% c(1,2,3)]
                                

dim(menopause.eset)
## Features  Samples
##    31456       99

## Basic prefiltering
## Reduce probes to those having an ENTREZ Id
entrezIds <- mget(featureNames(menopause.eset), envir=hgug4112aENTREZID)
haveEntrezIds <- names(entrezIds)[sapply(entrezIds,function(x) !is.na(x))]
numNoEntrezId <- length(featureNames(menopause.eset)) - length(haveEntrezIds)

menopause.eset <- menopause.eset[haveEntrezIds,]

#create new variable; group menopausal statuses = 1/2 together, 3 by itself
meno.status<-ifelse(menopause.eset$menopause < 3, 1, 2)

pData(menopause.eset)$meno.status<-meno.status
#make sure we have the correct number of people in each group
table(menopause.eset$meno.status)
# 1  2 
#76 23 

## Determine the interquartile range for each gene - will be used to filter results
## after applying the empirical Bayes step in Limma
## With limma

iqr.menopause.eset <- esApply(menopause.eset,1,IQR)



## Regression analysis using Limma

## Regression on age.decade
design.menopause <- model.matrix(~ menopause.eset$meno.status)
colnames(design.menopause) <- c("Intercept","menopause.limma")

menopause.limma.fit1 <- lmFit(menopause.eset,design.menopause)
menopause.limma.ebayes <- eBayes(menopause.limma.fit1)

## Apply IQR cutoff of median(iqr.menopause.eset)
menopause.limma.filtered <- menopause.limma.ebayes[iqr.menopause.eset>median(iqr.menopause.eset),]
##MD checking dimension
dim(menopause.limma.filtered)
#[1] 11969     2

#this is a summary of the results of the regression 
menopause.limma.summary <- topTable(menopause.limma.filtered,coef="menopause.limma",adjust.method="fdr",num=Inf)


pval.summary(menopause.limma.summary$P.Value)

## not a lot that is significane
#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1 
#    0      2     38    114    273    734  11969  


pval.summary(menopause.limma.summary$adj.P.Val)

#<1e-04 <0.001  <0.01 <0.025  <0.05   <0.1     <1 
#    0      0      0      0      0      0  11969  

## Compute qvalue if desired, however the correctness of using this method to
## Validity of computing q-value based on filtered p-values? - nothing is significant and we can stop
q.value.menopause <- qvalue(p=menopause.limma.summary$P.Value)
summary(q.value.menopause)

#        <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1    <1
#p-value      0      2    38    114   273  734 11969
#q-value      0      0     0      0     0    0 11969


qsummary(q.value.menopause, cuts = c(1e-04, 0.001, 0.01, 0.025, 0.05, 0.1,0.15, 0.2, 1))

###        <1e-04 <0.001 <0.01 <0.025 <0.05 <0.1 <0.15 <0.2    <1
#p-value      0      2    38    114   273  734  1269 1827 11969
#q-value      0      0     0      0     0    0     0    0 11969