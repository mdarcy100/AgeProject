collapseIDs<-function(x,allids=row.names(x),method="mean"){

       allids<-as.vector(allids)
       ids<- levels(as.factor(allids))
       x.col<- matrix(nrow=length(ids), ncol=dim(x)[2])

       if(length(ids)==dim(x)[1]){
                       dimnames(x)[[1]]<-allids
                       return(x)
       }

       for(i in 1:length(ids)){
               if(sum(allids==ids[i])>1){
                       indices <- allids==ids[i]
                       if(method=="mean"){
                               vals<-apply(x[indices,],2,mean)
                       }
                       if(method=="median"){
                               vals<-apply(x[indices,],2,median)
                       }
                       if(method=="stdev"){
                               temp<- x[indices,]
                               stdevs<- apply(temp,1,sd)
                               vals<- temp[match(max(stdevs),stdevs),]
                       }
                       x.col[i,] <- vals
               }else{
                       x.col[i,] <- x[allids==ids[i],]
               }
       }

       dimnames(x.col)<- list(ids,dimnames(x)[[2]])
       return(x.col)

}
castDirection<-function(x){
	if(x >= 0)
		return(1)
	else
		return(-1)	
}


strip_whitespace <- function(x)
{
    ## Strip leading and trailing whitespace.
    x <- sub("^[[:space:]]+", "", x)
    x <- sub("[[:space:]]+$", "", x)
    x
}

getReps <- function(expr.data){
    ## Find the sample namess containing rep[0-9], or 244[Kk] at the end
    rep.inds <- grep("(_rep[0-9]$|_244[kK]$)",colnames(expr.data))

    replicated.samps <- sapply(colnames(expr.data)[rep.inds],
                               function(x){
                                   unlist(strsplit(x,"_"))[1]
                               })

    unique.reps <- unique(replicated.samps)
    unique.reps <- unique.reps[unique.reps %in% colnames(expr.data)]

    reps <- sort(c(unique.reps,names(replicated.samps)))

    reps <- reps[reps %in% colnames(expr.data)]
    root.samps <- sapply(reps,function(x){unlist(strsplit(x,"_"))[1]})
    reps <- reps[-which(!(root.samps %in% root.samps[duplicated(root.samps)]))]
    list(all.reps=reps,unique.reps=unique.reps)

}

getRepsFast <- function(expr.data){
    all.samps <- colnames(expr.data)
    ## Find the sample namess containing rep[0-9], or 244[Kk] at the end
    rep.inds <- grep("(_rep[0-9]$|_244[kK]$)",all.samps)
    replicated.samps <- all.samps[rep.inds]
    rep.samp.names <- sapply(replicated.samps,
                                       function(x){
                                           unlist(strsplit(x,"_"))[1]
                                       })

    reps <- sort(c(unique(rep.samp.names),names(rep.samp.names)))
    reps <- reps[reps %in% all.samps]

    reps.ureps <- cbind(reps,
                        gsub("(_rep[0-9]$|_244[kK]|_rep[0-9]_244[kK]$)","",reps))

    reps.ureps <- reps.ureps[reps.ureps[,2] %in% names(table(reps.ureps[,2])[table(reps.ureps[,2])!=1]),]
    sampleId.reps <- reps.ureps[,2]
    list(all.reps=reps.ureps[,1], sampleId.reps=sampleId.reps)
}

lmCoeff <- function(eset, lm.vars, medianCenter=FALSE){

    if(medianCenter){
        exprs(eset) <-exprs(eset) -  apply(exprs(eset),1,
                                 function(x){median(x,na.rm=TRUE)})
    }

    cov.list <- vector("list",length = length(lm.vars))
    names(cov.list) <- lm.vars
    for(cov.name in lm.vars){
        cov.list[[cov.name]] <- pData(eset)[,cov.name]
    }
    lm.df <- data.frame(cov.list)

    gene.names <- featureNames(eset)
    genes.lm <- vector("list",length = length(gene.names))
    names(genes.lm) <- gene.names
    for(gene in gene.names){
        genes.lm[[gene]] <- summary(lm(exprs(eset)[gene,]~., data = lm.df))$coeff
    }
    genes.lm
}


lmCoeff.robust <- function(eset, lm.vars, medianCenter=FALSE){
    require(robustbase)
    if(medianCenter){
        exprs(eset) <-exprs(eset) -  apply(exprs(eset),1,
                                 function(x){median(x,na.rm=TRUE)})
    }

    cov.list <- vector("list",length = length(lm.vars))
    names(cov.list) <- lm.vars
    for(cov.name in lm.vars){
        cov.list[[cov.name]] <- pData(eset)[,cov.name]
    }
    lm.df <- data.frame(cov.list)

    gene.names <- featureNames(eset)
    genes.lm <- vector("list",length = length(gene.names))
    names(genes.lm) <- gene.names
    for(gene in gene.names){
        genes.lm[[gene]] <- summary(lmrob(exprs(eset)[gene,]~., data = lm.df))$coeff
    }
    genes.lm
}

getLmInfo <- function(lm.list,lm.vars){
    data.sum <- matrix(nrow = length(lm.list),ncol=(length(lm.vars)*2))
    coeff.names <- rownames(lm.list[[1]])
    stopifnot(all(lm.vars %in% coeff.names))

    rownames(data.sum) <- names(lm.list)
    cnames <- vector("character",length = ncol(data.sum))
    cnames[seq(1,length(cnames),2)] <- lm.vars
    cnames[seq(2,length(cnames),2)] <- paste("p_value",lm.vars,sep="_")
    colnames(data.sum) <- cnames

    for(gene in names(lm.list)){
        for(lm.var in lm.vars){
            data.sum[gene,] <- as.vector(t(lm.list[[gene]][lm.vars,c("Estimate","Pr(>|t|)")]))
        }
    }
    data.sum
}

makeAgeResidMat <- function(eset,data.sum,lm.var){
    beta.mat <- (rep(1,length(featureNames(eset))) %o% pData(eset)[,lm.var])*data.sum[,lm.var]
    resid.mat <- exprs(eset) - beta.mat
    exprs(eset) <- resid.mat
    eset
}


stripAGIjunk <- function(junk.names){
    junk.pre <- "AGI_HUM1_OLIGO_"
    no.junk.names <- unlist(sapply(junk.names,
                                   function(x){
                                       if(grepl(paste("^",junk.pre,sep=""),x)){
                                           gsub(paste("^",junk.pre,sep=""),"",x)
                                       }else{
                                           x
                                       }}))

    names(no.junk.names) <- NULL
    no.junk.names
}

heatmap.trunc <- function(x,clust.ids,trunc.val=2){
    ##sort.ind <- sort.int(clust.ids, index.return=TRUE)$ix

    x[x < -trunc.val] <- -trunc.val
    x[x >trunc.val] <- trunc.val

    ## Output in form for Cluster 3.0

    heatmap.2(x,
              dendrogram="row",
              trace = "none",
              Rowv = TRUE,
              Colv = FALSE,
              scale="none",
              col=greenred(9),
              ColSideColors = ifelse(clust.ids == 1, "blue","yellow")
              )
}


mod.findLargest <- function(gN, testStat, id.vec){
    tSsp <- split.default(testStat,id.vec)
    sapply(tSsp, function(x) names(which.max(x)))
}

pval.summary <- function(p.values, cuts = c(0,0.0001,0.001,0.01,0.025,0.05,.1,1)){
    pval.table <- cumsum(table(cut(p.values,cuts)))
    names(pval.table) <- paste('<',cuts[2:length(cuts)],sep='')
    pval.table
}
