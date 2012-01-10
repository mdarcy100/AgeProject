require(lattice)
require(pcurve)
require(rgl)

PCAsimp <- function(exp.data,cent=FALSE,scle=FALSE,
                    col.vec=NULL,pch.vec=NULL,rot.view=NULL,id.vec=NULL,key.vec=NULL){

  pc <- pca(t(exp.data),cent=cent,scle=scle)
  z <- pc$pcs[,3]
  x <- pc$pcs[,1]
  y <- pc$pcs[,2]
  if(is.character(key.vec)){
    key <- simpleKey(key.vec,points=TRUE,pch=19)
  }else{
    key <- NULL
  }
  
  cloud(z~x+y,
        col = col.vec,
        pch=pch.vec,
        xlab="PC1",ylab="PC2",zlab="PC3",
        scales = list(arrows=FALSE),
        screen = rot.view,
        border=TRUE,
        perspective=TRUE,
        distance=.4,
        zoom=.6,
        cex=2,
        colorkey=FALSE,
        key=key
        )
      }

PCArotate <- function(exp.data,cent=FALSE,scle=FALSE,
                      col.vec=NULL,text.vec=NULL,size=12){
  pc <- pca(t(exp.data),cent=FALSE,scle=FALSE)
  open3d()
  x <- pc$pcs[,1]
  y <- pc$pcs[,2]
  z <- pc$pcs[,3]
  plot3d(x, y, z,
         size=size,
         col = col.vec,)
}

PCAtext <- function(exp.data,text.vec,cent=FALSE,
                    scle=FALSE,col.vec="black",size=18,adj=0.5,type="p",radius=3,
                    cex=par3d("cex")
                    ){
  
  pc <- pca(t(exp.data),cent=FALSE,scle=FALSE)
  open3d()
  x <- pc$pcs[,1]
  y <- pc$pcs[,2]
  z <- pc$pcs[,3]
  plot3d(x, y, z,
         size=size,
         col = col.vec,
         type=type,
         radius=radius)
  text3d(x,y,z,text=text.vec,adj,cex=cex)
}


PCAloadvecs <- function(exp.data,cent=FALSE,scle=FALSE){

  pc <- pca(t(exp.data),cent=cent,scle=scle)
  z <- pc$pcs[,3]
  x <- pc$pcs[,1]
  y <- pc$pcs[,2]

  pc
}
