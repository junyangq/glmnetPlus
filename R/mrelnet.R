mrelnet=function(x,is.sparse,ix,jx,y,weights,offset,alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,jsd,intr,vnames,maxit,beta0){
  maxit=as.integer(maxit)
  weights=as.double(weights)
### compute the null deviance
  y=as.matrix(y)
  nr=ncol(y)
  if(nr>1){
    responseNames=dimnames(y)[[2]]
    if(is.null(responseNames))responseNames=paste("y",seq(nr),sep="")
  }

 storage.mode(y)="double"
   if(is.null(offset)){
    is.offset=FALSE}
  else{
    offset=as.matrix(offset)
    storage.mode(offset)="double"
    if(!all.equal(dim(offset),dim(y)))stop("Offset must match dimension of y")
    y=y-offset
    is.offset=TRUE
    }

  nulldev=0
  for(i in seq(nr)){
  ybar=weighted.mean(y[,i],weights)
  nulldev=nulldev+sum(weights* (y[,i]-ybar)^2)
}

  storage.mode(nr)="integer"

  if (!is.sparse && is.null(beta0)) {
    beta0 <- matrix(0, ncol(x), ncol(y))
  }

fit=if(is.sparse) {.Fortran("multspelnet",
        parm=alpha,nobs,nvars,nr,x,ix,jx,y,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,jsd,intr,maxit,
        lmu=integer(1),
        a0=double(nlam*nr),
        ca=double(nx*nlam*nr),
        ia=integer(nx),
        nin=integer(nlam),
        rsq=double(nlam),
        alm=double(nlam),
        nlp=integer(1),
        jerr=integer(1),PACKAGE="glmnetPlus"
        )
} else {
  dotCall64::.C64("multelnet",
       SIGNATURE = c("double", "integer", "integer", "integer", "double", "double", "double",
                     "integer", "double", "double", "integer", "integer", "integer", "double",
                     "double", "double", "integer", "integer", "integer", "integer", "double",
                     "integer", "double", "double", "integer", "integer", "double", "double",
                     "integer", "integer"),
       parm = alpha, nobs, nvars, nr, as.double(x), y, weights, 
       jd, vp, cl, ne, nx, nlam, flmin, 
       ulam, thresh, isd, jsd, intr, maxit, as.double(beta0),
       lmu = integer(1),
       a0 = double(nlam*nr),
       ca = double(nx*nlam*nr),
       ia = integer(nx),
       nin = integer(nlam),
       rsq = double(nlam),
       alm = double(nlam),
       nlp = integer(1),
       jerr = integer(1),
       INTENT = c(rep("r", 21), rep("w", 9)),
       PACKAGE = "glmnetPlus"
       )
}
# .Fortran("multelnet",
#           parm=alpha,nobs,nvars,nr,as.double(x),y,weights,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,jsd,intr,maxit,as.double(beta0),
#           lmu=integer(1),
#           a0=double(nlam*nr),
#           ca=double(nx*nlam*nr),
#           ia=integer(nx),
#           nin=integer(nlam),
#           rsq=double(nlam),
#           alm=double(nlam),
#           nlp=integer(1),
#           jerr=integer(1),PACKAGE="glmnetPlus"
#           )
if(fit$jerr!=0){
  errmsg=jerr(fit$jerr,maxit,pmax=nx,family="mrelnet")
  if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
  else warning(errmsg$msg,call.=FALSE)
}
  if(nr>1)
    outlist=getcoef.multinomial(fit,nvars,nx,vnames,nr,responseNames,center.intercept=FALSE)
  else
    outlist=getcoef(fit,nvars,nx,vnames)
  dev=fit$rsq[seq(fit$lmu)]
  outlist=c(outlist,list(dev.ratio=dev,nulldev=nulldev,npasses=fit$nlp,jerr=fit$jerr,offset=is.offset))
  class(outlist)="mrelnet"
  outlist
}
