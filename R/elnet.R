elnet=function(x,is.sparse,ix,jx,y,weights,offset,type.gaussian=c("covariance","naive"),alpha,nobs,nvars,jd,vp,cl,ne,nx,nlam,flmin,ulam,thresh,isd,intr,vnames,maxit,
               beta0,isg,plam,mem.save){
  maxit=as.integer(maxit)
  weights=as.double(weights)
  type.gaussian=match.arg(type.gaussian)

  ka=as.integer(switch(type.gaussian,
    covariance=1,
    naive=2,
    ))


 storage.mode(y)="double"
   if(is.null(offset)){
    is.offset=FALSE}
  else{
    storage.mode(offset)="double"
    is.offset=TRUE
    y=y-offset
  }
### compute the null deviance
  ybar=weighted.mean(y,weights)
  nulldev=sum(weights* (y-ybar)^2)
if(nulldev==0)stop("y is constant; gaussian glmnet fails at standardization step")

  if (!is.null(beta0) && (type.gaussian != "naive" || flmin < 1.0)) {
    stop("Initialization is not yet supported for automatic lambda sequence or non-naive update.\n")
  }

  if (isg && is.null(plam)) {
    warning("isg is on, but plam not provided - fall back to isg off")
    isg <- FALSE
  }

  if (!isg) plam <- -1.0

  if (!is.sparse && is.null(beta0)) {
    beta0 <- rep(0, ncol(x))
  }

  fit = if(is.sparse) {
    dotCall64::.C64("spelnet",
         SIGNATURE = c("integer", "double", "integer", "integer", "double", "integer", "integer", "double",
                       "double", "integer", "double", "double", "integer", "integer", "integer", "double",
                       "double", "double", "integer", "integer", "integer",
                       "integer", "double", "double", "integer", "integer",
                       "double", "double", "integer", "integer"),
         ka, parm = alpha, nobs, nvars, x, ix, jx, y,
         weights, jd, vp, cl, ne, nx, nlam, flmin,
         ulam, thresh, isd, intr, maxit,
         lmu = integer(1), a0 = double(nlam), ca = double(nx*nlam), ia = integer(nx), nin = integer(nlam),
         rsq = double(nlam), alm = double(nlam), nlp = integer(1), jerr = integer(1),
         INTENT = c(rep("r", 21), rep("w", 9)),
         PACKAGE = "glmnetPlus"
    )
  } else {
    is.tbl <- data.table::is.data.table(x)
    if (is.tbl) x <- as.matrix(x)
    dotCall64::.C64("elnet",
                    SIGNATURE = c("integer", "double", "integer", "integer", "double", "double",
                                  "double", "integer", "double", "double", "integer", "integer", "integer", "double",
                                  "double", "double", "integer", "integer", "integer", "double", "integer", "double",
                                  "integer", "double", "double", "integer", "integer",
                                  "double", "double", "integer", "integer", "double"),
                    ka, parm = alpha, nobs, nvars, x, y,
                    weights, jd, vp, cl, ne, nx, nlam, flmin,
                    ulam, thresh, isd, intr, maxit, beta0, isg, plam,
                    lmu = integer(1), a0 = double(nlam), ca = double(nx*nlam), ia = integer(nx), nin = integer(nlam),
                    rsq = double(nlam), alm = double(nlam), nlp = integer(1), jerr = integer(1), residuals = double(nobs*nlam),
                    INTENT = c(rep("rw", 4), ifelse(is.tbl||mem.save, "r", "rw"), rep("rw", 17), rep("w", 10)),
                    PACKAGE = "glmnetPlus"
                    )
  }

if(fit$jerr!=0){
  errmsg=jerr(fit$jerr,maxit,pmax=nx,family="gaussian")
  if(errmsg$fatal)stop(errmsg$msg,call.=FALSE)
  else warning(errmsg$msg,call.=FALSE)
}
  outlist=getcoef(fit,nvars,nx,vnames)
  dev=fit$rsq[seq(fit$lmu)]
  outlist=c(outlist,list(dev.ratio=dev,nulldev=nulldev,npasses=fit$nlp,jerr=fit$jerr,offset=is.offset))
  if ("residuals" %in% names(fit)) {
    residuals=matrix(fit$residuals, nrow=nobs)[, 1:fit$lmu, drop = FALSE]
    outlist[["residuals"]] <- residuals
  }
  class(outlist)="elnet"
  outlist
}
