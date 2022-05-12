MMest_multireg.formula <-
function(formula, data=NULL, ...){
  
  
  model.multiregresp<-function (data, type = "any"){
    
    if(attr(attr(data, "terms"), "response")){
      if(is.list(data) | is.data.frame(data)){
        v <- data[[1L]]
        if(is.data.frame(data) && is.vector(v)) v <- data[,1L,drop=FALSE]
        if(type == "numeric" && is.factor(v)){
          warning("using type=\"numeric\" with a factor response will be ignored")
        }else if(type == "numeric" | type == "double") 
          storage.mode(v) <- "double"
        else if(type != "any") 
          stop("invalid response type")
        if(is.matrix(v) && ncol(v) == 1L){ 
          if(is.data.frame(data)){v=data[,1L,drop=FALSE]}
          else{dim(v) <- NULL}}
        rows <- attr(data, "row.names")
        if(nrows <- length(rows)){
          if(length(v) == nrows) 
            names(v) <- rows
          else if(length(dd <- dim(v)) == 2L) 
            if(dd[1L] == nrows && !length((dn <- dimnames(v))[[1L]])) 
              dimnames(v) <- list(rows, dn[[2L]])
        }
        return(v)
      }else stop("invalid 'data' argument")
    }
    else return(NULL)
  }
  
  
  mt <- terms(formula, data = data)
  if(attr(mt, "response") == 0L) stop("response is missing in formula")
  mf <- match.call(expand.dots = FALSE)
  mf$... <- NULL
  mf[[1L]] <- as.name("model.frame")
  mf <- eval.parent(mf)
  miss <- attr(mf,"na.action")
  Y <- model.multiregresp(mf)
  Terms <- attr(mf, "terms")
  X <- model.matrix(Terms, mf)
  res <- MMest_multireg.default(X, Y, int = FALSE, ...)
  res$terms <- Terms
  cl <- match.call()
  cl[[1L]] <- as.name("MMest_multireg")
  res$call <- cl
  res$xlevels <- .getXlevels(mt, mf)
  if(!is.null(miss)) res$na.action <- miss
  return(res)
}
