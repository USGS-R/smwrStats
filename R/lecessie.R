# Compute the GOF test of le Cressie and van Houwelingen for logistic
#  regression models
#
# Coding history:
#    2009Jan25 DLLorenz Original Coding and begin of modifications
#    2011Aug22 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update for package
#    2012Aug28 DLLorenz Renamed to correct spelling
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2013Apr09 DLLorenz Added setGD to plot
#    2013Sep26 DLLorenz Sligth improvement to speed for large models
#

leCessie.test <- function(object, bandwidth, newterms) {
  ## Arguments:
  ##  object (a glm model object) the logistic regression model
  ##  bandwidth (numeric scalar) the bandwidth to use for the smoothing
  ##  newterms (formula) potential variables to add to the model. Should
  ##   be constructed like ~ . + X3
  ##
  ## Note there are very minor differences between this version and the SAS
  ## code by Saskia le Cessie.
  ## The errors begin with very small differences in the computed distances
  ## and propogate from there.
  ##
  ## Create the categorical distance function
  categorical <- function(x) {
    x <- as.matrix(x)
    retval <- apply(x, 2, function(y) {
      nrcats <- length(unique(y))
      ## allow for the case where the number of categories is 1
      if(nrcats == 1)
        disty <- rep(0, length(y)*(length(y)-1))
      else
        disty <- (dist(y, 'manhattan') != 0) * nrcats/(nrcats - 1)
    })
    retval <- rowSums(retval)
    return(retval)
  }
  ## Keep the original object first 
  object.orig <- object
  Dname <- as.character(object$call$formula)
  Dname <- paste(Dname[2], Dname[1], Dname[3]) # get it in the correct order
  ## Extract the components and compute the smoothed residuals and raw Q
  fits <- fitted(object)
  resids <- resid(object, 'response')
  N <- length(resids)
  if(!missing(newterms))
    object <- update(object, formula=newterms)
  obj.frame <- model.frame(object)
  if(class(obj.frame[[1]]) == 'matrix')
    stop("Cannot perform test where response is matrix")
  ## Compute distances
  distances <- lapply(obj.frame[,-1, drop=FALSE], function(x) {
    if(is.numeric(x))
      return(unclass(0.5*(dist(scale(x)))^2))
    if(class(x) == 'factor')
      return(unclass(categorical(x)))
    ## Must be a structure created from a function that creates a matrix
    return(unclass(0.5*dist(scale(x[[1]]))^2))
  })
  dist.mat <- matrix(0, N, N)
  dist.mat[lower.tri(dist.mat)] <- sqrt(rowSums(as.data.frame(distances)))
  dist.mat <- dist.mat + t(dist.mat)
  if(missing(bandwidth))
    bandwidth <- mean(dist.mat)
  R.raw <- pmax(1 - dist.mat/bandwidth, 0)
  smoothres <- t(resids) %*% R.raw
  Q.raw <- sum(smoothres * resids)
  ## Extract components and compute corrected estimated value of Q
  obj.mat <- model.matrix(object.orig)
  mu2 <- fits * (1 - fits)
  V <- diag(mu2)
  Vx <- diag(V)
  obj.hat <- Vx*obj.mat %*% solve(t(obj.mat) %*% (Vx * obj.mat)) %*% t(obj.mat)
  R.cor <- (diag(N) - obj.hat) # set up to re-use below
  R.cor <- R.cor %*% R.raw %*% R.cor
  E.Q <- sum(diag(R.cor) * mu2)
  ## Extract components and compute corrected std.err of Q
  mu4 <- mu2*(1-3*mu2)
  VarQ1 <- sum(diag(R.cor)^2*(mu4 - 3*mu2^2))
  R.tmp <- R.cor * rep(mu2, each=nrow(R.cor))
  VarQ2 <- 2*sum(diag(R.tmp %*% R.tmp))
  VarQ <- VarQ1 + VarQ2
  Test <- Q.raw * 2 * E.Q / VarQ
  df <- 2 * E.Q^2 / VarQ
  ## make the return value able to be printed with print.htest
  names(Test) <- "Chisq"
  names(df) <- "df"
  retval <- list(statistic=Test, parameters=df,
                 p.value=1 - pchisq(Test, df),
                 estimate=c(Q=Q.raw, 'E[Q]'=E.Q, 'se[Q]'=sqrt(VarQ)),
                 method="le Cessie-van Houwelingen GOF test",
                 data.name=Dname,
                 alternative="Some lack of fit\nnull hypothesis: No lack of fit",
                 ## end of htest material
                 object=object.orig, target.object=object,
                 bandwidth=bandwidth, max.distance=max(dist.mat),
                 smoothed.residuals=as.vector(smoothres),
                 distance.matrix=dist.mat,
                 hat=obj.hat)
  oldClass(retval) <- "lecessie"
  return(retval)
}

print.lecessie <- function(x, digits=4, ...) {
  ## Arguments:
  ##  x (lecessie object) the object to be printed
  ##  digits (integer scalar) the number of digits to print
  ##  ... (dots) not used, required for method function
  ##
  x.to.print <- x
  x.to.print$parameters <- round(x.to.print$parameters, digits) # fix this
  oldClass(x.to.print) <- "htest"
  print(x.to.print)
  cat("Distance between observations:\n")
  print(c(maximum=x$max.distance, bandwidth=x$bandwidth))
  cat("\n")
  invisible(x)
}

plot.lecessie <- function(x, which="All", set.up=TRUE, ...) {
  ## Arguments:
  ##  x (lecessie object) the object to be printed
  ##  which (character or numeric) which plots to plot
  ##  ... (dots) not used, required for method function
  ##
  ## Set up graphics page
  if(set.up) 
    setGD("LECESIE")
  target.mat <- model.matrix(x$target.object)
  target.names <- dimnames(target.mat)[[2L]]
  target.mat[, 1L] <- fitted(x$object)
  target.names[1L] <- "Fitted values"
  target.resid <- residuals(x$object, type="dev")
  if(tolower(which[1]) == "all")
    which <- seq(ncol(target.mat), 1L)
  else if(tolower(which[1]) == "first")
    which <- 1
  else if(is.character(which))
    which <- which(target.names %in% which)
  for(i in which) {
    xyPlot(target.mat[,i,drop=TRUE], x$smoothed.residuals,
           Plot=list(what="points", size=0.05), 
           ytitle="Smoothed residuals", xtitle=target.names[i],
           margin=c(NA, NA, 1.2, NA)) # leave room for title
    if(i == 1)
      addXY(target.mat[,i,drop=TRUE], target.resid,
            Plot=list(what="points", filled=FALSE, size=0.08))
  }
  invisible(x)
}
