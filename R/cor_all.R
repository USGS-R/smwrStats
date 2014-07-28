# Correlations between numeric columns
#
# Coding history:
#    2004Oct14 DLLorenz Initial coding for library.
#    2005Jul14 DLLorenz date fix
#    2006Feb03 DLLorenz Bug Fix when not all data were numeric
#    2008Apr21 DLLorenz Continuity Adjustment to Kendall's method
#    2011Jul28 DLLorenz Conversion to R
#    2011Oct25 DLLorenz Update for package
#    2012Mar08 DLLorenz Added summary function
#    2012Jun19 DLLorenz Tweeks in summary and other methods
#    2012Aug28 DLLorenz Change from diagPlot to plot
#    2012Dec21 DLLorenz Change gaussian to normal
#    2013Apr01 DLLorenz Bug fix in printed output
#    2013Apr09 DLLorenz Added setGD to plot
#    2013Apr31 DLLorenz Switched which to y-axis and bug fix
#    2014May27 DLLorenz Added log1p distribution

cor.all <- function(data, method="pearson", na.method="pairwise",
                    distribution="normal") {
  ## Arguments:
  ##  data (a data frame or matrix) the data on whihc to compute the cors
  ##  method (character scalar) the method to use for the cors
  ##  na.method (character scalar) how to handle missing values
  ##  distribution (character scalar) use log-transforms "lognormal"
  ##    or not "normal"
  ##
  ## Compute correlations between all numeric columns in data
  data.name <- deparse(substitute(data))
  data <- as.data.frame(data) # force to be a data. frame
  nums <- unlist(lapply(data, is.numeric))
  nams <- names(data)[nums]
  nums <- which(nums)
  N <- length(nams)
  if(N < 2)
    stop("Insufficient numeric data to compute correlations")
  ## Set up for handling NAs
  na.methods <- c("fail", "omit", "pairwise")
  na.meth <- pmatch(na.method, na.methods, nomatch=0)
  data <- data[, nums] # extract the data to be analyzed
  if(na.meth == 1) { # check if any NAs
    if(any(is.na(unlist(data))))
      stop("There are missing values in the numeric data")
  }
  else if(na.meth == 2)
    data <- na.omit(data)
  else if(na.meth == 0)
    stop("cor.all: Invalid na.method specified")
  ## Apply log-transform if requested
  distribution <- match.arg(distribution, c("normal", "lognormal", "log1p"))
  if(distribution == "lognormal") {
  	data <- as.data.frame(lapply(data, log))
  } else if(distribution == "log1p") 
  	data <- as.data.frame(lapply(data, log1p))
  cor.grid <- expand.grid(1:N, 1:N)
  cor.grid <- cor.grid[cor.grid[,1] < cor.grid[,2],]
  cor.stat <- matrix(1, ncol=N, nrow=N)
  cor.pval <- matrix(0, ncol=N, nrow=N)
  ## Construct matrix of counts of samples
  cor.count <-as.matrix( !is.na(data))
  cor.count <- t(cor.count) %*% cor.count
  ## Protect against warning messages for spearman and kendall options
  warn <- options("warn")
  options(warn=-1)
  method=match.arg(method, c("pearson", "spearman", "kendall"))
  for(i in 1:nrow(cor.grid)) {
    Xs <- cor.grid[i,1]
    Ys <- cor.grid[i,2]
    temp <- try(cor.test(data[[Xs]], data[[Ys]], alternative="two.sided",
                         method=method, continuity=TRUE))
    if(class(temp) != "Error") {
      cor.pval[Xs, Ys] <- cor.pval[Ys, Xs] <- temp$p.value
      cor.stat[Xs, Ys] <- cor.stat[Ys, Xs] <- temp$estimate
      meth <- temp$method
    } else {
      cor.pval[Xs, Ys] <- cor.pval[Ys, Xs] <- NA
      cor.stat[Xs, Ys] <- cor.stat[Ys, Xs] <- NA
    }
  }
  options(warn)
  dimnames(cor.stat) <- list(nams, nams)
  dimnames(cor.pval) <- list(nams, nams)
  retval <- list(estimates = cor.stat, p.values = cor.pval, counts = cor.count,
                 alternative = "two.sided", na.method = na.methods[na.meth],
                 method = meth, data.name = data.name, data=data,
                 call.method=method, distribution=distribution)
  oldClass(retval) <- "cor.all"
  return(retval)
}

print.cor.all <- function(x, digits = 4, lower=TRUE, ...) {
  ## Arguments:
  ##  x (cor.all object) the object to be printed
  ##  digits (integer scalar) the number of digits to print
  ##  lower (logical scalar) print only the lower tri?
  ##  ... (dots) not used, required for method function
  ##
  cat(paste("\n\t", x$method, "\n\n", sep = ""))
  cat("data: ", x$data.name, "\n\n")
  if(x$distribution == "lognormal") {
  	cat("All data were log-transformed\n")
  } else if(x$distribution == "log1p") 
  	cat("All data were log1p-transformed\n")
  est <- format(x$estimates, digits=digits)
  N <- ncol(est)
  val <- round(x$p.values, digits=digits)
  val0 <- (val == 0) & (abs(x$estimates) < 1)
  val <- format(val, width=digits+3, justify="right")
  if(any(val0)) {
    vallt <- paste("<0.", paste(rep(0, digits-1), collapse=''), "1", sep='')
    val[val0] <- vallt
  }
  val <- matrix(val, N, N) # Need to restore dim
  cnt <- format(x$counts)
  if(lower) {
    est <- est[-1, -N, drop=FALSE]
    val <- val[-1, -N, drop=FALSE]
    cnt <- cnt[-1, -N, drop=FALSE]
    N <- N - 1
    if(N > 1) {
      for(i in seq(1, N-1))
        for(j in seq(i+1, N)) {
          est[i, j] <- " "
          val[i, j] <- " "
          cnt[i, j] <- " "
        }
    }
  }
  outmat <- rbind(est[1,], val[1,], cnt[1,], ' ')
  if(N > 1)
    for(i in seq(2, N))
      outmat <- rbind(outmat, est[i,], val[i,], cnt[i,], ' ')
  ## Make row names descriptive
  dn <- dimnames(est)
  dn[[1]] <- rep(dn[[1]], each=4)
  dn[[1]] <- switch(x$call.method,
                    pearson=paste(dn[[1]], c("-cor", "-r!=0", "-N", "--"), sep=''),
                    spearman=paste(dn[[1]], c("-rho", "-r!=0", "-N", "--"), sep=''),
                    kendall=paste(dn[[1]], c("-tau", "-t!=0", "-N", "--"), sep=''))
  ## Blank line--blank row name
  dn[[1]] <- sub(".*--$", "", dn[[1]])
  dimnames(outmat) <- dn
  print(outmat, quote=FALSE)
  invisible(x)
}

plot.cor.all <- function(x, which="All", set.up=TRUE, ...) {
  ## Arguments:
  ##  x (a cor.all object) the object ot be plotted
  ##  which (character scalar) which plots to plot
  ##    'All' means full matrix, 'Lower' means lower tri implemented only 
  ##    when ncol is 5 or less, other is show correlations for variable
  ##  ... (dots) not used, required for method funtion
  ##
  ## Use splomPlot() with a neet little panel function for ncol(x$data <=4)
  ## and paired x-y plots n*(n-1)/2 otherwise. Use a 3.5x6.5 vertical figure layout.
  ## Set up graphics page
  if(set.up) 
    setGD("CORRELATIONS")
  ## Continue
  data <- x$data
  call.method <- x$call.method
  if(x$distribution == "lognormal") {
    names(data) <- paste("log", names(data), sep='.') # fix for plot only
    rownames(x$estimates) <- paste("log", rownames(x$estimates), sep='.')
    colnames(x$estimates) <- paste("log", colnames(x$estimates), sep='.')
    rownames(x$p.values) <- paste("log", rownames(x$p.values), sep='.')
    colnames(x$p.values) <- paste("log", colnames(x$p.values), sep='.')
  } else if(x$distribution == "log1p") {
  	names(data) <- paste("log1p", names(data), sep='.') # fix for plot only
  	rownames(x$estimates) <- paste("log1p", rownames(x$estimates), sep='.')
  	colnames(x$estimates) <- paste("log1p", colnames(x$estimates), sep='.')
  	rownames(x$p.values) <- paste("log1p", rownames(x$p.values), sep='.')
  	colnames(x$p.values) <- paste("log1p", colnames(x$p.values), sep='.')
  }
  N <- ncol(data)
  ckN <- as.integer(min(par("fin")) / 1.6249)
  if(which == "All") {
    Nok <- N <= ckN
    show.all <- TRUE
  } else if(which == "Lower") {
    Nok <- N <= (ckN + 1L)
    show.all <- FALSE
  } else { # Show which on the y-axis
    titles <- names(data)
    est <- c(pearson="cor", spearman="rho", kendall="tau")[call.method]
    if(x$distribution == "lognormal") {
    	which <- paste("log", which, sep='.')
    } else if(x$distribution == "log1p")
    	which <- paste("log1p", which, sep='.')
    for(j in titles) {
    	if(j != which) {
        AA.xy <- xyPlot(data[[j]], data[[which]], Plot=list(what="points"),
                        ylabels=5, xlabels=5,
                        xtitle=j, ytitle=which, margin=c(NA, NA, 1.6, 0.5))
        addTitle(paste(est, " = ", round(x$estimates[which,j], 3),
                       ", p.val = ", round(x$p.values[which,j], 3), sep=''))
        addSmooth(AA.xy$x, AA.xy$y, span=1.0, current=AA.xy)
      }
    }
    return(invisible(x))
  }
  ## Proceed with splom/paired plots
  if(Nok) {
    AA <- setSplom(num.variables=N, show.all=show.all, ymargin=3.5)
    splomPlot(data, AA, Panel=function(x, y, cur) {
      test <- cor.test(x, y, method=call.method)
      est <- names(test$est)
      tbl <- rbind(paste(est, " = ", round(test$est, 3), sep=''),
                   paste("p.val = ", round(test$p.value, 3), sep=''))
      addTable(tbl, "ul")
      return(addSmooth(x, y, span=1.0, current=cur))})
  } else { # more than 5 columns
    AA <- setLayout(width=3.5, height=c(3.35, 3.15), , num.rows=2, xtop=1.2)
    needpage <- FALSE
    titles <- names(data)
    est <- c(pearson="cor", spearman="rho", kendall="tau")[call.method]
    for(j in seq(N, 2L, -1L)) {
      for(i in seq(j - 1L, 1L, -1L)) {
        if(needpage)
          plot.new() # Sets up a new page
        needpage <- TRUE
        AA.m <- setGraph(1, AA)
        AA.m[3L] <- 1.6 # Need to allocate space for a title
        AA.xy <- xyPlot(data[[i]], data[[j]], Plot=list(what="points"),
                        ylabels=5, xlabels=5,
                        xtitle=titles[i], ytitle=titles[j], margin=AA.m)
        addTitle(paste(est, " = ", round(x$estimates[i,j], 3),
                       ", p.val = ", round(x$p.values[i,j], 3), sep=''))
        addSmooth(AA.xy$x, AA.xy$y, span=1.0, current=AA.xy)
        ## Now plot j,i
        AA.m <- setGraph(2L, AA)
        AA.xy <- xyPlot(data[[j]], data[[i]], Plot=list(what="points"),
                        ylabels=5, xlabels=5,
                        xtitle=titles[j], ytitle=titles[i], margin=AA.m)
        addSmooth(AA.xy$x, AA.xy$y, span=1.0, current=AA.xy)
      } # end of j
    } # End of i
  } # end of else (too many to splom)
  invisible(x)
}


summary.cor.all <- function(object, p.adjust="holm", variable=NULL, ...) {
  NN <- nrow(object$estimates)
  RN <- matrix(rownames(object$estimates), nrow=NN, ncol=NN)
  CN <- matrix(rep(colnames(object$estimates), each=NN), ncol=NN)
  if(is.null(variable))
    sel <- lower.tri(RN)
  else {
    pck <- which(variable == rownames(object$estimates))
    if(length(pck) != 1)
      stop(variable, " not found in list of variables")
    sel <- RN == variable
    sel[pck, pck] <- FALSE
  }
  retval <- data.frame(Var1=RN[sel], Var2=CN[sel], Cor=object$estimates[sel], 
                       Pval=p.adjust(object$p.values[sel], method=p.adjust),
                       Counts=object$counts[sel], stringsAsFactors=FALSE)
  if(p.adjust != "none")
    names(retval)[4] <- paste("Pval", p.adjust, sep='.')
  names(retval)[3] <- switch(object$call.method,
                             pearson="Cor",
                             spearman="Rho",
                             kendall="Tau")
  return(retval)
}
