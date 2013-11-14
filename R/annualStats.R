# retrieve dv data and compute summary stats
#
# Coding history:
#    2012Dec05 DLLorenz Original coding
#    2012Dec05          This version.

annualStats <- function(gage, dtype = "sw", year = "water",
                        begin.date = "", end.date = "", param = NULL,
                        complete=TRUE, statistic=mean, ...) {
  ## Arguments:
  ##  gage (character scalar) the USGS station number
  ##  dtype (character scalar) the type of data requested
  ##  year (character scalar) the type of year
  ##  begin.date (POSIX format string) the earliest date to retrieve data
  ##  end.date (POSIX format string) the latest date to retrieve data
  ##  param (character scalar) The parameter code to get if not default
  ##  complete (logical) require complete years (no missing values)
  ##   and remove provisional data
  ##  statistic, the function to compute the stat
  ##  dots additional arguments to the statistic
  ##
  ## Prelims and make sure dtype and year are valid
  stnam <- deparse(substitute(statistic))
  dtype <- match.arg(dtype, c("sw", "gw"))
  year <- match.arg(year, c("water", "calendar"))
  ## Get the param if necessary
  if(is.null(param) && dtype == "sw")
    param <- "00060"
  else if(is.null(param) && dtype == "gw")
    param <- "72019"
  ## Override the default behavior
  if (begin.date == "") 
    begin.date <- "1860-01-01"
  ## Get the data and set column names
  data <- suppressWarnings(readNWIS(gage=gage,
                                    dtype=paste(dtype, "dv", sep=""),
                                    begin.date=begin.date, end.date=end.date,
                                    param=param))
  col <- paste("X", param, sep="")
  colcd <- paste("X", param, "_cd", sep="")
  names(data)[4L] <- col
  names(data)[5L] <- colcd
  if(year == "water") { # make name an fix column later
    data$year <- waterYear(data$datetime)
    yrnam="WY"
  }
  else {
    data$year <- year(data$datetime)
    yrnam="CY"
  }
  ## If complete, first remove provisional data (colcd == "P")
  ##  and then check, remove incomplete years
  if(complete) {
    pick <- data[[colcd]] %cn% "P"
    data <- data[!pick, ]
    pick <- screenData(data$datetime, data[[col]], year=year, printit=FALSE)
    ## The screenData function returns a matrix of missing data counts for
    ## each month(column) for each year(row)
    pick <- rowSums(pick)
    pick <- pick[pick == 0L]
    if(length(pick) == 0L)
      stop("No valid data for annual summary statistic")
    data <- data[data$year %in% as.numeric(names(pick)), ]
  }
  ## Do the stats!
  retval <- tapply(data[[col]], data$year, FUN=statistic, ..., simplify=FALSE)
  retval <- do.call(rbind, retval)
  if(ncol(retval) == 1L) # otherwise assume that names are OK
    colnames(retval) <- stnam
  retval <- data.frame(agency_cd=data[1L, 1L], site_no=data[1L, 2L],
                       year=as.integer(rownames(retval)), param=param,
                       retval, stringsAsFactors=FALSE)
  ## Fix the year name
  names(retval)[3] <- yrnam
  if(!complete) { # tell 'em what they really got
    Nobs <- tapply(data[[col]], data$year, function(x) sum(!is.na(x)))
    retval$N <- as.vector(Nobs)
  }
  return(retval)
}
