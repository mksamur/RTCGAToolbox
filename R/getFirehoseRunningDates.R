#' Get standard data running dates.
#'
#' \code{getFirehoseRunningDates} returns the character vector for standard data release dates.
#'
#' @param last To list last n dates. (Default NULL)
#' @return A character vector for dates.
#' @examples
#' getFirehoseRunningDates()
#' getFirehoseRunningDates(last=2)
#' @export getFirehoseRunningDates
getFirehoseRunningDates <- function(last=NULL){
  check.integer <- function(N){
    !length(grep("[^[:digit:]]", as.character(N)))
  }
  if(is.null(last)){
    runDate <- read.table("http://www.canevolve.org/fmineRdate.txt",
                          header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    runDate <- as.character(runDate[,1])
    return(runDate)
  }
  else if(check.integer(last))
  {
    runDate <- read.table("http://www.canevolve.org/fmineRdate.txt",
                          header=FALSE, sep="\t", na.strings="NA", dec=".", strip.white=TRUE)
    runDate <- as.character(runDate[,1])
    if(last < length(runDate)){runDate <- runDate[1:last]}
    return(runDate)
  }
  else
  {
    stop('"last" must be integer')
  }
}