#' Get data analyze dates.
#'
#' \code{getFirehoseAnalyzeDates} returns the character vector for analyze release dates.
#'
#' @param last To list last n dates. (Default NULL)
#' @return A character vector for dates.
#' @examples
#' getFirehoseAnalyzeDates(last=2)
#' @export getFirehoseAnalyzeDates
getFirehoseAnalyzeDates <- function(last=NULL) {
    dates <-
    c("20160128", "20150821", "20150402", "20141017", "20140715",
    "20140416", "20140115", "20130923", "20130523", "20130421", "20130326",
    "20130222", "20130116", "20121221", "20121024", "20120913", "20120825",
    "20120725", "20120623", "20120525", "20120425", "20120321", "20120217",
    "20120124", "20111230", "20111128", "20111026", "20110921", "20110728",
    "20110525", "20110421", "20110327", "20110217", "20110114", "20101223")

    if (!is.null(last))
        head(dates, last)
    else
        dates
}
