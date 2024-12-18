#' Get standard data running dates.
#'
#' `getFirehoseRunningDates` returns the character vector for standard data release dates.
#'
#' @param last To list last n dates. (Default NULL)
#' @return A character vector for dates.
#' @examples
#' getFirehoseRunningDates()
#' getFirehoseRunningDates(last=2)
#' @export getFirehoseRunningDates
getFirehoseRunningDates <- function(last=NULL){
    dates <- c("20160128", "20151101", "20150821", "20150601", "20150402",
    "20150204", "20141206", "20141017", "20140902", "20140715", "20140614",
    "20140518", "20140416", "20140316", "20140215", "20140115", "20131210",
    "20131114", "20131010", "20130923", "20130809", "20130715", "20130623",
    "20130606", "20130523", "20130508", "20130421", "20130406", "20130326",
    "20130309", "20130222", "20130203", "20130116", "20121221", "20121206",
    "20121114", "20121102", "20121024", "20121020", "20121018", "20121004",
    "20120913", "20120825", "20120804", "20120725", "20120707", "20120623",
    "20120606", "20120525", "20120515", "20120425", "20120412", "20120321",
    "20120306", "20120217", "20120124", "20120110", "20111230", "20111206",
    "20111128", "20111115", "20111026")

    if (!is.null(last))
        head(dates, last)
    else
        dates
}
