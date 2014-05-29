\name{getFirehoseAnalyzeDates}
\alias{getFirehoseAnalyzeDates}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Firehose data analyze run dates
%%  ~~function to do ... ~~
}
\description{
This function return Firehose data analyze run dates to access recent and previous runs. This functon is useful for GISTIC processed copy number data. To call main funtions users must provide valid dates.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getFirehoseAnalyzeDates(last = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{last}{: Reports last x date(s) instead of all. 
%%     ~~Describe \code{last} here~~
}
}
\value{
This function returns a character vector for valid dates.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\author{
Mehmet Kemal Samur
%%  ~~who you are~~
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
