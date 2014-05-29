\name{getCNGECorrelation}
\alias{getCNGECorrelation}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Correlation between copy number data and gene expression data.
%%  ~~function to do ... ~~
}
\description{
This function takes a "FirehoseData" object as an input and calculates the coorelation between gene expression data and copy number data and returns "CorResult" list.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getCNGECorrelation(dataObject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataObject}{
  "FirehoseData" data object from "getFirehoseData" function.
%%     ~~Describe \code{dataObject} here~~
}
}
\value{
This function returns a list where each item will be "CorResult" object for one comparison.
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
\examples{
data(RTCGASample)
corRes = getCNGECorrelation(a2)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
