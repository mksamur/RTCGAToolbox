\name{getMutationRate}
\alias{getMutationRate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculates the mutation ratio among the samples.
%%  ~~function to do ... ~~
}
\description{
This function calculates the gene mutation frequency among the samples. The function takes a "FirehoseData" object as an input parameters and calculates the mutation rate by using the 'Mutation' field of the object. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getMutationRate(dataObject)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataObject}{: Must be "FirehoseData" object.
%%     ~~Describe \code{x} here~~
}
}
\value{
Returns a data frame where first column will be gene symbols and second column will be mutation ratio.
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
mutRate = getMutationRate(dataObject=a2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
