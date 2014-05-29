\name{getSurvival}
\alias{getSurvival}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Survival analysis based on gene expression value.
%%  ~~function to do ... ~~
}
\description{
Either 2 or 3 groups will be created by using expression data, depends on "numberofGroups" parameters. And difference between groups will be compared. A KM plot will be returned as a final product.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
getSurvival(dataObject, numberofGroups = 2, geneSymbols, sampleTimeCensor)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dataObject}{: A "FirehoseData" object includes the expression matrix.
%%     ~~Describe \code{dataObject} here~~
}
  \item{numberofGroups}{: Number of comparison group. If set as 2 median expression value will be used and lower than median samples will be compared higher than median expression samples. If set as 3, summary function output will be used for defining groups. First froup will be samples in the 1st Q. , second group will be samples between 1st Q and 3rd Q. and 3rd group memebers will be samples those expression value higher than 3rd Q. 
%%     ~~Describe \code{numberofGroups} here~~
}
  \item{geneSymbols}{: A gene symbols which users want to test.
%%     ~~Describe \code{geneSymbols} here~~
}
  \item{sampleTimeCensor}{: A data frame. First column must be sample barcodes, second column must be time and last column must be censor infromation. See "UserGuide" for details. 
%%     ~~Describe \code{sampleTimeCensor} here~~
}
}
\value{
Creates KM plots for the genes and calculates difference between groups.
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
clinicData = a2@Clinical
clinicData = clinicData[3:5,2:ncol(clinicData)]
clinicData[3,is.na(clinicData[3,])] = clinicData[2,is.na(clinicData[3,])]
survData <- data.frame(Samples=colnames(clinicData),
              Time=as.numeric(clinicData[3,]),
              Censor=as.numeric(clinicData[1,]))
getSurvival(dataObject=a2,geneSymbols=c("RRM2","FAM111B"),sampleTimeCensor=survData)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
