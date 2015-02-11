\name{extract}
\alias{extract}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Extract FirehoseData to Object
%%  ~~function to do ... ~~
}
\description{
This function will extract the dataset of choice from the "FirehoseData" object returned by the getFirehoseData function. 
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
extract(object, type)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{: A "FirehoseData" object from the getFirehoseData function.
%%     ~~Describe \code{x} here~~
}
  \item{type}{: The type of dataset from a list of possible datasets. Choices include: "RNAseq_Gene", "Clinic", "miRNASeq_Gene", "RNAseq2_Gene_Norm", "CNA_SNP", "CNV_SNP", "CNA_Seq", "CNA_CGH", "Methylation", "Mutation", "mRNA_Array", "miRNA_Array", "RPPA", "GISTIC_A", "GISTIC_T".
  }
}
\value{
This function returns an object based on the chosen dataset which can be a matrix or a list.
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\details{
The "GISTIC_A" type of dataset represents GISTIC data by all genes. "GISTIC_T"" represents data thresholded by genes. 
}
\author{
Levi Waldron
%%  ~~who you are~~
}
\examples{
\dontrun{
a5 <- extract(object=myS4obj,type="Clinic")}
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
