\name{sos}
\alias{sos}
\docType{data}
\title{SOS pathway time series data}
\description{
 In this package we offer the SOS data obtained from Uri Alon's lab (http://www.weizmann.ac.il/mcb/UriAlon/). SOS response is a general DNA repair system in bacteria which allows survival after DNA damage. This SOS dataset is taken from real experiment expression data in \emph{Escherichia coli}. It contains 8 genes under Experiment 3 (UV light intensities, 4:20 Jm^(-2)).
 }
\usage{data(sos)}
\format{
   sos.data is time series gene expression value data in ExpressionSet Class.
  assayData: Matrix with expression values of 8 genes in SOS pathway of  \emph{Escherichia coli}. These expression levels are observed at 50 time points. 
  phenoData: Sample data.frame includes label "time", which represents the value of time points.}
\references{
  M. Ronen, R. Rosenberg, B. I. Shraiman, and U. Alon. Assigning numbers to the arrows: parameterizing a gene regulation network by using accurate expression kinetics. Proceedings of the National Academy of Sciences of the United States of America, Aug. 2002. PMID: 12145321. 
}
\examples{
data(sos)
}
\keyword{Datasets}