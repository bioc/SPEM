\name{s_diff}
\alias{s_diff}
\alias{s_diff,ExpressionSet-method}
\title{
Calculate slopes from time points and time series matrix.
}
\description{
This function allows users calculate slopes from time points and time series data.
}
\usage{
\S4method{s_diff}{ExpressionSet}(TS_eSet)
}
\section{Methods}{
\describe{
\item{\code{signature(TS_eSet = "ExpressionSet")}}{
This method is created for function \code{s_diff}.
}
}}

\arguments{
  \item{TS_eSet}{
Time series data in ExpressionSet class. 
assayData: Matrix with n metabolite in row and m time points in column. 
phenoData: phenoData type. The sample data.frame should include the label "time", which represents the values of time points.
}
}

\value{
This function directly return a slope matrix. 
}

\author{
Yang, X-Y, Dent, Jennifer E. and Nardini, C.
}

\examples{

#########Load the SOS pathway data #######
data(sos)

#########Calculate results #######

Slope<-s_diff(sos)
}
\keyword{Slope}