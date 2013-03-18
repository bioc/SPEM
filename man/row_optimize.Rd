\name{row_optimize}
\alias{row_optimize}
\alias{row_optimize,ExpressionSet-method}
\title{
Calculate parameters for one row
}
\description{
This function calculates parameters for a single row in the expression data. If a large-size dataset will be calculated, this function is recommended. 
}
\usage{
\S4method{row_optimize}{ExpressionSet}(TS_eSet,S,beta, sparsity = 0.2, lbH = -3, ubH = 3, lbB = 0, ubB = 10)
}
\section{Methods}{
\describe{
\item{\code{signature(TS_eSet = "ExpressionSet")}}{
This method is created for the function \code{row_optimize}.
}
}}
\arguments{
  \item{TS_eSet}{
Time series data in ExpressionSet class
assayData: Matrix with n metabolite in row and m time points in column.
phenoData: Dataframe includes label "time", which represents the time points.
}
  \item{S}{
Slope of the row you want to calculated. You can either input a vector with length equal to the rows of assayData of TS_eSet, or use s_diff function in this package to calculate it.
}
  \item{beta}{
Initial beta.
}
  \item{sparsity}{
A threshold used to control the sparsity of reconstructed matrix. Values whose absolute value smaller than sparsity will be set to zero.
}
  \item{lbH}{
Lower boundary value of h.
}
  \item{ubH}{
Upper boundary value of h.
}
  \item{lbB}{
Lower boundary value of beta.
}
  \item{ubB}{
Upper boundary value of beta.
}
}
\details{
 In this SPEM package, we aim to reconstruct gene networks from time-series expression data using the S-system model. The input dataset should be as an ExpressionSet data container, describing, in assayData, expression data for n genes (rows) and m time points (columns), along with a vector of length m, which records the exact values of time points, thus showing the sample intervals in phenoData. SPEM will calculate the parameters alpha, g, beta and h of the S-system function set that best fits the dataset.

  
  In this function, user can calculate one row at a time. This function offers a parallel calculation option for users.
}
\value{
 This function return a vector bind with c(alpha, $g_i$, beta, $h_i$, Initial Beta, error).
}

\author{
Yang, X-Y, Dent, Jennifer E. and Nardini, C.
}



\examples{
#########Load the SOS pathway data #######
data(sos)

#########Set Slope and Initial Beta #######

Slope<- s_diff(sos) 
S<- Slope[1,] #S is the slope of the row you want to calculate. You can either input a vector yourself. 
beta<- runif(n=1,min=1,max=10)

#########Set parameters #######
sparsity<- 0.2
lbH<- -3
ubH<- 3
lbB<- 0
ubB<- 10

#########Calculate results #######

result_r<-row_optimize(sos,S,beta,sparsity,lbH,ubH,lbB,ubB)
}
\keyword{Reconstruction}