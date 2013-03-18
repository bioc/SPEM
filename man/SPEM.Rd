\name{SPEM}
\alias{SPEM}
\alias{SPEM,ExpressionSet-method}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
S-system parameter estimation method
}
\description{
This function calculates parameters of S-system from entire time series matrix.
}
\usage{
\S4method{SPEM}{ExpressionSet}(TS_eSet, n = 3, sparsity = 0.2, lbH = -3, ubH = 3, lbB = 0, ubB = 10)
}
\section{Methods}{
\describe{
\item{\code{signature(TS_eSet = "ExpressionSet")}}{
This method is created for function \code{SPEM}.
}
}}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{TS_eSet}{
Time series data in ExpressionSet class. 
assayData: Matrix with n metabolite in row and m time points in column. 
phenoData: phenoData type. The sample data.frame should include the label "time", which represents the values of time points.
}
  \item{n}{
Positive integer, SPEM will guess initial beta n times.
}
  \item{sparsity}{
A positive number. In order to force the interaction matrix to be sparse, interactions with absolute value smaller than "sparsity" will be set to zero.
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
 In this SPEM package, we aim to reconstruct gene networks from time-series expression data using the S-system model. The input dataset should be as an ExpressionSet data container, describing, in assayData, expression data for n genes (rows) and m time points (columns), along with a vector of length m, which records the exact values of time points, thus showing the sample intervals in phenoData. SPEM will calculate the parameters alpha, G, beta and H of the S-system function set that best fits the dataset.
}
\value{
  \item{alpha, G, beta, H}{Parameters of the reconstructed S-system.}
  \item{IniBeta}{Guess of the IniBeta value (Picked randomly by SPEM itself).}
  \item{error}{Regression error.} 
}


\author{
Yang, X-Y., Dent, Jennifer E. and Nardini, C.
}

\examples{

#########Generate Toy Model #######
#########
# If you want to calculate SOS dataset in this package, please read our vignette###
#Real dataset takes a long time to calculate. You may want to try function 'row_optimize' to compute it in parallel###


toy_expression_data<-matrix(data=abs(rnorm(12)),nrow=3,ncol=4, dimnames=list(paste("G",c(1:3),sep=''), paste("tp",c(0,2,4,6),sep="_")))
toy_timepoints_data<-data.frame(index=c(1:4), label=paste("tp",c(0,2,4,6),sep='_'), time=c(0,2,4,6),row.names=paste("tp",c(0,2,4,6),sep='_'))
toy_varMetadata<-data.frame(labelDescription=c("Index number","Label Detail", "Time points values"),row.names=c("index","label","time"))
toy_phenoData<-new("AnnotatedDataFrame", data=toy_timepoints_data,varMetadata=toy_varMetadata)
toy_ExpressionSet<-new("ExpressionSet", exprs=toy_expression_data,phenoData=toy_phenoData)

#########Set parameters #######
n<- 1 
sparsity<- 0.2
lbH<- -3
ubH<- 3
lbB<- 0
ubB<- 10
#########Calculate results #######

result<-SPEM(toy_ExpressionSet,n,sparsity,lbH,ubH,lbB,ubB)
}
\keyword{Bioinformatics, NetworkAnalysis, NetworkInference, Software}
