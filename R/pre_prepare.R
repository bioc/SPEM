setGeneric("pre_prepare",signature="TS_eSet",function(Slope,TS_eSet,n,sparsity,lbH,ubH,lbB,ubB)standardGeneric("pre_prepare"))
setMethod("pre_prepare","ExpressionSet",function(Slope,TS_eSet,n,sparsity,lbH,ubH,lbB,ubB)
{
	# TS_Set: Time series data, with n time points in row and m metabolite in column in ExpressionSet class
	# Slope: Slope file from the TS	
	# lbH: Low boundary value of h.
	# ubH: Up boundary value of h.
	# lbB: Low boundary value of Beta.
	# ubB: Up boundary value of Beta. 
	
	beta=runif(n,min=lbB,max=ubB)
	S=Slope
	cat(">>>calculating a certain gene.>>>\r\n")
	TS_eSet=TS_eSet
	op_candidate=sapply(X=beta,FUN=row_optimize,TS_eSet=TS_eSet,S=S,sparsity=sparsity,lbH=lbH,ubH=ubH,lbB=lbB,ubB=ubB)
	return(op_candidate[,order(op_candidate[nrow(op_candidate),])[1]])
})

