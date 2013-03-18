setGeneric("SPEM",signature="TS_eSet",function(TS_eSet,n=3,sparsity=0.2,lbH=-3,ubH=3,lbB=0,ubB=10)standardGeneric("SPEM"))
setMethod("SPEM","ExpressionSet",function(TS_eSet,n=3,sparsity=0.2,lbH=-3,ubH=3,lbB=0,ubB=10)
{
	# TS_eSet: Time series data in ExpressionSet Class
	# n: Times of initial beta guess.
	# sparsity: A threshold used to control the sparsity of reconstructed matrix. Values whose absolute value smaller than sparsity will be set to zero.
	# lbH: Low boundary value of h.
	# ubH: Up boundary value of h.
	# lbB: Low boundary value of Beta.
	# ubB: Up boundary value of Beta. 
	TS=t(exprs(TS_eSet))
	tp=as.numeric(as.vector(pData(TS_eSet)$time))
# Error checking.	
	if(is.null(TS)) # User must input time series expression data.
	{
		stop("Please input time series expression data as TS.")
	}
	if(is.null(tp))
	{
		stop("Please input time points as tp")
	}
	if(is.null(lbH)||is.null(ubH)||is.null(lbB)||is.null(ubB))# User should give a range of parameter H and beta
	{
		stop("Please input the range of your S system parameters")
	}
	if(!is.matrix(TS)) # TS must be a matrix
	{
		stop("error: TS should be a matrix!")
	}
	if(lbB>ubB) # Range setting
	{
		stop("error: ubB must bigger than lbB!")
	}
	if(lbH>ubH) # Range setting
	{
		stop("error: ubH must bigger than lbH!")
	}
	if(length(which(diff(tp)<=0))>0) # Time points setting
	{
		stop("error: time points must be a increasing vector!")
	}
	if(length(tp)!=nrow(TS)) # Time points setting
	{
		stop("error: time points should have the same length as number of columns in TS")
	}
	#if(n<0||is.integer(n)==FALSE)
	#{
	#	stop("n must be a positive integer.")
	#}
	if(sparsity<0)
	{
		stop("sparsity must be a positive value.")
	}

	cat(ncol(TS)," genes,",nrow(TS)," samples,",n," guesses of Initial Beta for each gene.\r\n",sep='')
	TS[which(TS==0)]=1e-7
	lbB=max(0,lbB)
	
	Slope=t(s_diff(TS_eSet))
	
	
	answer=apply(X=Slope,MARGIN=2,FUN=pre_prepare,TS_eSet=TS_eSet,n=n,sparsity=sparsity,lbH=lbH,ubH=ubH,lbB=lbB,ubB=ubB)
	
	answer=t(answer)
	alpha=answer[,1]
	g=answer[,2:(ncol(TS)+1)]
	colnames(g)=paste("g",as.character(1:ncol(TS)),sep='')
	rownames(g)=paste("g",as.character(1:ncol(TS)),sep='')
	beta=answer[,(ncol(TS)+2)]
	h=answer[,(ncol(TS)+3):(2*ncol(TS)+2)]
	colnames(h)=paste("h",as.character(1:ncol(TS)),sep='')
	rownames(h)=paste("h",as.character(1:ncol(TS)),sep='')
	IniBeta=answer[,(ncol(answer)-1)]
	error=answer[,ncol(answer)]
	
	#colnames(answer)=c("alpha",paste("g",as.character(1:ncol(TS)),sep=''),"beta",paste("h",as.character(1:ncol(TS)),sep=''),"error")
	#rownames(answer)=paste("Fuction",as.character(1:ncol(TS)),sep='')
	return(list(alpha=alpha,g=g,beta=beta,h=h,IniBeta=IniBeta,error=error))
})

