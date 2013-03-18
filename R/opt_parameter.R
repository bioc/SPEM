opt_parameter <-
function(TS,S,beta,meta,vt,lbH,ubH,lbB,ubB)
# TS: Time series data.
# S: The slope of one certain metabolite.
# vt: Initial guess of beta and h.
# lbH: Low boundary value of h.
# ubH: Up boundary value of h.
# lbB: Low boundary value of Beta.
# ubB: Up boundary value of Beta.
{
	h=c()
	g=c()
	
# Matrix preparing. See algorithm to get more detials.
	L=matrix(data=0,nrow=nrow(TS),ncol=ncol(TS)+1)
	L[,1]=1
	L[,2:(ncol(TS)+1)]=log(TS)
	C=solve(t(L)%*%L,tol=1e-18)%*%t(L)
	I=diag(nrow(TS))
	W1=L%*%C
	W=W1-I 
	
# Objective function.
	`log_betafunction_attrac_diff`<-function(x)
	{
		h=c()
		h[1:ncol(TS)]=0
		h[meta]=x[2:(length(meta)+1)]
		con=S+ssystem_set(0,rep(0,ncol(TS)),-(x[1]),h,TS)
		con[which(!(con>0))]=1e-10
		yd=log(con)
		bd=C%*%yd
		g=bd[2:length(bd)]
		alpha=exp(bd[1])
		res=ssystem_set(alpha,g,x[1],h,TS)
		t(res-S)%*%(res-S) # Return the error as object.
	}

# Constraint condition.
	nlin=function(x)
	{
		h=c()
		h[1:ncol(TS)]=0
		C=ssystem_set(0,rep(0,ncol(TS)),x[1],h,TS)
		return(C)
	}

	pars=c(beta,vt[2:length(vt)]) # Initial guess.
	UB=c(ubB,rep(ubH,length(meta))) # Range of value.
	LB=c(lbB,rep(lbH,length(meta))) # As above.

	ineqLB = rep(-Inf,length(S)) # Non-linear constrance.
	ineqUB = S # As above. 
	
	control = list(tol=1e-15,delta=1e-10)
	
	ret =try(suppressWarnings(solnp(pars=pars, fun=log_betafunction_attrac_diff, eqfun = NULL, eqB = NULL, ineqfun = nlin, ineqLB = ineqLB, ineqUB = ineqUB, UB = UB, LB = LB, control=control)),silent=TRUE)# Non-linear optimization function. We use the package Rsolnp to do it.

	if(class(ret)=="list")
	{
		beta=ret$par[1]
		h=rep(0,ncol(TS))
		h[meta]=ret$par[2:(length(meta)+1)]
		con=S+ssystem_set(0,rep(0,ncol(TS)),-beta,h,TS)
		yd=log(con)
		bd=C%*%yd #calculate alpha and g
		g=bd[2:length(bd)]
		g=as.vector(g)
		g[which(abs(g)<0.1)]=0
		alpha=exp(bd[1])
		res=ssystem_set(alpha,g,beta,h,TS)
		error=t(res-S)%*%(res-S)
	}
	else
	{
		alpha=NA
		beta=NA
		h[1:ncol(TS)]=NA
		g[1:ncol(TS)]=NA
		error=Inf
	}	
	if(is.na(error)==TRUE)
		error=Inf
	return(list(alpha=alpha,beta=beta,g=g,h=h,error=error))
}

