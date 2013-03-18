setGeneric("row_optimize",signature="TS_eSet",function(TS_eSet,S,beta,sparsity=0.2,lbH=-3,ubH=3,lbB=0,ubB=10)standardGeneric("row_optimize"))
setMethod("row_optimize","ExpressionSet",function(TS_eSet,S,beta,sparsity=0.2,lbH=-3,ubH=3,lbB=0,ubB=10)
{
	# TS_eSet: Time series data in ExpressionSet format
	# S: The slope of one certain metabolite.
	# beta: Initial guess of beta.
	# lbH: Low boundary value of h.
	# ubH: Up boundary value of h.
	# lbB: Low boundary value of Beta.
	# ubB: Up boundary value of Beta.
		
	
	TS=t(exprs(TS_eSet))
	if(is.null(lbH)||is.null(ubH)||is.null(lbB)||is.null(ubB))# User should give a range of parameter H and beta
	{
		stop("Please input the range of your S system parameters")
	}
	if(lbB>ubB) # Range setting
	{
		stop("error: ubB must bigger than lbB!")
	}
	if(lbH>ubH) # Range setting
	{
		stop("error: ubH must bigger than lbH!")
	}
	if(is.vector(S)==0)
	{
		stop("error:S must be a vector!")
	}
	if(length(S)!=nrow(TS))
	{
		stop("error: S should have the same length as number of columns in TS")
	}
	if(sparsity<0)
	{
		stop("sparsity must be a positive value.")
	}
	if(is.null(beta)){
		stop("please input initial beta")
	}
	if(is.null(S)){
		stop("please input Slope for your calculation")
	}
	
	cat("Initial Beta=",beta,"\r\n",sep='')
	TS[which(TS==0)]=1e-7
	lbB=max(0,lbB)
	meta=c(1:ncol(TS))
	IniBeta=beta
	initial_par=initial_beta_guess(beta,TS,S,meta)
	Vt=initial_par$Vt # Global parameter 
	vt=initial_par$vt # Local parameter
	if(length(which(abs(vt[2:length(vt)])<0.1111))==(length(vt)-1)) # local optimazition.
	{
		ct=1; # The counter. 
		meta2=meta
		if(vt[1]>ubB)
			ubB=vt[1]+5
		result1=list()
		vet_error=c()
		result1[[ct]]=opt_parameter(TS,S,beta,meta,vt,lbH,ubH,lbB,ubB)
		error_struc=length(which(c(result1[[ct]]$h,result1[[ct]]$g)!=0))*0.1
		vet_error[ct]=result1[[ct]]$error+error_struc
		aux1=which(abs(result1[[ct]]$h)<sparsity)
		if(length(aux1)==0)
			aux1=meta2
		while(length(aux1)<length(meta)&length(aux1)!=0) #while the h is too small,set h to zero
		{
			ct=ct+1
			meta=meta[-aux1]
			vt=initial_beta_guess(beta,TS,S,meta)$vt
			if(vt[1]>ubB)
				ubB=vt[1]+5
			result1[[ct]]=opt_parameter(TS,S,beta,meta,vt,lbH,ubH,lbB,ubB)
			aux1=which(abs(result1[[ct]]$h[meta])<sparsity)
			error_struc=length(which(c(result1[[ct]]$h,result1[[ct]]$g)!=0))*0.1
			vet_error[ct]=result1[[ct]]$error+error_struc
		}
		bb=which(vet_error==min(vet_error))
		aa=vet_error[ct]
		resultA=result1[[bb]]
		resultA$error=aa
		meta=meta2
		if(length(which(Vt[2:length(Vt)]>lbH))==(length(Vt)-1)&length(which(Vt[2:length(Vt)]<ubH))==(length(Vt)-1)) #global optimazition.
		{
			result2=list()
			vet_error=c()
			ct=1
			if(Vt[1]>ubB)
				ubB=Vt[1]+5
			result2[[ct]]=opt_parameter(TS,S,beta,meta,Vt,lbH,ubH,lbB,ubB)
			error_struc=length(which(c(result2[[ct]]$h,result2[[ct]]$g)!=0))*0.1
			vet_error[ct]=result2[[ct]]$error+error_struc
			aux1=which(abs(result2[[ct]]$h)<sparsity)
			if(length(aux1)==0)
				aux1=meta
			while(length(aux1)<length(meta)&length(aux1)!=0)
			{
				ct=ct+1
				meta=meta[-aux1]
				temp=initial_beta_guess(beta,TS,S,meta)
				Vt=temp$Vt
				if(Vt[1]>ubB)
					ubB=Vt[1]+5
				result2[[ct]]=opt_parameter(TS,S,beta,meta,Vt,lbH,ubH,lbB,ubB)
				aux1=which(abs(result2[[ct]]$h[meta])<sparsity)
				error_struc=length(which(c(result2[[ct]]$h,result2[[ct]]$g)!=0))*0.1
				vet_error[ct]=result2[[ct]]$error+error_struc
			}
			bb=which(vet_error==min(vet_error))
			aa=vet_error[ct]
			resultB=result2[[bb[1]]]
			resultB$error=aa
			meta=meta2
		}
		else
		{
			resultB=resultA
			resultB$error=Inf
		}
		if(resultA$error<resultB$error)
		{
			result=resultA
		}
		else
		{
			result=resultB
		}
	}
else
	{
		if(length(which(Vt[2:length(Vt)]>lbH))==(length(Vt)-1)&length(which(Vt[2:length(Vt)]<ubH))==(length(Vt)-1))# global optimization.
		{
			result2=list()
			vet_error=c()
			ct=1
			if(Vt[1]>ubB)
				ubB=Vt[1]+5
			result2[[ct]]=opt_parameter(TS,S,beta,meta,Vt,lbH,ubH,lbB,ubB)
			error_struc=length(which(c(result2[[ct]]$h,result2[[ct]]$g)!=0))*0.1
			vet_error[ct]=result2[[ct]]$error+error_struc
			aux1=which(abs(result2[[ct]]$h)<sparsity)
			if(length(aux1)==0)
				aux1=meta
			while(length(aux1)<length(meta)&length(aux1)!=0)
			{
				ct=ct+1
				meta=meta[-aux1]
				temp=initial_beta_guess(beta,TS,S,meta)
				Vt=temp$Vt
				if(Vt[1]>ubB)
					ubB=Vt[1]+5	
				result2[[ct]]=opt_parameter(TS,S,beta,meta,vt,lbH,ubH,lbB,ubB)
				aux1=which(abs(result2[[ct]]$h[meta2]<sparsity))
				error_struc=length(which(c(result2[[ct]]$h,result2[[ct]]$g)!=0))*0.1
				vet_error[ct]=result2[[ct]]$error+error_struc
			}
			bb=which(vet_error==min(vet_error))
			aa=vet_error[ct]
			result=result2[[bb]]
			result$error=aa
		}
		else
		{
			while(length(which(abs(vt[2:length(vt)])<0.1111))<(length(vt)-1))
			{
				beta=vt[1]+1
				temp=initial_beta_guess(beta,TS,S,meta)
				vt=temp$vt
			}
			ct=1
			meta=meta2
			if(vt[1]>ubB)
				ubB=vt[1]+5
			result1=list()
			vet_error=c()
			result1[[ct]]=opt_parameter(TS,S,beta,meta2,vt,lbH,ubH,lbB,ubB)
			error_struc=length(which(c(result1[[ct]]$h,result1[[ct]]$g)!=0))*0.1
			vet_error[ct]=result1[[ct]]$error+error_struc
			aux1=which(abs(result1[[ct]]$h)<sparsity)
			if(length(aux1)==0)
				aux1=meta2
			while(length(aux1)<length(meta2)&length(aux1)!=0)
			{
				ct=ct+1
				meta2=meta2[-aux1]
				temp=initial_beta_guess(beta,TS,S,meta2)
				vt=temp$vt
				if(vt[1]>ubB)
					ubB=vt[1]+5
				result1[[ct]]=opt_parameter(TS,S,beta,meta2,vt,lbH,ubH,lbB,ubB)
				aux1=which(abs(result1[[ct]]$h[meta2]<sparsity))
				error_struc=length(which(c(result1[[ct]]$h,result1[[ct]]$g)!=0))*0.1
				vet_error[ct]=result1[[ct]]$error+error_struc
			}
			bb=which(vet_error==min(vet_error))
			aa=vet_error[ct]
			result=result1[[bb]]
			result$error=aa
		}
	}
	##Set bound##
	alpha=result$alpha
	beta=result$beta
	g=result$g
	g[which(g>ubH)]=0#filter g that is out of bound
	g[which(g<lbH)]=0
	h=result$h
	##Compare very close g and h, we choose depend on alpha and beta###
	for(bound in 1:length(g)){
		if(round(g[bound],0)!=0&&round(h[bound],0)!=0){
			if(abs(round(g[bound],0))<abs(round(g[bound],0))&&((alpha<beta)[bound])){
				g[bound]=0
			}
			else{
				h[bound]=0
			}
		}
		else{
			if(abs(g[bound])<abs(h[bound])){
				g[bound]=0
			}
			else{
				h[bound]=0
			}
		}
	}
	result$g=g
	result$h=h
	return(c(alpha=result$alpha,g=result$g,beta=result$beta,h=result$h,InitialBeta=IniBeta,error=result$error))
}
)
