initial_h_guess <-
function(TS,S,beta,meta)
#TS: Time series data.
#S: The slope of one certain metabolite.
#beta: The beta value for initial guess of h.
#meta: The serial number of metabolites included in the calculation. e.g. meta=c(1,3) means only 1st and 3rd metabolites affect this function.
{
	thres=0.09 # The threshold for initial h. In this version. This threshold is a const.

	L=log(TS[which(S<0),meta]) #select the time points that the metabolite is down regulated.

# Global guess of h.
# To see whether the metabolite is mainly up regulated

	if(length(which(S<0))==0) # The betabolite keeps being up regulated.  	
	{
		L=log(TS[,meta])
		Vt=solve(t(L)%*%L,tol=1e-18)%*%t(L)%*%log(S+0.1/beta)
	}
	else if(length(which(S<0))==length(meta)) # The number down regulated time points is equal to the number of metabolites we concern.
	{
		L=log(TS[which(S<0),meta])
		Vt=solve(L,tol=1e-18)%*%log(2-S[which(S<0)]/beta)
	}
	else if(length(which(S<0))<length(meta))
	{
		Vt=c(rep(0.01,length(meta)))
	}
	else
	{
		L=log(TS[which(S<0),meta])
		Vt=solve(t(L)%*%L,tol=1e-18)%*%t(L)%*%log(2-S[which(S<0)]/beta)
	}
	
#Local guess of h.
#To calculate whether the down regulation is too steep. 

	LL=matrix(data=0,nrow=nrow(TS),ncol=1+length(meta))
	LL[,1]=1
	LL[,2:(1+length(meta))]=log(TS[,meta])
	if(length(which((S+exp(beta))<0))==0) #If there is no time point that is too steep.
	{
		vt=solve(t(LL)%*%LL,tol=1e-18)%*%t(LL)%*%log(S+exp(beta))
	}
	else # If there are still some steep points. 
	{
		vt=solve(t(LL)%*%LL,tol=1e-18)%*%t(LL)%*%log(S+exp(beta))
		con=which(abs(vt[2:length(vt)])<thres) # To see whether vt is under the requirement
		flag=0
		while(((length(which((S+exp(beta))<0))!=0)|(!vt[1]>0)|(length(con)==0))&&flag<20)# Set beta until there is no steep time point.
		{
			beta=beta+1
			vt=solve(t(LL)%*%LL,tol=1e-18)%*%t(LL)%*%log(S+exp(beta))
			con=which(abs(vt[2:length(vt)])<thres)
			flag=flag+1
		}
	}

	con=which(abs(vt[2:length(vt)])<thres)

	flag=0
	while(((!vt[1]>0)|length(con)!=length(meta))&&flag<20)# Let all the value in h get under the theshold
	{
		beta=beta+1
		vt=solve(t(LL)%*%LL,tol=1e-18)%*%t(LL)%*%log(S+exp(beta))
		con=which(abs(vt[2:length(vt)])<thres)
		flag=flag+1
	}

	Vt=c(beta,Vt)
	vt=as.vector(vt)
	return(list(Vt=Vt,vt=vt))
}

