ssystem_set <-
function(alpha,g,beta,h,TS)
#alpha,g,beta,h: The parameter of the S-system. See S-system function to get more details.
#TS: Time series data.
{

	G=matrix(data=rep(g,each=nrow(TS)),nrow=nrow(TS),ncol=length(g))#preparing g
	H=matrix(data=rep(h,each=nrow(TS)),nrow=nrow(TS),ncol=length(h))#preparing h

	ss=(alpha*apply(TS^G,1,FUN=prod)-beta*apply(TS^H,1,FUN=prod))#Calculating
	
	return(ss)
}

