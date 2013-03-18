initial_beta_guess <-
function(beta,TS,S,meta)
#TS: Time series data.
#S: The slope of one certain metabolite.
#beta: The beta value for initial guess of h.
#meta: The serial number of metabolites included in the calculation. e.g. meta=c(1,3) means only 1st and 3rd metabolites affect this function.
{
	initial_h=initial_h_guess(TS,S,beta,meta) #get the initial value of h.
	Vt=initial_h$Vt #global initial value of h.
	vt=initial_h$vt #local initial value of h.
	int=rep(0,ncol(TS))
	int[meta]=Vt[-1]#give the h to calculate the system value.
	global_s=S+ssystem_set(0,rep(0,ncol(TS)),-Vt[1],int,TS)
	con_g=which(!global_s>0)
	local_s=S+ssystem_set(0,rep(0,ncol(TS)),-beta,int,TS)
	con_l=which(!local_s>0)
	flag=0
	while(length(con_l)!=0&&flag<20) # local initial guess for beta.
	{
		beta=beta+1.1
		initial_h=initial_h_guess(TS,S,beta,meta)
		Vt=initial_h$Vt
		vt=initial_h$vt
		local_s=S+ssystem_set(0,rep(0,ncol(TS)),-beta,int,TS)
		con_l=which(!local_s>0)
		flag=flag+1
	}
	flag=0
	while(length(con_g)!=0&&flag<20) # local initial guess for beta.
	{
		beta=beta+1.1
		initial_h=initial_h_guess(TS,S,beta,meta)
		Vt=initial_h$Vt 
		vt=initial_h$vt 
		global_s=S+ssystem_set(0,rep(0,ncol(TS)),-beta,int,TS)
	    con_g=which(!global_s>0)
	    flag=flag+1
	}
	return(list(Vt=Vt,vt=vt))
}

