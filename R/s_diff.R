setGeneric("s_diff",signature="TS_eSet",function(TS_eSet)standardGeneric("s_diff"))
setMethod("s_diff","ExpressionSet",function(TS_eSet)
{   	
	TS=t(exprs(TS_eSet))
	tp=as.numeric(as.vector(pData(TS_eSet)$time))
	if(length(tp)!=nrow(TS)) # Time points setting
	{
		stop("error: time points should have the same length as number of columns in TS")
	}
	
	tp_2=seq(tp[1],tp[length(tp)],len=1000*length(tp))
	tp_2=c((tp_2[1]-diff(tp_2)[1]),tp_2,(tp_2[length(tp_2)]+diff(tp_2)[(length(tp_2)-1)]))
	xx=sort(unique(c(tp,tp_2)))
	pre_smooth=function(r) predict(smooth.spline(tp,r,spar=0),xx)$y
	dense=apply(TS,MARGIN=2,FUN=pre_smooth)
	cal=match(tp,xx)
	d=diff(xx)
	d_1=d[1:(length(d)-1)]
	d_2=d[2:length(d)]
	Slope=matrix(nrow=nrow(TS),ncol=ncol(TS))
	Slope[1:length(tp),]=(-d_2[cal-1]/d_1[cal-1]*dense[cal-1,]+d_1[cal-1]/d_2[cal-1]*dense[cal+1,])/(d_1[cal-1]+d_2[cal-1])+(1/d_1[cal-1]-1/d_2[cal-1])*dense[cal,]
	rownames(Slope)=rownames(TS)
	colnames(Slope)=colnames(TS)
	return(t(Slope))
}
)