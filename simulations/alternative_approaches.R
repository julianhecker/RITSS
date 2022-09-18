library(iSKAT)
library(sandwich)
library(data.table)
library(mgcv)

gesat=function(x, e, z, y)
{
	#GESAT(Z, Y, E, X=NULL, type="davies", lower=1e-20, upper=sqrt(nrow(Y))/log(nrow(Y)), nintervals=5, 
	#plotGCV=FALSE, plotfile=NA, scale.Z=TRUE, weights.Z=NULL, weights.V=NULL, out_type="C", 
	#impute.method = "fixed", is_check_genotype=TRUE, is_dosage=FALSE, missing_cutoff=0.15, SetID=NULL)
	
	#Z genotypes, Y phenotype, E environment, X covariates
	test=GESAT(Z=x, Y=as.matrix(y), E=as.matrix(e[,1]), X=as.matrix(cbind(z, e[,-1])), type="liu")
	return(test$pvalue)
}

# direct test (no robust statistic approach)
direct_test=function(ritss_obj)
{
    K=length(ritss_obj$sub_objs)
	stat=variance=0
	for(k in 1:K)
	{
	   stat=stat+sum(ritss_obj$sub_objs[[k]]$y_resid*ritss_obj$sub_objs[[k]]$score)
	   variance=variance+sum(ritss_obj$sub_objs[[k]]$y_resid**2*ritss_obj$sub_objs[[k]]$score**2)
	}
	z=stat/sqrt(variance)
	pval_overall=pnorm(-abs(z),0,1)*2
	return(pval_overall)
}

gam_test_sv <- function(y, x, e, z, e_index) 
{
	data=cbind(e, z)
	data=as.data.frame(data)
	colnames(data)=paste0("c",1:ncol(data))
	form="y~score+xi"
	for(i in 1:ncol(data))
	{
	  form=paste0(form, paste0("+s(c",i,", bs='cr', k=3)"))
	}
	pvs=rep(0, ncol(x))
	for(i in 1:ncol(x))
	{
		xi=x[,i]
		score=xi*e[,e_index]
		gam_fit=bam(as.formula(paste(form)), data=data)
		gam_fit_summary=summary(gam_fit)
		pv=unname(gam_fit_summary$p.pv[2])
		pvs[i]=pv
		
	}
	pvs=ifelse(pvs*ncol(x)>1, 1, pvs*ncol(x))
    return(pvs)
}