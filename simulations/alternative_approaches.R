library(iSKAT)
library(sandwich)
source("../R/bestsubsetregression.R")


gesat=function(x, e, z, y)
{
	test=GESAT(x, as.matrix(y), as.matrix(e[,1]), as.matrix(cbind(z, e[,-1])))
	return(test$pvalue)
}
linreg_inter_sv=function(x, e, z, y, inter_index)
{
	e_s=e[,inter_index]
	pval=rep(0, ncol(x))
	for(i in 1:ncol(x))
	{
		xe=x[,i]*e_s
		fit=lm(y~xe+x[,i]+e+z)
		sandwich_se <- diag(vcovHC(fit, type = "HC"))^0.5
		tmp=summary(fit)$coefficients[2,1]/sandwich_se[2]
		pval[i]=pnorm(-abs(tmp), 0, 1)*2
	}
	return(min(ifelse(pval*ncol(x)>=1,1,pval*ncol(x))))
}


linreg_2split=function(y, x, e, z, screening_function=screening_subprs, indices_env_factors=1)
{
  n_k=floor(length(y)/2)
  sub_sample_1=1:n_k
  sub_sample_2=(n_k+1):(length(y))
  y_1=y[sub_sample_1]; x_1=x[sub_sample_1,]; e_1=e[sub_sample_1,]; z_1=z[sub_sample_1,];
  y_2=y[sub_sample_2]; x_2=x[sub_sample_2,]; e_2=e[sub_sample_2,]; z_2=z[sub_sample_2,];
  ##############################################################
  
  screening_obj1=screening_function(y_1, x_1, e_1, z_1, indices_env_factors, FALSE)
  screening_obj2=screening_function(y_2, x_2, e_2, z_2, indices_env_factors, FALSE)
  
  number_scores1=length(screening_obj1$mapping_env_factors)
  score_env_mapping1=screening_obj1$mapping_env_factors
  
  number_scores2=length(screening_obj2$mapping_env_factors)
  score_env_mapping2=screening_obj2$mapping_env_factors
  
 
  scores_1=matrix(0, nrow=length(y_1), ncol=number_scores2)
  scores_2=matrix(0, nrow=length(y_2), ncol=number_scores1)
  
  for(k in 1:number_scores1)
  {
     betav=screening_obj1$betas[[k]]
	 env_index=score_env_mapping1[k]
	 scores_2[,k]=e_2[,env_index] * (as.matrix(x_2) %*% betav)
  }
  for(k in 1:number_scores2)
  {
     betav=screening_obj2$betas[[k]]
	 env_index=score_env_mapping2[k]
	 scores_1[,k]=e_1[,env_index] * (as.matrix(x_1) %*% betav)
  }
  score1=rowSums(scores_1); score2=rowSums(scores_2);
  
  fit1=lm(y_1~score1+x_1+e_1+z_1)
  fit2=lm(y_2~score2+x_2+e_2+z_2)
  
  sandwich_se1 <- diag(vcovHC(fit1, type = "HC"))^0.5
  sandwich_se2 <- diag(vcovHC(fit2, type = "HC"))^0.5
  
  tmp1=summary(fit1)$coefficients[2,1]/sandwich_se1[2]
  tmp2=summary(fit2)$coefficients[2,1]/sandwich_se2[2]
  
  pval_reg_score_robust1=pnorm(-abs(tmp1), 0, 1)*2
  pval_reg_score_robust2=pnorm(-abs(tmp2), 0, 1)*2
  
  pval=min(min(summary(fit1)$coefficients[2,4], summary(fit2)$coefficients[2,4])*2, 1)
  pval_robust=min(min(pval_reg_score_robust1, pval_reg_score_robust2)*2, 1)
  
  return(list("pval"=pval, "pval_robust"=pval_robust))
 
}

### took from RITSS R package code on 5/4/2022

screening_subprs=function(y, x, e, z, indices_env_factors, num_iterations=10, verbose=FALSE)
{
   n=nrow(x)
   if(length(indices_env_factors)!=1) return(NA);
   inter_index=indices_env_factors;
   if(inter_index<1 | inter_index > ncol(e)) return(NA);
   
   ############################
   # estimation of betas_x 
   
   tmp_data_mat=cbind(x, e, z)
   fit=lm(y~tmp_data_mat)
   betas_x=summary(fit)$coefficients[2:(ncol(x)+1),1]
   xw=x
   for(i in 1:ncol(x))
   {
     xw[,i]=xw[,i]*betas_x[i]
   }
   xwe=xw*e[,inter_index]
   if(verbose==TRUE) print("estimating beta-x done")
   ############################
   # regressing out x, e, and z
   
   yr=fit$residuals
   xe_r=xwe
   for(i in 1:ncol(xe_r))
   {
      fit=lm(xe_r[,i]~x[,i]+e+z)
	  xe_r[,i]=fit$residuals
   }
   if(verbose==TRUE) print("regressing out done")
   ############################
   indices_snps=1:ncol(xe_r)
   vars=apply(xe_r, 2, var)
   total_vars=sum(vars)
   vars=vars/total_vars
   tmp=as.data.frame(cbind(indices_snps, vars))
  
   
   ind_var=tmp$indices_snps[tmp$vars>1/nrow(tmp)] # keep only variants with reasonable variance
   if(length(ind_var)<10){ ind_var=tmp$indices_snps; }###-- update 4/11/2022
   best_subset_reg=best_subset_reg(yr, xe_r[,ind_var], length(ind_var), num_iterations=10, verbose)
   
   num_indices=1
   if(length(best_subset_reg$other)>1){num_indices=2}
   indices=matrix(0, nrow=num_indices, ncol=ncol(x))
   indices[1, ind_var[best_subset_reg$inter]]=1
   if(length(best_subset_reg$other)>1){indices[2, ind_var[best_subset_reg$other]]=1}
   if(verbose==TRUE) print("screening done")
   
   beta_1=rep(0, ncol(x));
   beta_2=rep(0, ncol(x));
   beta_1[indices[1,]==1]=betas_x[indices[1,]==1];
   if(num_indices==2) beta_2[indices[2,]==1]=betas_x[indices[2,]==1]
   mapping=rep(0, num_indices)
   mapping[1]=inter_index;
   if(num_indices==2) mapping[2]=inter_index;
   
   return(list("betas"=list(beta_1, beta_2), "mapping_env_factors"=mapping))
  
}
# best subset regression
best_subset_reg=function(yr, xe_r, S_max, num_iterations=10, verbose=FALSE)
{
   indices_snps=1:ncol(xe_r)
   seq_s=seq(10, S_max, by=5)
   folds=cut(seq(1,nrow(xe_r)),breaks=2,labels=FALSE)
   
   ind_regist1=matrix(0, nrow=length(seq_s),ncol= ncol(xe_r)+1)
   ind_regist2=matrix(0, nrow=length(seq_s),ncol= ncol(xe_r)+1)
   
   for(i in 1:2)
   {
		testIndexes <- which(folds==i,arr.ind=TRUE)
		test_x <- xe_r[testIndexes, ]; test_y= yr[testIndexes]; 
		train_x <- xe_r[-testIndexes, ]; train_y= yr[-testIndexes];
		
		L=get_L(train_y, train_x)
	    start_beta=rep(0, ncol(xe_r))
		ctr=0
	    for(s in seq_s)
	    {
			beta_tmp=gradient_algorithm(train_y, train_x, start_beta, L, s, num_iterations);
			tmp_inds=indices_snps[beta_tmp!=0]
			
			score=rowSums(test_x[, tmp_inds])
			fit=lm(test_y~score)
			z_sq=summary(fit)$coefficients[2,3]**2
			
			score_train=rowSums(train_x[, tmp_inds])
			fit=lm(train_y~score_train)
			z_sq_train=summary(fit)$coefficients[2,3]**2
			ind=rep(0, ncol(xe_r))
			ind[tmp_inds]=1;
			
			ctr=ctr+1
			if(i==1){ ind_regist1[ctr, tmp_inds]=1; ind_regist1[ctr, ncol(ind_regist1)]=z_sq;}
			if(i==2){ ind_regist2[ctr, tmp_inds]=1; ind_regist2[ctr, ncol(ind_regist1)]=z_sq;}
			
	    }
   }
   
   t1=max(ind_regist1[, ncol(ind_regist1)])
   t2=max(ind_regist2[, ncol(ind_regist2)])
   ind1=indices_snps[ind_regist1[ind_regist1[,ncol(ind_regist1)]==t1,1:ncol(xe_r)]==1]
   ind2=indices_snps[ind_regist2[ind_regist2[,ncol(ind_regist2)]==t2,1:ncol(xe_r)]==1]

   inter_ind=intersect(ind1, ind2)
   other_ind=unique(c(ind1, ind2))
   other_ind=other_ind[!(other_ind %in% inter_ind)]
   if(length(inter_ind)<5 | length(other_ind)<5){inter_ind=c(inter_ind, other_ind); other_ind=numeric(0)}
   
   return(list("inter"=inter_ind, "other"=other_ind))
}



