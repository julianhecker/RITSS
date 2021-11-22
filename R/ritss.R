#' RITSS
#'
#' @param y phenotype values.
#' @param x genotype data (can be expected counts).
#' @param e environmental factors.
#' @param z additional covariates. For example, genetic principal components.
#' @param ii1 indices for I_1.
#' @param ii2 indices for I_2.
#' @param ii3 indices for I_3. ii1, ii2, and ii3 should NOT overlap.
#' @param cut_off_p_value cut off p-value for the inclusion of scores. The first score is always tested.
#' @param screening_function function that performs the screening step.
#' @param indices_env_factors indices of environmental factors that are considered in the screening step.
#' @param verbose default is FALSE. If TRUE, RITSS outputs more details.
#'
#' @return list of objects including RITSS p-value, corresponding z-score. Also individual z-scores of the three substatistics, the three interaction scores, 
#' and information about the selected variants.
#' @export
#'
#' @examples
#' n=15000
#' m=50; d=5; p=5;
#' y=rnorm(n); x=matrix(rnorm(n*m), nrow=n, ncol=m); 
#' e=matrix(rnorm(n*d), nrow=n, ncol=d); z=matrix(rnorm(n*p), nrow=n, ncol=p);
#' ii1=1:5000; ii2=5001:10000; ii3=10001:15000;
#' ritss(y, x, e, z, ii1, ii2, ii3)

ritss=function(y, x, e, z, ii1, ii2, ii3, cut_off_p_value=0.05, screening_function=screening_subprs, indices_env_factors=1, verbose=FALSE)
{
    y_1=y[ii1]; y_2=y[ii2]; y_3=y[ii3];
	x_1=as.matrix(x[ii1,]); x_2=as.matrix(x[ii2,]); x_3=as.matrix(x[ii3,]);
	e_1=as.matrix(e[ii1,]); e_2=as.matrix(e[ii2,]); e_3=as.matrix(e[ii3,]);
	z_1=as.matrix(z[ii1,]); z_2=as.matrix(z[ii2,]); z_3=as.matrix(z[ii3,]);
	
	obj1=ritss_sub(x_1, e_1, z_1, y_1, x_2, e_2, z_2, y_2, x_3, e_3, z_3, y_3, cut_off_p_value, screening_function, indices_env_factors, verbose)
	obj2=ritss_sub(x_3, e_3, z_3, y_3, x_1, e_1, z_1, y_1, x_2, e_2, z_2, y_2, cut_off_p_value, screening_function, indices_env_factors, verbose)
	obj3=ritss_sub(x_2, e_2, z_2, y_2, x_3, e_3, z_3, y_3, x_1, e_1, z_1, y_1, cut_off_p_value, screening_function, indices_env_factors, verbose)

	stat=obj1$stat+obj2$stat+obj3$stat
	variance=obj1$var_est+obj2$var_est+obj3$var_est
	zsc=stat/sqrt(variance)
    pval_overall=pnorm(-abs(zsc),0,1)*2
	
	return(list("pval"=pval_overall, "z"=zsc, "z_1"=obj1$stat/sqrt(obj1$var_est), "z_2"=obj2$stat/sqrt(obj2$var_est), "z_3"=obj3$stat/sqrt(obj3$var_est), 
	"Ui_1"=obj3$Ui, "Ui_2"=obj2$Ui, "Ui_3"=obj1$Ui, "snps_in_score_i1"=obj3$snps_in_score, "snps_in_score_i2"=obj2$snps_in_score, "snps_in_score_i3"=obj1$snps_in_score))
	
}
# RITSS function for fixed i_1, i_2, and i_3
ritss_sub=function(x_1, e_1, z_1, y_1, x_2, e_2, z_2, y_2, x_3, e_3, z_3, y_3, cut_off_p_value=0.05, screening_function, indices_env_factors, verbose=FALSE)
{
  ##############################################################
  ### part I
  screening_obj=screening_function(y_1, x_1, e_1, z_1, indices_env_factors, verbose)
  if(verbose==TRUE) print("screening step done")
  number_scores=length(screening_obj$mapping_env_factors)
  score_env_mapping=screening_obj$mapping_env_factors
  included_env_factors=unique(score_env_mapping)
  
  ##############################################################
  ### part II

  n_k=floor(length(y_2)/2)
  sub_sample_1=1:n_k
  sub_sample_2=(n_k+1):(length(y_2))
  
 
  scores_i2=matrix(0, nrow=length(y_2), ncol=number_scores)
  scores_i3=matrix(0, nrow=length(y_3), ncol=number_scores)
  scores_x_i2=matrix(0, nrow=length(y_2), ncol=number_scores)
  scores_x_i3=matrix(0, nrow=length(y_3), ncol=number_scores)
  for(k in 1:number_scores)
  {
     betav=screening_obj$betas[[k]]
	 env_index=score_env_mapping[k]
	 scores_i2[,k]=e_2[,env_index] * (as.matrix(x_2) %*% betav)
	 scores_i3[,k]=e_3[,env_index] * (as.matrix(x_3) %*% betav)
	 scores_x_i2[,k]=as.matrix(x_2) %*% betav 
	 scores_x_i3[,k]=as.matrix(x_3) %*% betav 
  }
 
  lasso_obj=lasso_main_effects(as.matrix(scores_i2[sub_sample_1,]), x_2[sub_sample_1,], as.matrix(e_2[sub_sample_1,]), z_2[sub_sample_1,], y_2[sub_sample_1], x_3, e_3, z_3)
  
  pval_score=as.numeric(lasso_obj$pval_scores)
  beta_score=as.numeric(lasso_obj$beta_scores)
  if(verbose==TRUE) cat("I_2 filter score betas:", beta_score,"\n")
  if(verbose==TRUE) cat("I_2 filter score p-values:",pval_score,"\n")
  
  
  keep_scores=rep(0, number_scores); 
  keep_scores[pval_score<=cut_off_p_value]=1; keep_scores[1]=1; # always keep first score
  beta_score[keep_scores==0]=0
  beta_score[keep_scores==1]=1 
   
  
  
  ##############################################################
  ### part III
  
  variants_included=rep(0, ncol(x_3))
  
  ### interaction score computation
  for(k in 1:number_scores)
  {
     betav=screening_obj$betas[[k]]
	 if(keep_scores[k]==1) {variants_included[betav!=0]=1}
  }
  tmp_inds=1:ncol(x_3)
  snps_in_score=tmp_inds[variants_included==1]
  
  if(verbose==TRUE){cat("snps in score:", length(snps_in_score), "\n")}
  ##############################################################
  ### projections
  Ui_i3=scores_i3 %*% beta_score
  Ui_i3_prime=rep(0, nrow(x_3))
  score_tmp_indices=1:number_scores
  for(index in included_env_factors)
  {
     inds=score_tmp_indices[score_env_mapping==index]
	 if(verbose==TRUE){cat("env.factor: ", index,  inds, beta_score[inds], "\n")}
	 tmp_Ui_i2=as.matrix(scores_i2[,inds]) %*% as.vector(beta_score[inds])
	 tmp_Ui_i3=as.matrix(scores_i3[,inds]) %*% as.vector(beta_score[inds])
	 tmp_score_x_i2=as.matrix(scores_x_i2[,inds]) %*% as.vector(beta_score[inds])
     tmp_score_x_i3=as.matrix(scores_x_i3[,inds]) %*% as.vector(beta_score[inds])
	 trained=ACE_train(x_2[sub_sample_2,], as.matrix(e_2[sub_sample_2,]), z_2[sub_sample_2,], tmp_Ui_i2[sub_sample_2], index, tmp_score_x_i2[sub_sample_2])
     Ui_i3_prime=Ui_i3_prime+ACE_pred(x_3, e_3, z_3, tmp_Ui_i3, trained$rg, trained$rf, trained$rg2,  trained$rg2_indicators, trained$rf2, tmp_score_x_i3)
	 if(verbose==TRUE) cat("number of selected genetic variants in prediction of e: ", length(trained$x_inds), "\n")
  }
 
 
  y_resid=y_3-lasso_obj$predicted_y_main_effects
  
  ### test statistic computation
  S=sum(Ui_i3_prime*y_resid)
  var_est=sum(Ui_i3_prime^2*y_resid^2)
  if(verbose==TRUE){cat("S , sd(S): ", S, sqrt(var_est), "\n")}
  
  return(list("stat"=S, "var_est"=var_est, "betas"=screening_obj$betas, "snps_in_score"=snps_in_score, "Ui"=Ui_i3, "y_resid"=y_resid))
 
}
# screening step to identify interaction with subcomponent of genetic risk score
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
