#' RITSS
#'
#' @param y phenotype values.
#' @param x genotype data (can be expected counts).
#' @param e environmental factors.
#' @param z additional covariates. For example, genetic principal components.
#' @param indices partitioning of the n samples into K non-overlapping sub samples
#' @param cut_off_p_value cut off p-value for the inclusion of scores. The first score is always tested.
#' @param screening_function function that performs the screening step.
#' @param indices_env_factors indices of environmental factors that are considered in the screening step.
#' @param verbose default is FALSE. If TRUE, RITSS outputs more details.
#'
#' @return list of objects including RITSS p-value and corresponding z-score. 
#' @export
#'
#' @examples
#' n=15000
#' m=50; d=5; p=5;
#' y=rnorm(n); x=matrix(rnorm(n*m), nrow=n, ncol=m); 
#' e=matrix(rnorm(n*d), nrow=n, ncol=d); z=matrix(rnorm(n*p), nrow=n, ncol=p);
#' indices=create_splits(K=3, n=n, split_ratio=(1/3, 1/3, 1/3))
#' ritss(y=y, x=x, e=e, z=z, indices=indices)

ritss=function(y, x, e, z, indices, cut_off_p_value=0.05, screening_function=screening_subprs, indices_env_factors=1, verbose=FALSE)
{
	
	if(!is.numeric(y)) stop("y is not a numeric vector.")
	if(var(y)==0) stop("y has no variation.")
	if(length(unique(y))<=5) stop("current implementation is designed for quantitative traits. y seems to be discrete.")
	
    if(class(x)[1]!= "matrix") stop("x is not a matrix.")
	if(class(e)[1]!= "matrix") stop("e is not a matrix.")
	if(class(z)[1]!= "matrix") stop("z is not a matrix.")
	
	
    m=ncol(x); d=ncol(e); p=ncol(z); n=length(y);
	if(m<10) stop("please provide at least 10 genetic variants/SNPs.")

	if(n!=nrow(x) | n!=nrow(e) | n!=nrow(z)) stop("row dimension of x, e, or z does not match length of y.")
	if(!(indices_env_factors %in% 1:d)) stop("indices_env_factors invalid.")
	
	genetic_signals=rep(0, ncol(x))
    K=length(indices)
	sub_objs=list()
	stat=0; variance=0;
	
	for(k in 1:K)
	{
	     obj=ritss_sub(y=y, x=x, e=e, z=z, inds_test=indices[[k]]$inds_test, inds_screen=indices[[k]]$inds_screen, inds_main=indices[[k]]$inds_main,
		 inds_ace=indices[[k]]$inds_ace, cut_off_p_value=cut_off_p_value, screening_function=screening_function, indices_env_factors=indices_env_factors, verbose=verbose)
		 stat=stat + obj$stat
		 variance=variance + obj$var_est
		 genetic_signals[obj$snps_in_score]=genetic_signals[obj$snps_in_score]+1
		 sub_objs[[k]]=list(indices_test=indices[[k]]$inds_test, score=obj$Ui, score_prime=obj$Ui_prime, y_resid=obj$y_resid, 
		 y_resid_train=obj$y_resid_train, indices_main=indices[[k]]$inds_main, main_effects=obj$main_effects)
		 
	}
	###
    
	zsc=stat/sqrt(variance)
    pval_overall=pnorm(-abs(zsc),0,1)*2
	
	return(list("pval"=pval_overall, "z"=zsc, "genetic_signals"=genetic_signals, "sub_objs"=sub_objs))
	
}


ritss_sub=function(y, x, e, z, inds_test, inds_screen, inds_main, inds_ace, cut_off_p_value=0.05, screening_function=screening_subprs, indices_env_factors=1, verbose=FALSE)
{

  ##############################################################
  ### set up
  y_test=y[inds_test]; y_screen=y[inds_screen]; y_main=y[inds_main]; y_ace=y[inds_ace];
  x_test=as.matrix(x[inds_test,]); x_screen=as.matrix(x[inds_screen,]); x_main=as.matrix(x[inds_main,]); x_ace=as.matrix(x[inds_ace,]);
  e_test=as.matrix(e[inds_test,]); e_screen=as.matrix(e[inds_screen,]); e_main=as.matrix(e[inds_main,]); e_ace=as.matrix(e[inds_ace,]);
  z_test=as.matrix(z[inds_test,]); z_screen=as.matrix(z[inds_screen,]); z_main=as.matrix(z[inds_main,]); z_ace=as.matrix(z[inds_ace,]);
  n_var=ncol(x)
  n_test=length(inds_test)

  ##############################################################
  ### screening
  #screening_subprs=function(y, x, e, z, indices_env_factors, num_iterations=10, verbose=FALSE)
  if(verbose==TRUE) st=Sys.time()
  screening_obj=screening_function(y=y_screen, x=x_screen, e=e_screen, z=z_screen, indices_env_factors=indices_env_factors, verbose=verbose) ##!#
  if(verbose==TRUE){et=Sys.time(); cat("screening step done", et-st, "\n");} 
  number_scores=length(screening_obj$mapping_env_factors)
  score_env_mapping=screening_obj$mapping_env_factors
  included_env_factors=unique(score_env_mapping)
  
  ##############################################################
  ### main effects and ACE
  if(verbose==TRUE) st=Sys.time()
  scores_main=matrix(0, nrow=length(y_main), ncol=number_scores)
  scores_ace=matrix(0, nrow=length(y_ace), ncol=number_scores)
  scores_test=matrix(0, nrow=length(y_test), ncol=number_scores)
  scores_x_main=matrix(0, nrow=length(y_main), ncol=number_scores)
  scores_x_ace=matrix(0, nrow=length(y_ace), ncol=number_scores)
  scores_x_test=matrix(0, nrow=length(y_test), ncol=number_scores)
  
  for(k in 1:number_scores)
  {
     betav=screening_obj$betas[[k]]
	 env_index=score_env_mapping[k]
	 
	 scores_main[,k]=e_main[,env_index] * (as.matrix(x_main) %*% betav)
	 scores_ace[,k]=e_ace[,env_index] * (as.matrix(x_ace) %*% betav)
	 scores_test[,k]=e_test[,env_index] * (as.matrix(x_test) %*% betav)
	 
	 scores_x_main[,k]=as.matrix(x_main) %*% betav 
	 scores_x_ace[,k]=as.matrix(x_ace) %*% betav 
	 scores_x_test[,k]=as.matrix(x_test) %*% betav 
  }
  if(verbose==TRUE){et=Sys.time(); cat("set up step done", et-st, "\n");} 
  ### main effects
  #main_effects=function(scores, x, e, z, y, x_p, e_p, z_p)
  if(verbose==TRUE) st=Sys.time()
  lasso_obj=main_effects(scores=scores_main, x=x_main, e=e_main, z=z_main, y=y_main, scores_p=scores_test, x_p=x_test, e_p=e_test, z_p=z_test)
  if(verbose==TRUE){et=Sys.time(); cat("main effects step done", et-st, "\n");} 
  
  pval_score=as.numeric(lasso_obj$pval_scores)
  beta_score=as.numeric(lasso_obj$beta_scores)
  if(verbose==TRUE) cat("filter score betas:", beta_score,"\n")
  if(verbose==TRUE) cat("filter score p-values:",pval_score,"\n")
  
  keep_scores=rep(0, number_scores); 
  keep_scores[pval_score<=cut_off_p_value]=1; keep_scores[1]=1; # always keep first score
  beta_score[keep_scores==0]=0
  beta_score[keep_scores==1]=1
   
  variants_included=rep(0, n_var)
  
  ### interaction score computation
  for(k in 1:number_scores)
  {
     betav=screening_obj$betas[[k]]
	 if(keep_scores[k]==1) {variants_included[betav!=0]=1}
  }
  tmp_inds=1:n_var
  snps_in_score=tmp_inds[variants_included==1]
  
  if(verbose==TRUE){cat("snps in score:", length(snps_in_score), "\n")}
  
  ### ACE
  Ui_test=scores_test %*% beta_score
  Ui_test_prime=rep(0, n_test)
  if(verbose==TRUE) st=Sys.time()
  score_tmp_indices=1:number_scores
  for(index in included_env_factors)
  {
     inds=score_tmp_indices[score_env_mapping==index]
	 if(verbose==TRUE){cat("env.factor: ", index,  inds, beta_score[inds], "\n")}
	 
	 tmp_Ui_ace=as.matrix(scores_ace[,inds]) %*% as.vector(beta_score[inds])
	 tmp_Ui_test=as.matrix(scores_test[,inds]) %*% as.vector(beta_score[inds])
	 
	 tmp_score_x_ace=as.matrix(scores_x_ace[,inds]) %*% as.vector(beta_score[inds])
     tmp_score_x_test=as.matrix(scores_x_test[,inds]) %*% as.vector(beta_score[inds])
	 
	 ace=ACE(x=x_ace, e=e_ace, z=z_ace, ui=tmp_Ui_ace, index_e=index, score_x=tmp_score_x_ace,
	 x_p=x_test, e_p=e_test, z_p=z_test, ui_p=tmp_Ui_test, score_x_p=tmp_score_x_test)
     Ui_test_prime=Ui_test_prime+ace$ui_prime
	 if(verbose==TRUE) cat("number of selected genetic variants in prediction of e: ", length(ace$x_inds), "\n")
  }
 
  
  y_resid=y_test-lasso_obj$predicted_y_main_effects
  y_resid_train=y_main-lasso_obj$predicted_y_main_effects_train
  if(verbose==TRUE){et=Sys.time(); cat("ACE step done", et-st, "\n");}
  ##############################################################
  ### test statistic computation
  
  S=sum(Ui_test_prime*y_resid)
  var_est=sum(Ui_test_prime^2*y_resid^2)
  if(verbose==TRUE){cat("S , sd(S): ", S, sqrt(var_est), "\n")}
  
  return(list("stat"=S, "var_est"=var_est, "betas"=screening_obj$betas, "snps_in_score"=snps_in_score, "Ui"=Ui_test, "Ui_prime"=Ui_test_prime, "y_resid"=y_resid, "y_resid_train"=y_resid_train, 
  "main_effects"=lasso_obj$predicted_y_main_effects))
 
}
# screening step to identify interaction with subcomponent of genetic risk score
#' screening_subprs
#' @export
screening_subprs=function(y, x, e, z, indices_env_factors, num_iterations=10, verbose=FALSE)
{
   n=nrow(x)
   if(length(indices_env_factors)!=1) stop("error: only one environmental component allowed for this function.")
   inter_index=indices_env_factors;
   if(inter_index<1 | inter_index > ncol(e)) stop("error: index for environmental component invalid.")
   
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
   if(length(ind_var)<10){ ind_var=tmp$indices_snps; }
   
   best_subset_reg=best_subset_reg(yr=yr, xe_r=xe_r[,ind_var], S_max=length(ind_var))
   
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
best_subset_reg=function(yr, xe_r, S_max, num_iterations=10)
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
		
		C=get_C(train_y, train_x)
	    start_beta=rep(0, ncol(xe_r))
		ctr=0
	    for(size in seq_s)
	    {
			
			beta_tmp=gradient_algorithm(y=train_y, x=train_x, start_beta=start_beta, C=C, size=size, n_iter=num_iterations);
			tmp_inds=indices_snps[beta_tmp!=0]
			
			score=rowSums(test_x[, tmp_inds])
			fit=lm(test_y~score)
			z_sq=summary(fit)$coefficients[2,3]**2
			
			
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

# screening based on single variant GxE tests
#' screening_sv
#' @export
screening_sv=function(y, x, e, z, indices_env_factors, fdr_cutoff1=0.05, fdr_cutoff2=0.1, verbose=FALSE)
{
   n=nrow(x)
   if(length(indices_env_factors)!=1) stop("error: only one environmental component allowed for this function.")
   inter_index=indices_env_factors;
   if(inter_index<1 | inter_index > ncol(e)) stop("error: index for environmental component invalid.")
   
   ############################
  
   e_t=e[,inter_index]
   pval=rep(0, ncol(x))
   coef=rep(0, ncol(x))
   for(i in 1:ncol(x))
   {
		xe=x[,i]*e_t
		fit=lm(y~xe+x[,i]+e+z)
		pval[i]=summary(fit)$coefficients[2,4]
		coef[i]=summary(fit)$coefficients[2,1]
   }
   pval_adj=p.adjust(pval, method="BH")
   n1=sum(pval_adj<=fdr_cutoff1)
   n2=sum(pval_adj>fdr_cutoff1 & pval_adj<=fdr_cutoff2)
   mm1=pval_adj<=fdr_cutoff1
   mm2=pval_adj>fdr_cutoff1 & pval_adj<=fdr_cutoff2
   if(n1==0 & n2==0)
   {
      mm1=pval_adj==min(pval_adj)
	  n1=sum(mm1)
	  n2=0
   }
   if(n1==0 & n2>0)
   {
      mm1=pval_adj<=fdr_cutoff2
	  n1=sum(mm1)
	  n2=0
   }
   
   num_indices=1
   if(n1>0 & n2>0){num_indices=2}
   indices=matrix(0, nrow=num_indices, ncol=ncol(x))
   indices[1, mm1]=1
   if(num_indices==2){indices[2, mm2]=1}
   
   if(verbose==TRUE) print("screening done")
   
   beta_1=rep(0, ncol(x));
   beta_2=rep(0, ncol(x));
   beta_1[indices[1,]==1]=coef[indices[1,]==1];
   if(num_indices==2) beta_2[indices[2,]==1]=coef[indices[2,]==1]
   mapping=rep(0, num_indices)
   mapping[1]=inter_index;
   if(num_indices==2) mapping[2]=inter_index;
   
   return(list("betas"=list(beta_1, beta_2), "mapping_env_factors"=mapping))
  
}

# create sample splits according to prespecified fractions and K
#' create_splits
#' @export
create_splits=function(K, n, split_ratio)
{
    folds=sample(cut(seq(1,n),breaks=K,labels=FALSE))
	inds=1:n
	indices=list()
	for(k in 1:K)
	{
	   tmp=inds[folds==k]
	   tmpc=inds[folds!=k]
	   splits=sample(1:3, size=length(tmpc), replace=TRUE, prob=split_ratio)
	   tmp=list(inds_screen=tmpc[splits==1], inds_main=tmpc[splits==2], inds_ace=tmpc[splits==3], inds_test=tmp)
	   indices[[k]]=tmp
	}
	return(indices)
}
