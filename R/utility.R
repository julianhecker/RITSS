# estimate main effects while fitting candidate scores
lasso_main_effects=function(scores, x, e, z, y, x_p, e_p, z_p)
{
  xez_mat=create_lasso_matrix_scores_x_e_z_simple(scores, x, e, z, x_p, e_p, z_p)
  fit=lasso_fit_1se_refit(xez_mat$mat, y, xez_mat$inds_penalty, ncol(scores))
  y_pred=lasso_predict(xez_mat$mat, fit)
  mat=xez_mat$mat
  mat_p=xez_mat$mat_p
  for(i in 1:ncol(z))
  {
     tmp=x*z[,i]
	 tmp_p=x_p*z_p[,i]
	 for(j in 1:ncol(tmp))
	 {
			 fit=lm(y~tmp[,j]+y_pred)
			 pval=summary(fit)$coefficients[2,4]
			 if(pval<=0.01)
			 {
				mat=cbind(mat, tmp[,j])
				mat_p=cbind(mat_p, tmp_p[,j])
			 }
	 }
  }
  mat=as.matrix(mat)
  mat_p=as.matrix(mat_p)
  inds_penalty=(ncol(scores)+ncol(x)+ncol(e)+ncol(z)+1):ncol(mat)
  fit=lasso_fit_1se_refit(mat, y, inds_penalty, ncol(scores))
  predicted_y_main_effects=fit$intercept + as.matrix(mat_p) %*% fit$beta[-c(1:ncol(scores))]
  return(list("beta"=fit$beta_vec, "intercept"=fit$intercept, "pval_scores"=fit$pval_scores, "beta_scores"=fit$beta_scores, "predicted_y_main_effects"=predicted_y_main_effects))
  
}
# ACE algorithm: use trained data to predict objects in other part of data
ACE_pred=function(x, e, z, ui, rg, rf, rg2, rg2_indicators, rf2, score_x)
{
   xz_mat=create_lasso_matrix_e_given_x_and_z(x, z)
   rfp=lasso_predict(xz_mat, rg)
   uip=ui-rfp*score_x
   
   rfp=as.numeric(predict(rf, cbind(e, z))$predictions)
   uip=uip-rfp

   uip=uip-lasso_xz_pred(x, z, uip, rg2, rg2_indicators)
   
   rfp=as.numeric(predict(rf2, cbind(e, z))$predictions)
   uip=uip-rfp
   
   return(uip) 
}
# ACE algorithm: training
ACE_train=function(x, e, z, ui, index_e, score_x)
{

  xz_mat=create_lasso_matrix_e_given_x_and_z(x, z)
  rg=lasso_fit_1se(xz_mat, e[,index_e])
  rfp=lasso_predict(xz_mat, rg)
  uip=ui-rfp*score_x
  x_inds=1:ncol(x)
  x_inds=x_inds[rg$beta!=0]
  
 
  rf=grf::regression_forest(cbind(e,z), uip, num.trees = 100, seed=1, num.threads=4)
  rfp=as.numeric(predict(rf, cbind(e, z))$predictions)
  uip=uip-rfp
  
  lasso_obj=lasso_xz(x, z, uip)
  rg2=lasso_obj$rg
  uip=uip-lasso_obj$y_pred
  rg2_indicators=lasso_obj$indicators
  
  rf2=grf::regression_forest(cbind(e, z), uip, num.trees = 100, seed=1, num.threads=4)
  rfp=as.numeric(predict(rf2, cbind(e, z))$predictions)
  uip=uip-rfp
  
  return(list("rg"=rg, "rf"=rf, "rg2"=rg2, "rg2_indicators"=rg2_indicators, "rf2"=rf2, "x_inds"=x_inds))
}
# create matrix to predict based on trained LASSO
lasso_xz_pred=function(x, z, y, rg, indicators)
{
    xz_mat=as.matrix(cbind(x, z))
	for(i in 1:ncol(z))
	{
		 tmp=x*z[,i]
		 for(j in 1:ncol(tmp))
		 {
			 if(indicators[i,j]==1)
			 {
				xz_mat=cbind(xz_mat, tmp[,j])
			 }
		 }
	}
	xz_mat=as.matrix(xz_mat)
	y_pred=lasso_predict(xz_mat, rg)
	return(y_pred)
}
# LASSO training y~x+z
lasso_xz=function(x, z, y)
{
	xz_mat=as.matrix(cbind(x, z))
	rg=lasso_fit_1se(xz_mat, y)
	y_pred=lasso_predict(xz_mat, rg)
	indicators=matrix(0, nrow=ncol(z), ncol=ncol(x))
	for(i in 1:ncol(z))
	{
		 tmp=x*z[,i]
		 for(j in 1:ncol(tmp))
		 {
			 fit=lm(y~tmp[,j]+y_pred)
			 pval=summary(fit)$coefficients[2,4]
			 if(pval<=0.01)
			 {
				xz_mat=cbind(xz_mat, tmp[,j])
				indicators[i,j]=1
			 }
		 }
	}
	rg=lasso_fit_1se(as.matrix(xz_mat), y)
	y_pred=lasso_predict(xz_mat, rg)
	return(list("y_pred"=y_pred, "rg"=rg, "indicators"=indicators))
}

# LASSO fit 1 se solution
lasso_fit_1se=function(x, y)
{
   nf=10
   foldids=rep(seq(nf), length=length(y))
   
   cvlasso=glmnet::cv.glmnet(x, y, nfolds=nf, foldid=foldids)
   fit=glmnet::glmnet(x, y, lambda=cvlasso$lambda.1se)
   intercept=fit$a0
   beta_vec=as.numeric(fit$beta)
   
   return(list("beta"=beta_vec, "intercept"=intercept))
}
# predict using intercept and beta
lasso_predict=function(x, fit)
{
   return(fit$intercept+x%*%fit$beta)
}
# LASSO fit 1 se solution with refit
lasso_fit_1se_refit=function(x, y, inds_penalty, n_scores)
{
   nf=10
   foldids=rep(seq(nf), length=length(y))
   penalty = rep(0, ncol(x))
   penalty[inds_penalty]=1
   
   cvlasso=glmnet::cv.glmnet(x, y, nfolds=nf, foldid = foldids, penalty.factor=penalty)
   fit=glmnet::glmnet(x, y, penalty.factor=penalty, lambda=cvlasso$lambda.1se)
   intercept=fit$a0
   beta_vec=as.numeric(fit$beta)
   
   xp=x[,abs(beta_vec)!=0]
   fit=lm(y~xp)
   intercept=summary(fit)$coefficients[1,1]
   beta_vec[abs(beta_vec)!=0]=summary(fit)$coefficients[-1,1]
   
   beta_scores=beta_vec[1:n_scores]
   pval_scores=summary(fit)$coefficients[2:(n_scores+1),4]
   return(list("beta"=beta_vec, "intercept"=intercept, "pval_scores"=pval_scores, "beta_scores"=beta_scores))
}
# create matrix
create_lasso_matrix_e_given_x_and_z=function(x, z)
{
   mat=cbind(x, z)
   return(as.matrix(mat))
}
# create matrix
create_lasso_matrix_scores_x_e_z_simple=function(scores, x, e, z, x_p, e_p, z_p)
{
   mat=cbind(scores, x, e, z)
   ctr=ncol(scores)+ncol(x)+ncol(e)+ncol(z)
   mat_p=cbind(x_p, e_p, z_p)
   for(i in 1:ncol(e))
   {
      if(length(unique(e[,i]))>2){ mat=cbind(mat, (e[,i]-mean(e[,i]))**2); }
	  if(length(unique(e[,i]))>2) mat_p=cbind(mat_p, (e_p[,i]-mean(e_p[,i]))**2)
   }
   for(i in 1:ncol(z))
   {
      if(length(unique(z[,i]))>2){  mat=cbind(mat, (z[,i]-mean(z[,i]))**2);}
	  if(length(unique(z[,i]))>2)  mat_p=cbind(mat_p, (z_p[,i]-mean(z_p[,i]))**2)
   }
   for(i in 1:ncol(e))
   {
      mat=cbind(mat, z*e[,i])
	  mat_p=cbind(mat_p, z_p*e_p[,i])
   }
   mat=as.matrix(mat)
   mat_p=as.matrix(mat_p)
   inds_penalty=(ctr+1):ncol(mat)
   return(list("mat"=mat, "mat_p"=mat_p, "inds_penalty"=inds_penalty))
}

