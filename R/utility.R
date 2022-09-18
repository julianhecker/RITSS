ACE=function(x, e, z, ui, index_e, score_x, x_p, e_p, z_p, ui_p, score_x_p)
{
  ###################
  xz_mat=create_lasso_matrix_e_given_x_and_z(x=x, z=z)
  rg=lasso_fit_1se(x=xz_mat, y=e[,index_e])
  rfp=lasso_predict(x=xz_mat, fit=rg)
  ui=ui-rfp*score_x
  x_inds=1:ncol(x)
  x_inds=x_inds[rg$beta!=0]
  ###
  xz_mat_p=create_lasso_matrix_e_given_x_and_z(x=x_p, z=z_p)
  rfp=lasso_predict(x=xz_mat_p, fit=rg)
  ui_p=ui_p-rfp*score_x_p
  
  ####################
  data=cbind(e, z)
  colnames(data)=paste0("f",1:ncol(data))
  form="ui~"
  for(i in 1:(ncol(data)-1))
  {
	  if(length(unique(data[,i]))>2){form=paste0(form, paste0("s(f",i,",bs='cr',k=3)+"))}
	  if(length(unique(data[,i]))==2){form=paste0(form, paste0("f",i,"+"))}
  }
  i=ncol(data)
  if(length(unique(data[,i]))>2){form=paste0(form, paste0("s(f",i,",bs='cr',k=3)"))}
  if(length(unique(data[,i]))==2){form=paste0(form, paste0("f",i))}
  
  data=as.data.frame(data)
  gam_fit=mgcv::bam(as.formula(paste(form)), data=data)
  ui=gam_fit$residuals
  
  data=cbind(e_p, z_p)
  colnames(data)=paste0("f",1:ncol(data))
  data=as.data.frame(data) 
  prediction=predict(gam_fit, newdata=data)
  ui_p=ui_p-as.numeric(prediction)
  
  ####################
  lasso_obj=lasso_xz(x=x, z=z, y=ui)
  rg2=lasso_obj$rg
  ui=ui-lasso_obj$y_pred
  rg2_indicators=lasso_obj$indicators
  ###
  ui_p=ui_p-lasso_xz_pred(x=x_p, z=z_p, y=ui_p, rg=rg2, indicators=rg2_indicators)

  ####################
  data=cbind(e, z)
  colnames(data)=paste0("f",1:ncol(data))
  form="ui~"
  for(i in 1:(ncol(data)-1))
  {
	  if(length(unique(data[,i]))>2){form=paste0(form, paste0("s(f",i,",bs='cr',k=3)+"))}
	  if(length(unique(data[,i]))==2){form=paste0(form, paste0("f",i,"+"))}
  }
  i=ncol(data)
  if(length(unique(data[,i]))>2){form=paste0(form, paste0("s(f",i,",bs='cr',k=3)"))}
  if(length(unique(data[,i]))==2){form=paste0(form, paste0("f",i))}
  
  data=as.data.frame(data)
  gam_fit=mgcv::bam(as.formula(paste(form)), data=data)
  
  data=cbind(e_p, z_p)
  colnames(data)=paste0("f",1:ncol(data))
  data=as.data.frame(data) 
  prediction=predict(gam_fit, newdata=data)
  ui_p=ui_p-as.numeric(prediction)
   
  return(list("ui_prime"=ui_p, "x_inds"=x_inds)) 
}
main_effects=function(scores, x, e, z, y, scores_p, x_p, e_p, z_p, bool_xz=FALSE)
{
  n_scores=ncol(scores)
  data=cbind(scores, x)
  if(bool_xz==TRUE){ for(i in 1:ncol(z)){data=cbind(data, x*z[,i])}}
  colnames(data)=paste0("f",1:ncol(data))
  form="y~"
  for(i in 1:ncol(data))
  {
	  form=paste0(form, paste0("f",i,"+"))
  }
  colnames(data)[ncol(data)]=paste0("f",ncol(data))
  data=as.data.frame(data)
  ncol_tmp=ncol(data)
  data=cbind(data,e, z)
  ###!###
  data=as.data.frame(data)
  colnames(data)[(ncol_tmp+1):ncol(data)]=paste0("c",1:(ncol(data)-ncol_tmp))
  for(i in 1:(ncol(data)-1-ncol_tmp))
  {
	  
	  if(length(unique(data[,i+ncol_tmp]))>2){form=paste0(form, paste0("s(c",i,",bs='cr',k=3)+"))}
	  if(length(unique(data[,i+ncol_tmp]))==2){form=paste0(form, paste0("c",i,"+"))}
  }
  
  i=ncol(data)-ncol_tmp
  if(length(unique(data[,i+ncol_tmp]))>2){form=paste0(form, paste0("s(c",i,",bs='cr',k=3)"))}
  if(length(unique(data[,i+ncol_tmp]))==2){form=paste0(form, paste0("c",i))}
  
  gam_fit=mgcv::bam(as.formula(paste(form)), data=data)
  gam_fit_summary=summary(gam_fit)
 
  beta_scores=unname(gam_fit_summary$p.coeff[1+(1:n_scores)])
  pval_scores=unname(gam_fit_summary$p.pv[1+(1:n_scores)])
  
  data_init=data
  ###########################################
  data=cbind(scores_p, x_p)
  if(bool_xz==TRUE){ for(i in 1:ncol(z_p)){data=cbind(data, x_p*z_p[,i])}}
  colnames(data)=paste0("f",1:ncol(data))
  data=as.data.frame(data)
  ncol_tmp=ncol(data)
  data=cbind(data,e_p, z_p)
  ###!###
  data=as.data.frame(data)
  colnames(data)[(ncol_tmp+1):ncol(data)]=paste0("c",1:(ncol(data)-ncol_tmp))
 
  prediction=predict(gam_fit, newdata=data)
  prediction_train=predict(gam_fit, newdata=data_init)
  predicted_y_main_effects=as.numeric(prediction)-scores_p%*%as.numeric(beta_scores)
  predicted_y_main_effects_train=as.numeric(prediction_train)-scores%*%as.numeric(beta_scores)
  return(list("beta"=rep(0, ncol(x)), "intercept"=0, "pval_scores"=pval_scores, "beta_scores"=beta_scores, "predicted_y_main_effects"=predicted_y_main_effects, "predicted_y_main_effects_train"=predicted_y_main_effects_train))
  
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
	y_pred=lasso_predict(x=xz_mat, fit=rg)
	return(y_pred)
}
# LASSO training y~x+z
lasso_xz=function(x, z, y)
{
	xz_mat=as.matrix(cbind(x, z))
	rg=lasso_fit_1se(x=xz_mat, y=y)
	y_pred=lasso_predict(x=xz_mat, fit=rg)
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
	rg=lasso_fit_1se(x=as.matrix(xz_mat), y=y)
	y_pred=lasso_predict(x=xz_mat, fit=rg)
	return(list("y_pred"=y_pred, "rg"=rg, "indicators"=indicators))
}
# LASSO fit 1 se solution
lasso_fit_1se=function(x, y)
{
   nf=10
   foldids=rep(seq(nf), length=length(y))
   
   cvlasso=glmnet::cv.glmnet(x=x, y=y, nfolds=nf, foldid=foldids)
   fit=glmnet::glmnet(x=x, y=y, lambda=cvlasso$lambda.1se)
   intercept=fit$a0
   beta_vec=as.numeric(fit$beta)
   
   return(list("beta"=beta_vec, "intercept"=intercept))
}
# predict using intercept and beta
lasso_predict=function(x, fit)
{
   return(fit$intercept+x%*%fit$beta)
}

# create matrix
create_lasso_matrix_e_given_x_and_z=function(x, z)
{
   mat=cbind(x, z)
   return(as.matrix(mat))
}

