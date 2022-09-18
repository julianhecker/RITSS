######################################
###	Implementation of approximate best subset regression algorithm, 
### see section 3 of BY BERTSIMAS, KING, MAZUMDER Ann. Statist. 44(2): 813-852 (April 2016). DOI: 10.1214/15-AOS1388
##############################################################################


# compute gradient
gradient_lr=function(y, x, beta)
{
  tmp=-t(x)%*%(y-x%*%beta)
  return(tmp)
}
# recompute beta
recompute=function(y, x, beta, C, size)
{
  tmp=beta-1/C*gradient_lr(y=y, x=x, beta=beta)
  tmp2=sort(abs(tmp), decreasing = TRUE)[size+1]
  tmp[abs(tmp)<=tmp2]=0
  return(tmp)
}
# get C for normalization
get_C=function(y, x)
{
  mat_product=t(x)%*%x
  eigen_decomp=eigen(mat_product)
  lambda_max=eigen_decomp$values[1]
  C=lambda_max*1.2
  return(C)
}
# approximate best subset regression algorithm
gradient_algorithm=function(y, x, start_beta, C, size, n_iter)
{
  beta_tmp=start_beta
  for(i in 1:n_iter)
  {
    beta_new=recompute(y=y, x=x, beta=beta_tmp, C=C, size=size)
    beta_tmp=beta_new
  }
  return(beta_tmp)
}