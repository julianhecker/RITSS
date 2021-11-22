######################################
###	Implementation of approximate best subset regression algorithm, see section 3 of BY BERTSIMAS, KING, MAZUMDER Ann. Statist. 44(2): 813-852 (April 2016). DOI: 10.1214/15-AOS1388
##############################################################################


# compute gradient
gradient_lr=function(y, x, beta)
{
  tmp=-t(x)%*%(y-x%*%beta)
  return(tmp)
}
# recompute beta
recompute=function(y, x, beta, L, s)
{
  tmp=beta-1/L*gradient_lr(y, x, beta)
  tmp2=sort(abs(tmp), decreasing = TRUE)[s+1]
  tmp[abs(tmp)<=tmp2]=0
  return(tmp)
}
# get L for normalization
get_L=function(y, x)
{
  mat_product=t(x)%*%x
  eigen_decomp=eigen(mat_product)
  lambda_max=eigen_decomp$values[1]
  L=lambda_max*1.2
  return(L)
}
# approximate best subset regression algorithm
gradient_algorithm=function(y, x, start_beta, L, s, n_iter)
{
  beta_tmp=start_beta
  for(i in 1:n_iter)
  {
    beta_new=recompute(y, x, beta_tmp, L, s)
    beta_tmp=beta_new
  }
  return(beta_tmp)
}