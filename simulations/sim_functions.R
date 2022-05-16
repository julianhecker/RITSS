#### select the n=index variants with the smallest marginal association p-value or all genome-wide significant
select_x=function(y, x, e, z, index)
{
		fit=lm(y~x+e+z)
		pvals_x=summary(fit)$coefficients[2:(ncol(x)+1),4]
		ngws=sum(pvals_x<=5*10^-8)
		tmp=sort(pvals_x, decreasing = FALSE)[index]
		if(ngws>0.3*ncol(x) & ngws<=0.5*ncol(x)){tmp=5*10^-8}
		return(pvals_x<=tmp)
		
}
#### generate z
draw_z=function(min, max, n_z, n)
{
    z=matrix(runif(n * n_z, min=min, max=max), nrow=n, ncol=n_z)
	return(z)
}
#### generate x given z
draw_x_given_z=function(z, n_var, beta_PS, maf)
{
   n_z=ncol(z); n=nrow(z)
   seq=rep(1:n_z, n_var/n_z)
   mafs=matrix(maf, nrow=n, ncol=n_var)+beta_PS*z[,seq]
   x=matrix(0, nrow=n, ncol=n_var)
   for(i in 1:ncol(x)){x[,i]=rbinom(n, 2, mafs[,i])}
   return(x)
}
#### draw e given x and z
draw_e_given_x_and_z=function(z, x, n_e, beta_e_x, beta_e_z)
{
   n=nrow(x)
   e=matrix(rnorm(n*n_e, 0, 1), nrow=n, ncol=n_e)
   for(i in 1:n_e)
   {
		e[,i]=e[,i]+x %*% beta_e_x[,i] + z %*% beta_e_z[,i]
   }
   return(e)
}
#### normally distributed effects
draw_normal_effects=function(mean, sd, m, proportion, sign=FALSE)
{
   beta=rnorm(m, mean, sd)
   ind=rbinom(m, 1, proportion)
   beta[ind==0]=0
   if(sign==TRUE){ beta=abs(beta)}
   return(beta)
}
#### normally distributed effects (matrix)
draw_normal_effects_matrix=function(mean, sd, m, ncol, proportion, sign=FALSE)
{
   beta=matrix(rnorm(m * ncol, mean, sd), nrow=m, ncol=ncol)
   for(i in 1:ncol)
   {
       ind=rbinom(m, 1, proportion)
	   beta[ind==0,i]=0
	   if(sign==TRUE){ beta[,i]=abs(beta[,i])}
   }
   return(beta)
}
#### get environmental main effect
get_environmental_effect=function(e, exponent, beta_e)
{
   env_e=beta_e[1]*e[,1]^exponent
   
   n_e=ncol(e)
   if(n_e>=2)
   {
	   for(i in 2:n_e)
	   {
		 env_e=env_e+beta_e[i]*e[,i]
	   }
   }
   return(env_e)
}
#### get covariate main effect
get_covariate_effect=function(z, beta_z)
{
  env_z=z %*% beta_z
  return(env_z)
}
#### get genetic main effect
get_genetic_effect=function(x, beta_x)
{
  g=x %*% beta_x
  return(g)
}
#### get GxE effect (power studies) ##!#
get_gxe_effect=function(x, e, beta_xe, beta_x)
{
  ge=e[,1] * x %*% (beta_x * beta_xe)
  return(ge)
}##!#
#### get error
get_error=function(n, sd, beta_err_e, nne=FALSE, errors=numeric(0))
{
  eps=rnorm(n)*sd*(1+e[,1]*beta_err_e[1])
  if(nne==TRUE & length(errors)==0) {print("warning: NNE specified but no error provided, using normal errors."); errors=rnorm(n)}
  if(nne==TRUE) # assumes that length(errors)>=n
  {
     eps=errors[sample(1:length(errors), n, replace=F)]
	 eps=eps-mean(eps)
	 eps=eps/sd(eps)*sd
	 eps=eps+e[,1]*beta_err_e[1]*eps
  }
  return(eps)
}




