####################################
### load libraries and functions

library(sandwich)
library(RITSS)

source("sim_functions.R")
source("alternative_approaches.R") 

args = commandArgs(trailingOnly=TRUE)

#########################################

mean_xe=as.numeric(args[1]);
sd_xe=as.numeric(args[2]); 
prop_xe=as.numeric(args[3]);

########################################################################
### setup

# sd_x: std. err genotype effects
# sd_e: std. err env effects
# sd_err_e: std. err effect 1. env. factor on epsilon variance
# sd_e_x: std. err genotype effect on env. factors
# sd_e_z: std. err covariate effect on env. factors
# sd_eps: std. err epsilon
# sd_xe: std. err GxE effects
# prop_e_x: proportion non-zero genotype effects on env. factors
# prop_e_z: proportion non-zero covariate effects on env. factors
# prob_ex: proportion of non-zero GxE effects
# sign_e_x: put all genotype effects on env. factors in positive direction
# sign_e: put all env. factor effects in positive direction
# exponent: exponent of first env. factor in main effect

n_rep=1000
n=30000
n_var=100
n_e=5
n_z=5

sd_x=0.1
sd_e=0.4
sd_z=0.0
sd_err_e=0.0
sd_e_x=0.0
sd_e_z=0.0
sd_eps=1.0
beta_PS=0.0

prop_x=1
prop_e_x=0.2
prop_e_z=0.5

min=0.0
max=0.1
maf=0.3
exponent=1


seed=as.numeric(floor(mean_xe*20) + floor(sd_xe*1000) + floor(prop_xe*100000))


####################################
### simulation and analysis

set.seed(seed)

n_k=floor(n/3)
i1=1:n_k
i2=(n_k+1):(2*n_k)
i3=(2*n_k+1):n
results=matrix(0, nrow=n_rep, ncol=9) 
for(i in 1:n_rep)
{
        ### draw parameters and simulate phenotype
        beta_x=draw_normal_effects(0, sd_x, n_var, prop_x, FALSE)
        beta_e=draw_normal_effects(0, sd_e, n_e, 1, TRUE)
		beta_z=draw_normal_effects(0, sd_z, n_z, 1, FALSE)
		beta_xe=draw_normal_effects(mean_xe, sd_xe, n_var, prop_xe, TRUE) 
		beta_err_e=draw_normal_effects(0, sd_err_e, n_e, 1, TRUE)
		beta_e_x=draw_normal_effects_matrix(0, sd_e_x, n_var, n_e, prop_e_x, TRUE)
		beta_e_z=draw_normal_effects_matrix(0, sd_e_z, n_z, n_e, prop_e_z, FALSE)
		
		z=draw_z(min, max, n_z, n)
		x=draw_x_given_z(z, n_var, beta_PS, maf)
		e=draw_e_given_x_and_z(z, x, n_e, beta_e_x, beta_e_z)
        env_e=get_environmental_effect(e, exponent, beta_e)
		env_z=get_covariate_effect(z, beta_z) 
		genetic=get_genetic_effect(x, beta_x)
		gxe=get_gxe_effect(x, e, beta_xe, beta_x) 
		eps=get_error(n, sd_eps, beta_err_e, FALSE, errors=numeric(0))

		y=genetic + gxe + env_e + env_z + eps 
		y=as.numeric(y)
		
		
		y_1=y[i1]; y_2=y[i2]; y_3=y[i3];
		x_1=x[i1,]; x_2=x[i2,]; x_3=x[i3,];
		e_1=as.matrix(e[i1,]); e_2=as.matrix(e[i2,]); e_3=as.matrix(e[i3,])
		z_1=as.matrix(z[i1,]); z_2=as.matrix(z[i2,]); z_3=as.matrix(z[i3,]);
		
		####################################################
		### RITSS
		ritss_obj=ritss(y, x, e, z, i1, i2, i3, cut_off_p_value=0.05, verbose=FALSE)
		####################################################
		### alternative_approaches 
		lr=linreg_inter_sv(x, e, z, y, inter_index=1)
		lr2split=linreg_2split(y, x, e, z)
		gesat_p=gesat(x, e, z, y)
		###################################################	  
		results[i,1]=ritss_obj$z_1; results[i,2]=ritss_obj$z_2; results[i,3]=ritss_obj$z_3; 
		results[i,4]=ritss_obj$z; results[i,5]=ritss_obj$pval; 
		results[i,6]=lr; results[i,7]=lr2split$pval; results[i,8]=lr2split$pval_robust; results[i,9]=gesat_p; 

		cat(results[i,], "\n")
}
