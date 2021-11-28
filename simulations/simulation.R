####################################
### load libraries and functions

library(config)
library(sandwich)
library(RITSS)

source("sim_functions.R")
args = commandArgs(trailingOnly=TRUE)

#########################################
errors=NULL # add non-normal errors 
select="no"

#########################################
ps=args[1]
gec=args[2]
mem=args[3]
nne=args[4]
he=args[5]


########################################################################
### setup

# sd_x: std. err genotype effects
# sd_e: std. err env effects
# sd_err_e: std. err effect 1. env. factor on epsilon variance
# sd_e_x: std. err genotype effect on env. factors
# sd_e_z: std. err covariate effect on env. factors
# sd_eps: std. err epsilon
# prop_e_x: proportion non-zero genotype effects on env. factors
# prop_e_z: proportion non-zero covariate effects on env. factors
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

prop_e_x=0.2
prop_e_z=0.5

min=0.0
max=0.1
maf=0.3
exponent=1


seed=as.numeric(1)
if(ps=="yes"){beta_PS=1;sd_e_z=0.1;sd_z=0.1;seed=seed+2;}
if(gec=="yes"){sd_e_x=0.1;seed=seed+4;}
if(mem=="yes"){exponent=2;seed=seed+8;}
if(nne=="yes"){seed=seed+16;}
if(he=="yes"){sd_err_e=0.5;seed=seed+32}
if(select=="yes"){seed=seed+64;}



####################################
### simulation and analysis

set.seed(seed)

n_k=floor(n/3)
i1=1:n_k
i2=(n_k+1):(2*n_k)
i3=(2*n_k+1):n
results=matrix(0, nrow=n_rep, ncol=8)
for(i in 1:n_rep)
{
        ### draw parameters and simulate phenotype
        beta_x=draw_normal_effects(0, sd_x, n_var, 1, FALSE)
        beta_e=draw_normal_effects(0, sd_e, n_e, 1, TRUE)
		beta_z=draw_normal_effects(0, sd_z, n_z, 1, FALSE)
		beta_err_e=draw_normal_effects(0, sd_err_e, n_e, 1, TRUE)
		beta_e_x=draw_normal_effects_matrix(0, sd_e_x, n_var, n_e, prop_e_x, TRUE)
		beta_e_z=draw_normal_effects_matrix(0, sd_e_z, n_z, n_e, prop_e_z, FALSE)
		
		z=draw_z(min, max, n_z, n)
		x=draw_x_given_z(z, n_var, beta_PS, maf)
		e=draw_e_given_x_and_z(z, x, n_e, beta_e_x, beta_e_z)
        env_e=get_environmental_effect(e, exponent, beta_e)
		env_z=get_covariate_effect(z, beta_z) 
		genetic=get_genetic_effect(x, beta_x)
		eps=get_error(n, sd_eps, beta_err_e, nne, errors)

		y=genetic + env_e + env_z + eps
		y=as.numeric(y)
		###
		if(select==TRUE)
		{
		   index_select=n_var/2
		   inds=select_x(y, x, e, z, index_select)
		   x=x[,inds]
		}
		###
		
		y_1=y[i1]; y_2=y[i2]; y_3=y[i3];
		x_1=x[i1,]; x_2=x[i2,]; x_3=x[i3,];
		e_1=as.matrix(e[i1,]); e_2=as.matrix(e[i2,]); e_3=as.matrix(e[i3,])
		z_1=as.matrix(z[i1,]); z_2=as.matrix(z[i2,]); z_3=as.matrix(z[i3,]);
		
		####################################################
		### RITSS
		ritss_obj=ritss(y, x, e, z, i1, i2, i3, cut_off_p_value=0.05, verbose=FALSE)
		################################################
		
		### REG and REG-robust
		xe=ritss_obj$Ui_3
		fit=lm(y_3~xe+x_3+e_3+z_3)
		pval_reg_score=summary(fit)$coefficients[2,4]
		beta_reg_score=summary(fit)$coefficients[2,1]
		sandwich_se <- diag(vcovHC(fit, type = "HC"))^0.5
		tmp=beta_reg_score/sandwich_se[2]
		pval_reg_score_robust=pnorm(-abs(tmp), 0, 1)*2
		###################################################
		
		results[i,1]=ritss_obj$z_1; results[i,2]=ritss_obj$z_2; results[i,3]=ritss_obj$z_3; 
		results[i,4]=ritss_obj$z; results[i,5]=ritss_obj$pval; 
		results[i,6]=pval_reg_score; results[i,7]=beta_reg_score; results[i,8]=pval_reg_score_robust;

		cat(results[i,], "\n")
}
