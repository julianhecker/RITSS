####################################
### load libraries and functions

library(config)
library(sandwich)
library(MLmetrics)
library(RITSS)

config <- config::get(file = "../Rconfig.yml")

source("sim_functions.R")
source("alternative_approaches.R")

args = commandArgs(trailingOnly=TRUE)

#########################################
### load FEV1/FVC data for non-normal error simulations
ukb_pheno_data="../UKB_application/data/phenotype_and_covariate_data"
if(file.exists(ukb_pheno_data)){
	pheno=read.table(ukb_pheno_data,as.is=T, header=T, sep="\t")
	errors=pheno$fev1fvc
	errors=errors[!is.na(errors)]
}else{
	errors=numeric(0)
}

#########################################

scenario=as.numeric(args[1])
seed_num=as.numeric(args[2])
mean_xe=as.numeric(args[3]); # GxE power simulations
sd_xe=as.numeric(args[4]); # GxE power simulations
prop_xe=as.numeric(args[5]); # GxE power simulations
select=args[6] # SELECT:yes/no
nsims=as.numeric(args[7]) # number of simulations

K_param=as.numeric(args[8]) # number of splits K
c1=as.numeric(args[9]); c2=as.numeric(args[10]); c3=1-c1-c2; # determine splitting fractions (c_1, c_2, c_3)
# c_1 screening, c_2 main effects, c_3 ACE algorithm

########################################################################
### setup

n_rep=nsims
n=30000
n_var=100
n_e=5
n_z=2
#
sd_x=0.1
sd_e=0.1
sd_z=0.1
#
sd_err_e=0.0
sd_e_x=0.0
sd_e_z=0.1
sd_eps=1.0
#
beta_PS=1.0

prop_x=0.5 
prop_e_x=0.5
prop_e_z=0.5
#
min=0.0
max=0.1
maf=0.3


##########################################################################
seed=as.numeric(seed_num)

if(scenario==1){mem=0; nne=FALSE}
if(scenario==2){mem=1; sd_e=0.3; nne=FALSE} 
if(scenario==3){mem=1; sd_err_e=0.5; nne=TRUE}
if(scenario==4){sd_e_x=0.1; mem=1; nne=FALSE}
if(scenario==5){sd_e_x=0.1; mem=1; sd_err_e=0.2; nne=TRUE}


seed=seed+as.numeric(scenario)*1000
##########################################################################
### simulation and analysis

set.seed(seed)


results=matrix(0, nrow=n_rep, ncol=10) 
for(i in 1:n_rep)
{
        ### draw parameters and simulate phenotype
        beta_x=draw_normal_effects(mean=0, sd=sd_x, m=n_var, proportion=prop_x, sign=FALSE)
        beta_e=draw_normal_effects(mean=0, sd=sd_e, m=n_e, proportion=1, sign=TRUE)
		beta_z=draw_normal_effects(mean=0, sd=sd_z, m=n_z, proportion=1, sign=FALSE)
		
		beta_xe=draw_normal_effects_for_gxe(mean=mean_xe, sd=sd_xe, m=n_var, proportion=prop_xe, main_effects=beta_x, sign=TRUE)
		
		beta_err_e=draw_normal_effects(mean=0, sd=sd_err_e, m=n_e, proportion=1, sign=TRUE)
		beta_e_x=draw_normal_effects_matrix(mean=0, sd=sd_e_x, m=n_var, ncol=n_e, proportion=prop_e_x, sign=TRUE)
		beta_e_z=draw_normal_effects_matrix(mean=0, sd=sd_e_z, m=n_z, ncol=n_e, proportion=prop_e_z, sign=FALSE)
		
		z=draw_z(min=min, max=max, n_z=n_z, n=n)
		x=draw_x_given_z(z=z, n_var=n_var, beta_PS=beta_PS, maf=maf)
		e=draw_e_given_x_and_z(z=z, x=x, n_e=n_e, beta_e_x=beta_e_x, beta_e_z=beta_e_z) 
        env_e=get_environmental_effect(e=e, mem=mem, beta_e=beta_e)
		env_z=get_covariate_effect(z=z, beta_z=beta_z) 
		genetic=get_genetic_effect(x=x, beta_x=beta_x)
		gxe=get_gxe_effect(x=x, e=e, beta_xe=beta_xe, beta_x=beta_x) 
		eps=get_error(n=n, sd=sd_eps, beta_err_e=beta_err_e, nne=nne, errors=errors)

		y=genetic + gxe + env_e + env_z + eps
		y=as.numeric(y)
		
		###
		if(select=="yes"){
		   index_select=n_var/2
		   inds=select_x(y=y, x=x, e=e, z=z, index=index_select)
		   x=x[,inds]
		}
		###
		
		####################################################
		### RITSS
		
		indices=create_splits(K=K_param, n=n, split_ratio=c(c1, c2, c3))
		
		ritss_obj1=ritss(y=y, x=x, e=e, z=z, indices=indices, cut_off_p_value=config$cut_off_p_value, screening_function=screening_subprs, indices_env_factors=1, verbose=config$sim_verbose)
		ritss_obj2=ritss(y=y, x=x, e=e, z=z, indices=indices, cut_off_p_value=config$cut_off_p_value, screening_function=screening_sv, indices_env_factors=1, verbose=config$sim_verbose)

		
		results[i,1]=ritss_obj1$pval; 
		results[i,2]=ritss_obj2$pval;  
		
		### GESAT
		results[i,3]=gesat(x=x, e=e, z=z, y=y)
		
		### direct approach
		results[i,4]=direct_test(ritss_obj=ritss_obj1)
		results[i,5]=direct_test(ritss_obj=ritss_obj2)
		
		### GAM single variant
		results[i,6]=min(gam_test_sv(y=y, x=x, e=e, z=z, e_index=1))
	
		############################################################################################
		results[i,7]=results[i,8]=results[i,9]=results[i,10]=-9
		if(select=="no")
		{
			results[i,7]=sum(ritss_obj1$genetic_signals[abs(beta_xe)>0]>(K_param-1))/sum(abs(beta_xe>0))
			results[i,8]=sum(ritss_obj1$genetic_signals[abs(beta_xe)==0]>(K_param-1))/sum(abs(beta_xe==0))
			results[i,9]=sum(ritss_obj2$genetic_signals[abs(beta_xe)>0]>(K_param-1))/sum(abs(beta_xe>0))
			results[i,10]=sum(ritss_obj2$genetic_signals[abs(beta_xe)==0]>(K_param-1))/sum(abs(beta_xe==0))
		}
		
	
		cat(results[i,], "\n")
}


