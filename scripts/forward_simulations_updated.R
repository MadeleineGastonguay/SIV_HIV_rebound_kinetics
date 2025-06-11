# This script simulates HIV trajectories under various treatment combinations ####
# We use the population means and random effects estimated from human data and apply
# the fold-change treatment effects estiamted from SIV or SHIV data
# We sample 200 sets of population parameters (means, standard deviation of random effects, and treatment effects) from a multivariate
# normal distribution with a varcovar matrix output from monolix
### Note: although we consider correlation between population parameters determined from the human data and correlation between treatment effects
### determined from the SIV/SHIV data, we do not account for correlations between population parameters and treatment effects
# For each set of sampled population parameters, we sample the parameters for an individual from a normal distribution
# with the population mean and standard deviation


require(ggplot2)
require(deSolve)
require(MASS)
require(dplyr)
require(lamW)
require(scales)
require(rlist)
require(grid)
require(here)
require(patchwork)

source(here("scripts", "forward_simulation_functions.R"))



# Define fixed parameter values ####
# Some values depend on species

fixed <- list(
  N_E = 1e4, f = 0.9, dTar = 0.05, dI_macaque = 0.4, dI_human = 1, dP = 1e-3, 
  dE = 1, k_human = 10000, k_monkey = 50000, c = 23, final_time=200, 
  I1_human = 4 * 10 ^ (-6), I1_monkey = 4 * 10 ^ (-5),  C =.7, rho=10,
  detection_threshold=log10(50), v_one_infected = 0.0724, reactivation_rate = 1e-3, 
  LR_decay = 2.7/365, alpha_factor = 1e-8, fano_factor = 10, whole_body_conversion=3*10^4, washout_time = 3
)


N_E <- fixed$N_E #this is about what z approaches, in the dynamics between immune system and virus
#so we estimate at 1% of 10^6 CD4 cells per ml
f <- fixed$f
dTar <- fixed$dTar #target
dI <- fixed$dI_macaque #infected
dP <- fixed$dP #long lived
dE <- fixed$dE #short memory
k <- fixed$k_human #also k/c
c <- fixed$c
I1 <- fixed$I1_human
C <- fixed$C
final_time <- fixed$final_time
rho <- fixed$rho
washout_time <- fixed$washout_time



# Parse monolix output ####
results_folder <- here("results", "summary_single_fit")

parameter_file <- "SIV_SHIV_population_parameters.csv"
TLR7_varcov_file <- "study1_covarianceEstimatesSA.txt"
TLR7vac_varcov_file <- "study2_covarianceEstimatesSA.txt"
TLR7ab_varcov_file <- "study3_covarianceEstimatesSA.txt"

Human_params_file <- "human_population_parameters.csv"
Human_varcov_file <- "human_covarianceEstimatesSA.txt"

## Read in parameters 
parameters <- read_csv(here(results_folder, parameter_file))

## Results from TLR7 study
TLR7_pars <- parameters %>% filter(data == "study1") %>% select(-data) %>% column_to_rownames("parameter")

TLR7_varcov_mat <- read.csv(here(results_folder, TLR7_varcov_file), header = FALSE, row.names = 1)
colnames(TLR7_varcov_mat) <- rownames(TLR7_varcov_mat)

# Log-transform the population parameters 
for(row_name in rownames(TLR7_pars)){
  if(grepl("pop", row_name)){
    TLR7_pars[row_name, "value"] <- log(TLR7_pars[row_name, "value"])
  }
}

## Results from TLR7+Vac study
TLR7vac_pars <- parameters %>% filter(data == "study2") %>% select(-data) %>% column_to_rownames("parameter")

TLR7vac_varcov_mat <- read.csv(here(results_folder, TLR7vac_varcov_file), header = FALSE, row.names = 1)
colnames(TLR7vac_varcov_mat) <- rownames(TLR7vac_varcov_mat)

# Log-transform the population parameters 
for(row_name in rownames(TLR7vac_pars)){
  if(grepl("pop", row_name)){
    TLR7vac_pars[row_name, "value"] <- log(TLR7vac_pars[row_name, "value"])
  }
}

## Results from TLR7+Ab study
TLR7ab_pars <- parameters %>% filter(data == "study3") %>% select(-data) %>% column_to_rownames("parameter")

TLR7ab_varcov_mat <- read.csv(here(results_folder, TLR7ab_varcov_file), header = FALSE, row.names = 1)
colnames(TLR7ab_varcov_mat) <- rownames(TLR7ab_varcov_mat)

for(row_name in rownames(TLR7ab_pars)){
  if(grepl("pop", row_name)){
    TLR7ab_pars[row_name, "value"] <- log(TLR7ab_pars[row_name, "value"])
  }
}


## Results from human study
Human_pars <- read_csv(here(results_folder, Human_params_file)) %>% select(-data) %>% column_to_rownames("parameter")
Human_varcov_mat <- read.csv(here(results_folder, Human_varcov_file), header = FALSE, row.names = 1)
colnames(Human_varcov_mat) <- rownames(Human_varcov_mat)

# Log-transform the parameters (except treatment effects)
for(row_name in rownames(Human_pars)){
  if(grepl("pop", row_name)){
    Human_pars[row_name, "value"] <- log(Human_pars[row_name, "value"])
  }
}


# Construct covariate matrices for each study ####

#HUMAN MAIN TEXT FIGURE
#CHANGE I1, K, alternative scenarios
#Note that removal of study effect on NP baseline is handled later
#cov vec list is a vector with the covariates and treament_list is a list of the "treatment" set, one of "SIV", "SHIV", or "Both" 
#Cov vec list structure is (tlr7, vac, stud, ab)
cov_vec_list <- list(
  baseline = c(0, 0, 0, 0), # control
  tlr7_whitney = c(1, 0, 1, 0), # TLR7, study1
  tlr7_vac_barouch = c(1, 1, 0, 0), #TLR7 and Vac, study2
  antibody_tlr7=c(1, 0, 0, 1) # TLR7 and Ab
)
treatment_list <- list(baseline = "SIV", tlr7_whitney = "SIV", tlr7_vac_barouch = "SIV", antibody_tlr7 = "SHIV")
VL_cutoff <- log10(50)
base_organism <- "human"
#Human main
names_list <- c("No Treatment", "TLR7 (S1)", "TLR7 + Vaccine (S2)", "TLR7 + Antibody (S3)")



if(base_organism == "human"){
  I1 <- fixed$I1_human
  k <- fixed$k_human
  dI <- fixed$dI_human
}else{
  I1 <- fixed$I1_monkey
  k <- fixed$k_monkey
  dI <- fixed$dI_macaque
}


# Define parameters of each model (excluding random effects)
human_par_list <- rownames(Human_varcov_mat)
human_par_list <- human_par_list[!grepl("omega", human_par_list)]

TLR7_par_list <- rownames(TLR7_varcov_mat)
TLR7_par_list <- TLR7_par_list[!grepl("omega", TLR7_par_list)]

TLR7vac_par_list <- rownames(TLR7vac_varcov_mat)
TLR7vac_par_list <- TLR7vac_par_list[!grepl("omega", TLR7vac_par_list)]

TLR7ab_par_list <- rownames(TLR7ab_varcov_mat)
TLR7ab_par_list <- TLR7ab_par_list[!grepl("omega", TLR7ab_par_list)]




# Simulate many trajectories with varying population parameter values ####
num_its <- 200 # 200 trajectories
final_df <- data.frame()
quad_df <- data.frame()
final_df_data <- data.frame()
trajectory_record_final <- data.frame()

# Iterate through four scenarios: control, TLR7 in study 1, TLR7+Vac in study 2, and TLR7 + Ab in study 3
for(group_idx in 1:length(cov_vec_list)){ 
  group_name <- names(cov_vec_list)[group_idx]
  # Define covariates for scenario
  cov_vec_tmp <- cov_vec_list[[group_idx]]
  tmp_df <- data.frame()
  trajectory_record_tmp <- data.frame()
  for(i in 1:num_its){
    print(i)
    
    #### sample model params from multivariate normal distribution
    # TODO # MADDIE - WHY ARE OMEGAS NOT SAMPLED?
    fit_results_scaled_human <- sample_fitting_results_scaled(Human_pars[human_par_list, "value"], Human_varcov_mat[human_par_list, human_par_list])
    fit_results_scaled_tlr7 <- sample_fitting_results_scaled(TLR7_pars[TLR7_par_list, "value"], TLR7_varcov_mat[TLR7_par_list, TLR7_par_list])
    fit_results_scaled_tlr7vac <- sample_fitting_results_scaled(TLR7vac_pars[TLR7vac_par_list, "value"], TLR7vac_varcov_mat[TLR7vac_par_list, TLR7vac_par_list])
    fit_results_scaled_tlr7ab <- sample_fitting_results_scaled(TLR7ab_pars[TLR7ab_par_list, "value"], TLR7ab_varcov_mat[TLR7ab_par_list, TLR7ab_par_list])
    
    #### Format sampled treatment effects from SIV studies
    treatment_effects_tlr7 <- make_covariate_matrix_list(fit_results_scaled_tlr7)
    treatment_effects_tlr7vac <- make_covariate_matrix_list(fit_results_scaled_tlr7vac)
    treatment_effects_tlr7ab <- make_covariate_matrix_list(fit_results_scaled_tlr7ab)
    # treatment_effects_combined <- mapply(function(x, y) x + y, treatment_effects_tlr7vac, treatment_effects_tlr7ab, SIMPLIFY = F)
   
    treatment_effects_combined <- c()
    # if(treatment_list[group_idx] == "SIV"){
    #   treatment_effects_combined <- treatment_effects_tlr7vac
    # }
    # if(treatment_list[group_idx] == "SHIV"){
    #   treatment_effects_combined <- treatment_effects_tlr7ab
    # }
    # if(treatment_list[group_idx] == "Both"){
    #   treatment_effects_combined <- mapply(function(x, y) x + y, treatment_effects_tlr7vac, treatment_effects_tlr7ab, SIMPLIFY = F)
    # }
    if(group_name == "baseline"){
      treatment_effects_combined <- c()
    }else if(group_name == "tlr7_whitney"){
      treatment_effects_combined <- treatment_effects_tlr7
    }else if(group_name == "tlr7_vac_barouch"){
      treatment_effects_combined <- treatment_effects_tlr7vac
    }else if(group_name == "antibody_tlr7"){
      treatment_effects_combined <- treatment_effects_tlr7ab
    }else if(greoup_name == "both"){
      treatment_effects_combined <- mapply(function(x, y, z) x + y + z, treatment_effects_tlr7, treatment_effects_tlr7vac, treatment_effects_tlr7ab, SIMPLIFY = F)
    }
    
    omega_param_names <- c('omega_lambda', 'omega_betas', 'omega_ta', 'omega_m', 'omega_p', 'omega_NPs')
    
    if(base_organism == "human"){
      human_omega_vals <- Human_pars[omega_param_names, "value"]
      names(human_omega_vals) <- omega_param_names
      human_vals_tmp <- c(fit_results_scaled_human, human_omega_vals)
      # mean_info <- get_mean_params(human_vals_tmp)
      ## Add in means for fixed parameters:
      fixed_param <- rownames(Human_pars)[!rownames(Human_pars) %in% c(human_par_list, omega_param_names)]
      fixed_values <- Human_pars[fixed_param, "value"]
      names(fixed_values) <- fixed_param
      human_vals_tmp <- c(human_vals_tmp, fixed_values)
      mean_info <- get_mean_params(human_vals_tmp)
      if(group_name != "baseline") treatment_effects_combined$NP[3,3] <- 0 #no study effect on NP for humans
    }
    
    # if(base_organism == "macaque_tlr7vac"){
    #   TLR7vac_omega_vals <- TLR7vac_pars[omega_param_names, "value"]
    #   names(TLR7vac_omega_vals) <- omega_param_names
    #   TLR7vac_vals_tmp <- c(fit_results_scaled_tlr7vac, TLR7vac_omega_vals)
    #   mean_info <- get_mean_params(TLR7vac_vals_tmp)
    #   
    # }
    # 
    # if(base_organism == "macaque_tlr7ab"){
    #   TLR7ab_omega_vals <- TLR7ab_pars[omega_param_names, "value"]
    #   names(TLR7ab_omega_vals) <- omega_param_names
    #   TLR7ab_vals_tmp <- c(fit_results_scaled_tlr7ab, TLR7ab_omega_vals)
    #   mean_info <- get_mean_params(TLR7ab_vals_tmp)
    #   treatment_effects_combined$NP[1, 2] <- 0 #removed to prefer the other tlr7 effect
    # }
    # 
    # Sample parameters for an individual, given the population mean and standard deviation
    baseline <- sample_baseline_parameters_ls(1, mean_info)
    treatment_effects_tmp <- treatment_effects_combined
    

    #########################
    # Add treatment effects to parameters
    added_cvs <- add_covariate_effects(baseline, cov_vec_tmp, treatment_effects_tmp)
    
    added_cvs_log10 <- added_cvs / log(10)
    characteristics_results <- sim_traj_final(10^added_cvs_log10, return_full = T)
    characteristics_df <- characteristics_results[[1]]
    trajectory_df <- data.frame(characteristics_results[[2]])
    trajectory_df$id <- i
    characteristics_df$peak <- characteristics_df$peak/log(10)
    characteristics_df$final <- characteristics_df$final/log(10)
    trajectory_df$final <- characteristics_df$final
    characteristics_df$name <- group_name
    
    trajectory_record_tmp <- rbind(trajectory_record_tmp, trajectory_df)
    
    tmp_df <- rbind(tmp_df, characteristics_df)
  }
  #tmp_df <- tmp_df[complete.cases(tmp_df),]
  print(sum(complete.cases(tmp_df)))
  if(sum(complete.cases(tmp_df)) != num_its){
    print(tmp_df)
  }
  dens <- MASS::kde2d(tmp_df$peak, tmp_df$final, n = 10, lims = c(c(-3, 8), c(-3, 8)))
  densf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(dens$z))
  
  quad_df_tmp <- data.frame(quad_1 = round(sum(tmp_df$final >= VL_cutoff & tmp_df$peak >= VL_cutoff)/num_its, 2), 
                            quad_2 = round(sum(tmp_df$final < VL_cutoff & tmp_df$peak >= VL_cutoff)/num_its, 2), 
                            quad_3 = round(sum(tmp_df$final < VL_cutoff & tmp_df$peak < VL_cutoff)/num_its, 2),
                            quad_4 = round(sum(tmp_df$final >= VL_cutoff & tmp_df$peak < VL_cutoff)/num_its, 2),
                            name=group_name)
  quad_df <- rbind(quad_df, quad_df_tmp)
  
  #densf$z <- densf$z / max(densf$z) #normalize to same max
  densf$name <- group_name
  dx <- diff(dens$x[1:2])
  dy <- diff(dens$y[1:2])
  sz <- sort(dens$z)
  c1 <- cumsum(sz) * dx * dy
  densf$c1 <- c1
  densf$sz <- sz
  final_df <- rbind(final_df, densf)
  tmp_df$name <- group_name
  trajectory_record_tmp$name <- group_name
  
  final_df_data <- rbind(final_df_data, tmp_df)
  trajectory_record_final <- rbind(trajectory_record_final, trajectory_record_tmp)
}


prob <- c(.95, .90, .8, .5, .25)
plot(final_df_data$peak, final_df_data$final)
#contour(dens, levels=levels, labels=prob, add=T)

#dens <- MASS::kde2d(tmp_df$final, tmp_df$peak, n = 100)
#densf <- data.frame(expand.grid(x=dens$x, y=dens$y), z=as.vector(denz$z))


final_df_data$name_f <- factor(final_df_data$name, levels=names(cov_vec_list))
levels(final_df_data$name_f)
levels(final_df_data$name_f) <- names_list 

quad_df$name_f <- factor(quad_df$name, levels=names(cov_vec_list))
levels(quad_df$name_f)
levels(quad_df$name_f) <- names_list


trajectory_record_final$name_f <- factor(trajectory_record_final$name, levels=names(cov_vec_list))
levels(trajectory_record_final$name_f)
levels(trajectory_record_final$name_f) <- names_list


all_data_to_save <- list(trajectory_record_final = trajectory_record_final, final_df_data = final_df_data, quad_df=quad_df, cov_vec_list = cov_vec_list)

list.save(all_data_to_save, save_file_name)





############################ Read in data  #########################################

saved_data <- readRDS("human_main_wo_6_10_19.rds") #shiv_full.rds, human_main_full_wo.rds, "human_sup_full_wo.rds"
trajectory_record_final <- saved_data[[1]]
final_df_data <- saved_data[[2]]
quad_df <- saved_data[[3]]
cov_vec_list <- saved_data[[4]]

trajectory_record_final_cut <- trajectory_record_final[trajectory_record_final$id <= 25,]


#ggplot(data=final_df_data, aes(x=p_surv)) + geom_histogram() + facet_wrap(~name_f)
######


######### SIV ##############



#color by rebound time
m2_siv <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_wrap(~ name_f) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb+washout_time))) + scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
  theme_bw() +  
  xlab("Peak viral load (copies/mL)") + ylab("Set point viral load (copies/mL)") + coord_cartesian(xlim=c(-2.5, 6.5), ylim=c(-2.5, 6.5)) + 
  geom_vline(xintercept = VL_cutoff, linetype='dashed') + geom_hline(yintercept = VL_cutoff, linetype='dashed') + geom_text(x= 5, y = 6, aes(label = quad_1), data=quad_df, color='black') + 
  geom_text(x= 5, y = 0, aes(label = quad_2), data=quad_df, color='black') + geom_text(x= -1.5, y = 0, aes(label = quad_3), data=quad_df, color='black') + 
  scale_x_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  scale_y_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 
#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m2_siv)


ggsave("macaque_siv_predictions_scatter_6_10_19.pdf", plot=m2_siv, width=20, height=20, unit='cm')



m3_siv <- ggplot(trajectory_record_final_cut, aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap(~ name_f) + 
  geom_line(data = filter(trajectory_record_final_cut, pre_activated == "yes"), linetype="solid", size=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut, post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", trans = trim_tails(range = c(0,5)), breaks=c(0, 1, 2, 3, 4, 5), labels=c(expression(10^{0}), expression(10^{1}),expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m3_siv)


ggsave("macaque_siv_example_simulations_6_10_19.pdf", plot=m3_siv, width=20, height=20, unit='cm')


######### END SIV ############

######### SHIV #############




#color by rebound time
m2_shiv <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_wrap(~ name_f) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb+washout_time))) + scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
  theme_bw() +  
  xlab("Peak viral load (copies/mL)") + ylab("Set point viral load (copies/mL)") + coord_cartesian(xlim=c(-2.5, 6.5), ylim=c(-2.5, 6.5)) + 
  geom_vline(xintercept = VL_cutoff, linetype='dashed') + geom_hline(yintercept = VL_cutoff, linetype='dashed') + geom_text(x= 5, y = 6, aes(label = quad_1), data=quad_df, color='black') + 
  geom_text(x= 5, y = 0, aes(label = quad_2), data=quad_df, color='black') + geom_text(x= -1.5, y = 0, aes(label = quad_3), data=quad_df, color='black') + 
  scale_x_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  scale_y_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 
#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m2_shiv)


ggsave("macaque_shiv_predictions_scatter_6_10_19.pdf", plot=m2_shiv, width=20, height=20, unit='cm')


m3_shiv <- ggplot(trajectory_record_final_cut, aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap(~ name_f) + 
  geom_line(data = filter(trajectory_record_final_cut, pre_activated == "yes"), linetype="solid", size=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut, post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", trans = trim_tails(range = c(0,5)), breaks=c(0, 1, 2, 3, 4, 5), limits=c(0,5), labels=c(expression(10^{0}), expression(10^{1}),expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m3_shiv)


ggsave("macaque_shiv_example_simulations_6_10_19.pdf", plot=m3_shiv, width=20, height=20, unit='cm')


######### END SHIV ###########



######### HUMAN ##############




#color by rebound time
m2 <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_wrap(~ name_f, ncol=4) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb+washout_time))) + 
  scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, 
                       trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), 
                       labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
  theme_bw() +  
  xlab("Peak viral load (copies/mL)") + ylab("Set point viral load (copies/mL)") + coord_cartesian(xlim=c(-5, 7), ylim=c(-5, 7)) + 
  geom_vline(xintercept = VL_cutoff, linetype='dashed') + geom_hline(yintercept = VL_cutoff, linetype='dashed') + 
  geom_text(x= 5, y = 6, aes(label = quad_1), data=quad_df, color='black') + 
  geom_text(x= 5, y = 0, aes(label = quad_2), data=quad_df, color='black') + 
  geom_text(x= -1.5, y = 0, aes(label = quad_3), data=quad_df, color='black') + 
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6), labels=c(expression(10^{-4}),expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}),  expression(10^{6}))) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 
#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m2)

ggsave(here(out_folder, "human_predictions_scatter_4_14_25.pdf"), plot=m2,  width=25, heigh=10,  unit='cm')

m3 <- ggplot(trajectory_record_final_cut, aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap(~ name_f, ncol=4) + 
  geom_line(data = filter(trajectory_record_final_cut, pre_activated == "yes"), linetype="solid", size=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut, post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", 
                        # trans = trim_tails(range = c(-3,5))) +
                        breaks=c(-3,-1,1,3,5), labels=c(expression(10^{-3}), expression(10^{-1}), expression(10^{1}), expression(10^{3}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

show(m3)

ggsave(here(out_folder, "human_example_simulations_4_14_25.pdf"), plot=m3, width=25, height=10, unit='cm')


# Plot distribution of rebound times
m4 <- final_df_data %>% 
  filter(t_reb != Inf) %>%
  ggplot(aes((t_reb+washout_time))) + 
  geom_histogram(aes(y = after_stat(count/sum(count))), fill = "darkgray") + facet_wrap(~name_f, nrow =1) + 
  labs(x = "Rebound time", y = "Frequency") +
  theme_bw() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_log10(breaks = c((3), (7), (30), (90), (365)), 
                 labels = c("3 days", "1 week", "1 month", "3 months", "1 year"))  





# testplot <- rbind(m2, m3, m4)
# show(testplot)
# 
# gA <- ggplotGrob(m2)
# gB <- ggplotGrob(m3)
# gC <- ggplotGrob(m4)
# g <- rbind(gB, gA, gC, size="last")
# grid.newpage()
# grid.draw(g)

g <- m3 + m2  & plot_layout(ncol = 1)


ggsave(here(out_folder, "human_simulations_main_fixed_4_14_25.pdf"), width=25, height=18, units='cm', plot=g)



######HUMAN SUPPLEMENT ##########




#color by rebound time
m2 <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_wrap(~ name_f) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb+washout_time))) + scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
  theme_bw() +  
  xlab("Peak viral load (copies/mL)") + ylab("Set point viral load (copies/mL)") + coord_cartesian(xlim=c(-2.5, 6.5), ylim=c(-2.5, 6.5)) + 
  geom_vline(xintercept = VL_cutoff, linetype='dashed') + geom_hline(yintercept = VL_cutoff, linetype='dashed') + geom_text(x= 5, y = 6, aes(label = quad_1), data=quad_df, color='black') + 
  geom_text(x= 5, y = 0, aes(label = quad_2), data=quad_df, color='black') + geom_text(x= -1.5, y = 0, aes(label = quad_3), data=quad_df, color='black') + 
  scale_x_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  scale_y_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 
#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m2)

ggsave("human_predictions_scatter_6_10_19_supplemental.pdf", plot=m2, width=20, height=20, unit='cm')

m3 <- ggplot(trajectory_record_final_cut, aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap(~ name_f) + 
  geom_line(data = filter(trajectory_record_final_cut, pre_activated == "yes"), linetype="solid", size=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut, post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", trans = trim_tails(range = c(0,5)), breaks=c(0, 1, 2, 3, 4, 5), labels=c(expression(10^{0}), expression(10^{1}),expression(10^{2}), expression(10^{3}), expression(10^{4}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m3)




ggsave("human_example_simulations_6_10_19_supplemental.pdf", plot=m3, width=20, height=20, unit='cm')




######### END HUMAN ############



#Initial conditions
###testing

#what to plot? Peak and end of study viral load quadrants, 100 and 100 copies/mL
#just 3 possibilities 
#do all combinations of treatments, incl in humans
#vary each according to the 3 metrics



results <- read.table("fit_results/All_files/TLR7VAC_Indiv.txt", sep=",", stringsAsFactors = F, header=T)
params_test <- c(results[15,'lambda_SAEM',], 10^(-6) * results[15,'betas_SAEM',],  results[15,'ta_SAEM',],  results[15,'m_SAEM',],  results[15,'p_SAEM',],  results[15,'NPs_SAEM',] * 10^5)
params_test_log10 <- log10(params_test)
out_test <- sim_traj(params_test_log10)

plot(test2[,'time'], test2[,'V'], ylim=c(-4, 7))


parameters_log10 <- c(2.17099638, -5.62086757,  0.68748324,  2.52564814,  0.06567356,  5.12111436)

test <- sim_traj_v2(parameters_log10)
plot(test[,'time'], log10(test[,'V']))







########### PLOT HUMAN REACTIVATION TIMES INDIVIDUALLY #########
#Plot : time to cross VL 50 (w/out washout)
#Plot: distribution of r values
#Plot: V_ART for each individual

#waht are the parameter covariance estimates here?


human_individual_params <- paste(c(results_folder, "Humans_omega_pfixed_estimatedIndividualParameters.txt"), sep="", collapse="")
nnrti_file <- paste(c(results_folder, "chronic_rebound_nnrti.txt"), sep="", collapse="")

par_table <- read.table(human_individual_params, stringsAsFactors = F, header=T, sep=",")
nnrti_table <- read.table(nnrti_file, stringsAsFactors = F, header=T, sep="\t")
df1 <- data.frame(par_table)
df2 <- data.frame(nnrti_table)

df_merged <- merge(df1, df2, by="id")


ta_over_p <- c()
R0_record <- c()
stats_df <- data.frame()
for(i in 1:nrow(df_merged)){
  
  params_tmp <- c(df_merged[i,'lambda_SAEM',], 10^(-7) * df_merged[i,'betas_SAEM',],  df_merged[i,'ta_SAEM',],  df_merged[i,'m_SAEM',],  df_merged[i,'p_SAEM',],  df_merged[i,'NPs_SAEM',] * 10^5)
  ta_tmp <- params_tmp[3]
  
  characteristics_results <- sim_traj_final(params_tmp, return_full = T)
  characteristics_df <- characteristics_results[[1]]
  trajectory_df <- data.frame(characteristics_results[[2]])
  
  stats_df <- rbind(stats_df, characteristics_df)
  
  res_list <- calc_p_surv_R0(params_tmp)
  p_tmp <- res_list[1]
  R0_tmp <- res_list[2]
  ta_over_p <- c(ta_over_p, ta_tmp/p_tmp)
  R0_record <- c(R0_record, R0_tmp)
}

df_merged <- cbind(df_merged, stats_df)

df_merged$ta_over_p <- ta_over_p
df_merged$R0 <- R0_record

ggplot(data=df_merged, aes(x = ta_over_p, group = factor(nnrti), fill=factor(nnrti))) + geom_histogram(breaks = seq(from=0, to=200, by = 1), alpha = 0.4, position = "identity") + coord_cartesian(xlim=c(0, 30))

alison_plot <- ggplot(data=df_merged, aes(x=t_reb + washout_time, group=factor(nnrti), fill=factor(nnrti))) + geom_histogram(position="identity", alpha=0.4) + theme_bw()
show(alison_plot)


alison_plot2 <- ggplot(data=df_merged, aes(x=nnrti, y=t_reb + washout_time, group=factor(nnrti), fill=factor(nnrti))) + geom_boxplot() + theme_bw()
show(alison_plot2)


alison_plot3 <- ggplot(data=df_merged, aes(x=nnrti, y=r, group=factor(nnrti), fill=factor(nnrti))) + geom_boxplot() + theme_bw()
show(alison_plot3)


alison_plot4 <- ggplot(data=df_merged, aes(x=r, group=factor(nnrti), fill=factor(nnrti))) +  geom_histogram(position="identity", alpha=0.4) + theme_bw()
show(alison_plot4)


alison_plot5 <- ggplot(data=df_merged, aes(x=log10(V_ART), group=factor(nnrti), fill=factor(nnrti))) +  geom_histogram(position="identity", alpha=0.4) + theme_bw()
show(alison_plot5)


ggsave("nnrti_plots_5.pdf", width=6, height=6, plot=alison_plot5)


p1_nnrti <- ggplot(data=df_merged, aes(x = log10(ta_over_p), group = factor(nnrti), fill=factor(nnrti))) + geom_histogram(alpha = 0.4, position = "identity") + 
  theme_bw()

p3_nnrti <- ggplot(data=df_merged, aes(x = log10(R0), group = factor(nnrti), fill=factor(nnrti))) + geom_histogram(alpha = 0.4, position = "identity") + 
  theme_bw()



p2_nnrti <- ggplot(data=df_merged, aes(x = nnrti, y = log10(ta_over_p), color=factor(nnrti))) + geom_boxplot() + geom_hline(yintercept = log10(1/dI)) + 
  theme_bw()

p4_nnrti <- ggplot(data=df_merged, aes(x=log10(ta_over_p), y=R0, color=factor(nnrti))) + geom_point() + theme_bw()

show(p4_nnrti)

#ggsave("nnrti_histogram.pdf", plot=p1_nnrti, width=6, height=6)

#ggsave("nnrti_boxplot.pdf", plot=p2_nnrti, width=6, height=6)
ggsave("nnrti_R0_scatter.pdf", plot=p4_nnrti, width=6, height=6)

sums <- tapply(df_merged$ta_over_p, df_merged$nnrti, summary)
#difference in medians in 0.45

sums2 <- tapply(df_merged$R0, df_merged$nnrti, summary)

