# This script simulates HIV trajectories under various treatment combinations ####
# We use the population means and random effects estimated from human data and apply
# the fold-change treatment effects estimated from SIV or SHIV data
# We sample 200 sets of population parameters (means, standard deviation of random effects, and treatment effects) from a multivariate
# normal distribution with a varcovar matrix output from monolix
# We reject any samples with negative standard deviations, effectively sampling from a truncated normal distribution
### Note: although we consider correlation between population parameters determined from the human data and correlation between treatment effects
### determined from the SIV/SHIV data, we do not account for correlations between population parameters and treatment effects
# For each set of sampled population parameters, we sample the parameters for an individual from a normal distribution
# with the population mean and standard deviation

require(tidyverse)
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
select <- dplyr::select

source(here("scripts", "forward_simulation_functions.R"))



# Parse monolix output ####
results_folder <- here("results_eclipse", "summary_single_fit", "study4_comparison")

tlr7beta_varcov <- here("results_eclipse", "study4", "best_model", "FisherInformation", "covarianceEstimatesSA.txt")
tlr7lambda_varcov <- here("results_eclipse", "study4", "best_4.6", "FisherInformation", "covarianceEstimatesSA.txt")


## Format parameters for simulations ####
## Read in parameters 
tlr7beta_parameters <- read_csv(here("results_eclipse", "study4", "best_model", "populationParameters.txt")) %>% 
  column_to_rownames("parameter")
tlr7lambda_parameters <- read_csv(here("results_eclipse", "study4", "best_4.6", "populationParameters.txt")) %>% 
  column_to_rownames("parameter")


TLR7beta_varcov_mat <- read.csv(tlr7beta_varcov, header = FALSE, row.names = 1)
colnames(TLR7beta_varcov_mat) <- rownames(TLR7beta_varcov_mat)

# Log-transform the population parameters 
for(row_name in rownames(tlr7beta_parameters)){
  if(grepl("pop", row_name)){
    tlr7beta_parameters[row_name, "value"] <- log(tlr7beta_parameters[row_name, "value"])
  }
}


TLR7lambda_varcov_mat <- read.csv(tlr7lambda_varcov, header = FALSE, row.names = 1)
colnames(TLR7lambda_varcov_mat) <- rownames(TLR7lambda_varcov_mat)

# Log-transform the population parameters 
for(row_name in rownames(tlr7lambda_parameters)){
  if(grepl("pop", row_name)){
    tlr7lambda_parameters[row_name, "value"] <- log(tlr7lambda_parameters[row_name, "value"])
  }
}


TLR7beta_par_list <- rownames(TLR7beta_varcov_mat)
TLR7lambda_par_list <- rownames(TLR7lambda_varcov_mat)


# Construct covariate matrices for each study ####

#HUMAN MAIN TEXT FIGURE
#CHANGE I1, d1, alternative scenarios
#cov vec list is a vector with the covariates and treament_list is a list of the "treatment" set, one of "SIV", "SHIV", or "Both" 
#Cov vec list structure is (tlr7, vac, stud, ab)
cov_vec_list <- list(
  baseline = c(0, 0, 0, 0), # control
  tlr7_whitney_SIV = c(1, 0, 1, 0), # TLR7, study1
  tlr7_barouch_SIV = c(1, 0, 0, 0), #TLR7, study2
  vac_barouch_SIV = c(0, 1, 0, 0), #Vac, study2
  tlr7_vac_barouch_SIV = c(1, 1, 0, 0), #TLR7 and Vac, study2
  tlr7_barouch_SHIV=c(1, 0, 0, 0), # TLR7, study 3
  antibody_barouch_SHIV=c(0, 0, 0, 1), # Ab, study 3
  antibody_tlr7_barouch_SHIV=c(1, 0, 0, 1), # TLR7 and Ab, study 3
  antibody_vac_tlr7_walker_SHIV=c(1, 1, 0, 1), # TLR7 and Ab and Vac, study 4
  antibody_tlr7_walker_SHIV=c(1, 0, 0, 1), # TLR7 and Ab , study 4
  vac_tlr7_walker_SHIV=c(1, 1, 0, 0) # TLR7 and Vac , study 4
)

# treatment_list <- list(baseline = "SIV", tlr7_whitney = "SIV", tlr7_vac_barouch = "SIV", antibody_tlr7 = "SHIV", antibody_vac_tlr7 = "SHIV")
VL_cutoff <- log10(50)

model <- "TLR7beta" # either human, whitney_SIV, barouch_SIV, barouch_SHIV
model <- "TLR7lambda" # either human, whitney_SIV, barouch_SIV, barouch_SHIV

# #Human main
# names_list <- c("No Treatment", "TLR7 (S1)", "TLR7 + Vaccine (S2)", "TLR7 + Antibody (S3)", "TLR7 + Antibody + Vac (S4)")

# Define fixed parameter values ####
# Some values depend on species

fixed <- list(
  N_E = 1e4, f = 0.9, dTar = 0.05, dI_macaque = 0.4, dI_human = 1, dP = 1e-3, 
  dE = 1, k = 50000, c = 23, final_time=200, 
  I1_human = 4 * 10 ^ (-6), I1_monkey = 4 * 10 ^ (-5),  C =.7, rho=10,
  detection_threshold=log10(50), v_one_infected = 0.0724, reactivation_rate = 1e-3, 
  LR_decay = 2.7/365, alpha_factor = 1e-8, fano_factor = 10, whole_body_conversion=3*10^4, washout_time = 3,
  q = 1
)


N_E <- fixed$N_E #this is about what z approaches, in the dynamics between immune system and virus
#so we estimate at 1% of 10^6 CD4 cells per ml
f <- fixed$f
dTar <- fixed$dTar #target
dP <- fixed$dP #long lived
dE <- fixed$dE #short memory
# k <- fixed$k_human #also k/c
c <- fixed$c
C <- fixed$C
final_time <- fixed$final_time
rho <- fixed$rho
washout_time <- fixed$washout_time
k <- fixed$k
q <- fixed$q

I1 <- fixed$I1_monkey
dI <- fixed$dI_macaque


# Simulate many trajectories with varying population parameter values ####
num_its <- 50 #number of trajectories to simulate

virtual_trial <- function(model, num_its, cov_vec_list){
  
  # Initialize outputs ####
  final_df <- data.frame()
  quad_df <- data.frame()
  final_df_data <- data.frame()
  trajectory_record_final <- data.frame()
  parameter_df <- data.frame()
  
  
  # Sample population parameters for baseline organism: ####
  varcov_mat <- switch (model,
                        "TLR7beta" = TLR7beta_varcov_mat,
                        "TLR7lambda" = TLR7lambda_varcov_mat
  )
  par_df <- switch (model,
                    "TLR7beta" = tlr7beta_parameters,
                    "TLR7lambda" = tlr7lambda_parameters
  )
  par_list <- switch (model,
                      "TLR7beta" = TLR7beta_par_list,
                      "TLR7lambda" = TLR7lambda_par_list
  )
  
  base_organism_scaled <- matrix(nrow = num_its, ncol = nrow(par_df))
  baseline_individuals <- matrix(nrow = num_its, ncol = sum(grepl("pop", rownames(par_df))))
  for(i in 1:num_its){
    #### sample model parameters from multivariate normal distribution
    fit_results_scaled <- sample_fitting_results_scaled(par_df[par_list, "value"], varcov_mat[par_list, par_list])
    # Make sure omegas are positive:
    while(any(fit_results_scaled[grepl("omega",par_list)] < 0)){
      fit_results_scaled <- sample_fitting_results_scaled(par_df[par_list, "value"], varcov_mat[par_list, par_list])
    }
    
    ## Add in means for fixed parameters:
    fixed_param <- rownames(par_df)[!rownames(par_df) %in% par_list]
    fixed_values <- par_df[fixed_param, "value"]
    names(fixed_values) <- fixed_param
    value_tmp <- c(fit_results_scaled, fixed_values)
    
    # Sample parameters for an individual, given the population mean and standard deviation
    mean_info <- get_mean_params(value_tmp)
    baseline <- sample_baseline_parameters_ls(1, mean_info)
    
    base_organism_scaled[i,] <-  value_tmp
    baseline_individuals[i,] <- baseline
  }
  colnames(base_organism_scaled) <- names(c(fit_results_scaled, fixed_values))
  
  
  # Iterate through scenarios in cov_vec_list for treatment effects ####
  for(group_idx in 1:length(cov_vec_list)){ 
    group_name <- names(cov_vec_list)[group_idx]
    # Define covariates for scenario
    cov_vec_tmp <- cov_vec_list[[group_idx]]
    tmp_df <- data.frame()
    trajectory_record_tmp <- data.frame()
    parameters_temp <- data.frame()
    for(i in 1:num_its){
      print(i)
      
      
      # Sample treatment effects for corresponding study if it is different from the base organism:
      if(grepl(base_organism, group_name) | group_name == "baseline"){
        treatment_effects_tmp <- make_covariate_matrix_list(base_organism_scaled[i,])
      }else{
        
        treatment_study <- paste(rev(rev(str_split_1(group_name, "_"))[1:2]), collapse = "_")
        
        varcov_mat <- switch (treatment_study,
                              "whitney_SIV" = TLR7_varcov_mat,
                              "barouch_SIV" = TLR7vac_varcov_mat,
                              "barouch_SHIV" = TLR7ab_varcov_mat
        )
        par_df <- switch (treatment_study,
                          "whitney_SIV" = TLR7_pars,
                          "barouch_SIV" = TLR7vac_pars,
                          "barouch_SHIV" = TLR7ab_pars
        )
        par_list <- switch (treatment_study,
                            "whitney_SIV" = TLR7_par_list,
                            "barouch_SIV" = TLR7vac_par_list,
                            "barouch_SHIV" = TLR7ab_par_list
        )
        
        
        #### sample model parameters from multivariate normal distribution
        treatment_fit_results_scaled <- sample_fitting_results_scaled(par_df[par_list, "value"], varcov_mat[par_list, par_list])
        # Make sure omegas are positive:
        while(any(treatment_fit_results_scaled[grepl("omega",par_list)] < 0)){
          treatment_fit_results_scaled <- sample_fitting_results_scaled(par_df[par_list, "value"], varcov_mat[par_list, par_list])
        }
        
        treatment_effects_tmp <- make_covariate_matrix_list(treatment_fit_results_scaled)
      }
      
      
      
      #### Format sampled treatment effects from macaque studies
      ## The sampled effects above will be applied to the baseline parameters sampled from base organism
      ## Correlation amongst treatment effects is preserved, but correlation with population parameters is not possible to consider
      
      #########################
      # Add treatment effects to parameters
      added_cvs <- add_covariate_effects(baseline_individuals[i,], cov_vec_tmp, treatment_effects_tmp)
      
      added_cvs_log10 <- added_cvs / log(10)
      
      # cat("group: ", group_name)
      # cat("\nparam: ", 10^added_cvs_log10)
      
      characteristics_results <- sim_traj_final(10^added_cvs_log10, return_full = T, delay = T)
      characteristics_df <- characteristics_results[[1]]
      trajectory_df <- data.frame(characteristics_results[[2]])
      trajectory_df$id <- i
      characteristics_df$peak <- characteristics_df$peak/log(10)
      characteristics_df$final <- characteristics_df$final/log(10)
      trajectory_df$final <- characteristics_df$final
      characteristics_df$name <- group_name
      
      trajectory_record_tmp <- rbind(trajectory_record_tmp, trajectory_df)
      
      tmp_df <- rbind(tmp_df, characteristics_df)
      
      parameters_temp <- rbind(parameters_temp, added_cvs)
      names(parameters_temp) <- rownames(mean_info)
    }
    
    # Summarise results ####
    tmp_df2 <- tmp_df[!is.na(tmp_df$final),] # remove any simulations that didn't run because of NAN t_go
    n <- nrow(tmp_df2)
    quad_df_tmp <- data.frame(quad_1 = round(sum(tmp_df2$final >= VL_cutoff & tmp_df2$peak >= VL_cutoff)/n, 2), 
                              quad_2 = round(sum(tmp_df2$final < VL_cutoff & tmp_df2$peak >= VL_cutoff)/n, 2), 
                              quad_3 = round(sum(tmp_df2$final < VL_cutoff & tmp_df2$peak < VL_cutoff)/n, 2),
                              quad_4 = round(sum(tmp_df2$final >= VL_cutoff & tmp_df2$peak < VL_cutoff)/n, 2),
                              name=group_name)
    quad_df <- rbind(quad_df, quad_df_tmp)
    
    
    tmp_df$name <- group_name
    trajectory_record_tmp$name <- group_name
    parameters_temp$name <- group_name
    
    final_df_data <- rbind(final_df_data, tmp_df)
    trajectory_record_final <- rbind(trajectory_record_final, trajectory_record_tmp)
    parameter_df <- rbind(parameter_df, parameters_temp)
  }
  
  
  names_list <- c("No Treatment", "TLR7 (S1)", "TLR7 (S2)", "Vaccine (S2)", "TLR7 + Vaccine (S2)", "TLR7 (S3)", "Antibody (S3)", "TLR7 + Antibody (S3)")
  names(names_list) <- c("baseline", "tlr7_whitney_SIV", "tlr7_barouch_SIV", "vac_barouch_SIV", "tlr7_vac_barouch_SIV", "tlr7_barouch_SHIV", "antibody_barouch_SHIV", "antibody_tlr7_barouch_SHIV")
  names_list <- unname(names_list[names(cov_vec_list)])
  
  final_df_data$name_f <- factor(final_df_data$name, levels=names(cov_vec_list))
  levels(final_df_data$name_f) <- names_list
  
  quad_df$name_f <- factor(quad_df$name, levels=names(cov_vec_list))
  levels(quad_df$name_f) <- names_list
  
  trajectory_record_final$name_f <- factor(trajectory_record_final$name, levels=names(cov_vec_list))
  levels(trajectory_record_final$name_f) <- names_list
  
  parameter_df$name_f <- factor(parameter_df$name, levels=names(cov_vec_list))
  levels(parameter_df$name_f) <- names_list
  
  return(list(final_df_data, quad_df, trajectory_record_final, parameter_df))
}

tlr7beta_trials <- virtual_trial("TLR7beta", 100, cov_vec_list[c(1,9:11)])
tlr7lambda_trials <- virtual_trial("TLR7lambda", 100, cov_vec_list[c(1,9:11)])

combined_quad <- rbind(
  tlr7beta_trials[[2]] %>% mutate(mod = "TLR7:beta"), 
  tlr7lambda_trials[[2]] %>% mutate(mod = "TLR7:lambda")
) %>% 
  mutate(name = factor(name, levels = c("baseline", "antibody_tlr7_walker_SHIV", "vac_tlr7_walker_SHIV", "antibody_vac_tlr7_walker_SHIV"),
                       labels = c("Baseline", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac"))) 

rbind(
  tlr7beta_trials[[1]] %>% mutate(mod = "TLR7:beta"), 
  tlr7lambda_trials[[1]] %>% mutate(mod = "TLR7:lambda")
) %>% 
  mutate(name = factor(name, levels = c("baseline", "antibody_tlr7_walker_SHIV", "vac_tlr7_walker_SHIV", "antibody_vac_tlr7_walker_SHIV"),
                       labels = c("Baseline", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac"))) %>% 
  ggplot(aes(peak, final)) + 
  geom_point(data = . %>% filter(t_reb != Inf), aes(color = log10(t_reb)), size = 2) + 
  geom_point(data = . %>% filter(t_reb == Inf), color = "lightgray", size = 2) + 
  facet_grid(mod~name) + 
  geom_hline(yintercept = VL_cutoff, lty = "dashed", color = "gray") + 
  geom_vline(xintercept = VL_cutoff, lty = "dashed", color = "gray") + 
  geom_text(data = combined_quad, aes(5,5,label = quad_1)) + 
  geom_text(data = combined_quad, aes(5,-1,label = quad_2)) + 
  geom_text(data = combined_quad, aes(0,-1,label = quad_3)) + 
  scale_color_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, 
                        trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), 
                        labels=c("3 days", "1 week", "1 month", "3 months", "1 year")) + 
  labs(x = "Peak viral load", y = "Set point viral load")

rbind(
  tlr7beta_trials[[3]] %>% mutate(mod = "TLR7:beta"), 
  tlr7lambda_trials[[3]] %>% mutate(mod = "TLR7:lambda")
) %>% 
  mutate(name = factor(name, levels = c("baseline", "antibody_tlr7_walker_SHIV", "vac_tlr7_walker_SHIV", "antibody_vac_tlr7_walker_SHIV"),
                       labels = c("Baseline", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac"))) %>% 
  filter(id < 25) %>% 
  ggplot( aes(group = id, color=final,  x = time, y = V/log(10))) + facet_grid(mod~ name) + 
  geom_line(data = ~filter(.x, pre_activated == "yes"), linetype="solid", linewidth=0.3, color="darkgrey") +
  geom_line(data = ~filter(.x, post_activated == "yes"), linetype="solid") +
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



######### HUMAN MAIN FIGURE ##############

trajectory_record_final_cut <- trajectory_record_final[trajectory_record_final$id <= 20,]
treatment_subset <-  c("baseline", "tlr7_whitney_SIV", "tlr7_vac_barouch_SIV", "antibody_tlr7_barouch_SHIV")

#color by rebound time
m2 <- ggplot(final_df_data %>% filter(name %in% treatment_subset), aes(x = peak, y = final)) + facet_wrap(~ name_f, nrow=1) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data %>% filter(name %in% treatment_subset), t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data %>% filter(name %in% treatment_subset), t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb))) + 
  scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, 
                       trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), 
                       labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
  theme_bw() +  
  xlab("Peak viral load (copies/mL)") + ylab("Set point viral load (copies/mL)") + coord_cartesian(xlim=c(-5, 7), ylim=c(-5, 7)) + 
  geom_vline(xintercept = VL_cutoff, linetype='dashed') + geom_hline(yintercept = VL_cutoff, linetype='dashed') + 
  geom_text(x= 5, y = 6, aes(label = quad_1), data=quad_df%>% filter(name %in% treatment_subset), color='black') + 
  geom_text(x= 5, y = 0, aes(label = quad_2), data=quad_df%>% filter(name %in% treatment_subset), color='black') + 
  geom_text(x= -1.5, y = 0, aes(label = quad_3), data=quad_df%>% filter(name %in% treatment_subset), color='black') + 
  scale_x_continuous(breaks=c(-4, -2, 0, 2, 4, 6), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6), labels=c(expression(10^{-4}),expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}),  expression(10^{6}))) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 
#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m2)

m3 <- ggplot(trajectory_record_final_cut %>% filter(name %in% treatment_subset), aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap(~ name_f, nrow=1) + 
  geom_line(data = filter(trajectory_record_final_cut %>% filter(name %in% treatment_subset), pre_activated == "yes"), linetype="solid", linewidth=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut %>% filter(name %in% treatment_subset), post_activated == "yes"), linetype="solid") +
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


g <- m3 + m2  & plot_layout(ncol = 1)

ggsave(here("results_eclipse", "simulations", "human_simulations_main.pdf"), width=25, height=18, units='cm', plot=g)
ggsave(here( "paper_figures", "figure5.pdf"), width=25, height=15, units='cm', plot=g)



### Miscellaneous plots
# Plot distribution of parameters

parameter_dists <- parameter_df %>% pivot_longer(1:6, names_to = "parameter") %>% 
  ggplot(aes(name_f, exp(value))) + 
  geom_boxplot() + 
  facet_wrap(~parameter, scales = "free_y") + scale_y_log10()




######HUMAN SUPPLEMENT ##########

#color by rebound time
m2 <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_wrap(~ name_f) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb))) + scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
  theme_bw() +  
  xlab("Peak viral load (copies/mL)") + ylab("Set point viral load (copies/mL)") + coord_cartesian(xlim=c(-2.5, 6.5), ylim=c(-2.5, 6.5)) + 
  geom_vline(xintercept = VL_cutoff, linetype='dashed') + geom_hline(yintercept = VL_cutoff, linetype='dashed') + geom_text(x= 5, y = 6, aes(label = quad_1), data=quad_df, color='black') + 
  geom_text(x= 5, y = 0, aes(label = quad_2), data=quad_df, color='black') + geom_text(x= -1.5, y = 0, aes(label = quad_3), data=quad_df, color='black') + 
  scale_x_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  scale_y_continuous(breaks=c(-2, 0, 2, 4, 6), labels=c(expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}))) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

show(m2)

ggsave(here("paper_figures", "figure_S9.pdf"), plot=m2, width=20, height=20, unit='cm')

m3 <- ggplot(trajectory_record_final_cut , aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap(~ name_f) + 
  geom_line(data = filter(trajectory_record_final_cut , pre_activated == "yes"), linetype="solid", linewidth=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut , post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", 
                        breaks=c(-3,-1,1,3,5), labels=c(expression(10^{-3}), expression(10^{-1}), expression(10^{1}), expression(10^{3}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

ggsave(here("paper_figures", "figure_S8.pdf"), plot=m3, width=20, height=20, unit='cm')


######### END HUMAN ############



######### SIV ##############

final_df_data <- rbind(whitney_SIV_trials[[1]] %>% mutate(study = "S1") %>% 
                         filter(!name_f %in% c("TLR7 (S2)", "TLR7 + Vaccine (S2)", "TLR7 + Antibody (S3)", "TLR7 (S3)")), 
                       barouch_SIV_trials[[1]] %>% mutate(study = "S2") %>% 
                         filter(!name_f %in% c("TLR7 + Antibody (S3)", "TLR7 (S3)")))
quad_df <- rbind(whitney_SIV_trials[[2]] %>% mutate(study = "S1") %>% 
                   filter(!name_f %in% c("TLR7 (S2)", "TLR7 + Vaccine (S2)", "TLR7 + Antibody (S3)", "TLR7 (S3)")), 
                 barouch_SIV_trials[[2]] %>% mutate(study = "S2") %>% 
                   filter(!name_f %in% c("TLR7 + Antibody (S3)", "TLR7 (S3)")))
trajectory_record_final <- rbind(whitney_SIV_trials[[3]] %>% mutate(study = "S1") %>% 
                                   filter(!name_f %in% c("TLR7 (S2)", "TLR7 + Vaccine (S2)", "TLR7 + Antibody (S3)", "TLR7 (S3)")), 
                                 barouch_SIV_trials[[3]] %>% mutate(study = "S2") %>% 
                                   filter(!name_f %in% c("TLR7 + Antibody (S3)", "TLR7 (S3)")))
parameter_df <- rbind(whitney_SIV_trials[[4]] %>% mutate(study = "S1") %>% 
                        filter(!name_f %in% c("TLR7 (S2)", "TLR7 + Vaccine (S2)", "TLR7 + Antibody (S3)", "TLR7 (S3)")), 
                      barouch_SIV_trials[[4]] %>% mutate(study = "S2") %>% 
                        filter(!name_f %in% c("TLR7 + Antibody (S3)", "TLR7 (S3)")))

trajectory_record_final_cut <- trajectory_record_final[trajectory_record_final$id <= 20,]

#color by rebound time
m2_siv <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_grid(study~ name_f) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb))) + scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
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


# ggsave("macaque_siv_predictions_scatter_6_10_19.pdf", plot=m2_siv, width=20, height=20, unit='cm')
ggsave(here("paper_figures", "figure_S11.pdf"), plot=m2_siv, width=25, height=10, unit='cm')



m3_siv <- ggplot(trajectory_record_final_cut , aes(group = id, color=final,  x = time, y = V/log(10))) + facet_grid(study~ name_f) + 
  geom_line(data = filter(trajectory_record_final_cut , pre_activated == "yes"), linetype="solid", linewidth=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut , post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", 
                        breaks=c(-3,-1,1,3,5), labels=c(expression(10^{-3}), expression(10^{-1}), expression(10^{1}), expression(10^{3}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

show(m3_siv)

# ggsave("macaque_siv_example_simulations_6_10_19.pdf", plot=m3_siv, width=20, height=20, unit='cm')
ggsave(here("paper_figures", "figure_S10.pdf"), plot=m3_siv, width=25, height=10, unit='cm')


######### END SIV ############

######### SHIV #############

final_df_data <- barouch_SHV_trials[[1]] %>% filter(name_f != "TLR7 + Vaccine (S2)")
quad_df <- barouch_SHV_trials[[2]] %>% filter(name_f != "TLR7 + Vaccine (S2)")
trajectory_record_final <- barouch_SHV_trials[[3]] %>% filter(name_f != "TLR7 + Vaccine (S2)")
parameter_df <- barouch_SHV_trials[[4]] %>% filter(name_f != "TLR7 + Vaccine (S2)")

trajectory_record_final_cut <- trajectory_record_final[trajectory_record_final$id <= 20,]

#color by rebound time
m2_shiv <- ggplot(final_df_data, aes(x = peak, y = final)) + facet_wrap(~ name_f) + 
  #geom_point(alpha = 0.4, shape = 21, colour = "black", fill = NA, stroke = .15, size=1) +
  
  geom_point(data=filter(final_df_data, t_reb == Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", fill="lightgray") +
  geom_point(data=filter(final_df_data, t_reb != Inf), size=2.7, alpha=.55, shape=21, stroke=0.15, color="lightgray", aes(fill = log10(t_reb))) + scale_fill_distiller(name="Rebound\ntime", palette = "Spectral", direction=1, trans = trim_tails(range = c(log10(3),log10(365))), limits=c(log10(3), log10(365)), breaks=c(log10(3), log10(7), log10(30), log10(90), log10(365)), labels=c("3 days", "1 week", "1 month", "3 months", "1 year"))+
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

ggsave(here("paper_figures", "figure_S13.pdf"), plot=m2_shiv, width=20, height=20, unit='cm')
# ggsave("macaque_shiv_predictions_scatter_6_10_19.pdf", plot=m2_shiv, width=20, height=20, unit='cm')


m3_shiv <-  ggplot(trajectory_record_final_cut , aes(group = id, color=final,  x = time, y = V/log(10))) + facet_wrap( ~ name_f) + 
  geom_line(data = filter(trajectory_record_final_cut , pre_activated == "yes"), linetype="solid", linewidth=0.3, color="darkgrey") +
  geom_line(data = filter(trajectory_record_final_cut , post_activated == "yes"), linetype="solid") +
  scale_color_distiller(palette = "Spectral", direction=-1, name="Set point\nviral load", 
                        breaks=c(-3,-1,1,3,5), labels=c(expression(10^{-3}), expression(10^{-1}), expression(10^{1}), expression(10^{3}), expression(10^{5})))+
  theme_bw() +  
  xlab("Time after ART stop (days)") + ylab("Viral load (copies/mL)") + coord_cartesian(xlim=c(0, 120)) + #ylim=c(-2.5, 8)
  geom_hline(yintercept = log10(50), linetype="dashed") + 
  scale_y_continuous(breaks=c(-4, -2, 0, 2, 4, 6, 8), labels=c(expression(10^{-4}), expression(10^{-2}), expression(10^{0}), expression(10^{2}), expression(10^{4}), expression(10^{6}),  expression(10^{8}))) +
  #scale_x_continuous(breaks=c(0, 20, 40, 60, 80, 100, 120)) +
  theme(strip.background = element_blank(),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 


#scale_x_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6)) + scale_y_continuous(breaks=c(0, 2, 4, 6), labels=c(0, 2, 4, 6))

show(m3_shiv)

ggsave(here("paper_figures", "figure_S12.pdf"), plot=m3_shiv, width=20, height=20, unit='cm')
# ggsave("macaque_shiv_example_simulations_6_10_19.pdf", plot=m3_shiv, width=20, height=20, unit='cm')


######### END SHIV ###########







