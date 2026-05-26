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
library(ggh4x)
select <- dplyr::select

source(here("scripts", "forward_simulation_functions.R"))

theme_set(theme_classic(base_size = 12))

# Define ranges for sensitivity analysis ####

lambda_pop <- 50
beta_pop <- 5e-7
ta_pop <- 1
m_pop <- 1
p_pop <- 0.1
p_pop_strong <- 1
Np_pop <- 1e4

# num_its = 11
# lambdas <- exp(seq(log(lambda_range[1]), log(lambda_range[2]), length.out = num_its))
# lambdas <- c(1,  15, 35, 50, 80, 100, 200, 300, 500, 600, 800)
# betas <- (seq((beta_range[1]), (beta_range[2]), length.out = num_its))
# tas <- exp(seq(log(ta_range[1]), log(ta_range[2]), length.out = num_its))
# ms <- exp(seq(log(m_range[1]), log(m_range[2]), length.out = num_its))
# ps <- exp(seq(log(p_range[1]), log(p_range[2]), length.out = num_its))
# Nps <- exp(seq(log(Np_range[1]), log(Np_range[2]), length.out = num_its))

lambdas <- c(1e-7, 10^seq(0.5,log10(500), length.out = 8))
betas <- seq(0,20,by=2.5)*1e-7
tas <- 10^seq(log10(ta_range[1]), log10(ta_range[2]), length.out = 9)
ms <- 10^seq(-2,log10(300),length.out = 9)
ps <- 10^seq(-2,1,length.out = 9)
Nps <- 10^seq(2,6,length.out = 9)
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
dI <- fixed$dI_macaque #infected
dP <- fixed$dP #long lived
dE <- fixed$dE #short memory
c <- fixed$c
I1 <- fixed$I1_monkey
C <- fixed$C
final_time <- fixed$final_time
rho <- fixed$rho
washout_time <- fixed$washout_time
k <- fixed$k
q <- fixed$q



# Simulate many trajectories with varying population parameter values with a weak immune response ####
param_names <- c("lambda", "beta", "ta", "m", "p", "Np")

trajectories <- NULL
for(param in param_names){
  values <- eval(parse(text = paste0(param, "s")))
  parameters <- c(lambda_pop, beta_pop, ta_pop, m_pop, p_pop, Np_pop)
  idx <- which(param_names == param)
  for(i in 1:length(values)){
    parameters[idx] <- values[i]
    sim <- sim_traj_final(parameters, T, T)[[2]]
    trajectories <- rbind(trajectories, sim %>% mutate(parameter = param, parameter_value = values[i], run = i))
  }
}


strong_immune_trajectories <- NULL
for(param in param_names){
  values <- eval(parse(text = paste0(param, "s")))
  parameters <- c(lambda_pop, beta_pop, ta_pop, m_pop, p_pop_strong, Np_pop)
  for(i in 1:length(values)){
    idx <- which(param_names == param)
    parameters[idx] <- values[i]
    sim <- sim_traj_final(parameters, T, T)[[2]]
    strong_immune_trajectories <- rbind(strong_immune_trajectories, 
                                        sim %>% mutate(parameter = param, parameter_value = values[i], run = i))
  }
}

labels <- c(
  bquote("Target cell production rate ("~lambda~")"),
  bquote("Viral infectivity ("~beta~")"),
  bquote("Time to reactivation ("~t[a]~")"),
  bquote("Immune memory production rate ("~m~")"),
  bquote("Immune proliferation rate ("~p~")"),
  bquote("Immune avidity ("~N[p]~")")
)

weak_immune_response <- trajectories %>% 
  mutate(parameter = factor(parameter, levels = param_names, labels = labels)) %>% 
  ggplot(aes(time, V/log(10), group = interaction(run, parameter), color = run)) + 
  geom_line(data = . %>% filter(post_activated == "yes")) +
  # theme_classic() + 
  scale_color_viridis_c() + 
  facet_wrap(~parameter, labeller = "label_parsed", scales = "free")  +
  ggtitle("Weak immune response regime") +
  coord_cartesian(ylim = c(-1,8), xlim = c(0,200), expand = F) +
  scale_color_viridis_c(breaks = c(1,9), labels = c("low", "high")) +
  # theme_classic() +
  theme(legend.title = element_blank(), strip.background = element_blank(),
        plot.title = element_text(face = 2, hjust = 0)) +
  labs(x = "Time (days)", y = "Viral load (log10 copies/mL)")

weak_immune_response_non_immune_param <- trajectories %>% 
  filter(!parameter %in% c("m", "p", "Np")) %>%
  mutate(parameter = factor(parameter, levels = param_names, labels = labels)) %>% 
  ggplot(aes(time, V/log(10), group = interaction(run, parameter), color = run)) + 
  geom_line(data = . %>% filter(post_activated == "yes")) +
  # theme_classic() + 
  scale_color_viridis_c() + 
  facet_wrap(~parameter, labeller = "label_parsed", scales = "free")  +
  ggtitle("Weak immune response regime") +
  coord_cartesian(ylim = c(-1,8), xlim = c(0,200), expand = F) +
  scale_color_viridis_c(breaks = c(1,9), labels = c("low", "high")) +
  # theme_classic() +
  theme(legend.title = element_blank(), strip.background = element_blank(),
        plot.title = element_text(face = 2, hjust = 0)) +
  labs(x = "Time (days)", y = "Viral load (log10 copies/mL)")

strong_immune_response <- strong_immune_trajectories %>% 
  mutate(parameter = factor(parameter, levels = param_names, labels = labels)) %>% 
  ggplot(aes(time, V/log(10), group = interaction(run, parameter), color = run)) + 
  geom_line(data = . %>% filter(post_activated == "yes")) +
  # theme_classic() + 
  scale_color_viridis_c() + 
  facet_wrap(~parameter, labeller = "label_parsed", scales = "free")  +
  ggtitle("Strong immune response regime") +
  coord_cartesian(ylim = c(-1,8), xlim = c(0,200), expand = F) +
  scale_color_viridis_c(breaks = c(1,9), labels = c("low", "high")) +
  # theme_classic() +
  theme(legend.title = element_blank(), strip.background = element_blank(),
        plot.title = element_text(face = 2, hjust = 0)) +
  labs(x = "Time (days)", y = "Viral load (log10 copies/mL)")

strong_immune_response_immune_param <- strong_immune_trajectories %>% 
  filter(parameter %in% c("m", "p", "Np")) %>%
  mutate(parameter = factor(parameter, levels = param_names, labels = labels)) %>% 
  ggplot(aes(time, V/log(10), group = interaction(run, parameter), color = run)) + 
  geom_line(data = . %>% filter(post_activated == "yes")) +
  # theme_classic() + 
  scale_color_viridis_c() + 
  facet_wrap(~parameter, labeller = "label_parsed", scales = "free")  +
  ggtitle("Strong immune response regime") +
  coord_cartesian(ylim = c(-1,8), xlim = c(0,200), expand = F) +
  scale_color_viridis_c(breaks = c(1,9), labels = c("low", "high")) +
  # theme_classic() +
  theme(legend.title = element_blank(), strip.background = element_blank(),
        plot.title = element_text(face = 2, hjust = 0)) +
  labs(x = "Time (days)", y = "Viral load (log10 copies/mL)")


weak_immune_response + 
  strong_immune_response +
  plot_layout(ncol = 1, guides = "collect")


## Calculate relative sensitivity

# Start with base simulation for weak scenario
parameters <- c(lambda_pop, beta_pop, ta_pop, m_pop, p_pop, Np_pop)
sim_baseline <- sim_traj_final(parameters, T, T)[[2]]

# Increase each parameter by 10% and simulate:
epsilon <- 0.1

sensitivities <- NULL
for(param in param_names){
  idx <- which(param_names == param)
  parameters_perturbed <- parameters
  parameters_perturbed[idx] <- parameters_perturbed[idx]*(1+epsilon)
  sim_pert <- sim_traj_final(parameters_perturbed, T, T)[[2]]
  sensitivities <- rbind(sensitivities, sim_pert %>% mutate(parameter = param))
}

sensitivities_strong_immune <- NULL
parameters_strong <- c(lambda_pop, beta_pop, ta_pop, m_pop, p_pop_strong, Np_pop)
sim_baseline_strong <- sim_traj_final(parameters_strong, T, T)[[2]]
for(param in param_names){
  idx <- which(param_names == param)
  parameters_perturbed <- parameters_strong
  parameters_perturbed[idx] <- parameters_perturbed[idx]*(1+epsilon)
  sim_pert <- sim_traj_final(parameters_perturbed, T, T)[[2]]
  sensitivities_strong_immune <- rbind(sensitivities_strong_immune, sim_pert %>% mutate(parameter = param))
}


sensitivities %>% 
  ggplot(aes(time, V/log(10))) + geom_line(aes(color = "perturbed")) + 
  geom_line(data = sim_baseline) + 
  facet_wrap(~parameter)

sensitivities_strong_immune %>% 
  ggplot(aes(time, V/log(10))) + geom_line(aes(color = "perturbed")) + 
  geom_line(data = sim_baseline_strong) + 
  facet_wrap(~parameter)



weak_rel_sens_df <- sensitivities %>% 
  left_join(select(filter(sim_baseline, post_activated == "yes"), time, V_baseline = V), by = "time") %>% 
  filter(!is.na(V), !is.na(V_baseline)) %>% 
  mutate(immune = parameter %in% (c("m", "p", "Np")),
         parameter = factor(parameter, levels = c("lambda", "beta", "ta", "m", "p", "Np"), 
                            labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         relative_change = (exp(V) - exp(V_baseline))/exp(V_baseline) / epsilon) 

strong_rel_sens_df <- sensitivities_strong_immune %>% 
  left_join(select(filter(sim_baseline_strong, post_activated == "yes"), time, V_baseline = V), by = "time") %>% 
  filter(!is.na(V), !is.na(V_baseline)) %>% 
  mutate(immune = parameter %in% (c("m", "p", "Np")),
         parameter = factor(parameter, levels = c("lambda", "beta", "ta", "m", "p", "Np"), 
                            labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         relative_change = (exp(V) - exp(V_baseline))/exp(V_baseline) / epsilon) 

plot_rel_sensitivity <- function(df, immune_params){
  if(immune_params){
    cols <-  brewer_pal(palette = "Dark2")(6)[4:6]
  } else {
    cols <-  brewer_pal(palette = "Dark2")(6)[1:3]
  }
  
  df %>% 
    mutate(time = time - min(time)) %>% 
    filter(post_activated == "yes", immune == immune_params) %>% 
    ggplot(aes(time, relative_change, color = parameter)) + geom_line() + 
    scale_color_manual(values = cols, labels = scales::parse_format()) +
    xlim(c(0,200)) + ylim(c(-10,30)) +
    coord_cartesian(expand = F) +
    labs(x = "Time from reactivation (days)", y = "Relative sensitivity of\nviral load to parameter") + 
    guides(color = guide_legend(position = "inside")) + 
    theme(legend.position.inside = c(1,1), legend.justification = c(1,1), legend.title = element_blank(),
          legend.background = element_blank()) 
}

weak_non_immune_rel <- plot_rel_sensitivity(weak_rel_sens_df, immune_params = F)

weak_immune_rel <- plot_rel_sensitivity(weak_rel_sens_df, immune_params = T)

strong_non_immune_rel <- plot_rel_sensitivity(strong_rel_sens_df, immune_params = F)

strong_immune_rel <- plot_rel_sensitivity(strong_rel_sens_df, immune_params = T)



sup_sensitivity <- weak_immune_response + labs(tag = "a") + 
  free(wrap_elements(weak_non_immune_rel/weak_immune_rel)) + labs(tag = "b") +
  strong_immune_response + labs(tag = "c") +
  free(wrap_elements(strong_non_immune_rel/strong_immune_rel)) +  labs(tag = "d") + 
  plot_layout(widths = c(3, 1.25)) &
  # plot_annotation(tag_levels = c("a","a","a","a")) &
  theme(plot.tag = element_text(face = 2, vjust = -4, size = 20),
        strip.text = element_text(size = 12),
        plot.background = element_blank(),
        strip.clip = "off") 

ggsave(here("paper_figures", "supplemental_sensitivity_plot.pdf"), sup_sensitivity, width = 13.5, height = 10.5)

main_sensitivity <- weak_immune_response_non_immune_param + weak_non_immune_rel +
  strong_immune_response_immune_param + strong_immune_rel +
  plot_layout(widths = c(3, 1)) &
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(face = 2, size = 20),
        plot.title = element_blank(),
        strip.text = element_text(size = 12),
        plot.background = element_blank(),
        strip.clip = "off") 

ggsave(here("paper_figures", "main_sensitivity_plot.pdf"), main_sensitivity, width = 13.6, height = 6.5)
