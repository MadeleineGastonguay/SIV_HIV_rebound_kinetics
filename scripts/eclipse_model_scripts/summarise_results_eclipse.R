## This script collects and summarises results from fitting all four studies (2 SIV, 1 SHIV, and human data) separately

library(here)
library(tidyverse)

theme_set(theme_bw() + 
            theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
                  panel.border = element_rect(colour = "black"), aspect.ratio=1) )

# Define colors for plotting
colors <- c("0.0" = "black", "1.0" = "#00A3D3", "0.1" = "#FFD871", "1.1" ="#46C029")

source(here("scripts", "forward_simulation_functions.R"))

out_folder <- here("results_eclipse", "summary_single_fit")

if(!dir.exists(out_folder)) dir.create(out_folder)

# datasets <- c("study1", "study2", "study3", "study4")
datasets <- c("study1", "study2", "study3")


## Create a table with individual parameter estimates ####
all_indivParams <- NULL

for(dataset in datasets){
  # print(dataset)
  indivParams <- read_csv(here("results_eclipse", dataset, "best_model", "IndividualParameters", "estimatedIndividualParameters.txt")) %>% 
    mutate(data = dataset) 
  
  if(dataset != "study4"){
    if(!"ab" %in% names(indivParams)) indivParams <- indivParams %>% mutate(ab = 0)
    if(!"vac" %in% names(indivParams)) indivParams <- indivParams %>% mutate(vac = 0) 
    
    indivParams <- indivParams %>% mutate(intTLR7VAC = tlr7*vac,
                                          intTLR7AB = tlr7*ab,
                                          intTLR7ABVAC = 0)
    
    if("inttlr7ab" %in% names(indivParams)) indivParams <- select(indivParams, -inttlr7ab)
    
  }else{
    indivParams <- indivParams %>% mutate(ab = intTLR7AB, vac = intTLR7VAC, tlr7 = 1*(intTLR7AB==1 | intTLR7VAC == 1))
  }  
  all_indivParams <- rbind(all_indivParams, indivParams) %>% select(data, id, everything())
}

glimpse(all_indivParams)

write_csv(all_indivParams, here(out_folder, "SIV_SHIV_individual_parameters.csv"))

# ## Write out individual parameters for humans
# human_indivParam <- read_csv(here("results_eclipse", "human", "best_model", "IndividualParameters", "estimatedIndividualParameters.txt")) %>% 
#   mutate(data = "human") %>% 
#   select(data, id, everything())
# 
# write_csv(human_indivParam, here(out_folder, "human_individual_parameters.csv"))


## Create a table with individual predictions ####
all_indivFit <- data.frame()

for(dataset in datasets){
  indivFit <- read_delim(here("results_eclipse", dataset, "best_model", "predictions.txt")) %>% mutate(data = dataset)
  all_indivFit <- rbind(all_indivFit, indivFit) %>% select(data, id, everything())
}

write_csv(all_indivFit, here(out_folder, "SIV_SHIV_individual_predictions.csv"))



## Create a table with population parameter estimates ####

all_popparam <- data.frame()

for(dataset in datasets){
  popparam <- read_delim(here("results_eclipse", dataset, "best_model", "populationParameters.txt"), delim = ",") %>% 
    mutate(data = dataset) %>% 
    select(data, parameter, value, se_sa, rse_sa) 
  
  all_popparam <- rbind(all_popparam, popparam) 
}

write_csv(all_popparam, here(out_folder, "SIV_SHIV_population_parameters.csv"))

# human_popparam <- read_delim(here("results_eclipse", "human", "best_model", "populationParameters.txt"), delim = ",") %>% 
#   mutate(data = "human") %>% 
#   select(data, parameter, value, se_sa, rse_sa)
# 
# write_csv(human_popparam, here(out_folder, "human_population_parameters.csv"))

## Write out population parameters by group ####

get_group_parameters <- function(popparams, effects = c("tlr7", "vac", "ab")){
  # popparams <- as.data.frame(popparams)
  popparams <- popparams %>% filter(!grepl("omega", parameter))
  param <- popparams$value
  names(param) <- popparams$parameter
  
  covariates <- make_covariate_matrix_list(param)
  baseline <- log(param[grepl("pop", names(param))])
  names(baseline) <- str_remove(names(baseline), "_pop")
  
  # make sure covariates and baseline are in the same order
  covariates <- covariates[names(baseline)]
  
  
  control = exp(add_covariate_effects(baseline, c(0,0,0,0), covariates))
  tlr7 = exp(add_covariate_effects(baseline, c(1,0,0,0), covariates))
  
  
  if("vac" %in% effects){
    vac = exp(add_covariate_effects(baseline, c(0,1,0,0), covariates))
    tlr7vac = exp(add_covariate_effects(baseline, c(1,1,0,0), covariates))  
  }else{
    vac = NA
    tlr7vac = NA
    tlr7abvac = NA
  }
  
  if("ab" %in% effects){
    ab = exp(add_covariate_effects(baseline, c(0,0,0,1), covariates))
    tlr7ab = exp(add_covariate_effects(baseline, c(1,0,0,1), covariates))
  }else{
    ab = NA
    tlr7ab = NA
    tlr7abvac = NA
  }
  
  if(all(c("vac", "ab") %in% effects)){
    tlr7abvac = exp(add_covariate_effects(baseline, c(1,1,0,1), covariates))
    ab = NA
    tlr7 = NA
    vac = NA
  }
  
  
  out <- data.frame(
    parameter = names(baseline),
    control = control,
    tlr7 = tlr7,
    vac = vac, tlr7vac = tlr7vac,
    ab = ab, tlr7ab = tlr7ab,
    tlr7abvac = tlr7abvac
  )
  
  return(out)
  
}

popparams_by_group <- rbind(
  get_group_parameters(all_popparam %>% filter(data == "study1"), effects = "tlr7") %>% mutate(data = "study1"),
  get_group_parameters(all_popparam %>% filter(data == "study2"), effects = c("tlr7", "vac")) %>% mutate(data = "study2"),
  get_group_parameters(all_popparam %>% filter(data == "study3"), effects = c("tlr7", "ab")) %>% mutate(data = "study3")
  # get_group_parameters(all_popparam %>% filter(data == "study4"), effects = c("tlr7", "ab", "vac")) %>% mutate(data = "study4")
) %>% 
  select(data, parameter, everything()) %>% 
  remove_rownames()

write_csv(popparams_by_group, here(out_folder, "SIV_SHIV_population_parameters_by_group.csv"))



## Plot individual parameters by group:

temp_df <- all_indivParams %>% 
  select(data, id, contains("mode"), tlr7, ab, vac) 

plot_indiv_params <- function(indiv_df, popparams, param, log_y = FALSE, limits = NULL){
  
  indiv_df <- indiv_df %>% 
    select(data, id, paste(param, "mode", sep = "_"), tlr7, ab, vac)
  names(indiv_df)[names(indiv_df) == paste(param, "mode", sep = "_")] <- "estimate"
  popparams <- popparams %>% filter(parameter == param)
  
  if(all(indiv_df$ab == 0)){
    labels <- c(
      "0.0" = "Control",
      "1.0" = "TLR7",
      "0.1" = "Vac", 
      "1.1" = "TLR7 + Vac"
    )
    
    indiv_df <- indiv_df %>% mutate(treatment = factor(interaction(tlr7,vac), levels = names(labels), labels = labels) )
    
    popparams <- pivot_longer(popparams, !c(parameter, data), names_to = "treatment", values_to = "popmean") %>% 
      filter(!is.na(popmean)) %>% 
      mutate(treatment = fct_inorder(case_when(
        treatment == "control" ~ "Control",
        treatment == "tlr7" ~ "TLR7",
        treatment == "vac" ~ "Vac",
        treatment == "tlr7vac" ~ "TLR7 + Vac"
      )))
    
    
  }else if(all(indiv_df$vac == 0)){
    
    labels <- c(
      "0.0" = "Control",
      "1.0" = "TLR7",
      "0.1" = "Ab", 
      "1.1" = "TLR7 + Ab"
    )
    
    indiv_df <- indiv_df %>% mutate(treatment = factor(interaction(tlr7,ab), levels = names(labels), labels = labels) )
    
    popparams <- pivot_longer(popparams, !c(parameter, data), names_to = "treatment", values_to = "popmean") %>% 
      filter(!is.na(popmean)) %>% 
      mutate(treatment = fct_inorder(case_when(
        treatment == "control" ~ "Control",
        treatment == "tlr7" ~ "TLR7",
        treatment == "ab" ~ "Ab",
        treatment == "tlr7ab" ~ "TLR7 + Ab"
      )))
    
  }else{
    labels <- c("0.0" = "Control",
                "1.0" = "TLR7 + Vac", 
                "0.1" = "TLR7 + Ab",
                "1.1" = "TLR7 + Vac + Ab")
    
    indiv_df <- indiv_df %>% mutate(treatment = factor(interaction(vac, ab), levels = names(labels), labels = labels) )
    
    popparams <- pivot_longer(popparams, !c(parameter, data), names_to = "treatment", values_to = "popmean") %>% 
      filter(!is.na(popmean)) %>% 
      mutate(treatment = fct_inorder(case_when(
        treatment == "control" ~ "Control",
        treatment == "tlr7vac" ~ "TLR7 + Vac",
        treatment == "tlr7ab" ~ "TLR7 + Ab",
        treatment == "tlr7abvac" ~ "TLR7 + Vac + Ab"
      )))
    
  }
  
  scale <- ifelse(param == "NPs", 1e5, 1)
  
  colors <- c("0.0" = "black", "1.0" = "#00A3D3", "0.1" = "#FFD871", "1.1" ="#46C029")
  names(colors) <- labels
  
  y_label <- switch(param,
                    "lambda" = expression(lambda),
                    # "betas" = expression(beta),
                    "betas" = bquote(beta ~ "(×10"^-7*")"),
                    "ta" = expression(t[a]),
                    "NPs" = expression(N[p]),
                    "p" = expression(p),
                    "m" = expression(m)
  )
  
  plot <- indiv_df %>% 
    ggplot(aes(treatment, estimate*scale, color = treatment)) + 
    geom_violin(aes(treatment, estimate*scale, fill = treatment), alpha = 0.3, inherit.aes = F, color = NA) +
    geom_linerange(
      data = popparams, #%>% mutate(value = ifelse(parameter %in% log_scale, log10(value), value)), 
      aes(y = popmean*scale, xmin = as.numeric(treatment)-0.3, xmax = as.numeric(treatment)+0.3),
      lty = ifelse(unique(indiv_df$data) == "study1", "dashed", "solid"),
      linewidth = 1
    ) +
    ggbeeswarm::geom_beeswarm() + 
    labs(x = "", y = y_label) + 
    scale_color_manual(values=colors) + 
    scale_fill_manual(values=colors)
  
  if(log_y){
    plot <- plot + scale_y_log10(labels = scales::trans_format('log10',scales::math_format(10^.x)))
  }
  
  if(!is.null(limits)){
    plot <- plot + coord_cartesian(ylim = limits)
  }
  
  
  return(plot + theme_classic() + theme(legend.position = "none"))
}


study1_param <- plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                                  popparams_by_group %>% filter(data == "study1"), 
                                  param = "lambda", limits = c(0,600)) +
  ggtitle("Study 1: SIV-TLR7 Late ART") +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "betas", limits = c(0, 9)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "ta", log_y = T) +#, limits = c(10^-1,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "m", log_y = T) +#, limits = 10^c(-1,4)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "p", limits = c(0,4)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "NPs", log_y = T) &#, limits = 10^c(0,7.5)) &
  
  plot_layout(ncol = 1, axes = "collect")


study2_param <- plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                                  popparams_by_group %>% filter(data == "study2"), 
                                  param = "lambda", limits = c(0,80)) +
  ggtitle("Study 2: SIV-TLR7+VAC Early ART") +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "betas", limits = c(0, 20)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "ta", log_y = T, limits = c(10^-1,10^1)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "m", log_y = T, limits = 10^c(0,2)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "p", limits = c(0,1)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "NPs", log_y = T) &#, limits = 10^c(0,7.5)) &
  
  plot_layout(ncol = 1, axes = "collect")

(study1_param | study2_param) & plot_layout(widths = c(1,3), axes = "collect")


study3_param <- plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                                  popparams_by_group %>% filter(data == "study3"), 
                                  param = "lambda", limits = c(0,1500)) +
  ggtitle("Study 3: SHIV-TLR7+AB Early ART") +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "betas", limits = c(0,8)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "ta", log_y = T, limits = c(10^0,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "m", log_y = T) +#, limits = 10^c(-1,4)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "p", limits = c(0,1.5)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    # param = "NPs", log_y = T, limits = 10^c(2.5,7.5)) &
                    param = "NPs", log_y = T, limits = 10^c(3,5)) &
  
  plot_layout(ncol = 1, axes = "collect")

# study4_param <- plot_indiv_params(all_indivParams %>% filter(data == "study4"), 
#                                   popparams_by_group %>% filter(data == "study4"), 
#                                   param = "lambda", limits = c(0, 300)) +
#   ggtitle("Study 4: SHIV-TLR7+AB+Vac Early ART") +
#   plot_indiv_params(all_indivParams %>% filter(data == "study4"), 
#                     popparams_by_group %>% filter(data == "study4"), 
#                     param = "betas") + #, limits = c(1,5.5)) +
#   plot_indiv_params(all_indivParams %>% filter(data == "study4"), 
#                     popparams_by_group %>% filter(data == "study4"), 
#                     param = "ta", log_y = T, limits = c(10^-1,10^2)) +
#   plot_indiv_params(all_indivParams %>% filter(data == "study4"), 
#                     popparams_by_group %>% filter(data == "study4"), 
#                     param = "m", log_y = T, limits = 10^c(-1,4)) + 
#   plot_indiv_params(all_indivParams %>% filter(data == "study4"), 
#                     popparams_by_group %>% filter(data == "study4"), 
#                     param = "p", limits = c(0,4)) + 
#   plot_indiv_params(all_indivParams %>% filter(data == "study4"), 
#                     popparams_by_group %>% filter(data == "study4"), 
#                     # param = "NPs", log_y = T, limits = 10^c(2.5,7.5)) &
#                     param = "NPs", log_y = T, limits = 10^c(0,7.5)) &
#   
#   plot_layout(ncol = 1, axes = "collect")


# parameter_plots <- (study1_param | study2_param | study3_param | study4_param) & 
#   plot_layout(widths = c(1.5,3,3,3)) & theme(plot.title = element_text(hjust = 0.5))
parameter_plots <- (study1_param | study2_param | study3_param ) & 
  plot_layout(widths = c(1.5,3,3)) & theme(plot.title = element_text(hjust = 0.5))

ggsave(here(out_folder, "EBE_plots.pdf"), parameter_plots, width = 18, height = 12)


# Compare population parameters from all four studies:

compare_popparam <- all_popparam %>% 
  filter(grepl("pop", parameter)) %>% 
  rbind(
    all_popparam %>% 
      filter(parameter %in% c("lambda_pop", "betas_pop")) %>% 
      group_by(data) %>% 
      summarise(value = prod(value), parameter = "lambda*betas", se_sa = NA, rse_sa = NA)
  ) %>% 
  ggplot(aes(data, value, color = ifelse(data %in% c("study1", "study2"), "SIV", "SHIV"))) + 
  geom_point(size = 2) + 
  geom_errorbar(aes(ymin = value*(1-rse_sa/100), ymax = value*(1+rse_sa/100), 
                    lty = ifelse(data == "study1", "Late ART", "Early ART")), 
                width = 0.5) + 
  facet_wrap(~parameter, scales = "free_y") + 
  labs(lty = "ART start", color = "Virus", title = "Population parameters for control animals") 

ggsave(here(out_folder, "control_animal_parameters.png"), compare_popparam, width = 7, height = 5)

## Simulate typical trajectories for each treatment group ####

fixed <- list(
  N_E = 1e4, f = 0.9, dTar = 0.05, dI_macaque = 0.4, dI_human = 1, dP = 1e-3, 
  dE = 1, k = 50000, c = 23, final_time=200, 
  I1_human = 4 * 10 ^ (-6), I1_monkey = 4 * 10 ^ (-5),  C =.7, rho=10,
  detection_threshold=log10(50), v_one_infected = 0.0724, reactivation_rate = 1e-3, 
  LR_decay = 2.7/365, alpha_factor = 1e-8, fano_factor = 10, whole_body_conversion=3*10^4, washout_time = 3,
  tauQ <- 1 # delay
)


N_E <- fixed$N_E #this is about what z approaches, in the dynamics between immune system and virus
#so we estimate at 1% of 10^6 CD4 cells per ml
f <- fixed$f
dTar <- fixed$dTar #target
dI <- fixed$dI_macaque #infected
dP <- fixed$dP #long lived
dE <- fixed$dE #short memory
k <- fixed$k #also k/c
c <- fixed$c
I1 <- fixed$I1_human
C <- fixed$C
final_time <- fixed$final_time
rho <- fixed$rho
washout_time <- fixed$washout_time
tauQ <- fixed$tauQ # delay

I1 <- fixed$I1_monkey
dI <- fixed$dI_macaque

temp_df <- popparams_by_group %>% 
  pivot_longer(!c(data, parameter), names_to = "treatment") %>% 
  filter(!is.na(value)) %>% 
  pivot_wider(names_from = parameter)

typical_trajectories <- data.frame()
for(i in 1:nrow(temp_df)){
  params <- unlist(temp_df[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T, delay = T)[[2]]
  sim2 <- sim_traj_final(params, return_full = T, delay = T)[[2]]
  sim <- cbind(temp_df[i,1:2], sim)
  typical_trajectories <- rbind(typical_trajectories, sim )
}

head(typical_trajectories)

SIV_traj <- typical_trajectories %>% 
  filter(data %in% c("study1", "study2")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "vac" ~ "Vac",
    treatment == "tlr7vac" ~ "TLR7 + Vac"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment, lty = ifelse(data == "study1", "Study 1 (Late ART)", "Study 2 (Early ART)"))) + 
  geom_line(linewidth = 1, alpha = 0.7) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SIV TLR7 + Vac Studies", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
# legend.position = c(1,1), legend.justification = c(1, 1), legend.box = "horizontal") 

SIV_control_traj <- typical_trajectories %>% 
  filter(data %in% c("study1", "study2"), treatment == "control") %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "vac" ~ "Vac",
    treatment == "tlr7vac" ~ "TLR7 + Vac"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment, lty = ifelse(data == "study1", "Study 1 (Late ART)", "Study 2 (Early ART)"))) + 
  geom_line(linewidth = 1) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       color = NULL, linetype = NULL, title = "SIV TLR7 + Vac Studies") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) 
# guides(color = "none")

SHIV_traj <- typical_trajectories %>% 
  filter(data %in% c("study3", "study4")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "ab" ~ "Ab",
    treatment == "tlr7ab" ~ "TLR7 + Ab",
    treatment == "tlr7vac" ~ "TLR7 + Vac",
    treatment == "tlr7abvac" ~ "TLR7 + Vac + Ab"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment, lty = data)) + 
  geom_line(linewidth = 1) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) +
  # scale_color_manual(values = unname(colors)) + 
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SHIV TLR7 + Ab Study", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
# legend.position = c(1,1), legend.justification = c(1, 1), legend.box = "horizontal") 

SIV_traj_study1 <- typical_trajectories %>% 
  filter(data == "study1") %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "vac" ~ "Vac",
    treatment == "tlr7vac" ~ "TLR7 + Vac"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment, lty = ifelse(data == "study1", "Study 1 (Late ART)", "Study 2 (Early ART)"))) + 
  geom_line(linewidth = 1, alpha = 0.7) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SIV TLR7 + Vac Studies", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

trajectories <-  SIV_traj + SIV_control_traj  + theme(plot.title = element_blank()) + SHIV_traj + plot_layout(ncol = 1) 

ggsave(here(out_folder, "simulated_trajectories.pdf"), trajectories, width = 7, height = 10)

ggsave(here(out_folder, "study1_simulated_trajectories.pdf"), SIV_traj_study1, width = 7, height = 10/3)





## Compare TLR7 animals that don't rebound in Study 1 and Study 3 ####
TLR7_indiv_trajectories <- data.frame()
temp_df <- all_indivParams %>% filter(tlr7 == 1, ab == 0, vac == 0) %>% 
  select(data, id, contains("mode")) %>% 
  rename_with(function(x) str_remove(x, "_mode"), everything()) 
for(i in 1:nrow(temp_df)){
  params <- unlist(temp_df[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T)[[2]]
  sim <- cbind(temp_df[i,1:2], sim)
  TLR7_indiv_trajectories <- rbind(TLR7_indiv_trajectories, sim )
}

TLR7_indiv_fits <- TLR7_indiv_trajectories %>% 
  mutate(study = case_when(
    data == "study1" ~ "Study 1: SIV-TLR7 (Late ART)",
    data == "study2" ~ "Study 2: SIV-TLR7+VAC (Early ART)",
    data == "study3" ~ "Study 3: SHIV-TLR7+AB",
  )) %>% 
  ggplot(aes(time, V/log(10), group = id)) + 
  geom_line(data = . %>% filter(post_activated == "yes"), alpha = 0.5, color = colors[2]) + 
  geom_line(data = . %>% filter(post_activated == "no"),  color = "darkgray") + 
  geom_line(data = ~rbind(
    filter(.x, data == "study1" & id %in% c(15,19)),
    filter(.x, data == "study3" & id == 18)
  ) %>% filter(post_activated == "yes"), color = colors[2]) +
  facet_wrap(~study) + 
  geom_hline(data = all_indivFit %>%  mutate(study = case_when(
    data == "study1" ~ "Study 1: SIV-TLR7 (Late ART)",
    data == "study2" ~ "Study 2: SIV-TLR7+VAC (Early ART)",
    data == "study3" ~ "Study 3: SHIV-TLR7+AB",
  )) %>% group_by(data,study) %>% summarise(min = min(vl)),
  aes(yintercept = min), lty = "dashed", color = "darkgray") +
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)") 

TLR7_indiv_fits <- TLR7_indiv_fits + 
  geom_line(data = typical_trajectories %>% 
              filter(treatment == "tlr7", post_activated == "yes") %>% 
              mutate(study = case_when(
                data == "study1" ~ "Study 1: SIV-TLR7 (Late ART)",
                data == "study2" ~ "Study 2: SIV-TLR7+VAC (Early ART)",
                data == "study3" ~ "Study 3: SHIV-TLR7+AB",
              )) ,
            aes(time, V/log(10)), inherit.aes = F, linewidth = 1.5, color = colors[2]) + 
  coord_cartesian(xlim = c(0,180), ylim = c(-2.5, 6.5)) + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x),
    breaks = c(-2,0,2,4,6)
  ) +
  labs(title = "Comparison of Rebound Trajectories for TLR7 Animals")

ggsave(here(out_folder, "TLR7_individual_fits.pdf"), TLR7_indiv_fits, width = 10, height = 5)


SHIV_traj_study4 <- typical_trajectories %>% 
  filter(data %in% c("study4")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    # treatment == "tlr7" ~ "TLR7",
    # treatment == "ab" ~ "Ab",
    treatment == "tlr7ab" ~ "TLR7 + Ab",
    treatment == "tlr7vac" ~ "TLR7 + Vac",
    treatment == "tlr7abvac" ~ "TLR7 + Vac + Ab"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment)) + 
  geom_line(linewidth = 1) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) +
  scale_color_manual(values = unname(colors)) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SHIV TLR7 + Ab + Vac Study", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))

ggsave(here(out_folder, "study4_fits.pdf"), SHIV_traj_study4, width = 6, height = 3)


control_traj <- typical_trajectories %>% 
  # filter(data %in% c("study4")) %>% 
  filter(treatment == "control") %>%
  # mutate(treatment = fct_inorder(case_when(
  #   treatment == "control" ~ "Control",
  #   # treatment == "tlr7" ~ "TLR7",
  #   # treatment == "ab" ~ "Ab",
  #   treatment == "tlr7ab" ~ "TLR7 + Ab",
  #   treatment == "tlr7vac" ~ "TLR7 + Vac",
  #   treatment == "tlr7abvac" ~ "TLR7 + Vac + Ab"
  # ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = data,
             lty = ifelse(data %in% c("study1", "study2"), "SIV", "SHIV"))) + 
  geom_line(linewidth = 1) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) +
  scale_color_brewer(palette = "Set1") +
  # scale_color_manual(values = unname(colors)) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5))
