## This script collects and summarises results from fitting all four studies (2 SIV, 1 SHIV, and human data) separately

library(here)
library(tidyverse)
library(patchwork)

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

STUDY1_DATA <- read_delim(here("data", "derived", "STUDY1_DATA.txt"), delim = "\t")
STUDY2_DATA <- read_delim(here("data", "derived", "STUDY2_DATA.txt"), delim = "\t")
STUDY3_DATA <- read_delim(here("data", "source", "ABTLR7_DATA.txt"), delim = "\t")
HUMAN_DATA <- read_delim(here("data", "derived", "HUMAN_DATA.txt"), delim = "\t")

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

## Write out individual parameters for humans
human_indivParam <- read_csv(here("results_eclipse", "human", "best_model", "IndividualParameters", "estimatedIndividualParameters.txt")) %>%
  mutate(data = "human") %>%
  select(data, id, everything())

write_csv(human_indivParam, here(out_folder, "human_individual_parameters.csv"))


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

human_popparam <- read_delim(here("results_eclipse", "human", "best_model", "populationParameters.txt"), delim = ",") %>%
  mutate(data = "human") %>%
  select(data, parameter, value, se_sa, rse_sa)

write_csv(human_popparam, here(out_folder, "human_population_parameters.csv"))

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

plot_indiv_params <- function(indiv_df, popparams, param, log_y = FALSE, limits = NULL, lty = "solid"){
  
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
                    "lambda" = expression(italic(lambda)),
                    # "betas" = expression(beta),
                    "betas" = bquote(italic(beta ~ "(×10"^-7*")")),
                    "ta" = expression(italic(t[a])),
                    "NPs" = expression(italic(N[p])),
                    "p" = expression(italic(p)),
                    "m" = expression(italic(m))
  )
  
  plot <- indiv_df %>% 
    ggplot(aes(treatment, estimate*scale, color = treatment)) + 
    geom_violin(aes(treatment, estimate*scale, fill = treatment), alpha = 0.3, inherit.aes = F, color = NA) +
    geom_linerange(
      data = popparams, #%>% mutate(value = ifelse(parameter %in% log_scale, log10(value), value)), 
      aes(y = popmean*scale, xmin = as.numeric(treatment)-0.3, xmax = as.numeric(treatment)+0.3),
      linewidth = 1, lty = lty
    ) +
    ggbeeswarm::geom_beeswarm(size = 2.5) + 
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
                                  param = "lambda") +# , limits = c(0,600)) +
  ggtitle("Study 1: SIV-TLR7 Late ART") +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "betas") +#, limits = c(0, 9)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "ta", log_y = T) +#, limits = c(10^-1,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "m", log_y = T) +#, limits = 10^c(-1,4)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "p") +#, limits = c(0,4)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "NPs", log_y = T) &#, limits = 10^c(0,7.5)) &
  
  plot_layout(ncol = 1, axes = "collect")


study2_param <- plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                                  popparams_by_group %>% filter(data == "study2"), 
                                  param = "lambda", limits = c(0,100)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "betas", limits = c(10, 20)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "ta", log_y = T, limits = c(10^-1,10^1)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "NPs", log_y = T, limits = 10^c(2.5,7)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "p", limits = c(0,1))  + 
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "m", log_y = T, limits = 10^c(1,3)) &
  
  plot_layout(ncol = 1, axes = "collect") &
  theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12))

(study1_param | study2_param) & plot_layout(widths = c(1,3), axes = "collect")


study3_param <- plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                                  popparams_by_group %>% filter(data == "study3"), 
                                  param = "lambda", limits = c(0,400)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "betas", limits = c(0,20)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "ta", log_y = T, limits = c(10^0,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "NPs", log_y = T, limits = 10^c(3,5)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "p", limits = c(0,1.5)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    # param = "m", log_y = T, limits = 10^c(1,5)) &
                    param = "m", log_y = T, limits = 10^c(0,5)) &
  
  plot_layout(ncol = 1, axes = "collect") &
  theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12))

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

ggsave(here(out_folder, "study1_EBE_plots.pdf"), study1_param & 
         theme(plot.title = element_blank(), axis.title.y = element_text(size = 15), 
               axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12)),
       width = 4, height = 12)
ggsave(here(out_folder, "study2_EBE_plots.pdf"), study2_param, width = 6, height = 12)
ggsave(here(out_folder, "study3_EBE_plots.pdf"), study3_param, width = 6, height = 12)


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
dP <- fixed$dP #long lived
dE <- fixed$dE #short memory
k <- fixed$k #also k/c
c <- fixed$c
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
  geom_line(linewidth = 1) + 
  coord_cartesian(ylim = c(1,7), expand = FALSE, xlim = c(0,200)) + 
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
  coord_cartesian(ylim = c(1,7), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       color = NULL, linetype = NULL, title = "SIV TLR7 + Vac Studies") +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5)) +
  guides(color = "none")

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
  coord_cartesian(ylim = c(1,7), expand = FALSE, xlim = c(0,200)) +
  scale_color_manual(values = unname(colors)) +
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
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment)) + #, lty = ifelse(data == "study1", "Study 1 (Late ART)", "Study 2 (Early ART)"))) + 
  geom_line(linewidth = 2) + #, alpha = 0.7) + 
  coord_cartesian(ylim = c(0,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  # scale_linetype_manual(values = c("dotted", "solid")) +
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SIV TLR7 + Vac Studies", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank()) 

trajectories <-  SIV_traj + SIV_control_traj  + theme(plot.title = element_blank()) + SHIV_traj + plot_layout(ncol = 1) 

ggsave(here(out_folder, "simulated_trajectories.pdf"), trajectories, width = 7, height = 10)

ggsave(here(out_folder, "study1_simulated_trajectories.pdf"), SIV_traj_study1, width = 7, height = 10/3)

ggsave(here("presentation figures/study1_sim_traj.pdf"), SIV_traj_study1, width = 7, height = 4)


SIV_traj_study2 <- typical_trajectories %>% 
  filter(data %in% c("study2")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "vac" ~ "Vaccine",
    treatment == "tlr7vac" ~ "TLR7 + Vaccine"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment )) + 
  geom_line(linewidth = 2) + 
  geom_line(data = ~filter(., treatment == "Control"), linewidth = 2, alpha = 0.7, lty = "dotted") + 
  coord_cartesian(ylim = c(1,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SIV TLR7 + Vac Study", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())


ggsave("presentation figures/study2_sim_traj.pdf", SIV_traj_study2, width = 7, height = 4)

# Break it down by treatment for slides ####
SIV_traj_study2_c <- typical_trajectories %>% 
  filter(data %in% c("study2")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "vac" ~ "Vaccine",
    treatment == "tlr7vac" ~ "TLR7 + Vaccine"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment )) + 
  geom_line(data = ~filter(., treatment == "Control"), linewidth = 2) + 
  coord_cartesian(ylim = c(1,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors), drop = F) + 
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       color = NULL, lty = NULL, title = "Viral rebound") +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank()) + 
  guides(color = guide_legend(override.aes = list(lty = "solid")))

SIV_traj_study2_c_t <- SIV_traj_study2_c + 
  geom_line(data = ~filter(., treatment == "TLR7"), linewidth = 2, aes(color = "TLR7"), lty = "dashed") 

SIV_traj_study2_c_t_v <- SIV_traj_study2_c_t + 
  geom_line(data = ~filter(., treatment == "Vaccine"), linewidth = 2, aes(color = "Vaccine")) 

SIV_traj_study2_all <- SIV_traj_study2_c_t_v + 
  geom_line(data = ~filter(., treatment == "TLR7 + Vaccine"), linewidth = 2, aes(color = "TLR7 + Vaccine")) 


ggsave("presentation figures/study2_sim_traj_c.pdf", SIV_traj_study2_c, width = 7, height = 4)
ggsave("presentation figures/study2_sim_traj_c_t.pdf", SIV_traj_study2_c_t, width = 7, height = 4)
ggsave("presentation figures/study2_sim_traj_c_t_v.pdf", SIV_traj_study2_c_t_v, width = 7, height = 4)
ggsave("presentation figures/study2_sim_traj_all.pdf", SIV_traj_study2_all, width = 7, height = 4)

# add data to plot
study2_data <- read_delim(here("data", "derived", "STUDY2_DATA.txt"), delim = "\t") %>% 
  arrange(vac, tlr7) %>% 
  mutate(treatment = fct_inorder(case_when(
    tlr7 == 0 & vac == 0 ~ "Control",
    tlr7 == 1 & vac == 0 ~ "TLR7",
    tlr7 == 0 & vac == 1 ~ "Vaccine",
    tlr7 == 1 & vac == 1 ~ "TLR7 + Vaccine"
  )))

SIV_traj_study2_with_data <- typical_trajectories %>% 
  filter(data %in% c("study2")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "vac" ~ "Vaccine",
    treatment == "tlr7vac" ~ "TLR7 + Vaccine"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment )) + 
  geom_line(data = study2_data, alpha = 0.5, aes(time, vl, color = treatment, group = id), inherit.aes = F) +
  geom_line(linewidth = 1.5) + 
  coord_cartesian(ylim = c(1,7), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SIV TLR7 + Vac Study", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) + 
  facet_wrap(~treatment)

#####

SIV_traj_study3 <- typical_trajectories %>% 
  filter(data %in% c("study3")) %>% 
  mutate(treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "ab" ~ "Antibody",
    treatment == "tlr7ab" ~ "TLR7 + Antibody"
  ))) %>% 
  ggplot(aes(time, V/log(10), group = interaction(data, treatment), color = treatment )) + 
  geom_line(linewidth = 2) + 
  coord_cartesian(ylim = c(1,6), expand = FALSE, xlim = c(0,200)) + 
  scale_color_manual(values = unname(colors)) + 
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "Viral Rebound", color = NULL, lty = NULL) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5), panel.grid = element_blank())


ggsave("presentation figures/study3_sim_traj.pdf", SIV_traj_study3, width = 7, height = 4)

## Plot visual predictive check by group for study 1 ####

loadProject(here("results_eclipse", "study1", "best_model.mlxtran"))

plot_preferences <- getPlotPreferences("plotPredictionDistribution")

plot_pred_check <- function(TLR7 = 0){
  treatment <- as.character(interaction(TLR7))
  treatments <- c("control", "tlr7")
  names(treatments) <- c("0", "1")
  id <- which(treatment == names(treatments))
  treatment_name <- unname(treatments[id])
  print(treatment_name)
  
  plot_preferences$obs$color <- colors[id]
  plot_preferences$censObs$shape <- 1
  plot_preferences$obs$radius <- 2
  plot_preferences$censObs$radius <- 2
  plot_preferences$censObs$color <- colors[id]
  plot_preferences$bandsBase$fill <- colors[id]
  
  plotPredictionDistribution(
    stratify = list(list("tlr7", TLR7+1)),
    settings = list(obs = TRUE, cens = TRUE, median = FALSE),
    preferences = plot_preferences
  ) +
    geom_line(data = typical_trajectories %>% filter(treatment == treatment_name, data == "study1"), 
              aes(time, V/log(10)),
              linewidth = 1.5, inherit.aes = F, color = colors[id]) +
    scale_y_continuous(
      labels = scales::label_math(expr = 10^.x)
    ) + 
    labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
         color = NULL, lty = NULL) +
    theme_bw()  
}

study1_pred_dist <- plot_pred_check(0) + 
  ggtitle("Control") + 
  plot_pred_check(1) + 
  ggtitle("TLR7") + 
  plot_layout(ncol = 2) & 
  theme(legend.position = "none") &
  coord_cartesian(ylim = c(1,8), xlim = c(0,120)) &
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) & 
  geom_hline(yintercept = 1.7, lty = "dashed", color = "gray") 

ggsave("presentation figures/study1_pred_dist.pdf", study1_pred_dist, width = 9, height = 3.5)

## Plot visual predictive check by group for study 2 ####

loadProject(here("results_eclipse", "study2", "best_model.mlxtran"))

plot_preferences <- getPlotPreferences("plotPredictionDistribution")

plot_pred_check <- function(TLR7 = 0, VAC = 0){
  treatment <- as.character(interaction(TLR7, VAC))
  treatments <- c("control", "tlr7", "vac", "tlr7vac")
  names(treatments) <- c("0.0", "1.0", "0.1", "1.1")
  id <- which(treatment == names(treatments))
  treatment_name <- unname(treatments[id])
  print(treatment_name)
  
  plot_preferences$obs$color <- colors[id]
  plot_preferences$censObs$shape <- 1
  plot_preferences$obs$radius <- 2
  plot_preferences$censObs$radius <- 2
  plot_preferences$censObs$color <- colors[id]
  plot_preferences$bandsBase$fill <- colors[id]
  
  plotPredictionDistribution(
    stratify = list(filter = list(list("vac", VAC+1), list("tlr7", TLR7+1))),
    settings = list(obs = TRUE, cens = TRUE, median = FALSE),
    preferences = plot_preferences
  ) +
    geom_line(data = typical_trajectories %>% filter(treatment == treatment_name, data == "study2"), 
              aes(time, V/log(10)),
              linewidth = 1.5, inherit.aes = F, color = colors[id]) +
    scale_y_continuous(
      labels = scales::label_math(expr = 10^.x)
    ) + 
    labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
         color = NULL, lty = NULL) +
    theme_bw()  
}

study2_pred_dist <- plot_pred_check(0,0) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  ggtitle("Control") + 
  plot_pred_check(1,0) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank()) +
  ggtitle("TLR7") + 
  plot_pred_check(0,1) + 
  ggtitle("Vaccine") + 
  plot_pred_check(1,1) + 
  ggtitle("TLR7 + Vaccine") + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank()) +
  plot_layout(ncol = 2) & 
  theme(legend.position = "none") &
  coord_cartesian(ylim = c(1,6.5), xlim = c(0,170)) &
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) & 
  geom_hline(yintercept = log10(200), lty = "dashed", color = "gray") 

ggsave("presentation figures/study2_pred_dist.pdf", study2_pred_dist, width = 9, height = 7)


## Plot visual predictive check by group for study 3 ####

loadProject(here("results_eclipse", "study3", "best_model.mlxtran"))

plot_preferences <- getPlotPreferences("plotPredictionDistribution")

plot_pred_check <- function(TLR7 = 0, AB = 0){
  treatment <- as.character(interaction(TLR7, AB))
  treatments <- c("control", "tlr7", "ab", "tlr7ab")
  names(treatments) <- c("0.0", "1.0", "0.1", "1.1")
  id <- which(treatment == names(treatments))
  treatment_name <- unname(treatments[id])
  print(treatment_name)
  
  plot_preferences$obs$color <- colors[id]
  plot_preferences$censObs$shape <- 1
  plot_preferences$obs$radius <- 2
  plot_preferences$censObs$radius <- 2
  plot_preferences$censObs$color <- colors[id]
  plot_preferences$bandsBase$fill <- colors[id]
  
  plotPredictionDistribution(
    stratify = list(filter = list(list("ab", AB+1), list("tlr7", TLR7+1))),
    settings = list(obs = TRUE, cens = TRUE, median = FALSE),
    preferences = plot_preferences
  ) +
    geom_line(data = typical_trajectories %>% filter(treatment == treatment_name, data == "study3"), 
              aes(time, V/log(10)),
              linewidth = 1.5, inherit.aes = F, color = colors[id]) +
    scale_y_continuous(
      labels = scales::label_math(expr = 10^.x)
    ) + 
    labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
         color = NULL, lty = NULL) +
    theme_bw()  
}

study3_pred_dist <- plot_pred_check(0,0) + 
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_blank()) +
  ggtitle("Control") + 
  plot_pred_check(1,0) + 
  theme(axis.title = element_blank(), 
        axis.text = element_blank()) +
  ggtitle("TLR7") + 
  plot_pred_check(0,1) + 
  ggtitle("Ab") + 
  plot_pred_check(1,1) + 
  ggtitle("TLR7 + Ab") + 
  theme(axis.title.y = element_blank(), 
        axis.text.y = element_blank()) +
  plot_layout(ncol = 2) & 
  theme(legend.position = "none") &
  coord_cartesian(ylim = c(1,6.5), xlim = c(0,170)) &
  theme(panel.grid = element_blank(), plot.title = element_text(hjust = 0.5)) & 
  geom_hline(yintercept = log10(200), lty = "dashed", color = "gray") 

ggsave("presentation figures/study3_pred_dist.pdf", study3_pred_dist, width = 9, height = 7)

## Compare TLR7 animals that don't rebound in Study 1 and Study 3 ####
TLR7_indiv_trajectories <- data.frame()
temp_df <- all_indivParams %>% filter(tlr7 == 1, ab == 0, vac == 0) %>% 
  select(data, id, contains("mode")) %>% 
  rename_with(function(x) str_remove(x, "_mode"), everything()) 
for(i in 1:nrow(temp_df)){
  params <- unlist(temp_df[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T, delay = T)[[2]]
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

## Summary figure of drug effects:

# all_popparam %>% 
#   filter(grepl("beta_", parameter)) %>% 
#   select(-contains("_sa")) %>% 
#   mutate(
#     value = exp(value),
#     treatment = str_remove(parameter, "beta_")
#   ) %>% 


popparams_by_group %>% 
  mutate(tlr7_fc = tlr7/control,
         vac_fc = vac/control,
         ab_fc = ab/control,
         tlr7vac_fc = tlr7vac/control,
         tlr7ab_fc = tlr7ab/control) %>% 
  select(data, parameter, contains("fc")) %>% 
  pivot_longer(!c(data, parameter), names_to = "treatment") %>% 
  # filter(!is.na(value)) %>%
  # pivot_wider(names_from = parameter) %>% 
  mutate(treatment = str_remove(treatment, "_fc$"),
         value = ifelse(value >= 1, value, -1/value)) %>%
  # value = ifelse(value >= 1, value, 1/value)) %>%
  filter(value != 1) %>% 
  # filter(!is.na(value)) %>%
  ggplot(aes(parameter, treatment, fill = value)) + 
  geom_text(aes(label = abs(round(value)), color = value < 0), size = 5) +
  facet_wrap(~data, ncol = 1) +
  scale_color_manual(values = c("darkgreen", "darkred"))
# scale_fill_gradient2(midpoint = 0, limits = c(-370, 15))


## Plot individual trajectories for each study: ####


temp_df2 <- all_indivParams %>% select(data, id, contains("mode"), tlr7, ab, vac) %>% 
  mutate(treatment = case_when(tlr7 == 0 & ab == 0 & vac == 0 ~ "control",
                               tlr7 == 1 & ab == 0 & vac == 0 ~ "tlr7",
                               tlr7 == 0 & ab == 1 & vac == 0 ~ "ab",
                               tlr7 == 0 & ab == 0 & vac == 1 ~ "vac",
                               tlr7 == 1 & ab == 1 & vac == 0 ~ "tlr7ab",
                               tlr7 == 1 & ab == 0 & vac == 1 ~ "tlr7vac")) %>% 
  select(data, treatment, id, contains("mode")) %>% 
  rename_with(~gsub("_mode", "", .x), everything())

individual_trajectories <- data.frame()
for(i in 1:nrow(temp_df2)){
  print(i)
  params <- unlist(temp_df2[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T, delay = T)[[2]]
  sim <- cbind(temp_df2[i,1:3], sim)
  individual_trajectories <- rbind(individual_trajectories, sim )
}

# SIV individual fits:
study1_indiv_fits <- individual_trajectories %>% 
  filter(data == "study1") %>% 
  mutate(treatment = factor(treatment, levels = c("control", "tlr7", "vac", "tlr7vac"),
                            labels = c("Control", "TLR7", "Vac", "TLR7 + Vac"))) %>%
  arrange(data, treatment) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10))) + 
  geom_hline(yintercept = 1.7, lty = "dashed", color = "gray") +
  geom_point(data = STUDY1_DATA,
             aes(time, vl)) +
  geom_line(aes(color = treatment, group = interaction(data, id))) + 
  facet_wrap(~id, ncol = 8) + 
  scale_color_manual(values = unname(colors)) + 
  # xlim(c(0,120)) + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x),
    breaks = c(1,3,5,7), limits = c(1,7.5), expand = c(0,0)
  ) +
  scale_x_continuous(limits = c(0,120), expand = c(0,0)) + 
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold")) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "Study 1: SIV-TLR7 Individual Fits", color = "Treatment group")

study2_indiv_fits <- individual_trajectories %>% 
  filter(data == "study2") %>% 
  mutate(treatment = factor(treatment, levels = c("control", "tlr7", "vac", "tlr7vac"),
                            labels = c("Control", "TLR7", "Vac", "TLR7 + Vac"))) %>%
  arrange(data, treatment) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10))) + 
  geom_hline(yintercept = 1.7, lty = "dashed", color = "gray") +
  geom_point(data = STUDY2_DATA,
             aes(time, vl)) +
  geom_line(aes(color = treatment, group = interaction(data, id))) + 
  facet_wrap(~id, ncol = 9) + 
  scale_color_manual(values = unname(colors)) + 
  # xlim(c(0,120)) + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x),
    breaks = c(1,3,5,7), limits = c(1,7.5), expand = c(0,0)
  ) +
  scale_x_continuous(limits = c(0,180), expand = c(0,0), breaks = seq(0,180,by = 60)) + 
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold")) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "Study 2: SIV-TLR7-Vac Individual Fits", color = "Treatment group")

SIV_indiv_fits <- study1_indiv_fits + study2_indiv_fits + plot_layout(ncol = 1, heights = c(2,3))

ggsave(here("paper_figures", "figure_S5.pdf"), SIV_indiv_fits, width = 10, height = 10)

study3_indiv_fits <- individual_trajectories %>% 
  filter(data == "study3") %>%
  mutate(treatment = factor(treatment, levels = c("control", "tlr7", "ab", "tlr7ab"),
                            labels = c("Control", "TLR7", "Ab", "TLR7 + Ab"))) %>%
  arrange(data, treatment) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10))) + 
  geom_hline(yintercept = 2.3, lty = "dashed", color = "gray") +
  geom_point(data = STUDY3_DATA,
             aes(time, vl)) +
  geom_line(aes(color = treatment, group = interaction(data, id))) + 
  facet_wrap(~id, ncol = 8) + 
  scale_color_manual(values = unname(colors)) + 
  # xlim(c(0,120)) + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x),
    breaks = c(1,3,5,7), limits = c(1,7.5), expand = c(0,0)
  ) +
  scale_x_continuous(limits = c(0,180), expand = c(0,0), breaks = seq(0,180,by = 60)) + 
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold")) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "Study 3: SHIV-TLR7-Ab Individual Fits", color = "Treatment group")


ggsave(here("paper_figures", "figure_S6.pdf"), study3_indiv_fits, width = 10, height = 8)

## Human study

I1 <- fixed$I1_human
dI <- fixed$dI_human

temp_df2 <- human_indivParam %>% select(data, id, contains("mode"), nnrti) %>% 
  mutate(treatment = ifelse(nnrti == 0, "non-NNRTI", "NNRTI")) %>% 
  select(data, treatment, id, contains("mode")) %>% 
  rename_with(~gsub("_mode", "", .x), everything())

human_individual_trajectories <- data.frame()
for(i in 1:nrow(temp_df2)){
  print(i)
  params <- unlist(temp_df2[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T, delay = T)[[2]]
  sim <- cbind(temp_df2[i,1:3], sim)
  human_individual_trajectories <- rbind(human_individual_trajectories, sim )
}

human_indiv_fits <- human_individual_trajectories %>% 
  mutate(treatment = factor(treatment, levels = c("non-NNRTI", "NNRTI"))) %>%
  arrange(data, treatment) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10))) + 
  geom_hline(yintercept = 1.69897, lty = "dashed", color = "gray") +
  geom_point(data = HUMAN_DATA,
             aes(time, LOG10VL)) +
  geom_line(aes(color = treatment, group = interaction(data, id))) + 
  facet_wrap(~id, ncol = 9) + 
  scale_color_manual(values = c("red", "blue")) + 
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x),
    breaks = c(1,3,5,7), limits = c(1,7.5), expand = c(0,0)
  ) +
  scale_x_continuous(limits = c(0,180), expand = c(0,0), breaks = seq(0,180,by = 60)) + 
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold")) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "HIV-infected Individual Fits", color = "Treatment group")

ggsave(here("paper_figures", "figure_S7.pdf"), human_indiv_fits, width = 10, height = 10)



## Generate figure 4 for paper: ####


study1_param_paper <- plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                                        popparams_by_group %>% filter(data == "study1"), 
                                        param = "lambda", lty = "dashed", limits = c(0,800)) +
  ggtitle("Study 1: SIV-TLR7") +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "betas", lty = "dashed", limits = c(0, 20)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "ta", log_y = T, lty = "dashed", limits = c(10^-2,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "m", log_y = T, lty = "dashed", limits = 10^c(0,5)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "p", lty = "dashed", limits = c(0,5)) + 
  plot_indiv_params(all_indivParams %>% filter(data == "study1"), 
                    popparams_by_group %>% filter(data == "study1"), 
                    param = "NPs", log_y = T, lty = "dashed", limits = 10^c(2.5,7)) &
  
  plot_layout(ncol = 1, axes = "collect") &
  theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12))


study2_param_paper <- plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                                        popparams_by_group %>% filter(data == "study2"), 
                                        param = "lambda", limits = c(0,800)) +
  ggtitle("Study 2: SIV-TLR7-VAC") +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "betas", limits = c(0, 20)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "ta", log_y = T, limits = c(10^-2,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "m", log_y = T, limits = 10^c(0,5))+
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "p", limits = c(0,5))  + 
  plot_indiv_params(all_indivParams %>% filter(data == "study2"), 
                    popparams_by_group %>% filter(data == "study2"), 
                    param = "NPs", log_y = T, limits = 10^c(2.5,7)) &
  
  plot_layout(ncol = 1, axes = "collect") &
  theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12))



study3_param_paper <- plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                                        popparams_by_group %>% filter(data == "study3"), 
                                        param = "lambda", limits = c(0,800)) +
  ggtitle("Study 3: SHIV-TLR7-Ab") +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "betas", limits = c(0,20)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "ta", log_y = T, limits = c(10^-2,10^2)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "m", log_y = T, limits = 10^c(0,5)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "p", limits = c(0,5)) +
  plot_indiv_params(all_indivParams %>% filter(data == "study3"), 
                    popparams_by_group %>% filter(data == "study3"), 
                    param = "NPs", log_y = T, limits = 10^c(2.5,7)) &
  
  plot_layout(ncol = 1, axes = "collect") &
  theme(axis.title.y = element_text(size = 15), axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 12))

(study1_param_paper | study2_param_paper | study3_param_paper) & plot_layout(widths = c(1,2,2), axes = "collect") & theme(plot.title = element_text(hjust = 0.5))

trajectory_plots <- (SIV_traj + 
    ggtitle("STUDY1 + STUDY 2\n(SIV-TLR7 + SIV-TLR7-Vac)") + 
    guides(lty = guide_legend(position = "inside", order = 1), color = guide_legend(position = "inside", order = 2)) + 
    theme_classic() + 
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(), plot.title = element_text(hjust = 0.5, face = "bold"),
          legend.position.inside = c(0.6,0.85), legend.box = "horizontal", legend.background = element_blank())) +
  
  (SIV_control_traj + 
     theme_classic()  + 
     guides(lty = guide_legend(position = "inside")) + 
     theme(plot.title = element_blank(), legend.position.inside = c(0.6,0.85))) + 
  
  plot_spacer() + 
  
  (SHIV_traj + 
     ggtitle("STUDY 3: SHIV-TLR7-Ab") + 
     theme_classic()  + 
     guides(color = guide_legend(position = "inside"), lty = "none") + 
     theme(legend.position.inside = c(0.6,0.8),legend.background = element_blank(),
           plot.title = element_text(hjust = 0.5, face = "bold")))  +
  
  plot_layout(ncol = 1, heights= c(1,1,0.2,1)) 

figure4 <- (study1_param_paper & theme(plot.title = element_text(hjust = 0.5, face = "bold")) | 
              study2_param_paper & theme(plot.title = element_text(hjust = 0.5, face = "bold"), 
                                         axis.text.y = element_blank(), axis.title.y = element_blank(), 
                                         axis.line.y = element_blank(), axis.ticks.y = element_blank()) | 
              (plot_spacer()) | 
              (trajectory_plots) | 
              (plot_spacer()) | 
              study3_param_paper & theme(plot.title = element_text(hjust = 0.5, face = "bold"))) + 
  plot_layout(widths = c(1.25,3,0.25,2,0.25,3), axes = "collect") 

ggsave(here("paper_figures", "figure4.pdf"), figure4, width = 20, height = 12)
