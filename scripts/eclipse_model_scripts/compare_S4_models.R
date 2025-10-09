## This script collects and summarises results from fitting all four studies (2 SIV, 1 SHIV, and human data) separately

library(here)
library(tidyverse)
library(patchwork)
library(ggbeeswarm)
library(ggh4x)

theme_set(theme_bw() + 
            theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
                  panel.border = element_rect(colour = "black"), aspect.ratio=1) )

# Define colors for plotting
colors <- c("Control" = "black", "TLR7 + Antibody" = "#00A3D3", "TLR7 + Vaccine" = "#FFD871", "TLR7 + Antibody + Vaccine" ="#46C029")

source(here("scripts", "forward_simulation_functions.R"))

out_folder <- here("results_eclipse", "summary_single_fit", "study4_comparison")

if(!dir.exists(out_folder)) dir.create(out_folder)

datasets <- c("study4")

STUDY4_DATA <- read_delim(here("data", "source", "VACABTLR7_DATA.txt" ), delim = "\t")

ams <- c("best_model", "best_4.4", "best_4.6")

## Function to get treatment groups:
add_treatments <- function(df){
  df %>% mutate(Treatment = case_when(ab ==0 &vac ==0 &tlr7 ==0 ~ "Control",
                                      ab ==0 &vac ==0 &tlr7 ==1 ~ "TLR7",
                                      ab ==1 &vac ==0 &tlr7 ==0 ~ "Ab",
                                      ab ==0 &vac ==1 &tlr7 ==0 ~ "Vac",
                                      ab ==1 &vac ==0 &tlr7 ==1 ~ "TLR7 + Ab",
                                      ab ==0 &vac ==1 &tlr7 ==1 ~ "TLR7 + Vac",
                                      ab ==1 &vac ==1 &tlr7 ==1 ~ "TLR7 + Ab + Vac"),
                Treatment = factor(Treatment, levels = c("Control", "TLR7", "Vac", "Ab", "TLR7 + Vac", "TLR7 + Ab", "TLR7 + Ab + Vac")))
}

## Create a table with individual parameter estimates ####
all_indivParams <- all_indivFit <- all_popparam <- all_simParams <- NULL

for(am in ams){
  folder <- here("results_eclipse", "study4", am)
  
  indivParams <- read_csv(here(folder, "IndividualParameters", "estimatedIndividualParameters.txt")) %>% 
    mutate(am = am) 
  all_indivParams <- rbind(all_indivParams, indivParams) %>% select(am, id, everything())
  
  indivFit <- read_delim(here(folder, "predictions.txt")) %>% mutate(am = am) 
  all_indivFit <- rbind(all_indivFit, indivFit) %>% select(am, id, everything())
  
  popparam <- read_delim(here(folder, "populationParameters.txt"), delim = ",") %>% 
    mutate(am = am)  %>% 
    select(am, parameter, value, se_sa, rse_sa) 
  all_popparam <- rbind(all_popparam, popparam) 
  
  simParams <- read_csv(here(folder, "IndividualParameters", "simulatedIndividualParameters.txt")) %>% 
    mutate(am = am) 
  all_simParams <- rbind(all_simParams, simParams) %>% select(am, id, everything())
  
  
}

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
  get_group_parameters(all_popparam %>% filter(am == ams[1]), effects = c("tlr7", "ab", "vac")) %>% mutate(am = ams[1]),
  get_group_parameters(all_popparam %>% filter(am == ams[2]), effects = c("tlr7", "ab", "vac")) %>% mutate(am = ams[2]),
  get_group_parameters(all_popparam %>% filter(am == ams[3]), effects = c("tlr7", "ab", "vac")) %>% mutate(am = ams[3])
) %>% 
  select(am, parameter, everything()) %>% 
  remove_rownames()



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
  pivot_longer(!c(am, parameter), names_to = "treatment") %>% 
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

SHIV_traj <- typical_trajectories %>% 
  mutate(Treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "ab" ~ "Ab",
    treatment == "tlr7ab" ~ "TLR7 + Ab",
    treatment == "tlr7vac" ~ "TLR7 + Vac",
    treatment == "tlr7abvac" ~ "TLR7 + Ab + Vac"
  )),
  BICc = case_when(am == "best_model" ~ 904, 
                   am == "best_4.4"~ 927,
                   am == "best_4.6"~ 918),
  am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
              labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)"))) %>% 
  filter(post_activated == "yes") %>% 
  ggplot(aes(time, V/log(10), group = interaction(am, treatment), color = am)) + 
  geom_hline(yintercept = 1.1, lty = "dashed", color = "gray") +
  # geom_line(data = add_treatments(STUDY4_DATA), inherit.aes = F, aes(time, vl, group = id) , alpha = 0.5, linewidth = 0.25) +
  geom_point(data = add_treatments(STUDY4_DATA), inherit.aes = F, aes(time, vl, group = id) , alpha = 0.5) +
  geom_line(linewidth = 1) + 
  coord_cartesian(ylim = c(0,7), expand = FALSE, xlim = c(0,200)) +
  scale_color_viridis_d(direction = -1) +
  # scale_color_manual(values = unname(colors)) +
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "SHIV TLR7 + Ab + Vac Study", color = "Model", lty = NULL) +
  theme(plot.title = element_text(hjust = 0.5, face = 2), legend.position = "bottom") + 
  facet_wrap(~Treatment, dir = "v") 
# legend.position = c(1,1), legend.justification = c(1, 1), legend.box = "horizontal") 

ggsave(here(out_folder, "population_trajectories.png"), SHIV_traj, width = 7, height = 7)


## Plot individual trajectories for each study: ####


temp_df2 <- all_indivParams %>% select(am, id, contains("mode"), tlr7, ab, vac) %>% 
  mutate(treatment = case_when(tlr7 == 0 & ab == 0 & vac == 0 ~ "control",
                               tlr7 == 1 & ab == 0 & vac == 0 ~ "tlr7",
                               tlr7 == 0 & ab == 1 & vac == 0 ~ "ab",
                               tlr7 == 0 & ab == 0 & vac == 1 ~ "vac",
                               tlr7 == 1 & ab == 1 & vac == 0 ~ "tlr7ab",
                               tlr7 == 1 & ab == 0 & vac == 1 ~ "tlr7vac",
                               tlr7 == 1 & ab == 1 & vac == 1 ~ "tlr7abvac")) %>% 
  select(am, treatment, id, contains("mode")) %>% 
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

study4_indiv_fits <- individual_trajectories %>% 
  mutate(Treatment = fct_inorder(case_when(
    treatment == "control" ~ "Control",
    treatment == "tlr7" ~ "TLR7",
    treatment == "ab" ~ "Ab",
    treatment == "tlr7ab" ~ "TLR7 + Ab",
    treatment == "tlr7vac" ~ "TLR7 + Vac",
    treatment == "tlr7abvac" ~ "TLR7 + Ab + Vac"
  )),
  BICc = case_when(am == "best_model" ~ 904, 
                   am == "best_4.4"~ 927,
                   am == "best_4.6"~ 918),
  am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
              labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)"))) %>% 
  arrange(am, Treatment) %>% mutate(id = fct_inorder(as.character(id))) %>%  
  # filter(post_activated == "yes") %>% 
  ggplot(aes(time, V/log(10))) + 
  geom_hline(yintercept = 1.1, lty = "dashed", color = "gray") +
  geom_point(data = add_treatments(STUDY4_DATA),
             aes(time, vl, fill = Treatment), shape = 21, alpha = 0.5) +
  geom_line(data = . %>% filter(pre_activated == "yes"), aes(color = am, group = am), alpha = 0.75, linewidth = 0.5) + 
  geom_line(data = . %>% filter(pre_activated == "no"), aes(color = am, group = am), alpha = 0.75, linewidth = 1) + 
  facet_wrap(~id, ncol = 9) + 
  scale_color_viridis_d(direction = -1) +
  scale_fill_manual(values = unname(colors)) +
  scale_y_continuous(
    labels = scales::label_math(expr = 10^.x),
    breaks = c(-1,1,3,5,7),
    limits = c(-2,7), expand = c(0,0)
  ) +
  scale_x_continuous(limits = c(0,180), expand = c(0,0), breaks = seq(0,180,by = 60)) +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold")) + 
  labs(x = "Time after ART stop (days)", y = "Viral load (copies/ml)",
       title = "Study 4: SHIV-TLR7-Ab-Vac Individual Fits", color = "Model")

ggsave(here(out_folder, "individual_trajectories.png"), study4_indiv_fits, width = 11, height = 8)



## Compare individual parameter estimates

tmp_df <- rbind(
  popparams_by_group %>% 
    pivot_longer(!c(am, parameter)) %>% 
    filter(!is.na(value)),
  
  popparams_by_group %>% 
    pivot_longer(!c(am, parameter)) %>% 
    filter(!is.na(value), parameter %in% c("lambda", "betas")) %>%
    group_by(am, name) %>% 
    summarise(value = prod(value)*1e-7, parameter = "lambda*beta")
) %>% 
  mutate(ab = 1*grepl("ab", name),
         vac = 1*grepl("vac", name),
         tlr7 = 1*grepl("tlr7", name),
         am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
                     labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)")),
         value = case_when(parameter == "betas"~value*1e-7,
                           parameter == "NPs"~value*1e5,
                           .default = value),
         Treatment = factor(interaction(ab, vac), labels = c("Control", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac" )),
         Treatment = fct_rev(Treatment)) 

conditional_indiv_params <- all_simParams %>% 
  mutate(`lambda*beta` = lambda*betas*1e-7, betas=betas*1e-7, NPs = NPs*1e5) %>%
  pivot_longer(c(lambda, betas, ta, NPs, p, m, `lambda*beta`), names_to = "parameter") %>%
  # all_indivParams %>% mutate(`lambda*beta_mode` = lambda_mode*betas_mode) %>% 
  # pivot_longer(contains("mode"), names_to = "parameter") %>% 
  mutate(am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
                     labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)")),
         Treatment = factor(interaction(ab, vac), labels = c("Control", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac" )),
         Treatment = fct_rev(Treatment)) %>% 
  ggplot(aes(value, Treatment, fill = am, color = am)) + 
  geom_violin(scale = "width", alpha = 0.5, color = NA) + 
  geom_boxplot(fill = NA, width = 0.2, position = position_dodge(width = 0.9)) +
  # geom_beeswarm(dodge.width = 0.9, alpha = 0.75) +
  facet_wrap(~parameter, scales = "free_x") + 
  theme_minimal() + 
  scale_x_log10() + 
  # geom_point(data = tmp_df %>% mutate(parameter = paste0(parameter, "_mode")), position = position_dodge(width = 0.9), size = 3, shape = 17 ) +
  # geom_point(data = tmp_df , position = position_dodge(width = 0.9) ) +
  scale_color_viridis_d(direction = -1) + 
  scale_fill_viridis_d(direction = -1) + 
  labs(color = "Model", fill = "Model", y = "Treatment", x = "Parameter Estimate",
       title = "Distribution of individual parameter estimates sampled from the conditional distribution") + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = 2))

ggsave(here(out_folder, "individual_params_conditional.png"), conditional_indiv_params, width = 10, height = 8, bg = "white")

EBE_indiv_params <- all_indivParams %>% 
  mutate(`lambda*beta_mode` = lambda_mode*betas_mode*1e-7, betas_mode=betas_mode*1e-7, NPs_mode = NPs_mode*1e5) %>%
  pivot_longer(contains("mode"), names_to = "parameter") %>%
  mutate(am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
                     labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)")),
         Treatment = factor(interaction(ab, vac), labels = c("Control", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac" )),
         Treatment = fct_rev(Treatment)) %>% 
  ggplot(aes(value, Treatment, fill = am, color = am)) + 
  geom_violin(scale = "width", alpha = 0.5, color = NA) + 
  geom_beeswarm(dodge.width = 0.9, alpha = 0.75) +
  facet_wrap(~parameter, scales = "free_x") + 
  theme_minimal() + 
  scale_x_log10() + 
  geom_point(data = tmp_df %>% mutate(parameter = paste0(parameter, "_mode")), 
             position = position_dodge(width = 0.9),  shape = 108, size = 5 ) +
  # geom_point(data = tmp_df , position = position_dodge(width = 0.9) ) +
  scale_color_viridis_d(direction = -1) + 
  scale_fill_viridis_d(direction = -1) +
  labs(color = "Model", fill = "Model", y = "Treatment", x = "Parameter Estimate",
       title = "Distribution of EBEs") + 
  theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5, face = 2))

ggsave(here(out_folder, "individual_params_EBE.png"), EBE_indiv_params, width = 10, height = 8, bg = "white")

## Compare parameters by group:

pop_parameters_plot_groups <- rbind(
  popparams_by_group %>% 
    pivot_longer(!c(am, parameter), names_to = "Treatment") %>% 
    filter(!is.na(value)),
  
  popparams_by_group %>% 
    pivot_longer(!c(am, parameter), names_to = "Treatment") %>% 
    filter(!is.na(value), parameter %in% c("lambda", "betas")) %>% 
    group_by(am, Treatment) %>% summarise(parameter  = "lambda*beta", value = prod(value)*1e-7)
) %>% 
  mutate(am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
                     labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)")),
         Treatment = factor(Treatment, levels = c("control", "tlr7ab", "tlr7vac", "tlr7abvac"),
                            labels = c("Control" , "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac")),
         Treatment = fct_rev(Treatment),
         value = case_when(parameter == "NPs"~log10(value*1e5),
                           parameter == "betas"~value*1e-7,
                           .default = value),
         parameter = ifelse(parameter == "NPs", "log10(NP)", parameter),
         am = fct_rev(am)) %>% 
  ggplot(aes(value, Treatment, fill = am)) + 
  geom_col(position = position_dodge()) + 
  facet_wrap(~parameter, scales = "free_x") + 
  theme_minimal() + 
  scale_fill_viridis_d() +
  theme(legend.position = "bottom") + 
  labs(x = "Parameter Estimate", fill = "Model", title = "Parameter means for each treatment group") +
  guides(fill = guide_legend(reverse = TRUE)) + 
  theme(plot.title = element_text(face = 2, hjust = 0.5))


ggsave(here(out_folder, "pop_parameter_by_group.png"), pop_parameters_plot_groups, width = 8, height = 7, bg = "white")

## Compare parameter distributions

group_paramemters_sd <- left_join(
  popparams_by_group %>% 
    pivot_longer(!c(am, parameter), names_to = "Treatment", values_to = "mean") %>% 
    filter(!is.na(mean)),
  
  all_popparam %>% filter(grepl("omega", parameter)) %>% 
    separate(parameter, into = c("t", "parameter"), sep = "_") %>% 
    select(am, parameter, sd =value )
) 

density_data <- group_paramemters_sd %>%
  group_by(am, parameter, Treatment) %>%
  group_modify(~ {
    x_vals <- exp(seq((min(log(.x$mean) - 6*.x$sd)),
                      (max(log(.x$mean) + 6*.x$sd)),
                      length.out = 500))
    do.call(rbind, lapply(1:nrow(.x), function(i) {
      data.frame(
        # am = .x$am,
        # Treatment = .x$Treatment,
        x = x_vals,
        y = dnorm(log(x_vals), mean = log(.x$mean[i]), sd = .x$sd[i])
      )
    }))
  }) %>%
  ungroup()

distribution_plot <- density_data %>%
  mutate(x = case_when(
    parameter == "betas" ~ x*1e-7, 
    parameter == "NPs" ~ x*1e5,
    .default = x),
    parameter = factor(parameter, 
                       levels =  c("lambda","betas", "ta", "m", "p", "NPs"),
                       # labels =  c("lambda","beta~'(x'~10^-7~')'", "t['a']", "m", "p", "N['p']")),
                       labels =  c("lambda","beta", "t['a']", "m", "p", "N['p']")),
    Treatment = factor(Treatment, levels = c("control", "tlr7ab", "tlr7vac", "tlr7abvac"),
                       labels = c("Control", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac")),
    am = factor(am, levels = c("best_model", "best_4.6", "best_4.4"),
                labels = c("TLR7:Beta, Ab:Lambda (Best BICc)", "TLR7:Lambda, Ab:Lambda", "TLR7:Lambda, Ab:Beta (BH)"))) %>%
  ungroup %>% 
  group_by(parameter) %>% mutate(ynorm = y/max(y)) %>% 
  ggplot(aes(x = log10(x), y = ynorm, color = am)) +
  geom_line(linewidth = 1, alpha = 0.75) +
  facet_grid(Treatment~parameter, scales = "free_x", labeller = label_parsed) +
  theme_minimal() + 
  scale_x_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) +
  facetted_pos_scales(
    x = list(
      "lambda" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-1,4,by=2)),
      # "beta~'(x'~10^-7~')'" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(5, 10, 30, 50))),
      "beta" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-9,-1,by=2)),
      "t['a']" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(1,5,7,10))),
      "m" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(0,4,by=2)),
      "p" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,2,by=2))
    )
  ) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(),# panel.grid = element_blank(), 
        axis.line.x = element_line(),
        axis.text.x = element_text(size = 12),
        strip.text = element_text(size = 12),
        legend.position = "bottom") + 
  scale_color_viridis_d(direction = -1) + 
  labs(color = "Model", x = "Parameter Estimate",
       title = "Theoretical log-normal parameter distributions for each treatment group") + 
  theme(plot.title = element_text(hjust = 0.5, face = 2))

ggsave(here(out_folder, "parameter_distribution_group.png"), distribution_plot, width = 10, height = 7, bg = "white")

