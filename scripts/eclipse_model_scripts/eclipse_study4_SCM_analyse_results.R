# Load libraries
library(here)
library(ggh4x)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))
source(here("scripts", "forward_simulation_functions.R"))
theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("ignore", 3), rep("catcov", 3))

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 100, "p_pop" = 0.9)

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 10, height = 7)

write_csv(out$summary, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_summary.csv"))


## Barplot figure with results:
SCM_summary <- out$summary
SCM_summary_figure_study4 <- SCM_summary %>%
  filter(n_covar != 0) %>%
  mutate(id = 1:nrow(.)) %>%
  filter(id < min(id[n_covar > 7])) %>%
  group_by(n_covar) %>%
  mutate(shaded = min(BICc) + 2,
         n_covar_label = factor(as.character(n_covar), levels = 1:7,
                                labels = paste(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th"), "treatment effect added")),
         label = ifelse(BICc < threshold, param, ""),
         param = factor(param, levels = c("lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         label = factor(label, levels = c("","lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("","lambda", "beta", "t['a']", "m", "p", "N['p']"))) %>%
  mutate(param = fct_rev(param),
         treatment = factor(treatment, levels = rev(c("intTLR7VAC", "intTLR7AB", "intTLR7ABVAC")),
                            labels = rev(c("TLR7 + Vaccine", "TLR7 + Antibody", "TLR7 + Antibody + Vaccine")))
  ) %>%
  ggplot(aes(BICc - threshold, treatment, fill = param)) +
  geom_rect(data = ~filter(., BICc == min(BICc)),
            aes(xmin =BICc-threshold , xmax = shaded-threshold, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.5, inherit.aes = F) +
  geom_col(position = position_dodge(), color = NA) +
  geom_text(aes(BICc - threshold, label = label, color = param),
            parse = T, hjust = 1.1, position = position_dodge(width = 0.9), show.legend = F) +
  facet_wrap( ~ n_covar_label, ncol = 2, dir = "v") +
  # scale_y_discrete(labels = c("vac" = "Vaccine", "tlr7" = "TLR7", "inTLR7VAC" = "TLR7 * Vaccine")) +
  # scale_y_discrete(labels = rev(c( "Vaccine",  "TLR7", "TLR7 * Vaccine"))) +
  labs(x = "Change in BICc", y = "Treatment effect added to model", fill = "Parameter") +
  scale_fill_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0")), labels = function(x) parse(text = x)) +
  scale_color_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0"))) +
  guides(fill = guide_legend(reverse = T)) +
  # ggnewscale::new_scale_color() +
  # geom_col(data = . %>% filter(n_covar <= 3) %>% mutate(col = BICc == min(BICc)),
  #          aes(group = param, color = col), fill = NA, position = position_dodge(),  lty = "dashed", show.legend = F) +
  # scale_color_manual(values = c("#FFFFFF00", "black")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        strip.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

## Plot population trends at each step:
included_effects <- SCM_summary %>% 
  mutate(mod = 1:nrow(.)) %>% 
  filter(mod < min(mod[n_covar > 7])) %>% 
  group_by(n_covar) %>% 
  filter(BICc == min(BICc))

stepwise_preds <- stepwise_indivParam <- stepwise_popParam <-NULL
for(mod in included_effects$mod){
  num <- sprintf("%03d", mod)
  model_folder <- here(out_folder, "SCM", "ModelBuilding", str_interp("ModelBuilding_${num}"), str_interp("ModelBuilding_${num}"))
  predictions <- read_csv(here(model_folder, "predictions.txt")) %>% select(-contains("mode"))
  indivParam <- read_csv(here(model_folder, "IndividualParameters", "estimatedIndividualParameters.txt")) %>% select(-contains("mode"))
  popParam <- read_csv(here(model_folder, "populationParameters.txt")) %>% select(-contains("mode")) %>% select(parameter, value)
  stepwise_preds <- rbind(stepwise_preds, predictions %>% mutate(model = num))
  stepwise_indivParam <- rbind(stepwise_indivParam, indivParam %>% mutate(model = num))
  stepwise_popParam <- rbind(stepwise_popParam, popParam %>% mutate(model = num))
}

write_csv(stepwise_preds, here(out_folder, "SCM", "stepwise_predictions.csv"))
write_csv(stepwise_indivParam, here(out_folder, "SCM", "stepwise_individualParameters.csv"))
write_csv(stepwise_popParam, here(out_folder, "SCM", "stepwise_populationParameters.csv"))

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  ggplot(aes(time, popPred, group = interaction(id, model), color = treatment_effect)) + 
  geom_point(data = ~filter(.x, model == "001"), aes(time, vl), alpha = 0.5, color = "black") +
  geom_line(linewidth = 1) + 
  facet_wrap(ab ~ vac, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added")

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  filter(model != "001") %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(ab, vac) %>% 
  mutate(Treatment = fct_inorder(case_when(ab == 1 & vac == 1 ~ "TLR7 + Ab + Vac", 
                                           ab == 1 & vac == 0 ~ "TLR7 + Ab",
                                           ab == 0 & vac == 1 ~ "TLR7 + Vac",
                                           ab == 0 & vac == 0 ~ "Control"))) %>% 
  ggplot(aes(time, popPred, group = interaction(id, model), color = treatment_effect)) + 
  geom_point(aes(time, vl), alpha = 0.5, color = "black") +
  geom_line(linewidth = 1) + 
  facet_grid(treatment_effect ~ Treatment) + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added")

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_mean, group = interaction(id, model), color = treatment_effect)) + 
  geom_hline(yintercept = 1.04, lty = "dashed") +
  geom_point(data = ~filter(.x, model == "001"), aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_wrap(~id, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6))

## Simulate with R: ####

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
# k <- fixed$k_human #also k/c
c <- fixed$c
I1 <- fixed$I1_human
C <- fixed$C
final_time <- fixed$final_time
rho <- fixed$rho
washout_time <- fixed$washout_time
k <- fixed$k
q <- fixed$q
I1 <- fixed$I1_monkey
dI <- fixed$dI_macaque


temp_df2 <- stepwise_indivParam %>% select(model, id, contains("mean"), starts_with("int")) %>% 
  mutate(treatment = case_when(intTLR7AB == 0 & intTLR7VAC == 0 ~ "Control",
                               intTLR7AB == 1 & intTLR7VAC == 0 ~ "TLR7 + Ab",
                               intTLR7AB == 0 & intTLR7VAC == 1 ~ "TLR7 + Vac",
                               intTLR7AB == 1 & intTLR7VAC == 1 ~ "TLR7 + Vac + Ab")) %>% 
  select(model, treatment, id, contains("mean")) %>% 
  rename_with(~gsub("_mean", "", .x), everything())

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

save(individual_trajectories, file = here(out_folder, "SCM", "stepwise_indiv_sims.RData"))

individual_trajectories %>% rename(Treatment = treatment) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10), group = model, color = treatment_effect)) + 
  geom_hline(yintercept = 1.04, lty = "dashed") +
  geom_point(data = data, aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_wrap(~id, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6), xlim = c(0,168)) + 
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
       fill = "Treatment group") + 
  scale_fill_manual(values = colors, labels= c("Control", "TLR7+Antibody", "TLR7+Vaccine", "TLR7+Antibody+Vaccine"))

diff_id <- c(7132, 7133, 7137, 7138, 7142, 7145, 7147, 7149, 7152)

individual_trajectories %>% rename(Treatment = treatment) %>%
  filter(id %in% diff_id) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10), group = model, color = treatment_effect)) + 
  geom_hline(yintercept = 1.1, lty = "dashed") +
  geom_point(data = data %>% filter(id %in% diff_id), aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1) + 
  facet_wrap(~id, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6), xlim = c(0,168)) + 
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
       fill = "Treatment group") + 
  scale_fill_manual(values = colors, labels= c("TLR7+Antibody", "TLR7+Antibody+Vaccine"))

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

stepwise_popParam %>% 
  group_by(model) %>% 
  group_modify(~get_group_parameters(.x, effects = c("tlr7", "ab", "vac")) )

popparams_by_group <- rbind(
  get_group_parameters(stepwise_popParam %>% filter(model == "study1"), effects = "tlr7") %>% mutate(data = "study1"),
  get_group_parameters(all_popparam %>% filter(data == "study2"), effects = c("tlr7", "vac")) %>% mutate(data = "study2"),
  get_group_parameters(all_popparam %>% filter(data == "study3"), effects = c("tlr7", "ab")) %>% mutate(data = "study3")
  # get_group_parameters(all_popparam %>% filter(data == "study4"), effects = c("tlr7", "ab", "vac")) %>% mutate(data = "study4")
) %>% 
  select(data, parameter, everything()) %>% 
  remove_rownames()

## Population trajectories
temp_df <- stepwise_popParam %>% 
  mutate(treatment = case_when(intTLR7AB == 0 & intTLR7VAC == 0 ~ "Control",
                               intTLR7AB == 1 & intTLR7VAC == 0 ~ "TLR7 + Ab",
                               intTLR7AB == 0 & intTLR7VAC == 1 ~ "TLR7 + Vac",
                               intTLR7AB == 1 & intTLR7VAC == 1 ~ "TLR7 + Vac + Ab")) %>% 
  select(model, treatment, id, contains("mean")) %>% 
  rename_with(~gsub("_mean", "", .x), everything())

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


## Claculate difference between observation and prediction to see where we are over esitmating and estimating

stepwise_preds %>% 
  filter(id != 7122) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  mutate(error = ifelse(vl <= 1.1 & indivPred_mean <= 1.1, 0, indivPred_mean-vl)) %>% 
  ggplot(aes(time, error, color = treatment_effect)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_wrap(~id) +
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(y = "Difference between observed and predicted viral load", color = "Treatment effect added")

stepwise_preds %>% 
  filter(id %in% diff_id) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  mutate(error = ifelse(vl <= 1.1 & indivPred_mean <= 1.1, 0, indivPred_mean-vl)) %>% 
  ggplot(aes(time, error, color = treatment_effect)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(linewidth = 1) + 
  facet_wrap(~id) +
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(y = "Difference between observed and predicted viral load", color = "Treatment effect added")

## How do random effects change?

stepwise_popParam %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  filter(grepl("omega", parameter)) %>% 
  ggplot(aes(factor(n_covar), value, fill = treatment_effect)) + 
  geom_col() + 
  facet_wrap(~parameter) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added")

# how do treatment effects change?
stepwise_popParam %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  filter(grepl("beta_", parameter)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  ggplot(aes(factor(n_covar), value, fill = treatment_effect)) + 
  geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added")

# how do population parameters change?
stepwise_popParam %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  filter(grepl("pop", parameter)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  ggplot(aes(factor(n_covar), value, fill = treatment_effect)) + 
  geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added")

## Compared to prior model, how does next model change?

compare_df <- merge(
  stepwise_preds %>% 
    select(id, model, time, popPred) %>% 
    filter(model != "001") %>% 
    left_join(mutate(included_effects, model = sprintf("%03d", mod))),
  stepwise_preds %>% 
    select(id, model, time, popPred_ref = popPred) %>% 
    filter(model != "103") %>% 
    left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
    mutate(n_covar = n_covar + 1),
  by = c("n_covar", "id", "time")
  ) 

compare_df %>%
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, treatment.x, sep = ":")))) %>%
  ggplot(aes(time, popPred-popPred_ref, color = treatment_effect, group = interaction(id, model.x))) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_line(linewidth = 1) + 
  facet_wrap(ab~vac, labeller = "label_both") + 
  scale_color_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Time", y = "Change from model at previous step", color = "Treatment effect added")

compare_df %>%
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, treatment.x, sep = ":")))) %>%
  arrange(ab, vac) %>% 
  mutate(Treatment = fct_inorder(case_when(ab == 1 & vac == 1 ~ "TLR7 + Ab + Vac", 
                                           ab == 1 & vac == 0 ~ "TLR7 + Ab",
                                           ab == 0 & vac == 1 ~ "TLR7 + Vac",
                                           ab == 0 & vac == 0 ~ "Control"))) %>% 
  ggplot(aes(time, popPred-popPred_ref, color = treatment_effect, group = interaction(id, model.x))) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_line(linewidth = 1) + 
  facet_grid(treatment_effect~Treatment) + 
  scale_color_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Time", y = "Change from model at previous step", color = "Treatment effect added")

compare_df %>%
  filter(n_covar <=4 ) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, treatment.x, sep = ":")))) %>%
  ggplot(aes(time, popPred-popPred_ref, color = treatment_effect, group = interaction(id, model.x))) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_line(linewidth = 1) + 
  facet_wrap(ab~vac, labeller = "label_both") + 
  scale_color_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Time", y = "Change from model at previous step", color = "Treatment effect added")
