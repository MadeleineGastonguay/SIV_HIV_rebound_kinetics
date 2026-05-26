# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "phan_human"
project_name <- "no_cov_fixed.mlxtran"
out_folder <- here("results_eclipse", study)

if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_HIV_final.txt")

## Parameters to fix:
fixed <- c("lambda_pop" = 100, "p_pop" = 3)

## Model file
model_file <- here("model", "2regime_model_delay_HIV.txt")

## Data file
data_file <- here("data", "source", "phan_HIV_data.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", "catcov")

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

# Initial guess sampling over orders of magnitude ####

# ## Create project
# create_project(model_file, data_file, headerTypes, error_model,
#                exclude_random_effects, treatment_effects = NULL, normal_dist = NULL, fixed= fixed)
# 
# saveProject(here(out_folder, project_name))
# # loadProject(here(out_folder, project_name))

igs <- read_csv(here("results_eclipse", "human", "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
  mutate(estimate = signif(estimate, 1)) %>% 
  select(parameter, estimate)

alt_mod(
  out_folder, str_remove(project_name,"\\.mlxtran"),
  model_file, data_file, 
  igs,  treatment_effects = NULL, 
  headerTypes, error_model,
  exclude_random_effects,normal_dist = NULL, 
  fixed= fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 485.4712593   507.4712593   520.4298514   533.7136359     5.0000000     0.4942006 

alt_mod(
  out_folder, "Np_PTC_fixed",
  model_file, data_file, 
  igs,  treatment_effects = list("NPs" = "PTC"), 
  headerTypes, error_model,
  exclude_random_effects,normal_dist = NULL, 
  fixed= fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 466.2620265   490.2620265   504.3986724   517.6824570     5.0000000     0.3740311 

indiv_predictions <- read_csv(here(out_folder, "Np_PTC_fixed", "tables.txt"))
indiv_predictions %>% 
  ggplot(aes(time, log10(V_mode))) + 
  geom_line() + 
  geom_point(data= data, aes(time, vl)) + 
  facet_wrap(~id, scales = "free_x")

parameters <- read_csv(here(out_folder, "Np_PTC_fixed", "IndividualParameters", "estimatedIndividualParameters.txt"))

source("scripts/forward_simulation_functions.R")
sims <- parameters %>% 
  select(id, contains("mode")) %>% 
  rename_with(~str_remove(.x, "_mode")) %>% 
  mutate(betas = betas*1e-7, NPs = NPs*1e5) %>% 
  pmap_df(function(id, lambda, betas, ta, NPs, p, m){
    parameters <- c(lambda, betas, ta, m, p, NPs)
    return(mutate(sim_traj_final(parameters,T,T)[[2]], id = id))
  })

sims %>% 
  filter(post_activated == "yes") %>% 
  left_join(data %>% group_by(id, PTC) %>% summarise(max_time = max(time))) %>% 
  filter(time < max_time) %>% 
  ggplot(aes(time, V/log(10), color = factor(PTC))) + 
  geom_hline(yintercept = log10(50), lty = "dashed") +
  geom_line() + 
  geom_point(data = data, aes(time, vl, shape = factor(censor))) + 
  facet_wrap(~id, scales = "free_x", ncol = 6) + 
  scale_shape_manual(values = c(19,21)) +
  theme(legend.position = "none") + 
  scale_color_brewer(palette = "Set1") + 
  coord_cartesian(ylim = c(0,6))

filtered_params <- parameters %>% 
  select(id, contains("mode")) %>% 
  rename_with(~str_remove(.x, "_mode")) %>% 
  mutate(betas = betas*1e-7, NPs = NPs*1e5)

sims %>% 
  filter(post_activated == "yes") %>% 
  left_join(data %>% distinct(id, PTC)) %>% 
  ggplot(aes(time, (E-I)/log(10), color = factor(PTC), group=id)) + 
  geom_line() + 
  scale_color_brewer(palette = "Set1") 

sims %>% 
  mutate(phan_I = exp(V)*23/3220) %>% 
  ggplot(aes(exp(I), phan_I )) + geom_point()

sims %>% 
  filter(post_activated == "yes") %>% 
  left_join(data %>% distinct(id, PTC)) %>% 
  left_join(filtered_params) %>% 
  mutate(expansion = p*exp(P)/(NPs + exp(V))) %>% 
  ggplot(aes(time, log10(expansion), color = factor(PTC), group=id)) + 
  geom_line() + 
  scale_color_brewer(palette = "Set1") 

sims %>% 
  group_by(id) %>% filter(time == max(time)) %>% 
  left_join(data %>% distinct(id, PTC)) %>% 
  left_join(filtered_params) %>% 
  ggplot(aes(NPs, V/log(10))) + 
  geom_point() + 
  scale_x_log10()
