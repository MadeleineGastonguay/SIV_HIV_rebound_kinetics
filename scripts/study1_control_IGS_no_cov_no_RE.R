# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
project_name <- "IGS_no_cov_no_RE.mlxtran"
out_folder <- here("results", "study1_control")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_final.txt")

## Data file
data_file <- here("data", "derived","STUDY1_CONTROL_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens")

## Random effects to exclude
exclude_random_effects <- c("lambda", "p", "betas", "NPs", "ta", "m") # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

# Initial guess sampling over orders of magnitude ####

## Define initial value ranges
initialValue_ranges_narrow <- read_csv(here(out_folder, "initialValue_ranges_narrow.csv")) %>% 
  column_to_rownames("parameter")


## Create project
create_project(model_file, data_file, headerTypes, error_model,
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL)

# saveProject(here(out_folder, project_name))
loadProject(here(out_folder, project_name))


## initial guess sampling from narrow range
# Use initial condition of 1 for all random effects
parameter_results <- initial_guess_sampling(
  initialValue_ranges_narrow, 
  nits = 25,
  record_shrinkage = F,
  record_indiv_params = F,
  set_variance = 1
)

parameter_results2 <- initial_guess_sampling(
  initialValue_ranges_narrow, 
  nits = 25,
  record_shrinkage = F,
  record_indiv_params = F,
  set_variance = 1
)

tabestimates <- rbind(tabestimates, parameter_results2$tabestimates %>% mutate(run = run + 25))
tabiters <- rbind(tabiters, parameter_results2$tabiters %>% mutate(run = run + 25))
phase <- rbind(phase, parameter_results2$phase)
initialvalues <- rbind(initialvalues, parameter_results2$initialvalues %>% mutate(run = run + 25))


## Save results
for(name in names(parameter_results)){
  assign(name, parameter_results[[name]])
}

# save(tabestimates, tabiters, phase, initialvalues, 
#      file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))


plot_bicc_v_estimate(tabestimates %>% filter(bicc > 309), remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates  %>% filter(bicc > 309), remove_bicc_outliers = T, #remove_estimate_outliers = T,
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates %>% filter(bicc > 309), remove_bicc_outliers = T, remove_estimate_outliers = T,
                  out_folder = out_folder, project_name = project_name)

plot_vs_rses(tabestimates %>% filter(bicc > 309), y = "estimate")

estimate_boxplots(tabestimates %>% filter(bicc > 309), out_folder = out_folder, project_name = project_name)


pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  print(parameter_plots(
    parameter, remove_bicc_outliers = F, remove_estimate_outliers = F,
    log10_tabestimates = NULL, tabestimates = tabestimates %>% filter(run != 17),
    tabiters = tabiters %>% filter(run != 17), phase = phase[-17,], initialvalues = initialvalues %>% filter(run != 17)
  ))
}
dev.off()

predictions <- read_csv("results/study1_control/IGS_no_cov_no_RE/predictions.txt")
predictions %>% distinct(time, popPred) %>% 
  ggplot(aes(time, popPred)) + 
  geom_line(linewidth = 1) + facet_wrap(~id) + 
  geom_point(aes(time, vl), data = predictions)

# 
# # Run profile likelihood for lambda and p ####
# igs <- tabestimates %>% filter(bicc == min(bicc)) %>% mutate(estimate = signif(estimate, 1)) 
# 
# lambda_profile <- profile_likelihood("lambda_pop", c(40,400), igs, nits = 60)
# 
# # write_csv(lambda_profile, here(out_folder, "lambda_profile_likelihood.csv"))
# lambda_profile <- read_csv( here(out_folder, "lambda_profile_likelihood.csv"))
# 
# lambda_profile %>% 
#   distinct(lambda_pop, bicc) %>% 
#   ggplot(aes(lambda_pop, bicc)) + 
#   geom_point()
# 
# # Best lambda is 230
# lambda_profile %>% filter(bicc == min(bicc))
# 
# 
# p_profile <- profile_likelihood("p_pop", c(0.5,5), igs, nits = 60)
# # write_csv(p_profile, here(out_folder, "p_profile_likelihood.csv"))
# p_profile <- read_csv( here(out_folder, "p_profile_likelihood.csv"))
# 
# p_profile %>% 
#   distinct(p_pop, bicc) %>% 
#   ggplot(aes(p_pop, bicc)) + 
#   geom_point()
# 
# # Best p is 2.5
# p_profile %>% filter(bicc == min(bicc))
# 
# 
