# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
project_name <- "IGS_no_cov_removed_RE_p.mlxtran"
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
exclude_random_effects <- c("lambda", "p") # Either NULL or a character vector of parameters with no random effects

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

saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))


## initial guess sampling from narrow range
# Use initial condition of 1 for all random effects
parameter_results <- initial_guess_sampling(
  initialValue_ranges_narrow, 
  nits = 50,
  record_shrinkage = T,
  record_indiv_params = T,
  set_variance = 1
)

# parameter_results2 <- initial_guess_sampling(
#   initialValue_ranges_narrow, 
#   nits = 40,
#   record_shrinkage = T,
#   record_indiv_params = T,
#   set_variance = 1
# )

## Save results
for(name in names(parameter_results)){
  assign(name, parameter_results[[name]])
}

# tabestimates <- rbind(tabestimates, parameter_results2$tabestimates %>% mutate(run = run + 10))
# tabiters <- rbind(tabiters, parameter_results2$tabiters %>% mutate(run = run + 10))
# phase <- rbind(phase, parameter_results2$phase)
# initialvalues <- rbind(initialvalues, parameter_results2$initialvalues %>% mutate(run = run + 10))
# shrinkage <- rbind(shrinkage, parameter_results2$shrinkage %>% mutate(run = run + 10))
# indivParams <- rbind(indivParams, parameter_results2$indivParams %>% mutate(run = run + 10))
# indivFit <- rbind(indivFit, parameter_results2$indivFit %>% mutate(run = run + 10))

# save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit,
#      file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

# load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams)
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
    ggplot(aes(time, indivPred_mode, color = log10(!!as.symbol(param)), group = interaction(id, run))) + 
      geom_line(alpha = 0.5) + 
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1)
  )
}
dev.off()

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "population_trajectories_by_param.pdf"))
indivFit %>% 
  left_join(tabestimates %>% distinct(run, bicc)) %>%
  ggplot(aes(time, popPred, color = bicc, group = interaction(run))) + 
  geom_line(alpha = 0.5) + 
  # facet_wrap(~id) +
  scale_color_viridis_c(direction = -1) 

for(param in names(indivParams)[2:7]){
  print(
    indivFit %>% 
      left_join(tabestimates %>% filter(parameter == str_interp("${param}_pop"))) %>%
      ggplot(aes(time, popPred, color = estimate, group = interaction( run))) + 
      geom_line(alpha = 0.5) + 
      # facet_wrap(~id) +
      scale_color_viridis_c(direction = -1) + 
      labs(color = param)
  )
}
dev.off()


plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates, remove_bicc_outliers = T, #remove_estimate_outliers = T,
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, remove_bicc_outliers = T, #remove_estimate_outliers = T,
                  out_folder = out_folder, project_name = project_name)

plot_vs_rses(tabestimates)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)


pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  print(parameter_plots(
    parameter, remove_bicc_outliers = F, remove_estimate_outliers = F,
    log10_tabestimates = NULL, tabestimates = tabestimates,
    tabiters = tabiters, phase = phase, initialvalues = initialvalues
  ))
}
dev.off()

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
