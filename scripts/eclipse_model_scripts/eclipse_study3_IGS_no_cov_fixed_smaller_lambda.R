# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# If FALSE, results will be loaded from previous runs of this script (after running project setup)
RERUN <- TRUE

# Project setup ####
## Results folder
study <- "study3"
project_name <- "IGS_no_cov_fixed_param_smaller_lambda.mlxtran"
out_folder <- here("results_eclipse", study)
if(!dir.exists(out_folder)) dir.create(out_folder)

## Parameters to fix:
fixed <- c("lambda_pop" = 200, "p_pop" = 1)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "source", "ABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", ncol(data)- 4))

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]


# Run Analysis ####

if(RERUN){
  ## Create Project
  create_project(model_file, data_file, headerTypes, error_model,
                 exclude_random_effects, treatment_effects = NULL, 
                 normal_dist = NULL, fixed = fixed)
  
  saveProject(here(out_folder, project_name))
  
  
  # Initial guess sampling from a narrower range of values ####
  
  ## Calculate narrower ranges for IGS:
  initialValue_ranges_narrow <- read_csv(here(out_folder, "initialValue_ranges_narrow.csv")) %>% column_to_rownames("parameter")

  
  ## initial guess sampling from narrow range
  # Use initial condition of 1 for all random effects
  parameter_results <- initial_guess_sampling(
    initialValue_ranges_narrow, 
    nits = 50,
    record_shrinkage = T,
    record_indiv_params = T,
    set_variance = 1,
    fixed = fixed,
    record_table = T
  )
  
  ## Save results
  for(name in names(parameter_results)){
    assign(name, parameter_results[[name]])
  }
  
  save(list = names(parameter_results),
       file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
  
}else{
  # Load Project
  loadProject(here(out_folder, project_name))  
  # Load results from IGS over narrow range
  load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
}


# Visualize results from narrow IGS ####

plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

## Assess identifiability
correlate_estimates(tabestimates, remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)), remove_bicc_outliers = T)

compare_estimates(tabestimates, remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  if(parameter %in% names(fixed)) next
  print(parameter_plots(
    parameter, remove_bicc_outliers = F, remove_estimate_outliers = F,
    log10_tabestimates = NULL, tabestimates = tabestimates,
    tabiters = tabiters, phase = phase, initialvalues = initialvalues
  ))
}
dev.off()

## Plot individual fits
indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams)
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
      ggplot(aes(time, indivPred_mode, color = (!!as.symbol(param)), group = interaction(id, run))) + 
      geom_point(data = . %>% filter(run == min(run)), alpha = 0.5, inherit.aes = F, aes(time, vl)) +
      geom_line(alpha = 0.5) + 
      geom_line(data = . %>% merge(distinct(tabestimates, run, bicc) %>% filter(bicc == min(bicc))),
                color = "black", lty = "dashed") +
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1)
  )
}
dev.off()

## Plot population fits
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

# values of P are too large to be physiologically relevant ....
# Could reduce NE to reduce m
table_predictions %>% 
  left_join(distinct(tabestimates, bicc, run)) %>% 
  ggplot(aes(time, P_mode, group = factor(run), color = bicc)) + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

## Reasonable RSE for the best estimates... high for m
tabestimates %>% filter(bicc == min(bicc))

# Save best estimates to use as initial conditions for SCM
best_estimates <- tabestimates %>% filter(bicc == min(bicc))
write_csv(best_estimates, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "best_estimates.csv"))
