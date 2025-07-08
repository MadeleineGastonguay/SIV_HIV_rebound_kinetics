# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "human"
project_name <- "IGS_no_cov_fixed_param.mlxtran"
out_folder <- here("results_eclipse", study)

## Model file
model_file <- here("model", "2regime_model_HIV_final.txt")

## Parameters to fix:
fixed <- c("lambda_pop" = 100, "p_pop" = 3)

## Model file
model_file <- here("model", "2regime_model_delay_HIV.txt")

## Data file
data_file <- here("data", "derived", "HUMAN_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "ignore", "observation", "cens", "catcov")

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

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
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL, fixed= fixed)

saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))


# Initial guess sampling from a narrower range of values ####

## initial guess sampling from narrow range
# Use initial condition of 1 for all random effects
parameter_results <- initial_guess_sampling(
  initialValue_ranges_narrow, 
  nits = 100,
  record_shrinkage = T,
  record_indiv_params = T,
  set_variance = 1,
  fixed = fixed
)

## Save results
for(name in names(parameter_results)){
  assign(name, parameter_results[[name]])
}

save(list = names(parameter_results),
     file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

# load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))


indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams) %>% left_join(distinct(tabestimates, run, bicc))
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
    ggplot(aes(time, indivPred_mode, color = (!!as.symbol(param)), group = interaction(id, run))) + 
      geom_point(data = . %>% filter(run == 1), aes(y = LOG10VL), alpha = 0.5, color= "black") +
      geom_line(alpha = 0.5) + 
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1) + 
      geom_line(data = . %>% filter(bicc == min(bicc)), lty = "dashed", color = "black")
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
      ggplot(aes(time, popPred, color = (estimate), group = interaction( run))) + 
      geom_line(alpha = 0.5) + 
      # facet_wrap(~id) +
      scale_color_viridis_c(direction = -1) + 
      labs(color = param) + 
      geom_line(data = . %>% filter(bicc == min(bicc)), lty = "dashed", color = "black")
  )
}
dev.off()


plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates, remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates,
                  out_folder = out_folder, project_name = project_name)

plot_vs_rses(tabestimates, project_name = project_name, out_folder = out_folder)

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

#  RSE looks good
tabestimates  %>% filter(bicc == min(bicc))


# Plot shrinkage to evaluate random effects
# Shrinkage is large for beta and p 
shrinkage %>% 
  ggplot(aes(parameter, conditionalMode)) + 
  geom_boxplot()


## Write best estimates:

best_estimates <- tabestimates %>%  filter(bicc == min(bicc))
write_csv(best_estimates, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "best_estimates.csv"))

## Save the model with best BICc as best model for human study because SCM found no significant treatment effects:

igs <- initialvalues %>% filter(run == unique(best_estimates$run)) %>% rename(estimate = initialValue)

loadProject(here(out_folder, "IGS_no_cov_fixed_param.mlxtran"))
saveProject(here("results_eclipse", study, "best_model.mlxtran"))

set_base(igs, fixed = fixed)
runScenario()
saveProject(here("results_eclipse", study, "best_model.mlxtran"))
generate_report(here("results_eclipse", study, "best_model.mlxtran"))


