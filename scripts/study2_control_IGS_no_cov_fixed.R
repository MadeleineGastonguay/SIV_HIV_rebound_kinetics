# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
project_name <- "IGS_no_cov_fixed_param.mlxtran"
out_folder <- here("results", "study2_control")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_final.txt")

## Data file
data_file <- here("data", "derived","STUDY2_CONTROL_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens")

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

# Parameters to fix
# fixed <- c("lambda_pop" = 70, "p_pop" = 0.7)
fixed <- c("lambda_pop" = 100, "p_pop" = 1)

## Create project
create_project(model_file, data_file, headerTypes, error_model,
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL, fixed= fixed)

saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))


# Initial guess sampling over orders of magnitude ####

## Calculate narrower ranges for IGS:
initialValue_ranges_narrow <- read_csv(here(out_folder, "initialValue_ranges_narrow.csv")) %>% 
  column_to_rownames("parameter") %>% select(-best)

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

save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit,
     file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

# load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

best_estimates <- tabestimates %>% filter(bicc == min(bicc))
write_csv(best_estimates, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "best_estimates.csv"))

# # Run model with initial guesses set from best fit:
# igs <- initialvalues %>% filter(run == unique(filter(tabestimates, bicc == min(bicc))$run))
# set_base(igs %>% rename(estimate = initialValue))
# runScenario()
# saveProject(here(out_folder, project_name))

indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams) %>% left_join(distinct(tabestimates, bicc, run))
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
    ggplot(aes(time, indivPred_mode, color = (!!as.symbol(param)), group = interaction(id, run))) + 
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
  scale_color_viridis_c(direction = -1) +
  geom_line(data = . %>% filter(bicc == min(bicc)), lty = "dashed", color = "black")

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

correlate_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = F,
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = F,
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)


calculate_BICw(tabestimates) %>%  filter(cumprob <= 0.9) %>% pull(run) -> runs
distinct(tabestimates, run, bicc) %>% ggplot(aes("BICc", bicc, color = run %in% runs)) + ggbeeswarm::geom_beeswarm()

plot_bicc_v_estimate(tabestimates, color = runs)
plot_bicc_v_estimate(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

correlate_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, remove_estimate_outliers = F,
                    out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

compare_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, remove_estimate_outliers = F,
                  out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

estimate_boxplots(tabestimates  %>% filter(run %in% runs), out_folder = out_folder, project_name = project_name , file_addendum = "_best_runs")

plot_vs_rses(tabestimates, out_folder, project_name)

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



# Plot shrinkage to evaluate random effects
# Shrinkage is large for beta and p 
shrinkage %>% 
  ggplot(aes(parameter, conditionalMode)) + 
  geom_boxplot()
