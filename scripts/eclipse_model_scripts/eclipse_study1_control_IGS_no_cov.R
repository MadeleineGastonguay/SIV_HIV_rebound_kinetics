# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# If FALSE, results will be loaded from previous runs of this script (after running project setup)
RERUN <- FALSE

# Project setup ####
## Results folder
project_name <- "IGS_no_cov.mlxtran"
out_folder <- here("results_eclipse", "study1_control")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "derived","STUDY1_CONTROL_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens")

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant")
names(error_model) <- names(data)[headerTypes == "observation"]

# Initial guess sampling over orders of magnitude ####

## Define initial value ranges
initialValue_ranges <- tibble(
  lambda = c(10,1000), 
  betas = c(0.001,100), # Will be multiplied by 10^-7 in the model 
  m = c(0.01, 100),
  p = c(0.01,10), 
  NPs = c(10^2, 10^7)*10^-5, # Will be multiplied by 10^5 in the model
  ta = c(10^-2,10)
) 

# Run Analysis ####

if(RERUN){
  
  
  ## Create project
  create_project(model_file, data_file, headerTypes, error_model,
                 exclude_random_effects, treatment_effects = NULL, normal_dist = NULL)
  
  saveProject(here(out_folder, project_name))
  
  ## Run convergence over orders of magnitude
  initialValue_ranges_log10 <- initialValue_ranges %>% 
    t() %>% as.data.frame %>% 
    setNames(c("min", "max")) %>% 
    rownames_to_column("parameter") %>% 
    mutate(parameter = paste0(parameter, "_pop")) 
  
  # Add ranges from 1 to 5 for random effect standard deviations
  parameters <- getPopulationParameterInformation()$name
  initialValue_ranges_log10 <- rbind(
    initialValue_ranges_log10,
    tibble(parameter = parameters[!grepl("pop", parameters)],
           min = 1, max = 5)
  ) %>% 
    column_to_rownames("parameter")
  
  
  # sample new initial estimates
  log10_parameter_results <- initial_guess_sampling(
    initialValue_ranges_log10, 
    log10 = parameters[grepl("pop", parameters)], 
    nits = 100,
    record_shrinkage = T,
    record_indiv_params = T
  )
  
  ## Save results
  for(name in names(log10_parameter_results)){
    assign(name, log10_parameter_results[[name]])
  }
  
  save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit,
       file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "log10_IGS.RData"))
  
  log10_tabestimates <- tabestimates
  log10_initialvalues <- initialvalues
  
  runs <- calculate_BICw(log10_tabestimates) %>% filter(cumprob <= 0.95) %>% pull(run)
  
  
  
  # Best runs have unrealistically large m, going to remove them
  exclude_runs <- log10_tabestimates %>% filter(bicc < 194) %>% pull(run) %>% unique
  
  new_runs <- calculate_BICw(log10_tabestimates %>% filter(!run %in% exclude_runs)) %>% filter(cumprob <= 0.95) %>% pull(run)
  
  
  ## Calculate narrower ranges for IGS:
  initialValue_ranges_narrow <- get_narrow_ranges(log10_tabestimates %>% filter(!run %in% exclude_runs), p_threshold = 0.95) %>% select(-best)
  
  write_csv(initialValue_ranges_narrow, here(out_folder, "initialValue_ranges_narrow.csv"))
  
  initialValue_ranges_narrow <- initialValue_ranges_narrow %>% 
    column_to_rownames("parameter")
  
  ## initial guess sampling from narrow range
  # Use initial condition of 1 for all random effects
  parameter_results <- initial_guess_sampling(
    initialValue_ranges_narrow, 
    nits = 50,
    record_shrinkage = T,
    record_indiv_params = T,
    set_variance = 1
  )
  
  parameter_results2 <- initial_guess_sampling(
    initialValue_ranges_narrow, 
    nits = 50,
    record_shrinkage = T,
    record_indiv_params = T,
    set_variance = 1
  )
  
  ## Save results
  for(name in names(parameter_results)){
    assign(name, parameter_results[[name]])
  }
  
  # combine elements of parmater_results2 with parameter_results and update the run number in results 2 to add 50 
  tabestimates <- rbind(tabestimates, parameter_results2$tabestimates %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  tabiters <- rbind(tabiters, parameter_results2$tabiters %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  phase <- rbind(phase, parameter_results2$phase)
  initialvalues <- rbind(initialvalues, parameter_results2$initialvalues %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  shrinkage <- rbind(shrinkage, parameter_results2$shrinkage %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  indivParams <- rbind(indivParams, parameter_results2$indivParams %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  indivFit <- rbind(indivFit, parameter_results2$indivFit %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  
  
  save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit,
       file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
  
}else{
  # Load Project
  loadProject(here(out_folder, project_name))  
  # Load results from IGS over orders of magnitude
  load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "log10_IGS.RData"))
  log10_tabestimates <- tabestimates
  log10_initialvalues <- initialvalues
  # Load results from IGS over narrow range
  load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
}

# Visualize results from wide IGS ####
plot_parameter_ranges(initialValue_ranges, log10_tabestimates %>% filter(!run %in% exclude_runs), p_threshold = 0.95)

distinct(log10_tabestimates, run, bicc) %>% 
  ggplot(aes("BICc", bicc, color = run %in% runs)) + 
  geom_beeswarm()

plot_bicc_v_estimate(log10_tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = T, color = runs,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_log10")

distinct(log10_tabestimates, run, bicc) %>% 
  ggplot(aes("BICc", bicc, color = run %in% new_runs)) + 
  geom_beeswarm()

plot_bicc_v_estimate(log10_tabestimates, remove_bicc_outliers = T,  color = new_runs,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_log10_low_m")

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter), run %in% new_runs), remove_bicc_outliers = T)
compare_estimates(tabestimates %>% filter(run %in% new_runs), remove_bicc_outliers = T)


# Visualize results from narrow IGS ####
## Look at parameter ranges from each step
plot_parameter_ranges(initialValue_ranges, log10_tabestimates, tabestimates, p_threshold = 0.95)


plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)


## Assess identifiability
correlate_estimates(tabestimates , remove_bicc_outliers = T, remove_estimate_outliers = T,
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = T,
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  print(parameter_plots(
    parameter, remove_bicc_outliers = F, remove_estimate_outliers = F,
    log10_tabestimates = log10_tabestimates, tabestimates = tabestimates,
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
      ggplot(aes(time, indivPred_mode, color = log10(!!as.symbol(param)), group = interaction(id, run))) + 
      geom_line(alpha = 0.5) + 
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1)
  )
}
dev.off()

## Plot population fits
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "population_trajectories_by_param_filtered_runs.pdf"))
indivFit %>% 
  filter(run %in% runs) %>% 
  left_join(tabestimates %>% distinct(run, bicc)) %>%
  ggplot(aes(time, popPred, color = bicc, group = interaction(run))) + 
  geom_line(alpha = 0.5) + 
  # facet_wrap(~id) +
  scale_color_viridis_c(direction = -1) 

for(param in names(indivParams)[2:7]){
  print(
    indivFit %>% 
      filter(run %in% runs) %>% 
      left_join(tabestimates %>% filter(parameter == str_interp("${param}_pop"))) %>%
      ggplot(aes(time, popPred, color = log10(estimate), group = interaction( run))) + 
      geom_line(alpha = 0.5) + 
      # facet_wrap(~id) +
      scale_color_viridis_c(direction = -1) + 
      labs(color = param)
  )
}
dev.off()


# Run model with initial guesses set from best fit:
igs <- initialvalues %>% filter(run == unique(filter(tabestimates, bicc == min(bicc))$run))
set_base(igs %>% rename(estimate = initialValue))
runScenario()
saveProject(here(out_folder, project_name))


## Decide on parameter bounds for profile likelihood #####

tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate))

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(parameter == "lambda_pop") %>% 
  ggplot(aes(estimate, bicc, color = a)) + 
  geom_point() + 
  scale_x_log10() +
  geom_vline(xintercept = 50) +
  geom_vline(xintercept = 150) +
  coord_cartesian(xlim = c(10, 1e3))

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(parameter == "m_pop") %>% 
  ggplot(aes(estimate, bicc, color = a)) + 
  geom_point() + 
  scale_x_log10() +
  geom_vline(xintercept = 10) +
  geom_vline(xintercept = 1000) +
  coord_cartesian(xlim = c(1, 1e3))

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(parameter == "p_pop") %>% 
  ggplot(aes(estimate, bicc, color = a)) + 
  geom_point() +
  # scale_x_log10() +
  coord_cartesian(xlim = c(1e-2, 50)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 2)

#####
igs <- data.frame(
  parameter = unique(tabestimates$parameter),
  initialValue = 1,
  method = "MLE"
) %>% 
  mutate(initialValue = case_when(
    parameter == "betas_pop" ~ 15,
    parameter == "lambda_pop" ~ 50,
    parameter == "m_pop" ~ 10,
    parameter == "NPs_pop" ~ 1,
    parameter == "p_pop" ~ 1,
    parameter == "ta_pop" ~ 0.2,
    .default = initialValue
  ))


rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(grepl("pop", parameter), bicc < 208, bicc > 195,
         !(parameter == "m_pop" & estimate < 1e-1)) %>% 
  ggplot(aes(estimate, bicc, color = a )) + 
  geom_point(alpha = 0.5) + 
  scale_x_log10() + 
  geom_smooth(se = FALSE) +
  geom_smooth(color = "black")  + 
  facet_wrap(~parameter, scales = "free_x") +
  geom_vline(data = igs %>% filter(grepl("pop", parameter)),
             aes(xintercept = initialValue), linetype = "dashed", color = "black")


## Run profile likelihood #####
if(RERUN){
  lambda_profile <- profile_likelihood("lambda_pop", c(20,1000), igs, nits = 99, record_predictions = T)
  write_csv(lambda_profile$tabestimates, here(out_folder, "lambda_profile_likelihood.csv"))
  write_csv(lambda_profile$predictions, here(out_folder, "lambda_profile_likelihood_table.csv"))
  
  p_profile <- profile_likelihood("p_pop", c(0.1,5), igs, nits = 50, record_predictions = T)
  write_csv(p_profile$tabestimates, here(out_folder, "p_profile_likelihood.csv"))
  write_csv(p_profile$predictions, here(out_folder, "p_profile_likelihood_table.csv"))
  
  m_profile <- profile_likelihood("m_pop", c(10, 1000), igs, nits = 100, record_predictions = T)
  write_csv(m_profile$tabestimates, here(out_folder, "m_profile_likelihood.csv"))
  write_csv(m_profile$predictions, here(out_folder, "m_profile_likelihood_table.csv"))
}else{
  lambda_profile <- read_csv( here(out_folder, "lambda_profile_likelihood.csv"))
  p_profile <- read_csv( here(out_folder, "p_profile_likelihood.csv"))
  m_profile <- read_csv( here(out_folder, "m_profile_likelihood.csv"))
}