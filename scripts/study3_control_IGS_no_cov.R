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
study <- "study3_control"
project_name <- "IGS_no_cov.mlxtran"
out_folder <- here("results", study)
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_final.txt")

## Data file
data_file <- here("data", "derived", "STUDY3_CONTROL_DATA.txt")
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
  ## Create Project
  create_project(model_file, data_file, headerTypes, error_model,
                 exclude_random_effects, treatment_effects = NULL, normal_dist = NULL)
  
  saveProject(here(out_folder, project_name))
  
  # Initial guess sampling over orders of magnitude ####
  initialValue_ranges_log10 <- initialValue_ranges %>% 
    t() %>% as.data.frame %>% 
    setNames(c("min", "max")) %>% 
    rownames_to_column("parameter") %>% 
    mutate(parameter = paste0(parameter, "_pop")) 
  
  ## Add ranges from 1 to 5 for random effect standard deviations
  parameters <- getPopulationParameterInformation()$name
  initialValue_ranges_log10 <- rbind(
    initialValue_ranges_log10,
    tibble(parameter = parameters[!grepl("pop", parameters)],
           min = 1, max = 5)
  ) %>% 
    column_to_rownames("parameter")
  
  
  ## sample new initial estimates
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
  
  # Initial guess sampling from a narrower range of values ####
  
  ## Calculate narrower ranges for IGS:
  initialValue_ranges_narrow <- get_narrow_ranges(tabestimates, p_threshold = 0.95) %>% 
    column_to_rownames("parameter") %>% select(-best)
  
  plot_parameter_ranges(initialValue_ranges, log10_tabestimates, NULL, 0.95)
  
  write_csv(initialValue_ranges_narrow %>% rownames_to_column("parameter"), here(out_folder, "initialValue_ranges_narrow.csv"))
  
  ## initial guess sampling from narrow range
  # Use initial condition of 1 for all random effects
  parameter_results <- initial_guess_sampling(
    initialValue_ranges_narrow, 
    nits = 100,
    record_shrinkage = T,
    record_indiv_params = T,
    set_variance = 1
  )
  
  ## Save results
  for(name in names(parameter_results)){
    assign(name, parameter_results[[name]])
  }

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
## Look at parameter ranges from each step
plot_parameter_ranges(initialValue_ranges, log10_tabestimates, p_threshold = 0.95)

## Color by "best" runs according to BICc
runs <- calculate_BICw(log10_tabestimates) %>% filter(cumprob <= 0.95) %>% pull(run)
plot_bicc_v_estimate(log10_tabestimates, remove_bicc_outliers = T, color = runs,
                     out_folder = out_folder, project_name = project_name, narrow = F)



# Visualize results from narrow IGS ####
## Look at parameter ranges from each step
plot_parameter_ranges(initialValue_ranges, log10_tabestimates, tabestimates, p_threshold = 0.95)

plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

# filter to best runs
runs <- calculate_BICw(tabestimates) %>% filter(cumprob <= 0.9) %>% pull(run)
distinct(tabestimates, run, bicc) %>% ggplot(aes("BICc", bicc, color = run %in% runs)) + ggbeeswarm::geom_beeswarm()
plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F, color=runs)
                     
## Assess identifiability
correlate_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")
compare_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")
estimate_boxplots(tabestimates  %>% filter(run %in% runs), 
                  out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

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
    ggplot(aes(time, indivPred_mode, color = !!as.symbol(param), group = interaction(id, run))) + 
      geom_point(data = . %>% filter(run == 1), alpha = 0.5, inherit.aes = F, aes(time, vl)) +
      geom_line(alpha = 0.5) + 
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



# Run profile likelihood for lambda, p, and m ####
if(RERUN){
  igs <- tabestimates %>% filter(bicc == min(bicc)) %>% mutate(estimate = signif(estimate, 1)) 
  
  lambda_profile <- profile_likelihood("lambda_pop", c(40,400), igs, nits = 60)
  write_csv(lambda_profile, here(out_folder, "lambda_profile_likelihood.csv"))
  
  p_profile <- profile_likelihood("p_pop", c(0.5,5), igs, nits = 60)
  write_csv(p_profile, here(out_folder, "p_profile_likelihood.csv"))
  
  m_profile <- profile_likelihood("m_pop", c(10, 160), igs, nits = 60)
  write_csv(m_profile, here(out_folder, "m_profile_likelihood.csv"))

}else{
  lambda_profile <- read_csv( here(out_folder, "lambda_profile_likelihood.csv"))
  p_profile <- read_csv( here(out_folder, "p_profile_likelihood.csv"))
  m_profile <- read_csv( here(out_folder, "m_profile_likelihood.csv"))
}

# Visualize Profile Likelihood results ####
## Lambda
lp <- lambda_profile %>% 
  distinct(lambda_pop, bicc) %>%
  ggplot(aes(lambda_pop, bicc)) + 
  geom_point() 

# Best lambda is 125
lambda_profile %>% filter(bicc == min(bicc))


## p
pp <- p_profile %>% 
  distinct(p_pop, bicc) %>%
  ggplot(aes(p_pop, bicc)) + 
  geom_point()

# Best p is 1
p_profile %>% filter(bicc == min(bicc))


## m
# The likelihood is very flat
mp <- m_profile %>%
  distinct(m_pop, bicc) %>%
  ggplot(aes(m_pop, bicc)) + 
  geom_point() 

# Best m is 43
m_profile %>% filter(bicc == min(bicc))

lp + pp + mp + 
  plot_layout(guides = "collect", axis_titles = "collect_y") & ylim(c(254, 263)) &
  plot_annotation(title = "Profile Likelihood Results From Study 3 Animals") & 
  theme(plot.title = element_text(hjust = 0.5)) &
  geom_smooth()
