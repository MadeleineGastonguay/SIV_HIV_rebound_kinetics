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
study <- "study4_control"
project_name <- "IGS_no_cov.mlxtran"
out_folder <- here("results_eclipse", study)
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "derived", "STUDY4_CONTROL_DATA.txt")
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
  
  write_csv(initialValue_ranges_narrow %>% rownames_to_column("parameter"), here(out_folder, "initialValue_ranges_narrow.csv"))
  
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
  
  tabestimates <- rbind(tabestimates, parameter_results2$tabestimates %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  tabiters <- rbind(tabiters, parameter_results2$tabiters %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  phase <- rbind(phase, parameter_results2$phase)
  initialvalues <- rbind(initialvalues, parameter_results2$initialvalues %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  shrinkage <- rbind(shrinkage, parameter_results2$shrinkage %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  indivParams <- rbind(indivParams, parameter_results2$indivParams %>% mutate(run = run + max(parameter_results$tabestimates$run)))
  indivFit <- rbind(indivFit, parameter_results2$indivFit %>% mutate(run = run + max(parameter_results$tabestimates$run)))

  save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit,
       file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
  
  # Save best estimate 
  write_csv(tabestimates %>% filter(bicc == min(bicc)),
            here(out_folder, str_remove(project_name, "\\.mlxtran$"), "best_estimates.csv"))
  
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

## Color by "best" runs according to BICc
runs <- calculate_BICw(log10_tabestimates) %>% filter(cumprob <= 0.95) %>% pull(run)
distinct(log10_tabestimates, run, bicc) %>% 
  ggplot(aes("BICc", bicc, color = run %in% runs)) + 
  geom_beeswarm()

plot_bicc_v_estimate(log10_tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = T,
                     color = runs,
                     out_folder = out_folder, project_name = project_name, narrow = F,
                     file_addendum = "_log10")

## Look at parameter ranges from each step
plot_parameter_ranges(initialValue_ranges, log10_tabestimates, p_threshold = 0.95)

## Assess identifiability
correlate_estimates(log10_tabestimates %>% filter(grepl("pop", parameter), run %in% runs), remove_bicc_outliers = T)
compare_estimates(log10_tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T)



# Visualize results from narrow IGS ####
## Look at parameter ranges from each step
plot_parameter_ranges(initialValue_ranges, log10_tabestimates, tabestimates, p_threshold = 0.95)

plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = F, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

tabestimates %>% filter(parameter %in% c("lambda_pop", "betas_pop")) %>% 
  group_by(run, bicc) %>% summarise(prod = prod(estimate)) %>% 
  ggplot(aes(prod, bicc)) + 
  geom_point()

## Assess identifiability
correlate_estimates(tabestimates, remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name)
compare_estimates(tabestimates, remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name)
estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  print(parameter_plots(
    parameter, remove_bicc_outliers = T, remove_estimate_outliers = F,
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
    ggplot(aes(time, indivPred_mode, color = (!!as.symbol(param)), group = interaction(id, run))) + 
      geom_point(data = . %>% filter(run == min(run)), alpha = 0.5, inherit.aes = F, aes(time, vl)) +
      geom_line(alpha = 0.5) + 
      geom_line(data = . %>% left_join(distinct(tabestimates, run, bicc) %>% filter(bicc == min(bicc))),
                lty = "dashed", color = "black") +
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

## Make sure best estimates are reasonable in terms of P:
best <- tabestimates %>% filter(bicc == min(bicc)) %>% pull(run) %>% unique
initialvalues %>% filter(run == best) %>% select(-run) %>% rename(name = parameter) %>% setPopulationParameterInformation()
runScenario()

variables <- read_csv(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "tables.txt"))
variables %>% 
  ggplot(aes(time, P_mode)) + 
  geom_line() + 
  facet_wrap(~id, scales = "free_y") + 
  scale_y_log10() + 
  geom_line(data = indivFit %>% filter(run == best), aes(y = 10^indivPred_mode, color = "vl"))

# Decide on profile likelihood bounds
tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate))

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(parameter == "lambda_pop", bicc < 380) %>% 
  ggplot(aes(estimate, bicc, color = a)) + 
  geom_point() + 
  scale_x_log10() + 
  geom_vline(xintercept = 200) +
  geom_vline(xintercept = 2000) + 
  coord_cartesian(xlim = c(1e2, 1e4))

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(parameter == "m_pop", bicc < 380) %>% 
  ggplot(aes(estimate, bicc, color = a)) + 
  geom_point() + 
  scale_x_log10() + 
  geom_vline(xintercept = 700) +
  geom_vline(xintercept = 3000)

#####
igs <- data.frame(
  parameter = unique(tabestimates$parameter),
  initialValue = 1,
  method = "MLE"
) %>% 
  mutate(initialValue = case_when(
    parameter == "betas_pop" ~ 30,
    parameter == "lambda_pop" ~ 700,
    parameter == "m_pop" ~ 1500,
    parameter == "NPs_pop" ~ 0.2,
    parameter == "p_pop" ~ 0.7,
    parameter == "ta_pop" ~ 8,
    .default = initialValue
  ))

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(grepl("pop", parameter), bicc < 375,
         !(parameter == "betas_pop" & estimate < 1)) %>%
  ggplot(aes(estimate, bicc, color = a )) + 
  geom_point(alpha = 0.5) + 
  scale_x_log10() + 
  geom_smooth(se = FALSE) +
  geom_smooth(color = "black")  + 
  facet_wrap(~parameter, scales = "free_x") +
  geom_vline(data = igs %>% filter(grepl("pop", parameter)),
             aes(xintercept = initialValue), lty = "dashed", color = "black")

rbind(log10_tabestimates %>% mutate(a = "log10"),
      tabestimates %>% mutate(a = "narrow")) %>% 
  filter(grepl("pop", parameter)) %>% 
  group_by(parameter, a) %>% 
  filter(bicc == min(bicc))



# Run profile likelihood for lambda, p, and m ####
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

# Visualize Profile Likelihood results ####
## Lambda
lp <- lambda_profile %>% 
  filter(parameter == "m_pop") %>% 
  # distinct(lambda_pop, bicc) %>% 
  ggplot(aes(lambda_pop, bicc, color = estimate > 70)) + 
  geom_point() + 
  labs(color = "m_pop > 70") 

# Best lambda is 235
lambda_profile %>% filter(bicc == min(bicc))

# Filter out big m range
lambda_profile %>% 
  filter(estimate <= 50, parameter == "m_pop") %>% 
  # distinct(lambda_pop, bicc) %>%
  ggplot(aes(lambda_pop, bicc)) + 
  geom_point()

# Best lambda is 400
lambda_profile %>% filter(estimate <= 50, parameter == "m_pop") %>%  filter(bicc == min(bicc))

## p
pp <- p_profile %>% 
  filter(parameter == "m_pop") %>% 
  # distinct(p_pop, bicc) %>%
  ggplot(aes(p_pop, bicc, color = estimate > 70)) + 
  geom_point() +
  labs(color = "m_pop > 70")

# Best p is 3.7
p_profile %>% filter(bicc == min(bicc))

# Filter out big m range
p_profile %>% 
  filter(estimate <= 60, parameter == "m_pop") %>%
  distinct(p_pop, bicc) %>%
  ggplot(aes(p_pop, bicc)) + 
  geom_point() 

# Best p is 3
p_profile %>% filter(estimate <= 60, parameter == "m_pop") %>% filter(bicc == min(bicc))


## m
# The likelihood is very flat
mp <- m_profile %>%
  distinct(m_pop, bicc) %>%
  ggplot(aes(m_pop, bicc, color = m_pop < 70)) + 
  geom_point() + 
  annotate(geom = "rect", xmin = -Inf, xmax = Inf, ymin = min(m_profile$bicc), ymax = min(m_profile$bicc) + 2, alpha = 0.25) + 
  labs(color = "m_pop > 70")

lp + geom_smooth(data = . %>% filter(estimate < 70), show.legend = F) +
  pp + geom_smooth(data = . %>% filter(estimate < 70), show.legend = F) +
  mp + #geom_smooth() +
  plot_layout(guides = "collect", axis_titles = "collect_y") & ylim(c(430,465)) &
  plot_annotation(title = "Profile Likelihood Results From Study 1 Animals") & 
  theme(plot.title = element_text(hjust = 0.5))
