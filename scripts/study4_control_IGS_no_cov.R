# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
project_name <- "IGS_no_cov.mlxtran"
out_folder <- here("results", "study4_control")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_final.txt")

## Data file
data_file <- here("data", "derived","STUDY4_CONTROL_DATA.txt")
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


## Create project
create_project(model_file, data_file, headerTypes, error_model,
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL)

# saveProject(here(out_folder, project_name))
loadProject(here(out_folder, project_name))

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
# load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "log10_IGS.RData"))

log10_tabestimates <- tabestimates
log10_initialvalues <- initialvalues

runs <- calculate_BICw(log10_tabestimates) %>% filter(cumprob <= 0.95) %>% pull(run)
plot_bicc_v_estimate(log10_tabestimates, remove_bicc_outliers = T, color = runs)
distinct(log10_tabestimates, run, bicc) %>% ggplot(aes("BICc", bicc, color = run %in% runs)) + ggbeeswarm::geom_beeswarm()

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)), remove_bicc_outliers = T)

compare_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = T)


# Initial guess sampling from a narrower range of values ####
plot_parameter_ranges(initialValue_ranges, log10_tabestimates, p_threshold = 0.95)

## Calculate narrower ranges for IGS:
initialValue_ranges_narrow <- get_narrow_ranges(tabestimates, p_threshold = 0.95) %>% select(-best)

write_csv(initialValue_ranges_narrow, here(out_folder, "initialValue_ranges_narrow.csv"))

initialValue_ranges_narrow <- initialValue_ranges_narrow %>% 
  column_to_rownames("parameter")

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

# load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

# Run model with initial guesses set from best fit:
igs <- initialvalues %>% filter(run == unique(filter(tabestimates, bicc == min(bicc))$run))
set_base(igs %>% rename(estimate = initialValue))
runScenario()
saveProject(here(out_folder, project_name))

## Find runs that converge to best BICc:
runs <- calculate_BICw(tabestimates) %>% filter(cumprob <= .95) %>% pull(run)
ggplot(distinct(tabestimates, run ,bicc), aes("", bicc)) + 
  ggbeeswarm::geom_beeswarm(aes(color = run %in% runs))

plot_bicc_v_estimate(tabestimates, color = runs)

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


plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

plot_bicc_v_estimate(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_filtered_runs")

correlate_estimates(tabestimates , remove_bicc_outliers = T, remove_estimate_outliers = T,
                    out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, remove_estimate_outliers = F,
                    out_folder = out_folder, project_name = project_name, file_addendum = "_filtered_runs")

compare_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = T,
                  out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, remove_estimate_outliers = T,
                  out_folder = out_folder, project_name = project_name, file_addendum = "_filtered_runs")

estimate_boxplots(tabestimates %>% filter(run %in% runs), 
                  out_folder = out_folder, project_name = project_name, file_addendum = "_filtered_runs")

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


plot_parameter_ranges(initialValue_ranges, log10_tabestimates, tabestimates, p_threshold = 0.7)

# Run profile likelihood for lambda and p ####
igs <- tabestimates %>% filter(bicc == min(bicc)) %>% mutate(estimate = signif(estimate, 1)) 

lambda_profile <- profile_likelihood("lambda_pop", c(5,200), igs, nits = 60)

write_csv(lambda_profile, here(out_folder, "lambda_profile_likelihood.csv"))
# lambda_profile <- read_csv( here(out_folder, "lambda_profile_likelihood.csv"))

lp <- lambda_profile %>% 
  filter(parameter == "ta_pop") %>%
  ggplot(aes(lambda_pop, bicc)) + 
  geom_point()

lambda_profile %>% 
  ggplot(aes(lambda_pop, estimate, color = bicc)) + 
  geom_point() + 
  facet_wrap(~parameter, scales = "free_y") + 
  color_bicc()

# Best lambda is 60
lambda_profile %>% filter(bicc == min(bicc))


# Profile likelihood for p ###

p_profile <- profile_likelihood("p_pop", c(0.4,2), igs, nits = 50)
write_csv(p_profile, here(out_folder, "p_profile_likelihood.csv"))
# p_profile <- read_csv( here(out_folder, "p_profile_likelihood.csv"))

pp <- p_profile %>% 
  filter(parameter == "ta_pop") %>%
  # distinct(p_pop, bicc) %>%
  ggplot(aes(p_pop, bicc)) + 
  geom_point() 

# Best p is 0.7
p_profile %>% filter(bicc == min(bicc))

p_profile %>% 
  ggplot(aes(p_pop, estimate, color = bicc)) + 
  geom_point() + 
  facet_wrap(~parameter, scales = "free_y") + 
  color_bicc()

m_profile <- profile_likelihood("m_pop", c(20, 80), igs, nits = 60)
m_profile2 <- profile_likelihood("m_pop", c(1, 20), igs, nits = 30)
m_profile <- rbind(m_profile, m_profile2)
write_csv(m_profile, here(out_folder, "m_profile_likelihood.csv"))

mp <- m_profile %>% 
  distinct(m_pop, bicc) %>%
  ggplot(aes(m_pop, bicc)) + 
  geom_point() 

# best m is 47
m_profile %>% filter(bicc == min(bicc))

lp + pp + mp + 
  plot_layout( axis_titles = "collect_y", guides = "collect") & ylim(c(366,373)) &
  plot_annotation(title = "Profile Likelihood Results From Study 1 Control Animals") & 
  theme(plot.title = element_text(hjust = 0.5)) &
  labs(y = "BICc") 
