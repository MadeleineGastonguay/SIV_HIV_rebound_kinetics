# Load libraries
library(here)
library(ggh4x)
library(RColorBrewer)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
project_name <- "IGS_no_cov_fixed_param_removeRE_p.mlxtran"
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
exclude_random_effects <- c("p") # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

# Parameters to fix
fixed <- c("lambda_pop" = 200, "p_pop" = 3)

## Create project
create_project(model_file, data_file, headerTypes, error_model,
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL, fixed= fixed)

saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))


# Initial guess sampling over orders of magnitude ####

## Calculate narrower ranges for IGS:
initialValue_ranges_narrow <- read_csv(here(out_folder, "initialValue_ranges_narrow.csv")) %>% 
  filter(!parameter %in% paste("omega", exclude_random_effects, sep = "_"))

initialValue_ranges_narrow <- initialValue_ranges_narrow %>% 
  column_to_rownames("parameter")

## initial guess sampling from narrow range
# Use initial condition of 1 for all random effects
parameter_results <- initial_guess_sampling(
  initialValue_ranges_narrow, 
  nits = 10,
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

# # Run model with initial guesses set from best fit:
# igs <- initialvalues %>% filter(run == unique(filter(tabestimates, bicc == min(bicc))$run))
# set_base(igs %>% rename(estimate = initialValue))
# runScenario()
# saveProject(here(out_folder, project_name))

indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams)
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
      ggplot(aes(time, indivPred_mode, color = !!as.symbol(param), group = interaction(id, run))) + 
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

correlate_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = F,
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = F,
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)


# calculate_BICw(tabestimates) %>% ggplot(aes(bicc, cumprob)) + geom_point()
# calculate_BICw(tabestimates) %>%  filter(cumprob <= 0.7) %>% pull(run) -> runs
runs <- tabestimates %>% filter(bicc <= 195) %>% pull(run) %>% unique
best_low_m <- calculate_BICw(tabestimates %>% filter(!run %in% runs)) %>% filter(cumprob <= 0.95) %>% pull(run)

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "trajectories_by_regime.pdf"))
tabestimates %>% 
  mutate(regime = case_when(run %in% runs ~ 'Best high m', run %in% best_low_m ~ "Best low m", .default = "Sub-optimal BICc")) %>% 
  ggplot(aes(estimate, bicc, color = regime)) + 
  geom_point() + 
  facet_wrap(~parameter, scales = "free_x") + 
  scale_color_brewer(palette = "Set1")

indivFit %>% 
  left_join(tabestimates %>% distinct(run, bicc)) %>%
  mutate(regime = case_when(run %in% runs ~ 'Best high m', run %in% best_low_m ~ "Best low m", .default = "Sub-optimal BICc")) %>% 
  ggplot(aes(time, popPred, color = bicc, group = interaction(run))) + 
  geom_point(data = . %>% group_by(regime) %>% filter(run == min(run)) %>% ungroup,
             inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  scale_color_viridis_c(direction = -1) +
  facet_wrap2(~regime, strip = strip_themed(background_x = elem_list_rect(fill = grDevices::adjustcolor(brewer.pal(3, "Set1"), alpha.f = 0.5))) ) 

indivDF %>% 
  mutate(regime = case_when(run %in% runs ~ 'Best high m', run %in% best_low_m ~ "Best low m", .default = "Sub-optimal BICc")) %>% 
  ggplot(aes(time, indivPred_mode, color = regime, group = interaction(id, run))) + 
  geom_point(data = . %>% filter(run == 1),
             inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.1) + 
  facet_wrap(~id) + 
  scale_color_brewer(palette = "Set1") + 
  guides(color = guide_legend(override.aes = list(alpha = 1)))
dev.off()

distinct(tabestimates, run, bicc) %>% ggplot(aes("BICc", bicc, color = run %in% runs)) + ggbeeswarm::geom_beeswarm()

plot_bicc_v_estimate(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

plot_bicc_v_estimate(tabestimates %>% filter(run %in% best_low_m), remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_best_low_m")

correlate_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, remove_estimate_outliers = F,
                    out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

correlate_estimates(tabestimates %>% filter(run %in% best_low_m), remove_bicc_outliers = T, remove_estimate_outliers = F,
                    out_folder = out_folder, project_name = project_name, file_addendum = "_best_low_m")

compare_estimates(tabestimates %>% filter(run %in% runs), remove_bicc_outliers = T, remove_estimate_outliers = F,
                  out_folder = out_folder, project_name = project_name, file_addendum = "_best_runs")

compare_estimates(tabestimates %>% filter(run %in% best_low_m), remove_bicc_outliers = T, remove_estimate_outliers = F,
                  out_folder = out_folder, project_name = project_name, file_addendum = "_best_low_m")

estimate_boxplots(tabestimates  %>% filter(run %in% runs), out_folder = out_folder, project_name = project_name , file_addendum = "_best_runs")

estimate_boxplots(tabestimates  %>% filter(run %in% best_low_m), out_folder = out_folder, project_name = project_name , file_addendum = "_best_low_m")

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



# Plot shrinkage to evaluate random effects
# Shrinkage is large for beta and p 
shrinkage %>% 
  ggplot(aes(parameter, as.numeric(conditionalMode))) + 
  geom_boxplot()


best_estimates <- tabestimates %>% filter(bicc == min(bicc))
write_csv(best_estimates, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "best_estimates.csv"))
