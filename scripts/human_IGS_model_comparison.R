# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "human"
out_folder <- here("results", study, "IGS_model_comparison")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_HIV_final_twash.txt")

## Data file
data_file <- here("data", "derived", "HUMAN_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "ignore", "observation", "cens", "catcov")


## Random effects to exclude
exclude_random_effects <- "tw" # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 200, "p_pop" = 3, "tw_pop" = 3)

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]


## Sample initial values for all parameters #####
# We will sample population parameters from the narrow range determined from first IGS 
# We will set all standard deviations to 1 and treatment effects to 0
initialValue_ranges_narrow <- read_csv(here("results", study, "initialValue_ranges_narrow.csv")) %>% 
  column_to_rownames("parameter")

parameters <- rownames(initialValue_ranges_narrow)
parameters <- c(parameters, "tw_pop")

initialvalues <- NULL
for(i in 1:100){
  popini <- sapply(parameters, function(j){
    if(j %in% names(fixed)){
      fixed[[j]]
    }else if(grepl("pop", j)){
      runif(1, initialValue_ranges_narrow[j,1], initialValue_ranges_narrow[j,2])  
    }else{
      1
    }
  })
  
  # set sampled values as new initial estimates
  newpopparams <- data.frame(parameter = parameters, initialValue = popini) %>% 
    mutate(method = ifelse(parameter %in% names(fixed), "FIXED", "MLE")) %>% 
    remove_rownames()
  initialvalues <- rbind(initialvalues, newpopparams %>% mutate(run = i))
}


if(!file.exists(here(out_folder, "initial_values.csv"))) write_csv(initialvalues, here(out_folder, "initial_values.csv"))

predefined_initialvalues <- read_csv(here(out_folder, "initial_values.csv"))

predefined_initialvalues %>% 
  ggplot(aes(parameter, initialValue)) +
  geom_boxplot() + 
  facet_wrap(~parameter, scales = "free")

# Run IGS #####

# am1 #####
# NNRTI effect on twash
alt_mod <- "am1"
loadProject(here("results", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_result_am1 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = 100, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, record_condition_number = T
)


for(name in names(parameter_result_am1)){
  assign(name, parameter_result_am1[[name]])
}

save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit, condition_number,
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

# load(here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  if(parameter %in% names(fixed)) next
  print(parameter)
  print(
    parameter_plots(
      parameter,
      remove_bicc_outliers = F, 
      remove_estimate_outliers = F,
      NULL, tabestimates, tabiters,
      phase, initialvalues
    )
  )
}
dev.off()

plot_vs_rses(tabestimates, out_folder, project_name)

plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

plot_bicc_v_estimate(tabestimates)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates,
                  out_folder = out_folder, project_name = project_name)


# AM2 #####
# Model without any covariates
alt_mod <- "am2"
loadProject(here("results", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_result_am2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = 100, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, record_condition_number = T
)


for(name in names(parameter_result_am2)){
  assign(name, parameter_result_am2[[name]])
}

save(tabestimates, tabiters, phase, initialvalues, shrinkage, indivParams, indivFit, condition_number,
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

# load(here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  if(parameter %in% names(fixed)) next
  print(parameter)
  print(
    parameter_plots(
      parameter,
      remove_bicc_outliers = F, 
      remove_estimate_outliers = F,
      NULL, tabestimates, tabiters,
      phase, initialvalues
    )
  )
}
dev.off()

plot_vs_rses(tabestimates, out_folder, project_name)

plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name,
  
)


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)



# Compare results ####

ams <- c("am1", "am2")

all_tabestimates <- NULL
all_initialvalues <- NULL
all_pvals <- NULL
all_indivParams <- NULL
all_indivFit <- NULL
all_condition_number <- NULL 
for(am in ams){
  load(here(out_folder, str_interp("IGS_${am}"), "IGS.RData")) 
  all_tabestimates <- rbind(all_tabestimates, tabestimates %>% mutate(am = am))
  all_initialvalues <- rbind(all_initialvalues, initialvalues %>% mutate(am = am))
  # if(am == "am2") pvals <- data.frame()
  # all_pvals <- rbind(all_pvals, pvals %>% mutate(am = am))
  all_indivParams <- rbind(all_indivParams, indivParams %>% mutate(am = am))
  all_indivFit <- rbind(all_indivFit, indivFit %>% mutate(am = am))
  all_condition_number <- rbind(all_condition_number, condition_number %>% mutate(am = am))
}

all_tabestimates <- all_tabestimates %>% 
  mutate(am = factor(am, levels = ams)) %>% 
  group_by(run) %>% 
  mutate(best = am[which.min(bicc)]) %>% 
  ungroup

# Look at results #####

all_tabestimates %>% select(parameter, am, estimate, best, run)  %>% 
  left_join(all_initialvalues) %>%
  filter(grepl("pop", parameter)) %>% 
  # mutate(am2 = factor(am, levels = ams[-1], labels = c("Base model", "Alternative Np", "Alternative ta", "Alternative m"))) %>% 
  ggplot(aes(initialValue, estimate, color = am, shape = am == best)) + 
  geom_point(alpha = 0.75, size = 2) + 
  facet_wrap(~parameter, scales = "free") + 
  labs(x = "Initial guess", y = "Estimate", color = "Model", shape = "Best Fit") 

# Compare model performances ####

## out of 100 runs, the number of times each model was the best according to BICc
all_tabestimates %>% 
  distinct(run, best) %>% 
  count(best)

# best      n
# am1      20
# am2      80


## out of 100 runs, the number of times each model was less than 2 BICc from the best model
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  rename(best_model = best) %>% 
  group_by(run) %>% 
  mutate(best = min(bicc),
         equal_to_best = abs(bicc - best) <= 2) %>% 
  ungroup %>% 
  filter(equal_to_best) %>% 
  count(am, equal_to_best) 

# am    equal_to_best     n
# am1   TRUE             43
# am2   TRUE             93

## out of 100 runs, AM1 (twash effect) significantly out-performed the model without covariates only 7 times
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  pivot_wider(names_from = am, values_from = bicc) %>% 
  count(am2 > (am1 + 2))

# `am2 > (am1 + 2)`     n
#  FALSE                93
#  TRUE                  7

all_tabestimates_bicc <- all_tabestimates %>% distinct(run, best, am, ofv, bicc)

ncovar <- all_tabestimates %>% filter(grepl("beta_", parameter)) %>% distinct(parameter, am) %>% count(am) 

all_tabestimates_bicc %>% 
  group_by(run) %>% 
  mutate(equal_to_best = bicc <= min(bicc) + 2) %>% 
  ggplot(aes(am, bicc)) + 
  geom_violin() + 
  ggbeeswarm::geom_beeswarm(alpha = 0.75, size = 2,
                            aes(color = ifelse(am == best, "Best Model", ifelse(equal_to_best, "Equal to best model", "Not the best")))) + 
  scale_color_manual(values = c("blue", "darkorchid", "darkgray")) +
  labs(x = "Model", y = "BICc", color = "") +
  theme(legend.position = "bottom") +
  coord_flip() 


# How many times does AM2 outperform AM1?
plot_df <- all_tabestimates %>% 
  distinct(run, best, bicc, am) %>% 
  pivot_wider(names_from = am, values_from = bicc) %>%
  mutate(color = case_when(am2 > (am1 + 2) ~ "Alternative model better", am1 > (am2 + 2) ~ "SCM model better", TRUE ~ "Equal performance"))

plot_df %>%
  ggplot(aes(am2, am1)) + 
  
  geom_abline(alpha = 0.7) +
  geom_abline(intercept = -2, lty = "dashed", alpha = 0.7) +
  geom_abline(intercept = 2, lty = "dashed", alpha = 0.7)  +
  
  # Shaded regions using geom_ribbon
  geom_ribbon(aes(ymin = -Inf, ymax = am2 - 2), fill = "darkred", alpha = 0.1) +  # Below bottom dashed line
  annotate("text", x = 1435, y = 1422, label = str_interp("Alternative\nmodel better\n(${count(plot_df, color)$n[1]})"), size = 4, color = "darkred", fontface = "bold") +
  
  geom_ribbon(aes(ymin = am2 - 2, ymax = am2 + 2), fill = "darkgreen", alpha = 0.1) +  # Between dashed lines
  annotate("text", x = 1436, y = 1436, label = str_interp("Equal\nperformance\n(${count(plot_df, color)$n[2]})"), size = 4, color = "darkgreen", fontface = "bold") +
  
  geom_ribbon(aes(ymin = am2 + 2, ymax = Inf), fill = "darkgray", alpha = 0.1) +   # Above top dashed line
  annotate("text", x = 1427, y = 1442, label = str_interp("SCM model better\n(${count(plot_df, color)$n[3]})"), size = 4, color = "darkgray", fontface = "bold") +
  
  geom_point(aes(color = color), alpha = 0.7, size = 2) +
  
  labs(x = "No covariates", y = "NNRTI on twash", shape = "") +
  scale_color_manual(values = c("darkred", "darkgreen", "darkgray")) + 
  scale_shape_manual(values = c(15, 17)) +
  guides(color = "none") 

# Check condition numbers

all_condition_number %>% 
  left_join(distinct(tabestimates, run, bicc)) %>% 
  ggplot(aes(condition_number, am, color = bicc)) + geom_violin(scale = "width") + geom_beeswarm(alpha = 0.5)  + 
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = NA) +
  scale_x_log10() + 
  geom_vline(xintercept = 100, linetype = 2) +
  color_bicc() + 
  geom_point(data = . %>% group_by(am) %>% filter(bicc == min(bicc)), size = 2, shape = 17, alpha = 1, color = "black")




## Save best results from best model! ####

best_run <- all_tabestimates %>% filter(am == "am2") %>% filter(bicc == min(bicc)) %>% pull(run) %>% unique

igs <- all_initialvalues %>% filter(am == "am2", run == best_run) %>% select(-am, -run)

loadProject(here("results", "human", "SCM", "am2.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)
saveProject(here("results", "human", "best_model.mlxtran"))

runScenario()

saveProject(here("results", "human", "best_model.mlxtran"))

generate_report(here("results", "human", "best_model.mlxtran"))
