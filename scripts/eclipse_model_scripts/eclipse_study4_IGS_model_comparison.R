# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "study4"
out_folder <- here("results_eclipse", study, "IGS_model_comparison")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", ncol(data) - 4))


## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 100, "p_pop" = 0.9)

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]


## Sample initial values for all parameters #####
# We will sample population parameters from the narrow range determined from first IGS 
# We will set all standard deviations to 1 and treatment effects to 0
initialValue_ranges_narrow <- read_csv(here("results_eclipse", study, "initialValue_ranges_narrow.csv")) %>% 
  column_to_rownames("parameter")

parameters <- rownames(initialValue_ranges_narrow)
nits <- 50

if(!file.exists(here(out_folder, "initial_values.csv"))){
  initialvalues <- NULL
  for(i in 1:nits){
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
  
  write_csv(initialvalues, here(out_folder, "initial_values.csv"))
}

predefined_initialvalues <- read_csv(here(out_folder, "initial_values.csv"))

predefined_initialvalues %>% 
  ggplot(aes(parameter, initialValue)) +
  geom_boxplot() + 
  facet_wrap(~parameter, scales = "free")

# Run IGS #####

# Final model from SCM #####
alt_mod <- "final_model"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_final_model <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_final_model)){
  assign(name, parameter_results_final_model[[name]])
}

save(list = names(parameter_results_final_model),
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


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)





# AM1 #####
alt_mod <- "am1"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am1 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am1)){
  assign(name, parameter_results_am1[[name]])
}

save(list = names(parameter_results_am1),
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


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)





# AM3 #####
alt_mod <- "am3"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am3 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3)){
  assign(name, parameter_results_am3[[name]])
}

save(list = names(parameter_results_am3),
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


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)


# AM3.1 #####
alt_mod <- "am3.1"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am3.1 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3.1)){
  assign(name, parameter_results_am3.1[[name]])
}

save(list = names(parameter_results_am3.1),
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


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)





# AM3 without fixed parameters #####
alt_mod <- "am3"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))

project_name <- str_interp("IGS_${alt_mod}_no_fixed.mlxtran")
saveProject(here(out_folder, project_name))
popparam <- getPopulationParameterInformation()
popparam$method <- "MLE"
setPopulationParameterInformation(popparam)
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am3_no_fixed <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = NULL, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3_no_fixed)){
  assign(name, parameter_results_am3_no_fixed[[name]])
}

save(list = names(parameter_results_am3_no_fixed),
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


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)


# AM3.1 without fixed parameters #####
alt_mod <- "am3.1"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))

project_name <- str_interp("IGS_${alt_mod}_no_fixed.mlxtran")
saveProject(here(out_folder, project_name))
popparam <- getPopulationParameterInformation()
popparam$method <- "MLE"
setPopulationParameterInformation(popparam)
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am3.1_no_fixed <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = NULL, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3.1_no_fixed)){
  assign(name, parameter_results_am3.1_no_fixed[[name]])
}

save(list = names(parameter_results_am3.1_no_fixed),
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


correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)




# Compare results ####

ams <- c("final_model", "am1", "am3", "am3.1")#, "am3_no_fixed", "am3.1_no_fixed")

ams <- c("am3_no_fixed", "am3.1_no_fixed")

all_tabestimates <- NULL
all_initialvalues <- NULL
all_pvals <- NULL
all_indivParams <- NULL
all_indivFit <- NULL
all_condition_numbers <- NULL
all_table_predictions <- NULL
all_shrinkage <- NULL
for(am in ams){
  load(here(out_folder, str_interp("IGS_${am}"), "IGS.RData")) 
  all_tabestimates <- rbind(all_tabestimates, tabestimates %>% mutate(am = am))
  all_initialvalues <- rbind(all_initialvalues, initialvalues %>% mutate(am = am))
  all_pvals <- rbind(all_pvals, pvals %>% mutate(am = am))
  all_indivParams <- rbind(all_indivParams, indivParams %>% mutate(am = am))
  all_indivFit <- rbind(all_indivFit, indivFit %>% mutate(am = am))
  all_condition_numbers <- rbind(all_condition_numbers, condition_number %>% mutate(am = am))
  all_table_predictions <- rbind(all_table_predictions, table_predictions %>% mutate(am = am))
  all_shrinkage <- rbind(all_shrinkage, shrinkage %>% mutate(am = am))
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

## out of 50 runs, the number of times each model was the best according to BICc
all_tabestimates %>% 
  distinct(run, best) %>% 
  count(best)

# best            n
# <fct>       <int>
#   1 final_model    50


## out of 100 runs, the number of times each model was less than 2 BICc from the best model
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  rename(best_model = best) %>% 
  group_by(run) %>% 
  mutate(best = min(bicc),
         equal_to_best = abs(bicc - best) <= 2) %>% 
  ungroup %>% 
  filter(equal_to_best) %>% 
  count(am, equal_to_best) 

#  am          equal_to_best     n
# <fct>       <lgl>         <int>
#   1 final_model TRUE             50

# ## out of 100 runs, AM3 significantly out-performed only AM6.3 8 times
# all_tabestimates %>% distinct(run, bicc, am, best) %>% 
#   pivot_wider(names_from = am, values_from = bicc) %>% 
#   count(am6.3 > (am3 + 2))
# 

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
# annotate(geom = "text", x = c(1:nrow(ncovar)), y = 470, label = paste((ncovar$n), "covariates"), hjust = 1)

# Look at significance of treatment effects ####

# IN am1, tlr7*Ab on p is never signficiant and tlr7*vac on p is rarely signficiant
# In bigger models, ta effect is more significant than in small models
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  ggplot(aes(am, log10(padj), color = padj >= 0.05)) + 
  ggbeeswarm::geom_beeswarm() + 
  facet_wrap(~parameter, ncol = 3) +
  geom_hline(yintercept = log10(0.05)) +
  scale_color_manual(values = c("black", "gray")) + 
  theme(legend.position = "none")

# TLR7*AB on beta is signifcant ~50% but TLR7*AB on lambda is signficant >80%
# TLR7*VAC on beta is significant ~80% even though it doesnt make sense
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  filter(!parameter %in% c("beta_NPs_intTLR7AB_1", "beta_NPs_intTLR7VAC_1")) %>% 
  group_by(parameter, am) %>% 
  count(significant = padj < 0.05) %>%
  mutate(percent = n/50*100) %>% 
  filter(significant) %>% view


# Look at size of random effects for each model ####
all_tabestimates %>% filter(grepl("omega", parameter)) %>% 
  mutate(parameter = factor(parameter, levels = paste0("omega_", c("lambda", "betas", "ta", "NPs", "p", "m")))) %>% 
  ggplot(aes(am, exp(estimate))) + 
  geom_violin() +
  # geom_boxplot(width = 0.1) +
  ggbeeswarm::geom_beeswarm(aes(color = bicc), alpha = 0.7) + 
  facet_wrap(~parameter, scales = "free_x") + 
  coord_flip() + 
  annotate(geom="point", x = 1, y = 1, shape = NA)  +
  labs(y = "Random effect estimate (fold change)", x = "Model") + 
  color_bicc()


all_condition_numbers %>%
  ggplot(aes(am, condition_number)) +
  geom_violin() + 
  geom_beeswarm(dodge.width = 0.9) + 
  geom_beeswarm(data =. %>% left_join(distinct(all_tabestimates, am, run, bicc)) %>% 
                  group_by(am) %>% filter(bicc == min(bicc)),
                dodge.width = 0.9, shape = 17, size = 4, color = "red") + 
  scale_y_log10() +
  geom_hline(yintercept = 100, linetype = "dashed") +
  coord_flip()

all_tabestimates %>% 
  left_join(all_condition_numbers) %>% 
  group_by(am, run) %>% 
  # Remove am1 becuase it is overfit and am6+ since they have tlr7*vac on ta
  # filter(!am %in% c("am1", "am6", "am6.3", "am6.1")) %>% 
  mutate(n_covar = sum(grepl("beta_", parameter)), 
         lambda_effect = any(grepl("beta_lambda", parameter))) %>% 
  ggplot(aes(condition_number, bicc, color = am, shape = factor(n_covar))) + 
  geom_point() + 
  geom_vline(xintercept = 100, lty = "dashed") + 
  scale_x_log10() + 
  facet_wrap(~lambda_effect)
