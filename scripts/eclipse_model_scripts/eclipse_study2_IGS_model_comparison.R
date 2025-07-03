# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "study2"
out_folder <- here("results_eclipse", study, "IGS_model_comparison")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "derived",str_interp("${toupper(study)}_DATA.txt"))
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", ncol(data) - 4))


## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 50, "p_pop" = 0.7)

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

# AM1 #####
alt_mod <- "am1"
loadProject(here("results_eclipse", study, "SCM_new_fixed", str_interp("${alt_mod}.mlxtran")))
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
# AM1 with vaccine effect on lambda
alt_mod <- "am3"
loadProject(here("results_eclipse", study, "SCM_new_fixed", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_result_am3 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_result_am3)){
  assign(name, parameter_result_am3[[name]])
}

save(list = names(parameter_result_am3),
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

# AM4 #####
# AM1 with tlr7 on Np instead of tlr7*vac
alt_mod <- "am4"
loadProject(here("results_eclipse", study, "SCM_new_fixed", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am4 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am4)){
  assign(name, parameter_results_am4[[name]])
}

save(list = names(parameter_results_am4),
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

# load( here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

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



plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)




# Compare results ####

ams <- c("am1", "am3", "am4")

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

# best      n
# <fct> <int>
#   1 am1       5
# 2 am3      45

# am4 (with tlr7 effect on Np instead of tlr7*vac) is never the best model


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
# <fct> <lgl>         <int>
#   1 am1   TRUE             24
# 2 am3   TRUE             50

# Am1 and am3 give equally good fit 24/50 times


## out of 50 runs, AM1 never significantly out-performed AM3
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  pivot_wider(names_from = am, values_from = bicc) %>% 
  count(am3 > (am1 + 2))

# `am3 > (am1 + 2)`     n
# <lgl>             <int>
#   1 FALSE                50

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


# Look at significance of treatment effects ####

# tlr7*vac on Np is almost always significant, and is significant more often than tlr7 on Np
# tlr7 on Np is never significant at the 0.01 level
# vaccine on lambda is not always significant at the 0.05 or 0.01 level
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  ggplot(aes(am, log10(padj))) + 
  ggbeeswarm::geom_beeswarm(color = "darkgray") + 
  facet_wrap(~parameter, ncol = 3) +
  geom_hline(aes(yintercept = log10(0.05), color = "p-value = 0.05")) +
  geom_hline(aes(yintercept = log10(0.01), color = "p-value = 0.01")) +
  scale_color_bmj(name = "")
  

# vaccine on lambda is significant at the 0.05 level in 74% of cases and significant at the 0.01 level in 44%
# tlr7 on Np is significant at the 0.05 level in 84% of cases and at the 0.01 level in 0%
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  filter(!parameter %in% c("beta_NPs_vac_1", "beta_NPs_intTLR7VAC_1", "beta_ta_vac_1")) %>% 
  group_by(parameter, am) %>% 
  count(significant = padj < 0.05) %>% 
  filter(significant) 


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

# condition numbers tend to be bigger for am3 (with lambda), but there are still some less than 100
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
  # Remove am4 because the treatment effect is not significant
  filter(!am %in% c("am4")) %>% 
  mutate(n_covar = sum(grepl("beta_", parameter)), 
         lambda_effect = any(grepl("beta_lambda", parameter))) %>% 
  ggplot(aes(condition_number, bicc, color = am)) + 
  geom_point() + 
  geom_vline(xintercept = 100, lty = "dashed") + 
  scale_x_log10() 

# Compare BICcs from best models (they are similar):
all_tabestimates %>% 
  distinct(run, am, bicc) %>% 
  group_by(am) %>% filter(bicc == min(bicc))

# run am     bicc
# <dbl> <fct> <dbl>
#   1    28 am1    518.
# 2     2 am3    516.
# 3    14 am4    521.


## Save best model ####

best_run <- all_tabestimates %>% filter(am == "am1") %>% filter(bicc == min(bicc))
igs <- all_initialvalues %>% filter(am == "am1", run == unique(best_run$run))

loadProject(here(out_folder, "IGS_am1.mlxtran"))
saveProject(here("results_eclipse", study, "best_model.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)
runScenario()
saveProject(here("results_eclipse", study, "best_model.mlxtran"))

generate_report(here("results_eclipse", study, "best_model.mlxtran"))
