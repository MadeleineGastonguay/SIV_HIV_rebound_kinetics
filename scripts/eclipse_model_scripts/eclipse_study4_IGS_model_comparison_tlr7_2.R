# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "study4"
out_folder <- here("results_eclipse", study, "IGS_model_comparison_tlr7_2")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", 3), rep("ignore", 3))


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

# AM1 #####
alt_mod <- "am2"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am2)){
  assign(name, parameter_results_am2[[name]])
}

save(list = names(parameter_results_am2),
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
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
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


# AM3.2 #####
alt_mod <- "am3.2"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am3.2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3.2)){
  assign(name, parameter_results_am3.2[[name]])
}

save(list = names(parameter_results_am3.2),
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

# AM4.4 #####
alt_mod <- "am4.4"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am4.4 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am4.4)){
  assign(name, parameter_results_am4.4[[name]])
}

save(list = names(parameter_results_am4.4),
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

# AM4.6 #####
alt_mod <- "am4.6"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am4.6 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am4.6)){
  assign(name, parameter_results_am4.6[[name]])
}

save(list = names(parameter_results_am4.6),
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


# AM3.2 without fixed parameters #####
alt_mod <- "am3.2"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}_no_fixed.mlxtran")
saveProject(here(out_folder, project_name))
temp <- getPopulationParameterInformation()
temp$method <- "MLE"
setPopulationParameterInformation(temp)
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

predefined_initialvalues_no_fixed <- predefined_initialvalues %>% mutate(method = "MLE")

# Run convergence for narrower range of values
parameter_results_am3.2_no_fixed <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = NULL, 
  predefined_initialvalues = predefined_initialvalues_no_fixed, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3.2_no_fixed)){
  assign(name, parameter_results_am3.2_no_fixed[[name]])
}

save(list = names(parameter_results_am3.2_no_fixed),
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


# AM4.4 without fixed parameters #####
alt_mod <- "am4.4"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}_no_fixed.mlxtran")
saveProject(here(out_folder, project_name))
temp <- getPopulationParameterInformation()
temp$method <- "MLE"
setPopulationParameterInformation(temp)
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

predefined_initialvalues_no_fixed <- predefined_initialvalues %>% mutate(method = "MLE")

# Run convergence for narrower range of values
parameter_results_am4.4_no_fixed <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = NULL, 
  predefined_initialvalues = predefined_initialvalues_no_fixed, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am4.4_no_fixed)){
  assign(name, parameter_results_am4.4_no_fixed[[name]])
}

save(list = names(parameter_results_am4.4_no_fixed),
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



# AM4.6 without fixed parameters #####
alt_mod <- "am4.6"
loadProject(here("results_eclipse", study, "SCM_tlr7_2", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}_no_fixed.mlxtran")
saveProject(here(out_folder, project_name))
temp <- getPopulationParameterInformation()
temp$method <- "MLE"
setPopulationParameterInformation(temp)
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

predefined_initialvalues_no_fixed <- predefined_initialvalues %>% mutate(method = "MLE")

# Run convergence for narrower range of values
parameter_results_am4.6_no_fixed <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = NULL, 
  predefined_initialvalues = predefined_initialvalues_no_fixed, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am4.6_no_fixed)){
  assign(name, parameter_results_am4.6_no_fixed[[name]])
}

save(list = names(parameter_results_am4.6_no_fixed),
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

# ams <- c("am2", "am3.1", "am3.2", "am4.4", "am4.6")
ams <- c( "am3.2_no_fixed", "am4.4_no_fixed", "am4.6_no_fixed")

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

## out of 100 runs, the number of times each model was the best according to BICc
all_tabestimates %>% 
  distinct(run, best) %>% 
  count(best)

# best      n
#   1 am2       6
# 2 am3.1    35
# 3 am3.2     9


## out of 50 runs, the number of times each model was less than 2 BICc from the best model
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  rename(best_model = best) %>% 
  group_by(run) %>% 
  mutate(best = min(bicc),
         equal_to_best = abs(bicc - best) <= 2) %>% 
  ungroup %>% 
  filter(equal_to_best) %>% 
  count(am, equal_to_best) 

# am    equal_to_best     n
#   1 am2   TRUE             15
# 2 am3.1 TRUE             38
# 3 am3.2 TRUE             13

## out of 100 runs, AM3.2 significantly out-performed AM3.1 12 times
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  pivot_wider(names_from = am, values_from = bicc) %>% 
  count(am3.1 > (am3.2 + 2))

# `am3.1 > (am3.2 + 2)`     n
#   1 FALSE                    38
# 2 TRUE                     12

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

# Look at significance of treatment effects ###
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  ggplot(aes(am, log10(padj), color = padj >= 0.05)) + 
  ggbeeswarm::geom_beeswarm() + 
  facet_wrap(~parameter, ncol = 3) +
  geom_hline(yintercept = log10(0.05)) +
  scale_color_manual(values = c("black", "gray")) + 
  theme(legend.position = "none")

# Ab on m is significant in only 24% of runs for am 3.1
# TLR7 on lambda is significant in only 26% of runs for am4.4 vs 84% for TLR7 on beta in am3.2
# ab on beta is signficaint in only 28% of runs for am4.4 vs 92% for ab on lambda in am3.2
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  group_by(parameter, am) %>% 
  count(significant = padj < 0.05) %>% 
  filter(significant) %>% 
  mutate(p = n/50*100) %>% view


# Look at size of random effects for each model ####
# 6.1 and 5 have bigger RE on m (they dont have an m or lambda or p effet)

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
  mutate(n_covar = sum(grepl("beta_", parameter)), 
         lambda_effect = any(grepl("beta_lambda", parameter))) %>% 
  ggplot(aes(condition_number, bicc, color = am, shape = factor(n_covar))) + 
  geom_point() + 
  geom_vline(xintercept = 100, lty = "dashed") + 
  scale_x_log10() + 
  facet_wrap(~lambda_effect)



## Compare Simulations for AM3.2 and AM4.4 (lamdba and beta effects swapped)

# AM3.2 can better capture peak viremia for 7160, 7220
# For several fully censored animals, AM3.2 assumes immune control but AM4.4 assumes delayed rebound followed by immune control

all_table_predictions %>% 
  left_join(all_tabestimates_bicc) %>% 
  filter(am %in% c("am3.2", "am4.4", "am4.6")) %>%
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  # pivot_longer(contains("mode"), names_to = "variable") %>%
  ggplot(aes(time, log10(V_mode), color = am, group = interaction(run, am))) + 
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line() + 
  facet_wrap(~id) 

# AM 3.2 predicts higher set point viremia in vac only and lower in ab only and ab + vac
# AM3.2 predicts faster exponential growht in vac only and lsower in control
# Am3.2 predicts higher peak viremia in vac only and lower in ab only
all_indivFit %>% 
  filter(am %in% c("am3.2", "am4.4", "am4.6")) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  left_join(all_tabestimates_bicc) %>% 
  # pivot_longer(contains("mode"), names_to = "variable") %>%
  ggplot(aes(time, popPred, color = am, group = interaction(run, am))) + 
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  geom_line(data = ~group_by(.x, am) %>% filter(bicc == min(bicc)), color = "black") + 
  facet_wrap(ab ~ vac, labeller = "label_both")

all_indivFit %>% 
  filter(am %in% c("am3.2", "am4.4", "am4.6")) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  left_join(all_tabestimates_bicc) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  # pivot_longer(contains("mode"), names_to = "variable") %>%
  ggplot(aes(time, popPred, color = am, group = interaction(run, am))) + 
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(linewidth = 1) + 
  # geom_line(data = ~group_by(.x, am) %>% filter(bicc == min(bicc)), color = "black") + 
  facet_wrap(ab ~ vac, labeller = "label_both")

# Compare group level parameters
group_paramemters <- left_join(
  all_tabestimates %>% 
    filter(am %in% c("am3.2", "am4.4", "am4.6"), grepl("pop", parameter)) %>% 
    left_join(all_condition_numbers) %>% 
    filter(condition_number < 100) %>% 
    group_by(am) %>% filter(bicc == min(bicc)) %>% 
    select(parameter, estimate, am) %>% 
    rbind(data.frame(am = rep(c("am3.2", "am4.4", "am4.6"), each = 2),
                     parameter = rep(c("lambda_pop", "p_pop"), 3),
                     estimate = unname(rep(fixed, 3)))) %>% 
    separate(parameter, sep = "_", into = c("parameter", "t")) %>% 
    select(-t) %>% rename(pop = estimate),
  
  all_tabestimates %>% 
    filter(am %in% c("am3.2", "am4.4", "am4.6"), grepl("beta_", parameter)) %>% 
    left_join(all_condition_numbers) %>% 
    filter(condition_number < 100) %>% 
    group_by(am) %>% filter(bicc == min(bicc)) %>% 
    select(parameter, estimate, am) %>% 
    separate(parameter, sep = "_", into = c("t", "parameter", "treatment", "t2")) %>% 
    select(-t, -t2) %>% 
    pivot_wider(names_from = treatment, values_from = estimate, values_fill = 0) 
) %>% 
  select(am, parameter, pop, tlr7, ab, vac) %>% 
  replace_na(list(ab=0,tlr7=0,vac = 0)) %>% 
  group_by(am, parameter) %>% 
  summarise(Control = pop,
            TLR7Ab = pop*exp(tlr7+ab),
            TLR7Vac = pop*exp(tlr7+vac),
            TLR7AbVac = pop*exp(tlr7+vac+ab)
  ) 

sds <- all_tabestimates %>% 
  filter(am %in% c("am3.2", "am4.4", "am4.6"), grepl("omega", parameter)) %>% 
  left_join(all_condition_numbers) %>% 
  filter(condition_number < 100) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  select(parameter, estimate, am) %>% 
  separate(parameter, sep = "_", into = c("t", "parameter")) %>% 
  select(-t) %>% rename(sd = estimate)

group_paramemters %>% 
  pivot_longer(!c(am, parameter), names_to = "Treatment") %>% 
  ggplot(aes(Treatment, value, fill = am)) + 
  geom_col(position = position_dodge()) + 
  geom_col(data = ~group_by(.x, am, Treatment) %>% filter(parameter %in% c("betas", "lambda")) %>% 
             summarise(parameter = "lambda*beta", value = prod(value)), position = position_dodge()) +
  facet_wrap(~parameter, scales = "free_y") 


pop_params <- group_paramemters %>% 
  pivot_longer(!c(am, parameter), names_to = "Treatment") %>%  
  left_join(sds)  %>% 
  rename(mean = value)

density_data <- pop_params %>%
  group_by(am, parameter, Treatment) %>%
  group_modify(~ {
    x_vals <- exp(seq((min(log(.x$mean) - 4*.x$sd)),
                      (max(log(.x$mean) + 4*.x$sd)),
                      length.out = 500))
    do.call(rbind, lapply(1:nrow(.x), function(i) {
      data.frame(
        # am = .x$am,
        # Treatment = .x$Treatment,
        x = x_vals,
        y = dnorm(log(x_vals), mean = log(.x$mean[i]), sd = .x$sd[i])
      )
    }))
  }) %>%
  ungroup()

distribution_plot <- density_data %>%
  mutate(x = case_when(#parameter == "betas" ~ x*1e-7, 
    parameter == "NPs" ~ x*1e5,
    .default = x),
    parameter = factor(parameter, 
                       levels =  c("lambda","betas", "ta", "m", "p", "NPs"),
                       labels =  c("lambda","beta~'(x'~10^-7~')'", "t['a']", "m", "p", "N['p']")),
    Treatment = factor(Treatment, levels = c("Control", "TLR7Ab", "TLR7Vac", "TLR7AbVac"),
                       labels = c("Control", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Ab + Vac"))) %>%
  ungroup %>% 
  group_by(parameter) %>% mutate(ynorm = y/max(y)) %>% 
  ggplot(aes(x = log10(x), y = ynorm, color = am)) +
  geom_line(linewidth = 1) +
  facet_grid(Treatment~parameter, scales = "free_x", labeller = label_parsed, switch = "x") +
  # labs(title = "Normal Distributions by Parameter",
  #      x = "", y = "Density") +
  theme_minimal() + 
  # scale_x_log10() +
  scale_x_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) +
  # facetted_pos_scales(
  #   x = list(
  #     "lambda" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(10, 30, 100, 300))),
  #     "beta" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(5, 10, 30, 50)))
  #   )
  # ) +
  # facetted_pos_scales(
  #   x = list(
  #     "lambda" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(10, 30, 100, 300)), limits = log10(c(10,350))),
  #     "beta~'(x'~10^-7~')'" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(5, 10, 30, 50))),
  #     "t['a']" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,2)),
  #     "m" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,4,by=2), limits = c(-2.5,4.5)),
  #     "p" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-1,1))
  #   )
  # ) +
  theme(axis.text.y = element_blank(), axis.title = element_blank(),
        axis.ticks.y = element_blank(),# panel.grid = element_blank(), 
        axis.line.x = element_line(),
        axis.text.x = element_text(size = 15),
        strip.placement = "outside", strip.text = element_text(size = 20))



## Save best model ####

## Select best model as best BICc for am3.2 with condition # < 100

all_tabestimates_bicc %>% filter(am == "am3.2") %>% 
  left_join(all_condition_numbers) %>% 
  ggplot(aes(condition_number, bicc)) + geom_text(aes(label = run)) + 
  scale_x_log10()

best_run <- all_tabestimates_bicc %>% filter(am == "am3.2") %>% 
  left_join(all_condition_numbers) %>% 
  filter(condition_number < 100) %>% 
  filter(bicc == min(bicc)) %>% 
  pull(run) 

igs <- all_initialvalues %>% filter(am == "am3.2", run == best_run)

loadProject(here(out_folder, "IGS_am3.2.mlxtran"))
saveProject(here("results_eclipse", study, "best_model.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)
runScenario()
saveProject(here("results_eclipse", study, "best_model.mlxtran"))
generate_report(here("results_eclipse", study, "best_model.mlxtran"))


# Save best model for AM4.4 ####

all_tabestimates_bicc %>% filter(am == "am4.4") %>% 
  left_join(all_condition_numbers) %>% 
  ggplot(aes(condition_number, bicc)) + geom_text(aes(label = run)) + 
  scale_x_log10()

best_run <- all_tabestimates_bicc %>% filter(am == "am4.4") %>% 
  left_join(all_condition_numbers) %>% 
  filter(condition_number < 100) %>% 
  filter(bicc == min(bicc)) %>% 
  pull(run) 

igs <- all_initialvalues %>% filter(am == "am4.4", run == best_run)

loadProject(here(out_folder, "IGS_am4.4.mlxtran"))
saveProject(here("results_eclipse", study, "best_4.4.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)
runScenario()
saveProject(here("results_eclipse", study, "best_4.4.mlxtran"))
generate_report(here("results_eclipse", study, "best_4.4.mlxtran"))


# Save best model for AM4.6 ####

all_tabestimates_bicc %>% filter(am == "am4.6") %>% 
  left_join(all_condition_numbers) %>% 
  ggplot(aes(condition_number, bicc)) + geom_text(aes(label = run)) + 
  scale_x_log10()

best_run <- all_tabestimates_bicc %>% filter(am == "am4.6") %>% 
  left_join(all_condition_numbers) %>% 
  filter(condition_number < 100) %>% 
  filter(bicc == min(bicc)) %>% 
  pull(run) 

igs <- all_initialvalues %>% filter(am == "am4.6", run == best_run)

loadProject(here(out_folder, "IGS_am4.6.mlxtran"))
saveProject(here("results_eclipse", study, "best_4.6.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)
runScenario()
saveProject(here("results_eclipse", study, "best_4.6.mlxtran"))
generate_report(here("results_eclipse", study, "best_4.6.mlxtran"))
