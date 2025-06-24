# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "study1"
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
fixed <- c("lambda_pop" = 100, "p_pop" = 1)

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

initialvalues <- NULL
nits <- 50
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


if(!file.exists(here(out_folder, "initial_values.csv"))) write_csv(initialvalues, here(out_folder, "initial_values.csv"))
if(!file.exists(here(out_folder, "initial_values2.csv"))) write_csv(initialvalues, here(out_folder, "initial_values2.csv"))

predefined_initialvalues <- read_csv(here(out_folder, "initial_values.csv"))
predefined_initialvalues2 <- read_csv(here(out_folder, "initial_values2.csv"))

rbind(predefined_initialvalues, predefined_initialvalues2) %>% 
  ggplot(aes(parameter, initialValue)) +
  geom_boxplot() + 
  facet_wrap(~parameter, scales = "free")

# Run IGS #####

# Final SCM model #####
# Final model from SCM:
# TLR7 effect on ta, lambda, and p
alt_mod <- "final_model"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_result_final_model <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = TRUE, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE, 
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_result_final_model)){
  assign(name, parameter_result_final_model[[name]])
}

save(list = names(parameter_result_final_model),
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


# There are two regimes for this model
plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

low_m_regime <- tabestimates %>% filter(parameter == "m_pop", estimate < 10) %>% pull(run)

plot_vs_rses(tabestimates, out_folder, project_name)

correlate_estimates(
  tabestimates %>% filter(run %in% low_m_regime), 
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_low_m"
)

correlate_estimates(
  tabestimates %>% filter(run %in% low_m_regime, grepl("pop", parameter))
)

compare_estimates(
  tabestimates %>% filter(run %in% low_m_regime), 
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_low_m"
)

compare_estimates(
  tabestimates %>% filter(!run %in% low_m_regime), 
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_big_m"
)

pdf(here(out_folder, gsub(".mlxtran", "", project_name), "table_predictions.pdf"), width = 10, height = 10)

table_predictions %>% 
  ggplot(aes(time, V_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_point(data = data, aes(time, 10^vl), inherit.aes = F, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, P_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, E_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, I_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, T_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

dev.off()


# AM1 #####
# AM1 from SCM 
# TLR7 effect on ta and lambda
alt_mod <- "am1"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_result_am1 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)

# Run convergence for narrower range of values
parameter_result_am1.2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues2, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)

for(name in names(parameter_result_am1)){
  assign(name, parameter_result_am1[[name]])
}

save(list = names(parameter_result_am1),
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

# load(here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

max_run <- max(tabestimates$run)

# combine elements of parameter_results_am3.2 with parameter_results_am3.2 and update the run number in results 2 to add 50 
tabestimates <- rbind(tabestimates, parameter_result_am1.2$tabestimates %>% mutate(run = run + max_run))
tabiters <- rbind(tabiters, parameter_result_am1.2$tabiters %>% mutate(run = run + max_run))
phase <- rbind(phase, parameter_result_am1.2$phase)
initialvalues <- rbind(initialvalues, parameter_result_am1.2$initialvalues %>% mutate(run = run + max_run))
shrinkage <- rbind(shrinkage, parameter_result_am1.2$shrinkage %>% mutate(run = run + max_run))
indivParams <- rbind(indivParams, parameter_result_am1.2$indivParams %>% mutate(run = run + max_run))
indivFit <- rbind(indivFit, parameter_result_am1.2$indivFit %>% mutate(run = run + max_run))
table_predictions <- rbind(table_predictions, parameter_result_am1.2$table_predictions %>% mutate(run = run + max_run))
pvals <- rbind(pvals, parameter_result_am1.2$pvals %>% mutate(run = run + max_run))
condition_number <- rbind(condition_number, parameter_result_am1.2$condition_number %>% mutate(run = run + max_run))

save(list = names(parameter_result_am1.2),
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS_100its.RData"))



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


# There is clearly two parameter regimes: high m / low ta and low ta / high m
# The low ta regime has less measurement noise and smaller random effects for m
plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

low_m_regime <- tabestimates %>% filter(parameter == "m_pop", estimate < 20) %>% pull(run)

plot_vs_rses(tabestimates %>% filter(run %in% low_m_regime), out_folder, project_name, 
             file_addendum = "_low_m")

correlate_estimates(
  tabestimates %>% filter(run %in% low_m_regime), 
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_low_m"
)

correlate_estimates(
  tabestimates %>% filter(run %in% low_m_regime, grepl("pop", parameter))
)

compare_estimates(
  tabestimates %>% filter(run %in% low_m_regime),
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_low_m"
)

compare_estimates(
  tabestimates %>% filter(!run %in% low_m_regime),
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_big_m"
)

estimate_boxplots(tabestimates %>% filter(run %in% low_m_regime), 
                  out_folder = out_folder, project_name = project_name,
                  file_addendum = "_low_m")

pdf(here(out_folder, gsub(".mlxtran", "", project_name), "parameter_plots_low_m.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  if(parameter %in% names(fixed)) next
  print(parameter)
  print(
    parameter_plots(
      parameter,
      remove_bicc_outliers = F, 
      remove_estimate_outliers = F,
      NULL, tabestimates %>% filter(run %in% low_m_regime), tabiters %>% filter(run %in% low_m_regime),
      phase[low_m_regime,], initialvalues %>% filter(run %in% low_m_regime)
    )
  )
}
dev.off()

pdf(here(out_folder, gsub(".mlxtran", "", project_name), "table_predictions.pdf"), width = 10, height = 10)

table_predictions %>% 
  ggplot(aes(time, V_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_point(data = data, aes(time, 10^vl), inherit.aes = F, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, P_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, E_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, I_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, T_mode, color = run %in% low_m_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

dev.off()


# AM3 #####
# AM3 from SCM 
# TLR7 effect on ta, lambda, and m
alt_mod <- "am3"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am3 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)

parameter_results_am3.2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues2, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_am3)){
  assign(name, parameter_results_am3[[name]])
}


save(list = names(parameter_results_am3),
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

# load( here(out_folder, gsub(".mlxtran", "", project_name), "IGS.RData"))

max_run <- max(tabestimates$run)

# combine elements of parameter_results_am3.2 with parameter_results_am3.2 and update the run number in results 2 to add 50 
tabestimates <- rbind(tabestimates, parameter_results_am3.2$tabestimates %>% mutate(run = run + max_run))
tabiters <- rbind(tabiters, parameter_results_am3.2$tabiters %>% mutate(run = run + max_run))
phase <- rbind(phase, parameter_results_am3.2$phase)
initialvalues <- rbind(initialvalues, parameter_results_am3.2$initialvalues %>% mutate(run = run + max_run))
shrinkage <- rbind(shrinkage, parameter_results_am3.2$shrinkage %>% mutate(run = run + max_run))
indivParams <- rbind(indivParams, parameter_results_am3.2$indivParams %>% mutate(run = run + max_run))
indivFit <- rbind(indivFit, parameter_results_am3.2$indivFit %>% mutate(run = run + max_run))
table_predictions <- rbind(table_predictions, parameter_results_am3.2$table_predictions %>% mutate(run = run + max_run))
pvals <- rbind(pvals, parameter_results_am3.2$pvals %>% mutate(run = run + max_run))
condition_number <- rbind(condition_number, parameter_results_am3.2$condition_number %>% mutate(run = run + max_run))

save(list = names(parameter_results_am3.2),
     file = here(out_folder, gsub(".mlxtran", "", project_name), "IGS_100its.RData"))

load(here(out_folder, gsub(".mlxtran", "", project_name), "IGS_100its.RData"))


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


# two regimes
plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

low_ta_regime <- tabestimates %>% filter(parameter == "ta_pop", estimate <= 0.4) %>% pull(run)

correlate_estimates(
  tabestimates %>% filter(run %in% low_ta_regime), 
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_low_ta")

# No lingering correlations
correlate_estimates(tabestimates %>% filter(run %in% low_ta_regime, grepl("pop", parameter)))

compare_estimates(tabestimates %>% filter(run %in% low_ta_regime), 
                  out_folder = out_folder, project_name = project_name,
                  file_addendum = "_low_ta")

estimate_boxplots(tabestimates %>% filter(run %in% low_ta_regime), 
                  out_folder = out_folder, project_name = project_name,
                  file_addendum = "_low_ta")


pdf(here(out_folder, gsub(".mlxtran", "", project_name), "table_predictions.pdf"), width = 10, height = 10)
table_predictions %>% 
  ggplot(aes(time, V_mode, color = run %in% low_ta_regime, group = factor(run)))  + 
  geom_point(data = data, aes(time, 10^vl), inherit.aes = F, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, P_mode, color = run %in% low_ta_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, E_mode, color = run %in% low_ta_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, I_mode, color = run %in% low_ta_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

table_predictions %>% 
  ggplot(aes(time, T_mode, color = run %in% low_ta_regime, group = factor(run)))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  labs(color = "low m, low ta regime")

dev.off()



# Compare results ####

ams <- c("final_model","am1", "am3")

all_tabestimates <- NULL
all_initialvalues <- NULL
all_pvals <- NULL
all_indivParams <- NULL
all_indivFit <- NULL
all_condition_numbers <- NULL
all_table_predictions <- NULL
all_shrinkage <- NULL
for(am in ams){
  file <- here(out_folder, str_interp("IGS_${am}"), "IGS_100its.RData")
  if(file.exists(file)){
    load(file)
  }else{
    load(here(out_folder, str_interp("IGS_${am}"), "IGS.RData")) 
  }
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
  ggplot(aes(initialValue, estimate, color = am, shape = am == best)) + 
  geom_point(alpha = 0.75, size = 2) + 
  facet_wrap(~parameter, scales = "free") + 
  labs(x = "Initial guess", y = "Estimate", color = "Model", shape = "Best Fit") 

# Compare model performances ####

## out of 100 runs, the number of times each model was the best according to BICc
all_tabestimates %>% 
  distinct(run, best) %>% 
  count(best)

# best            n
# final_model    37
# am1             3
# am3            10

## filter to low ta, low m, lower measurement error regime:
regime1 <- all_tabestimates %>% filter(parameter == "a", estimate < 0.32) %>%
  mutate(t = paste(am, run, sep = ".")) %>% pull(t)

all_tabestimates %>% 
  ggplot(aes(estimate, bicc, color = am, shape = paste(am, run, sep = ".") %in% regime1)) + 
  geom_point() + facet_wrap(~parameter, scales = "free_x")



## out of 100 runs, the number of times each model was less than 2 BICc from the best model
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  filter(paste(am, run, sep = ".") %in% regime1) %>% 
  rename(best_model = best) %>% 
  group_by(run) %>% 
  mutate(best = min(bicc),
         equal_to_best = abs(bicc - best) <= 2) %>% 
  ungroup %>% 
  filter(equal_to_best) %>% 
  count(am, equal_to_best) 

# am          equal_to_best     n
# 1 final_model TRUE              9
# 2 am1         TRUE             11
# 3 am3         TRUE             25


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
  coord_flip() + 
  annotate(geom = "text", x = c(1:nrow(ncovar)), y = 440, label = paste((ncovar$n), "covariates"), hjust = 1) + 
  facet_wrap(~ (paste(am, run, sep = ".") %in% regime1), scales = "free_x")


# How frequently is the tlr7 effect on m and p signficant?
all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr"),
         regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  ggplot(aes(am, log10(padj), color = regime1)) + 
  ggbeeswarm::geom_beeswarm(dodge.width = 0.5) + 
  facet_wrap(~parameter, ncol = 3) +
  geom_hline(yintercept = log10(0.05)) + 
  labs(color = "Low m, low ta regime")
# scale_color_manual(values = c("black", "gray")) + 
# theme(legend.position = "none")

# tlr7 effect on p in final SCM model is never significant in the desired regime
# tlr7 effect on m in am 3 is signficant less than 50% of the time in the desired regime
all_pvals %>% 
  filter(am != "final_model") %>% 
  mutate(padj = p.adjust(p.value, method = "fdr"),
         regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  left_join(distinct(., am, run, regime1) %>% add_count(am, regime1, name = "n_r1")) %>% 
  filter(padj < 0.05) %>% 
  count(regime1,parameter, am , n_r1)  %>% 
  mutate(p = n/n_r1*100)

# Look at size of random effects for each model ####
# There's no model with clearly large random effects

all_tabestimates %>% filter(grepl("omega", parameter)) %>% 
  mutate(parameter = factor(parameter, levels = paste0("omega_", c("lambda", "betas", "ta", "NPs", "p", "m")))) %>% 
  ggplot(aes(am, exp(estimate))) + 
  geom_violin() +
  # geom_boxplot(width = 0.1) +
  ggbeeswarm::geom_beeswarm(aes(color = bicc), alpha = 0.7) + 
  facet_wrap(paste(am, run, sep = ".") %in% regime1 ~parameter,  nrow = 2) + 
  coord_flip() + 
  annotate(geom="point", x = 1, y = 1, shape = NA)  +
  labs(y = "Random effect estimate (fold change)", x = "Model") + 
  color_bicc() +
  scale_y_log10()

all_tabestimates %>% filter(grepl("omega", parameter)) %>% 
  mutate(parameter = factor(parameter, levels = paste0("omega_", c("lambda", "betas", "ta", "NPs", "p", "m")))) %>% 
  ggplot(aes(am, (estimate))) + 
  geom_violin() +
  # geom_boxplot(width = 0.1) +
  ggbeeswarm::geom_beeswarm(aes(color = bicc), alpha = 0.7) + 
  facet_wrap(paste(am, run, sep = ".") %in% regime1 ~parameter,  nrow = 2) + 
  coord_flip() + 
  annotate(geom="point", x = 1, y = 1, shape = NA)  +
  labs(y = "Random effect estimate", x = "Model") + 
  color_bicc()


all_condition_numbers %>%
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  ggplot(aes(am, condition_number, color = regime1)) +
  geom_violin() + 
  geom_beeswarm(dodge.width = 0.9) + 
  geom_beeswarm(data =. %>% left_join(distinct(all_tabestimates, am, run, bicc)) %>% 
                  group_by(am, regime1) %>% filter(bicc == min(bicc)),
                dodge.width = 0.9, shape = 17, size = 4) + 
  scale_y_log10() +
  geom_hline(yintercept = 100, linetype = "dashed") +
  coord_flip()

all_indivFit %>%
  left_join(distinct(data, id, tlr7)) %>% 
  filter(am != "final_model") %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  ggplot(aes(time, popPred, color = am, group = interaction(am, run))) +
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  facet_grid(regime1~tlr7, labeller = "label_both")

all_indivFit %>%
  left_join(distinct(data, id, tlr7)) %>% 
  left_join(distinct(all_tabestimates, am, bicc, run)) %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter(am != "final_model", regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  ggplot(aes(time, indivPred_mode, color = am, group = interaction(am, run))) +
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.5, linewidth = 1) + 
  facet_wrap(~id)


all_indivParams %>%
  left_join(distinct(all_tabestimates, am, bicc, run)) %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter(am != "final_model", regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  select(-run, -bicc, -regime1) %>% 
  pivot_longer(!c(id, am, tlr7), names_to = "parameter") %>% 
  pivot_wider(names_from = am) %>% 
  ggplot(aes(am1, am3, color = factor(tlr7))) + 
  geom_point() + 
  geom_abline() + 
  facet_wrap(~parameter, scales = "free")  + 
  scale_x_log10() +
  scale_y_log10()

pop_param <- all_tabestimates %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter(regime1, am != "final_model") %>% 
  group_by(am) %>% filter(bicc == min(bicc)) 

pop_param <- rbind(
  pop_param %>% 
    filter(grepl("pop", parameter)) %>% 
    mutate(parameter = gsub("_pop", "", parameter), tlr7 = 0) %>% 
    select(parameter, value = estimate, tlr7) ,
  
  pop_param %>% 
    filter(!grepl("omega_", parameter), parameter != "a") %>% 
    rowwise() %>% 
    mutate(p = ifelse(grepl("pop", parameter), gsub("_pop", "", parameter),
                      unlist(str_split(parameter, "_"))[2]), 
           estimate = ifelse(grepl("pop", parameter), estimate, exp(estimate))) %>% 
    select(p, parameter, estimate) %>% 
    group_by(am, p) %>% summarise(value = prod(estimate, na.rm = TRUE)) %>%
    mutate(value = ifelse(p == "lambda", value*100, value), tlr7 = 1) %>% 
    rename(parameter = p)
)

## There is some evidence for a tlr7 effect on m, but is it strong enough to include in the model???
all_indivParams %>%
  left_join(distinct(all_tabestimates, am, bicc, run)) %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter(am != "final_model", regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  select(-run, -bicc, -regime1) %>% 
  pivot_longer(!c(id, am, tlr7), names_to = "parameter")  %>% 
  ggplot(aes(factor(tlr7), value, color = am)) + 
  geom_violin(scale = "width") + geom_beeswarm(dodge.width = 0.9) + 
  facet_wrap(~parameter, scales = "free_y") + 
  scale_y_log10() +
  geom_point(data = pop_param, shape = 17, size = 4, position = position_dodge(width = 0.9)) 


# Ridiculously large RSE for m in am3
all_tabestimates %>% 
  filter(am != "final_model", paste(am, run, sep = ".") %in% regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) 

# Best BICcs are not noticeably different
all_tabestimates %>% 
  filter(am != "final_model", paste(am, run, sep = ".") %in% regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc))  %>% distinct(am, bicc)
