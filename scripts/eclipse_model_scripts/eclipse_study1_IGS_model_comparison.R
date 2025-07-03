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
fixed <- c("lambda_pop" = 200, "p_pop" = 3)

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
nits <- 100

if(!file.exists(here(out_folder, "initial_values.csv"))){
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

# Final SCM model #####
# Final model from SCM:
# TLR7 effect on ta, beta, and p
alt_mod <- "final_model"
loadProject(here("results_eclipse", study, "SCM_new_fixed", str_interp("${alt_mod}.mlxtran")))
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

low_m_regime <- tabestimates %>% filter(parameter == "m_pop", estimate < 100) %>% pull(run)
# 40 estimates in this regime
length(low_m_regime)

plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  color = low_m_regime
)

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

estimate_boxplots(tabestimates %>% filter(run %in% low_m_regime),
                  out_folder = out_folder, project_name = project_name,
                  file_addendum = "_low_m")

estimate_boxplots(tabestimates %>% filter(!run %in% low_m_regime),
                  out_folder = out_folder, project_name = project_name,
                  file_addendum = "_big_m")

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


# AM2 #####
# AM2 from SCM 
# TLR7 effect on ta and lambda
alt_mod <- "am2"
loadProject(here("results_eclipse", study, "SCM_new_fixed", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_am2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
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


# There is clearly two parameter regimes: high m / low ta and low ta / high m
# The low ta regime has less measurement noise and smaller random effects for m
plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

low_m_regime <- tabestimates %>% filter(parameter == "m_pop", estimate < 200) %>% pull(run)

plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  color = low_m_regime
)

plot_vs_rses(tabestimates %>% filter(run %in% low_m_regime), out_folder, project_name, 
             file_addendum = "_low_m")

correlate_estimates(
  tabestimates %>% filter(run %in% low_m_regime), 
  out_folder = out_folder, project_name = project_name,
  file_addendum = "_low_m"
)


# No lingering correlations in the low m regime
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


estimate_boxplots(tabestimates %>% filter(!run %in% low_m_regime), 
                  out_folder = out_folder, project_name = project_name,
                  file_addendum = "_big_m")

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



# Compare results ####

ams <- c("final_model","am2")

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
# <fct>       <int>
#   1 final_model    52
# 2 am2            48

# The final model was the best more frequently because it landed in the immune control regime more often

## filter to low ta, low m, lower measurement error regime:
regime1 <- all_tabestimates %>% filter(parameter == "omega_ta", estimate > 1) %>%
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
# <fct>       <lgl>         <int>
#   1 final_model TRUE             17
# 2 am2         TRUE             85


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


# How frequently is the tlr7 effect on p signficant?
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
all_pvals %>% 
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
  # filter(am != "final_model") %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  ggplot(aes(time, popPred, color = am, group = interaction(am, run))) +
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.5) + 
  facet_grid(regime1~tlr7, labeller = "label_both")

all_indivFit %>%
  left_join(distinct(data, id, tlr7)) %>% 
  left_join(distinct(all_tabestimates, am, bicc, run)) %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter( regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  ggplot(aes(time, indivPred_mode, color = am, group = interaction(am, run))) +
  geom_point(data = data, inherit.aes = F, aes(time, vl), alpha = 0.5) +
  geom_line(alpha = 0.5, linewidth = 1) + 
  facet_wrap(~id)


all_indivParams %>%
  left_join(distinct(all_tabestimates, am, bicc, run)) %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter( regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  select(-run, -bicc, -regime1) %>% 
  pivot_longer(!c(id, am, tlr7), names_to = "parameter") %>% 
  pivot_wider(names_from = am) %>% 
  ggplot(aes(am2, final_model, color = factor(tlr7))) + 
  geom_point() + 
  geom_abline() + 
  facet_wrap(~parameter, scales = "free")  + 
  scale_x_log10() +
  scale_y_log10()

pop_param <- all_tabestimates %>% 
  mutate(regime1 = paste(am, run, sep = ".") %in% regime1) %>% 
  filter(regime1) %>% 
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
  filter( regime1) %>% 
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
  filter( paste(am, run, sep = ".") %in% regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) 

# Best BICcs are not noticeably different
all_tabestimates %>% 
  filter( paste(am, run, sep = ".") %in% regime1) %>% 
  group_by(am) %>% filter(bicc == min(bicc))  %>% distinct(am, bicc)


## Save best model ####

# Run with minimum BICc has massive RSE for m and ta... find the next best
best_run <- all_tabestimates %>% 
  filter(am == "am2", paste(am, run, sep = ".") %in% regime1) %>% 
  filter(bicc == min(bicc))

best_run <- all_tabestimates %>% 
  filter(am == "am2", paste(am, run, sep = ".") %in% regime1) %>% 
  filter(bicc != min(bicc)) %>% 
  filter(bicc == min(bicc))

# P-vlaue for tlr7 effect on lambda is not signficant... 
all_pvals %>% filter(am == "am2", run == unique(best_run$run))
all_condition_numbers %>% filter(am == "am2", run == unique(best_run$run))

temp_pval <- all_pvals %>% 
  filter(am == "am2") %>% 
  select(-statistics) %>% 
  pivot_wider(names_from = parameter, values_from = p.value) %>% 
  mutate(signif = case_when(beta_lambda_tlr7_1 < 0.05 & beta_ta_tlr7_1 < 0.05 ~ "both",
                            beta_lambda_tlr7_1 < 0.05 & beta_ta_tlr7_1 > 0.05 ~ "lambda",
                            beta_lambda_tlr7_1 > 0.05 & beta_ta_tlr7_1 < 0.05 ~ "ta",
                            .default = "neither"))

all_tabestimates %>% 
  filter(am == "am2", paste(am, run, sep = ".") %in% regime1) %>% 
  distinct(am, run, bicc) %>% 
  left_join(all_condition_numbers) %>% 
  left_join(temp_pval) %>% 
  ggplot(aes(signif, bicc, color = condition_number <= 100)) + geom_beeswarm()

all_tabestimates %>% 
  filter(am == "am2", paste(am, run, sep = ".") %in% regime1) %>% 
  left_join(all_condition_numbers) %>% 
  left_join(temp_pval) %>% 
  # filter(grepl("pop", parameter) | grepl("beta_", parameter)) %>% 
  group_by(signif) %>% filter(bicc == min(bicc)) %>% 
  ggplot(aes(signif, estimate)) + geom_col() + 
  geom_errorbar(data = . %>% filter(rses < 1e3),
                aes(ymin = estimate*(1-rses/100), ymax = estimate*(1+rses/100), color = rses < 100), width = 0.1) +
  facet_wrap(~parameter, scales = "free_y") 


# Use estimate that has best BICc when both effects are significant (also has smaller RSE)
best_run <- all_tabestimates %>% 
  filter(am == "am2", paste(am, run, sep = ".") %in% regime1) %>% 
  left_join(select(temp_pval, am, run, signif)) %>% 
  filter(signif == "both") %>% 
  select(-signif) %>% 
  filter(bicc == min(bicc))

igs <- all_initialvalues %>% filter(am == "am2", run == unique(best_run$run))

loadProject(here(out_folder, "IGS_am2.mlxtran"))
saveProject(here("results_eclipse", study, "best_model.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)
runScenario()
saveProject(here("results_eclipse", study, "best_model.mlxtran"))
generate_report(here("results_eclipse", study, "best_model.mlxtran"))

