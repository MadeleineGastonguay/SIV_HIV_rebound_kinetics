# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "study3"
out_folder <- here("results_eclipse", study, "IGS_model_comparison")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "source", "ABTLR7_DATA.txt")
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

predefined_initialvalues <- read_csv(here(out_folder, "initial_values.csv"))

predefined_initialvalues %>% 
  ggplot(aes(parameter, initialValue)) +
  geom_boxplot() + 
  facet_wrap(~parameter, scales = "free")

# Run IGS #####

# Final SCM model #####
# Final model from SCM:
# a on beta and NP, tlr7 on m
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

plot_vs_rses(tabestimates, out_folder, project_name)

plot_bicc_v_estimate(
  tabestimates, log10_x = F,
  out_folder = out_folder, project_name = project_name
)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

# Still a correlation between beta and ta and ta and m, but there isn't too much variation in the parameters
correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

table_predictions <- left_join(table_predictions, distinct(tabestimates, run, bicc))
pdf(here(out_folder, gsub(".mlxtran", "", project_name), "table_predictions.pdf"), width = 10, height = 10)

table_predictions %>% 
  ggplot(aes(time, V_mode,  group = factor(run), color = bicc))  + 
  geom_point(data = data, aes(time, 10^vl), inherit.aes = F, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

table_predictions %>% 
  ggplot(aes(time, P_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

table_predictions %>% 
  ggplot(aes(time, E_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

table_predictions %>% 
  ggplot(aes(time, I_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

table_predictions %>% 
  ggplot(aes(time, T_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

dev.off()

# AM2 #####
# AM2 from SCM 
# swap ab effect on p for Np 
alt_mod <- "am2"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_result_am2 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_result_am2)){
  assign(name, parameter_result_am2[[name]])
}

save(list = names(parameter_result_am2),
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

# Still a strong correlation between beta and ta and ta and m
correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)

table_predictions <- left_join(table_predictions, distinct(tabestimates, run, bicc))
pdf(here(out_folder, gsub(".mlxtran", "", project_name), "table_predictions.pdf"), width = 10, height = 10)

table_predictions %>% 
  ggplot(aes(time, V_mode,  group = factor(run), color = bicc))  + 
  geom_point(data = data, aes(time, 10^vl), inherit.aes = F, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

table_predictions %>% 
  ggplot(aes(time, P_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

table_predictions %>% 
  ggplot(aes(time, E_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

table_predictions %>% 
  ggplot(aes(time, I_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

table_predictions %>% 
  ggplot(aes(time, T_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

dev.off()


# BH1.1 #####
# ab effect on ta, Np, beta and tlr7 effect on m
alt_mod <- "bh1.1"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_bh1 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_bh1)){
  assign(name, parameter_results_bh1[[name]])
}

save(list = names(parameter_results_bh1),
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

# still a strong correlation between beta and ta and m and ta
correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates,
                  out_folder = out_folder, project_name = project_name)


table_predictions <- left_join(table_predictions, distinct(tabestimates, run, bicc))
pdf(here(out_folder, gsub(".mlxtran", "", project_name), "table_predictions.pdf"), width = 10, height = 10)

table_predictions %>% 
  ggplot(aes(time, V_mode,  group = factor(run), color = bicc))  + 
  geom_point(data = data, aes(time, 10^vl), inherit.aes = F, alpha = 0.5) +
  geom_line(alpha = 0.5) +
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

table_predictions %>% 
  ggplot(aes(time, P_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc()

table_predictions %>% 
  ggplot(aes(time, E_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

table_predictions %>% 
  ggplot(aes(time, I_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

table_predictions %>% 
  ggplot(aes(time, T_mode,  group = factor(run), color = bicc))  + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() + 
  color_bicc() 

dev.off()



# BH3 #####
# ab effect on Np, beta and tlr7 effect on m, interaction on Np
alt_mod <- "bh3"
loadProject(here("results_eclipse", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_bh3 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = nits, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, 
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = TRUE, record_table = TRUE
)


for(name in names(parameter_results_bh3)){
  assign(name, parameter_results_bh3[[name]])
}

save(list = names(parameter_results_bh3),
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

# Still a correlation between ta an beta and m and ta
correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)




# Compare results ####

ams <- c("final_model","am2", "bh1.1", "bh3")

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
# final_model    15
# am2            35

# The biological hypotheses are never the best

## out of 50 runs, the number of times each model was less than 2 BICc from the best model
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  rename(best_model = best) %>% 
  group_by(run) %>% 
  mutate(best = min(bicc),
         equal_to_best = abs(bicc - best) <= 2) %>% 
  ungroup %>% 
  filter(equal_to_best) %>% 
  count(am, equal_to_best) 

# am          equal_to_best     n
# final_model TRUE             40
# am2         TRUE             49
# bh1.1       TRUE              1
# bh3         TRUE              1


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
  annotate(geom = "text", x = c(1:nrow(ncovar)), y = 770, label = paste((ncovar$n), "covariates"), hjust = 1)


# How frequently are the treatment effects significant?
# tlr7*ab effect is never significant on Np (BH3)
# ab effect is never significant on ta (BH1.1)

all_pvals %>% 
  mutate(padj = p.adjust(p.value, method = "fdr")) %>% 
  ggplot(aes(am, log10(padj), color = padj < 0.05)) + 
  ggbeeswarm::geom_beeswarm() + 
  facet_wrap(~parameter, ncol = 3) +
  geom_hline(yintercept = log10(0.05)) + 
  labs(color = "Low m, low ta regime") +
  scale_color_manual(values = c("black", "gray")) +
  theme(legend.position = "none")


# Look at size of random effects for each model ####
# There's no model with clearly large random effects

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


# Compare condition numbers:
all_condition_numbers %>%
  ggplot(aes(am, condition_number)) +
  geom_violin() + 
  geom_beeswarm() + 
  geom_beeswarm(data = . %>% left_join(distinct(all_tabestimates, am, run, bicc)) %>% 
                  group_by(am) %>% filter(bicc == min(bicc)),
                shape = 17, size = 4, color = "red") + 
  scale_y_log10() +
  geom_hline(yintercept = 100, linetype = "dashed") +
  coord_flip()


## look at RSE for best fits of two similar models:
# The best model BICs are indistinguishible
all_tabestimates %>% filter(am %in% c("final_model", "am2")) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  arrange(parameter, am) %>% filter(abs(rses) > 50) %>% view


## Compare EBEs for best models:

pop_param <- all_tabestimates %>% 
  filter(am %in% c("final_model", "am2")) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) 

pop_param <- rbind(
  pop_param %>% 
    filter(grepl("pop", parameter)) %>% 
    mutate(parameter = gsub("_pop", "", parameter), tlr7 = 0, ab = 0) %>% 
    select(parameter, value = estimate, tlr7, ab) ,
  
  pop_param %>% 
    filter(grepl("pop", parameter) | grepl("tlr7_1", parameter)) %>% 
    rowwise() %>% 
    mutate(p = ifelse(grepl("pop", parameter), gsub("_pop", "", parameter),
                      unlist(str_split(parameter, "_"))[2]), 
           estimate = ifelse(grepl("pop", parameter), estimate, exp(estimate))) %>% 
    select(p, parameter, estimate) %>% 
    group_by(am, p) %>% summarise(value = prod(estimate, na.rm = TRUE)) %>%
    mutate(value = ifelse(p == "lambda", value*100, value), tlr7 = 1, ab = 0) %>% 
    rename(parameter = p),
  
  pop_param %>% 
    filter(grepl("pop", parameter) | grepl("ab_1", parameter)) %>% 
    rowwise() %>% 
    mutate(p = ifelse(grepl("pop", parameter), gsub("_pop", "", parameter),
                      unlist(str_split(parameter, "_"))[2]), 
           estimate = ifelse(grepl("pop", parameter), estimate, exp(estimate))) %>% 
    select(p, parameter, estimate) %>% 
    group_by(am, p) %>% summarise(value = prod(estimate, na.rm = TRUE)) %>%
    mutate(value = ifelse(p == "lambda", value*100, value), tlr7 = 0, ab = 1) %>% 
    rename(parameter = p),
  
  pop_param %>% 
    filter(grepl("pop", parameter) | grepl("beta_", parameter)) %>% 
    rowwise() %>% 
    mutate(p = ifelse(grepl("pop", parameter), gsub("_pop", "", parameter),
                      unlist(str_split(parameter, "_"))[2]), 
           estimate = ifelse(grepl("pop", parameter), estimate, exp(estimate))) %>% 
    select(p, parameter, estimate) %>% 
    group_by(am, p) %>% summarise(value = prod(estimate, na.rm = TRUE)) %>%
    mutate(value = ifelse(p == "lambda", value*100, value), tlr7 = 1, ab = 1) %>% 
    rename(parameter = p)
)

## There is some evidence for a tlr7 effect on m, but is it strong enough to include in the model???
all_indivParams %>%
  filter(am %in% c("final_model", "am2")) %>% 
  left_join(distinct(all_tabestimates, am, bicc, run)) %>% 
  left_join(distinct(data, id, tlr7, ab) %>% 
              mutate(id = as.character(id))) %>% 
  group_by(am) %>% filter(bicc == min(bicc)) %>% 
  select(-run, -bicc) %>% 
  pivot_longer(!c(id, am, tlr7, ab), names_to = "parameter")  %>% 
  ggplot(aes(interaction(tlr7, ab), value, color = am)) + 
  geom_violin(scale = "width") + geom_beeswarm(dodge.width = 0.9) + 
  facet_wrap(~parameter, scales = "free_y") + 
  scale_y_log10() +
  geom_point(data = pop_param, shape = 17, size = 4, position = position_dodge(width = 0.9))


