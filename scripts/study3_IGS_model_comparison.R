# Load libraries
library(here)
library(ggsci)

source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
## Results folder
study <- "study3"
out_folder <- here("results", study, "IGS_model_comparison")
if(!dir.exists(out_folder)) dir.create(out_folder)

## Model file
model_file <- here("model", "2regime_model_final.txt")

## Data file
data_file <- here("data", "source", "ABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov",3))


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
initialValue_ranges_narrow <- read_csv(here("results", study, "initialValue_ranges_narrow.csv")) %>% 
  column_to_rownames("parameter")

parameters <- rownames(initialValue_ranges_narrow)

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

# AM1 #####
# TLR7 on m, Ab on beta and NP
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
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = T
)


for(name in names(parameter_result_am1)){
  assign(name, parameter_result_am1[[name]])
}

save(tabestimates, tabiters, phase, initialvalues, shrinkage, pvals, indivParams, indivFit, condition_number,
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

# plot_bicc_v_estimate(tabestimates)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(tabestimates, out_folder, project_name)


## Plot individual and population trajectories

indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams) %>% 
  left_join(distinct(tabestimates, run ,bicc)) 

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
      ggplot(aes(time, indivPred_mode, color = !!as.symbol(param), group = interaction(id, run))) + 
      geom_point(data=. %>% filter(run == 1), alpha = 0.5, aes(time, vl), inherit.aes = F) +
      geom_line(alpha = 0.25) + 
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1) + 
      geom_line(data =. %>% filter(bicc == min(bicc)), color = "black", lty = "dashed") 
  )
}
dev.off()

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "population_trajectories_by_param.pdf"), width = 10, height = 7)

for(param in names(indivParams)[2:7]){
  if(paste(param, "pop", sep = "_") %in% names(fixed)) next
  
  print(
    indivFit %>% 
      left_join(tabestimates %>% filter(parameter == str_interp("${param}_pop"))) %>%
      left_join(distinct(data, id, tlr7, ab) )  %>% 
      ggplot(aes(time, popPred, color = estimate, group = interaction( run))) + 
      geom_line(alpha = 0.25) + 
      facet_grid(ab~tlr7, labeller="label_both") +
      scale_color_viridis_c(direction = -1, option = "A") + 
      labs(color = paste0(param, "_pop")) + 
      geom_line(data =. %>% filter(bicc == min(bicc)), lty = "dashed", color = "black", alpha = 0.7) +
      
      tabestimates %>% 
      filter(parameter == str_interp("${param}_pop")) %>% 
      ggplot(aes(param, estimate)) + 
      geom_violin(scale = "width") + 
      geom_beeswarm(aes(color = bicc), alpha = 0.5) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), shape = 17, size = 2, aes(color = bicc)) +
      color_bicc() +
      
      plot_layout(ncol = 2, widths = c(2,0.5))
    
  )
  
}
dev.off()



# AM2 #####
# TLR7 on m, Ab on beta and p
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
  record_pvals = TRUE, record_indiv_params = TRUE,
  record_condition_number = T
)


for(name in names(parameter_result_am2)){
  assign(name, parameter_result_am2[[name]])
}

save(tabestimates, tabiters, phase, initialvalues, shrinkage, pvals, indivParams, indivFit, condition_number,
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

# plot_bicc_v_estimate(tabestimates)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name,
                    remove_bicc_outliers = T)

compare_estimates(tabestimates, remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers = T)$tabestimates, 
                  out_folder, project_name)

## PLot individual trajectories

indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams) %>% 
  left_join(distinct(tabestimates, run ,bicc)) 

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
      ggplot(aes(time, indivPred_mode, color = !!as.symbol(param), group = interaction(id, run))) + 
      geom_point(data=. %>% filter(run == 1), alpha = 0.5, aes(time, vl), inherit.aes = F) +
      geom_line(alpha = 0.25) + 
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1) + 
      geom_line(data =. %>% filter(bicc == min(bicc)), color = "black", lty = "dashed") 
  )
}
dev.off()

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "population_trajectories_by_param.pdf"), width = 10, height = 7)

for(param in names(indivParams)[2:7]){
  if(paste(param, "pop", sep = "_") %in% names(fixed)) next
  
  print(
    indivFit %>% 
      left_join(tabestimates %>% filter(parameter == str_interp("${param}_pop"))) %>%
      left_join(distinct(data, id, tlr7, ab) )  %>% 
      ggplot(aes(time, popPred, color = estimate, group = interaction( run))) + 
      geom_line(alpha = 0.25) + 
      facet_grid(ab~tlr7, labeller="label_both") +
      scale_color_viridis_c(direction = -1, option = "A") + 
      labs(color = paste0(param, "_pop")) + 
      geom_line(data =. %>% filter(bicc == min(bicc)), lty = "dashed", color = "black", alpha = 0.7) +
      
      tabestimates %>% 
      filter(parameter == str_interp("${param}_pop")) %>% 
      ggplot(aes(param, estimate)) + 
      geom_violin(scale = "width") + 
      geom_beeswarm(aes(color = bicc), alpha = 0.5) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), shape = 17, size = 2, aes(color = bicc)) +
      color_bicc() +
      
      plot_layout(ncol = 2, widths = c(2,0.5))
    
  )
  
}
dev.off()


# BH3 #####
# BH3: TLR7 on m, Ab on Betas, Ab and TLR7*AB on Np
alt_mod <- "bh3"
loadProject(here("results", study, "SCM", str_interp("${alt_mod}.mlxtran")))
project_name <- str_interp("IGS_${alt_mod}.mlxtran")
saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run convergence for narrower range of values
parameter_results_bh1 <- initial_guess_sampling(
  initialValue_ranges_narrow, log10 = NULL, nits = 100, 
  set_variance = 1, record_shrinkage = T, fixed = fixed, 
  predefined_initialvalues = predefined_initialvalues, record_pvals = TRUE,
  record_indiv_params = T, record_condition_number = T
)


for(name in names(parameter_results_bh1)){
  assign(name, parameter_results_bh1[[name]])
}

save(tabestimates, tabiters, phase, initialvalues, shrinkage, pvals, indivParams, indivFit, condition_number,
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

# plot_bicc_v_estimate(tabestimates)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name,
                    remove_bicc_outliers = T)

compare_estimates(tabestimates, remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name)

estimate_boxplots(filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers = T)$tabestimates, 
                  out_folder, project_name)

## PLot individual trajectories

indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams) %>% 
  left_join(distinct(tabestimates, run ,bicc)) 

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
      ggplot(aes(time, indivPred_mode, color = !!as.symbol(param), group = interaction(id, run))) + 
      geom_point(data=. %>% filter(run == 1), alpha = 0.5, aes(time, vl), inherit.aes = F) +
      geom_line(alpha = 0.25) + 
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1) + 
      geom_line(data =. %>% filter(bicc == min(bicc)), color = "black", lty = "dashed") 
  )
}
dev.off()

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "population_trajectories_by_param.pdf"), width = 10, height = 7)

for(param in names(indivParams)[2:7]){
  if(paste(param, "pop", sep = "_") %in% names(fixed)) next
  
  print(
    indivFit %>% 
      left_join(tabestimates %>% filter(parameter == str_interp("${param}_pop"))) %>%
      left_join(distinct(data, id, tlr7, ab) )  %>% 
      ggplot(aes(time, popPred, color = estimate, group = interaction( run))) + 
      geom_line(alpha = 0.25) + 
      facet_grid(ab~tlr7, labeller="label_both") +
      scale_color_viridis_c(direction = -1, option = "A") + 
      labs(color = paste0(param, "_pop")) + 
      geom_line(data =. %>% filter(bicc == min(bicc)), lty = "dashed", color = "black", alpha = 0.7) +
      
      tabestimates %>% 
      filter(parameter == str_interp("${param}_pop")) %>% 
      ggplot(aes(param, estimate)) + 
      geom_violin(scale = "width") + 
      geom_beeswarm(aes(color = bicc), alpha = 0.5) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), shape = 17, size = 2, aes(color = bicc)) +
      color_bicc() +
      
      plot_layout(ncol = 2, widths = c(2,0.5))
    
  )
  
}
dev.off()



# Compare results ####

ams <- c("am1", "am2", "bh3")

all_tabestimates <- NULL
all_initialvalues <- NULL
all_pvals <- NULL
all_indivParams <- NULL
all_indivFit <- NULL
all_condition_number <- NULL
all_shrinkage <- NULL
for(am in ams){
  load(here(out_folder, str_interp("IGS_${am}"), "IGS.RData"))
  all_tabestimates <- rbind(all_tabestimates, tabestimates %>% mutate(am = am))
  all_initialvalues <- rbind(all_initialvalues, initialvalues %>% mutate(am = am))
  all_pvals <- rbind(all_pvals, pvals %>% mutate(am = am))
  all_indivParams <- rbind(all_indivParams, indivParams %>% mutate(am = am))
  all_indivFit <- rbind(all_indivFit, indivFit %>% mutate(am = am))
  all_condition_number <- rbind(all_condition_number, condition_number %>% mutate(am = am))
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
# 1 am1      48
# 2 am2      51
# 3 bh3       1


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
# am1   TRUE             83
# am2   TRUE             91
# bh3   TRUE              4

## out of 100 runs, the model with Ab effect on p outperformed the model with Np 17 times
all_tabestimates %>% distinct(run, bicc, am, best) %>%
  pivot_wider(names_from = am, values_from = bicc) %>%
  count(am1 > (am2 + 2))

# `am1 > (am2 + 2)`     n
# FALSE                83
# TRUE                 17

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



# The interaction effect on Np in BH3 is never significant ####
all_pvals %>%
  mutate(padj = p.adjust(p.value, method = "fdr")) %>%
  ggplot(aes(am, log10(padj), color = padj >= 0.05)) +
  ggbeeswarm::geom_beeswarm() +
  facet_wrap(~parameter, ncol = 3) +
  geom_hline(yintercept = log10(0.05)) +
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


# Compare individual parameters bewteen best models
all_indivParams %>%
  pivot_longer(!c(id, run, am), names_to = "parameter", values_to = "EBE") %>%
  pivot_wider(names_from = am, values_from = EBE) %>%
  left_join(distinct(data, id, tlr7, ab) %>% mutate(id = as.character(id))) %>%
  ggplot(aes(am1, am2, color = interaction(tlr7,ab), fill = interaction(tlr7,ab))) +
  geom_abline() +
  geom_point(size = 2, alpha = 0.5) +
  # geom_point(data = bests, shape = 24, size = 2, color = "black") +
  facet_wrap(~parameter, scales = "free") +
  scale_x_log10() + scale_y_log10()

# Check condition numbers

all_condition_number %>% 
  left_join(distinct(tabestimates, run, bicc)) %>% 
  ggplot(aes(condition_number, am, color = bicc)) + geom_violin(scale = "width") + geom_beeswarm(alpha = 0.5)  + 
  geom_boxplot(outlier.shape = NA, width = 0.1, fill = NA) +
  scale_x_log10() + 
  geom_vline(xintercept = 100, linetype = 2) +
  color_bicc() + 
  geom_point(data = . %>% group_by(am) %>% filter(bicc == min(bicc)), size = 2, shape = 17, alpha = 1, color = "black")



# Plot shrinkage

all_shrinkage %>% 
  ggplot(aes(parameter, conditionalMode, color = am)) + 
  geom_boxplot() + 
  ylim(c(0,100))

# Biggest shrinkage for beta and p, followed by lambda and Np

# Check RSE for AM1

all_tabestimates %>% filter(am == "am1") %>% 
  group_by(am) %>% filter(bicc == min(bicc))

all_tabestimates %>% filter(am == "am2") %>% 
  group_by(am) %>% filter(bicc == min(bicc))

# Big RSE for m_pop and omega_betas





## Save best results ###

best_run <- all_tabestimates %>% filter(am == "am1") %>% filter(bicc == min(bicc)) %>% pull(run) %>% unique

igs <- all_initialvalues %>% filter(am == "am1", run == best_run) %>% select(-am, -run)

loadProject(here("results", "study3", "IGS_model_comparison", "IGS_am1.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)

saveProject(here("results", "study3", "best_model.mlxtran"))

runScenario()

saveProject(here("results", "study3", "best_model.mlxtran"))

generate_report(here("results", "study3", "best_model.mlxtran"))
