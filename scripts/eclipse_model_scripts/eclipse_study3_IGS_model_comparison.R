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

plot_bicc_v_estimate(tabestimates)

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)



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

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)







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

correlate_estimates(tabestimates, 
                    out_folder = out_folder, project_name = project_name)

correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)))

compare_estimates(tabestimates, 
                  out_folder = out_folder, project_name = project_name)






# Compare results ####

ams <- c("final_model","am2", "bh1.1", "bh3")

all_tabestimates <- NULL
all_initialvalues <- NULL
all_pvals <- NULL
all_indivParams <- NULL
all_indivFit <- NULL
for(am in ams){
  load(here(out_folder, str_interp("IGS_${am}"), "IGS.RData")) 
  all_tabestimates <- rbind(all_tabestimates, tabestimates %>% mutate(am = am))
  all_initialvalues <- rbind(all_initialvalues, initialvalues %>% mutate(am = am))
  all_pvals <- rbind(all_pvals, pvals %>% mutate(am = am))
  all_indivParams <- rbind(all_indivParams, indivParams %>% mutate(am = am))
  all_indivFit <- rbind(all_indivFit, indivFit %>% mutate(am = am))
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
# am1      75
# am2      25


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
# am1   TRUE             80
# am2   TRUE             51

## out of 100 runs, AM2 significantly out-performed AM1 20 times
all_tabestimates %>% distinct(run, bicc, am, best) %>% 
  pivot_wider(names_from = am, values_from = bicc) %>% 
  count(am1 > (am2 + 2))

# `am1 > (am2 + 2)`     n
# FALSE                80
# TRUE                 20

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
  annotate(geom = "text", x = c(1:2), y = 470, label = paste(rev(ncovar$n), "covariates"), hjust = 1)



# How many times does AM2 outperform AM1?
plot_df <- all_tabestimates %>% 
  distinct(run, best, bicc, am) %>% 
  pivot_wider(names_from = am, values_from = bicc) %>%
  mutate(color = case_when(am1 > (am2 + 2) ~ "Alternative model better", am2 > (am1 + 2) ~ "SCM model better", TRUE ~ "Equal performance"))

# There are two clusters where AM2 performs better than AM1:
AM2_better_c1 <- plot_df %>% filter(am2 < 430 & am1 > 430 & am1 > (am2 + 2) ) %>% pull(run) %>% unique 
AM2_better_c2 <- plot_df %>% filter(am2 > 430 & am1 > 438 & am1 > (am2 + 2 )) %>% pull(run) %>% unique 

plot_df %>%
  mutate(shape = ifelse(run %in% AM2_better_c1, "Larger m_pop in AM2", 
                        ifelse(run %in% AM2_better_c2, "Larger TLR7:lambda in AM2", "AM1 equal or better"))) %>%
  ggplot(aes(am1, am2)) + 
  
  geom_abline(alpha = 0.7) +
  geom_abline(intercept = -2, lty = "dashed", alpha = 0.7) +
  geom_abline(intercept = 2, lty = "dashed", alpha = 0.7)  +
  
  # Shaded regions using geom_ribbon
  geom_ribbon(aes(ymin = -Inf, ymax = am1 - 2), fill = "darkred", alpha = 0.1) +  # Below bottom dashed line
  annotate("text", x = 427, y = 420, label = str_interp("Alternative\nmodel better\n(${count(plot_df, color)$n[1]})"), size = 4, color = "darkred", fontface = "bold") +
  
  geom_ribbon(aes(ymin = am1 - 2, ymax = am1 + 2), fill = "darkgreen", alpha = 0.1) +  # Between dashed lines
  annotate("text", x = 430, y = 430, label = str_interp("Equal\nperformance\n(${count(plot_df, color)$n[2]})"), size = 4, color = "darkgreen", fontface = "bold") +
  
  geom_ribbon(aes(ymin = am1 + 2, ymax = Inf), fill = "darkgray", alpha = 0.1) +   # Above top dashed line
  annotate("text", x = 425, y = 450, label = str_interp("SCM model better\n(${count(plot_df, color)$n[3]})"), size = 4, color = "darkgray", fontface = "bold") +
  
  geom_point(data = . %>% filter(color != "Alternative model better"), aes(color = color), alpha = 0.7, size = 2) +
  geom_point(data = . %>% filter(color == "Alternative model better"), aes(color = color, shape = shape), alpha = 0.7, size = 2) +
  labs(x = "TLR7 on ta and lambda", y = "TLR7 on ta, lambda, and Np", shape = "") +
  scale_color_manual(values = c("darkred", "darkgreen", "darkgray")) + 
  scale_shape_manual(values = c(15, 17)) +
  guides(color = "none") 

# AM2 outperforms AM1 for two different reasons:
# If both models assume delayed reactivation is the mecnaism for control in aviremic monkeys,
# AM2 does better when AM1 estimates a smaller TLR7 effect on lambda (cluster 1)
# If AM2 assumes immune control is the mechanism for control in aviremic monkeys,
# it out performs AM1 when AM1 estimates a smaller m_pop (meaning it infers dealyed rebound as the mechanism)
all_tabestimates %>% 
  select(run, parameter, estimate,  am) %>% 
  pivot_wider(names_from = am, values_from = estimate) %>% 
  mutate(color = ifelse(run %in% AM2_better_c1, "AM2 better because of m_pop", 
                        ifelse(run %in% AM2_better_c2, "AM2 better because of TLR7:lambda", "AM1 equal or better"))) %>%
  ggplot(aes(am1, am2, color = color)) + 
  geom_point() + 
  geom_abline() +
  facet_wrap(~parameter, scales = "free")

all_indivFit %>% filter(run %in% AM2_better_c1) %>% 
  ggplot(aes(time, indivPred_mode, color = am, group = interaction(id, run, am))) + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id)

all_indivFit %>% filter(run %in% AM2_better_c2) %>% 
  ggplot(aes(time, indivPred_mode, color = am, group = interaction(id, run, am))) + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id)

immune_control_am1 <- all_indivFit %>% 
  filter(am == "am1", id == 15, time < 75, indivPred_mode > 0.7) %>% 
  pull(run) %>% unique
immune_control_am2 <- all_indivFit %>% 
  filter(am == "am2", id == 15, time < 75, indivPred_mode > 0.7) %>% 
  pull(run) %>% unique


# When both models infer mechanism of control to be delayed reactivation, the model with a TLR7 effect 
# on Np signicantly outperforms the model without one in only 8% of runs
plot_df %>% filter((!run %in% immune_control_am1) & (!run %in% immune_control_am2)) %>% 
  count(color)

# If both models infer the mechanism of control is a delay in rebound, AM2 outperforms AM1 
# when AM1 estimates a smaller TLR7 effect on lambda
all_tabestimates %>% 
  left_join(distinct(plot_df, run, color)) %>% 
  mutate(both_delayed = (!run %in% immune_control_am1) & (!run %in% immune_control_am2)) %>% 
  filter(both_delayed) %>% 
  select(parameter, estimate, run, am, color) %>%
  pivot_wider(names_from = am, values_from = estimate) %>% 
  ggplot(aes(am1, am2, color = color)) +
  geom_abline() +
  geom_point(alpha = 0.7, size = 2) + 
  facet_wrap(~parameter, scales = "free") 
  

# The TLR7 effect on Np is never significant ####
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
  mutate(am = factor(am, levels = c("am1", "am2"), labels = c("no TLR7 effect on Np", "TLR7 effect on Np")),
    parameter = factor(parameter, levels = paste0("omega_", c("lambda", "betas", "ta", "NPs", "p", "m")))) %>% 
  ggplot(aes(am, exp(estimate))) + 
  geom_violin() +
  # geom_boxplot(width = 0.1) +
  ggbeeswarm::geom_beeswarm(aes(color = bicc), alpha = 0.7) + 
  facet_wrap(~parameter, scales = "free_x") + 
  coord_flip() + 
  annotate(geom="point", x = 1, y = 1, shape = NA)  +
  labs(y = "Random effect estimate (fold change)", x = "Model") + 
  color_bicc()



all_indivParams %>% 
  pivot_longer(!c(id, run, am), names_to = "parameter", values_to = "EBE") %>% 
  pivot_wider(names_from = am, values_from = EBE) %>% 
  left_join(distinct(data, id, tlr7, vac) %>% mutate(id = as.character(id))) %>% 
  ggplot(aes(am1, am3, color = interaction(tlr7,vac), fill = interaction(tlr7,vac))) + 
  geom_abline() + 
  geom_point(size = 2, alpha = 0.5) + 
  geom_point(data = bests, shape = 24, size = 2, color = "black") +
  facet_wrap(~parameter, scales = "free") + 
  scale_x_log10() + scale_y_log10()


plot_indiv_traj <- function(indivDF, data, ID, param, option = "A", 
                            boxplot_coord = c(NA, NA, 0, 3), boxplot_position = "inset",
                            log_boxplot = F, log10_color = F){
  plot_df <- indivDF %>% filter(id == ID) %>% 
    rename(color = str_interp("indivParam_${param}")) 
  data_df <- data %>% filter(id == ID)
  
  if(param == "NPs") plot_df$color <- plot_df$color*1e5
  
  plot <- plot_df %>% 
    ggplot(aes(time, indivPred_mode, group = interaction(id, run), color = color)) 
  
  if(log10_color) plot <- plot + geom_line(alpha = 0.25, aes(color = log10(color)))
  
  plot <- plot + 
    geom_line(data = . %>% filter(bicc == min(bicc)), color = "black", linewidth = 0.5) + 
    scale_color_viridis_c(direction = -1, option = option) + 
    geom_point(data = data_df, aes(time, vl), inherit.aes = F) + 
    labs(x = "Time", y = "Predicted Viral Load", color = ifelse(log10_color, str_interp("log10(${param})") ,param), 
         title = str_interp("Predictions for animal ${ID}"))
  
  boxplot <- plot_df %>% 
    distinct(id, run, color, bicc) %>% 
    ggplot(aes("", color)) + 
    geom_boxplot(outlier.shape = NA) + 
    ggbeeswarm::geom_beeswarm(aes(color = bicc), alpha = 0.7) + 
    scale_color_viridis_c(direction = -1) + 
    theme_minimal() +
    theme(axis.title.x = element_blank(), 
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_rect(color = NA),
          legend.key.size = unit(0.5, "lines"),
          legend.text = element_text(size = 5),
          legend.title = element_text(size = 6)) +
    labs(y = param, color = "BICc")
  
  if(log_boxplot){
    boxplot <- boxplot + scale_y_log10()
  }
  
  
  if(boxplot_position == "inset"){
    boxplot <- ggplotGrob(boxplot)
    
    xmax <- data_df %>% pull(time) %>% max()
    xmax <- 0.9*xmax
    xmin <- ifelse(is.na(boxplot_coord[1]), xmax*0.6, boxplot_coord[1])
    xmax <- ifelse(is.na(boxplot_coord[2]), xmax, boxplot_coord[2])
    
    ymax <- ifelse(max(data_df$vl) < boxplot_coord[4], 0.9*max(data_df$vl), boxplot_coord[4])
    
    plot + annotation_custom(boxplot, xmin = xmin, xmax = xmax, 
                             ymin = boxplot_coord[3], ymax = ymax)
  }else if(boxplot_position == "right"){
    plot + boxplot + plot_layout(ncol = 2, widths = c(2,0.5)) 
  }
}

# Save best model results ###

# Using AM1, best results from regime with lower m

best_run <- all_tabestimates %>% filter(am == "am1", !run %in% am1_immune_regime) %>% filter(bicc == min(bicc)) %>% pull(run) %>% unique

igs <- all_initialvalues %>% filter(am == "am1", run == best_run) %>% select(-am, -run)

loadProject(here("results_eclipse", "study1", "IGS_model_comparison", "IGS_am1.mlxtran"))

set_base(igs %>% rename(estimate = initialValue), fixed = fixed)

saveProject(here("results_eclipse", "study1", "best_model.mlxtran"))

runScenario()

saveProject(here("results_eclipse", "study1", "best_model.mlxtran"))


generate_report(here("results_eclipse", "study1", "best_model.mlxtran"))
