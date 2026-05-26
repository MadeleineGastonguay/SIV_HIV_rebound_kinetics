# Load libraries
library(here)
library(RColorBrewer)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# If FALSE, results will be loaded from previous runs of this script (after running project setup)
RERUN <- FALSE

# Project setup ####
## Results folder
study <- "study1"
project_name <- "IGS_no_cov_fixed_param.mlxtran"
out_folder <- here("results_eclipse", study)
if(!dir.exists(out_folder)) dir.create(out_folder)

## Parameters to fix:
fixed <- c("lambda_pop" = 200, "p_pop" = 3)

## Model file
model_file <- here("model", "2regime_model_delay.txt")

## Data file
data_file <- here("data", "derived", str_interp("${toupper(study)}_DATA.txt"))
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", "catcov")

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]


# Run Analysis ####

if(RERUN){
  ## Create Project
  create_project(model_file, data_file, headerTypes, error_model,
                 exclude_random_effects, treatment_effects = NULL,
                 normal_dist = NULL, fixed = fixed)
  
  saveProject(here(out_folder, project_name))
  
  
  # Initial guess sampling from a narrower range of values ####
  
  ## Calculate narrower ranges for IGS:
  initialValue_ranges_narrow <- read_csv(here(out_folder, "initialValue_ranges_narrow.csv")) %>% column_to_rownames("parameter")

  
  ## initial guess sampling from narrow range
  # Use initial condition of 1 for all random effects
  parameter_results <- initial_guess_sampling(
    initialValue_ranges_narrow, 
    nits = 100,
    record_shrinkage = T,
    record_indiv_params = T,
    set_variance = 1,
    fixed = fixed,
    record_table = T
  )
  
  ## Save results
  for(name in names(parameter_results)){
    assign(name, parameter_results[[name]])
  }
  
  save(list = names(parameter_results),
       file = here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
  
}else{
  # Load Project
  loadProject(here(out_folder, project_name))  
  # Load results from IGS over narrow range
  load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))
}


# Visualize results from narrow IGS ####

plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name)

regime1 <- tabestimates %>% 
  filter(parameter == "omega_m", estimate < 3.5) %>% 
  pull(run)

# 62% of runs fall in the regime of latent reservoir reduction (with worse BICc)
length(regime1)

plot_bicc_v_estimate(tabestimates, remove_bicc_outliers = T, log10_x = F, color = regime1)

plot_bicc_v_estimate(tabestimates %>% filter(run %in% regime1), remove_bicc_outliers = T, log10_x = F,
                     out_folder = out_folder, project_name = project_name, file_addendum = "_low_m")

## Create figure showing difference in regimes: ####
dual_regime<-tabestimates %>% 
  mutate(regime = ifelse(run %in% regime1, "latent reservoir", "immune control")) %>% 
  group_by(regime) %>% 
  filter(bicc == min(bicc)) %>% 
  ungroup()

df <- merge(
  dual_regime %>% 
    filter(grepl("pop", parameter)) %>% 
    mutate(parameter = gsub("_pop", "", parameter)) %>% 
    select( parameter, mean = estimate, regime) %>% 
    rbind(data.frame(parameter = rep(c("lambda", "p"), each = 2), 
                     mean = rep(unname(fixed), each = 2 ),
                     regime = c("latent reservoir", "immune control"))),
  
  dual_regime %>% 
    filter(grepl("omega", parameter)) %>% 
    mutate(parameter = gsub("omega_", "", parameter)) %>% 
    select(parameter, sd = estimate, regime)
) 

density_data <- df %>%
  group_by(parameter) %>%
  group_modify(~ {
    x_vals <- exp(seq((min(log(.x$mean) - 4*.x$sd)),
                      (max(log(.x$mean) + 4*.x$sd)),
                      length.out = 500))
    do.call(rbind, lapply(1:nrow(.x), function(i) {
      data.frame(
        x = x_vals,
        y = dnorm(log(x_vals), mean = log(.x$mean[i]), sd = .x$sd[i]),
        regime = as.factor(.x$regime[i])
      )
    }))
  }) %>%
  ungroup()

make_dual_regime_distribution_plot <- function(param, limits = NULL, breaks = NULL){
  dual_regime_distribution_plot <- density_data %>%
    left_join(dual_regime %>% distinct(run, bicc, regime) %>% mutate(run = factor(run))) %>%
    filter(parameter == param) %>% 
    mutate(x = case_when(#parameter == "betas" ~ x*1e-7, 
      parameter == "NPs" ~ x*1e5,
      .default = x),
      parameter = factor(parameter, 
                         levels =  c("lambda","betas", "ta", "m", "p", "NPs"),
                         labels =  c("lambda","beta~'(x'~10^-7~')'", "t['a']", "m", "p", "N['p']"))) %>%
    ggplot(aes(x = log10(x), y = y, color = regime, group = regime)) +
    geom_line(linewidth = 1) +
    facet_wrap(~parameter, scales = "free", ncol = 1, labeller = label_parsed, strip.position = "bottom") +
    theme_minimal() + 
    guides(color = guide_legend(position = "inside")) +
    theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
          axis.ticks.y = element_blank(), panel.grid = element_blank(), axis.line.x = element_line(),
          strip.text = element_text(size = 15), strip.placement = "outside",
          legend.position.inside = c(0,1), legend.justification = c(0,1), legend.text = element_text(size = 12),
          legend.title = element_text(size = 12)) + 
    scale_color_brewer(name = "Regime", palette = "Set1", labels = c("Immune\ncontrol", "Reduced\nreservoir")) +
    labs(x = "")
  
  # Adjust limits:
  if(param %in% c("lambda", "betas", "p")){
    dual_regime_distribution_plot <- dual_regime_distribution_plot + 
      scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(breaks), limits = limits)
  }else{
    dual_regime_distribution_plot <- dual_regime_distribution_plot + 
      scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = log10(breaks), limits = limits)
  }
  
  dual_regime_distribution_plot
  
}

make_dual_regime_EBE_plot <- function(param, limits = NULL, breaks = NULL, log10 = FALSE, ...){
  df <- indivParams %>% 
    merge(distinct(dual_regime, run, regime)) %>% 
    pivot_longer(c(lambda, betas, ta, NPs, p, m), names_to = "parameter", values_to = "EBE") %>% 
    filter(parameter == param) %>% 
    # mutate(EBE = ifelse(parameter == "NPs", EBE*1e5, ifelse(parameter == "betas", EBE*1e-7, EBE))) %>% 
    mutate(EBE = ifelse(parameter == "NPs", EBE*1e5, EBE)) %>% 
    select(id, parameter, EBE, regime) %>% 
    pivot_wider(names_from = regime, values_from = EBE) %>% 
    left_join(distinct(mutate(data, id = as.character(id)), id, tlr7)) %>% 
    mutate(parameter = factor(parameter, levels = c("lambda", "betas", "ta", "m", "p", "NPs"),
                              labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")))
  
  if(log10){
    df <- df %>% mutate(`latent reservoir` = log10(`latent reservoir`), `immune control` = log10(`immune control`))
  }
  
  # Axis labels:
  param_label_x <- switch(param, 
                          "lambda" = expression(atop(lambda*" EBE in reduced",  " reservoir regime")),
                          "betas" = expression(atop(beta*" EBE in reduced",  " reservoir regime (x"*10^-7*")")),
                          "ta" = expression(atop(t['a']*" EBE in reduced",  " reservoir regime")),
                          "m" = expression(atop("m"*" EBE in reduced",  " reservoir regime")),
                          "p" = expression(atop("p"*" EBE in reduced",  " reservoir regime")),
                          "NPs" = expression(atop(N['p']*" EBE in reduced",  " reservoir regime"))
  )
  
  param_label_y <- switch(param, 
                          "lambda" = expression(atop(lambda*" EBE in immune",  "control regime")),
                          "betas" = expression(atop(beta*" EBE in immune",  "control regime (x"*10^-7*")")),
                          "ta" = expression(atop(t['a']*" EBE in immune",  "control regime")),
                          "m" = expression(atop("m"*" EBE in immune",  "control regime")),
                          "p" = expression(atop("p"*" EBE in immune",  "control regime")),
                          "NPs" = expression(atop(N['p']*" EBE in immune",  "control regime"))
  )
  
  dual_regime_EBE_plot <- df %>% ggplot(aes((`latent reservoir`), (`immune control`), color = factor(tlr7))) + 
    geom_abline( lty = "dashed", alpha = 0.5) + 
    geom_point(size = 3, alpha = 0.5) + 
    ggrepel::geom_text_repel(data = . %>% filter(id %in% c(15, 19)),
                             aes(label = id), min.segment.length = 0,
                             show.legend = F, box.padding = 0.6) +
    scale_color_manual(values = c("0" = "black", "1" = "#00A3D3"), labels = c("Control", "TLR7")) + 
    theme_minimal() +
    theme(legend.title = element_blank(), strip.text = element_text(size = 15)) + 
    # labs(x = bquote(str_interp("${param_label} EBE in reduced reservoir regime")), y = bquote(str_interp("${param_label} EBE in immune control regime")))
    labs(x = param_label_x, y = param_label_y)
  
  if(!log10){
    dual_regime_EBE_plot + 
      scale_x_continuous(limits = (limits), ...) + 
      scale_y_continuous(limits = (limits), ...) 
  }else{
    if(log10(diff(limits)) > 2){
      dual_regime_EBE_plot + 
        scale_x_continuous(limits = log10(limits), labels = scales::label_math(expr = 10^.x), breaks = log10(breaks)) + 
        scale_y_continuous(limits = log10(limits), labels = scales::label_math(expr = 10^.x), breaks = log10(breaks))  
    }else{
      dual_regime_EBE_plot + 
        scale_x_continuous(limits = log10(limits), labels = function(breaks) 10^breaks, breaks = log10(breaks)) + 
        scale_y_continuous(limits = log10(limits), labels = function(breaks) 10^breaks, breaks = log10(breaks))
    }
    
  }
  
}

dual_regime_EBE_plot <- make_dual_regime_EBE_plot("lambda", limits = c(90,400), log10 = F) + 
  guides(color = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0,1), legend.justification = c(-0.1,1)) +
  make_dual_regime_EBE_plot("betas", limits = c(5,10), log10 = F) + 
  theme(legend.position = "none") +
  make_dual_regime_EBE_plot("ta", limits = c(0.1,100), log10 = T, breaks = c(0.1, 1, 10, 100)) + 
  theme(legend.position = "none") +
  make_dual_regime_EBE_plot("m", limits = 10^c(-2,6), log10 = T, breaks = 10^seq(-2,6,by=2)) + 
  theme(legend.position = "none") +
  make_dual_regime_EBE_plot("p", limits = c(2,4), log10 = F) + 
  theme(legend.position = "none") +
  make_dual_regime_EBE_plot("NPs", limits = 10^c(4.75,7.25), log10 = T, breaks = 10^seq(5,7))  + 
  theme(legend.position = "none") +
  plot_layout(ncol = 1) &
  theme(#axis.text.x = element_text(size = 12, color = brewer.pal(2,"Set1")[2]),
    #axis.text.y = element_text(size = 12, color = brewer.pal(2,"Set1")[1]),
    axis.title.x = element_text(size = 12, color = brewer.pal(2,"Set1")[2]),
    axis.title.y = element_text(size = 12, color = brewer.pal(2,"Set1")[1]))


dual_regime_distribution_plot <- make_dual_regime_distribution_plot("lambda", breaks = (c(30, 100, 300, 1000))) + 
  guides(color = guide_legend(position = "inside")) +
  theme(legend.position.inside = c(0,1), legend.justification = c(0,1)) +
  make_dual_regime_distribution_plot("betas", breaks = (c(5, 7, 10))) +
  theme(legend.position = "none") +
  make_dual_regime_distribution_plot("ta", breaks = 10^seq(-2,2)) +
  theme(legend.position = "none") +
  make_dual_regime_distribution_plot("m", breaks = 10^seq(-4,10,by=4)) +
  theme(legend.position = "none") +
  make_dual_regime_distribution_plot("p", breaks = (c(1, 3, 10))) +
  theme(legend.position = "none") +
  make_dual_regime_distribution_plot("NPs", breaks = 10^seq(3,9,by=2)) +
  theme(legend.position = "none") +
  plot_layout(ncol = 1) 




non_rebound_fit <- indivFit %>% 
  filter(id %in% c(15, 19)) %>% 
  merge(distinct(dual_regime, run, regime)) %>% 
  ggplot(aes(time, indivPred_mode)) + 
  geom_hline(yintercept = 1.7, lty = "dashed", color = "gray") +
  geom_point(aes(y = vl)) +
  geom_line(linewidth = 1, aes(color = regime)) +
  facet_wrap(~id) + 
  scale_color_brewer(palette = "Set1", labels = c("Immune control regime", "Reduced reservoir regime")) +
  labs(x = "Time after ART stop (days)", y = "Viral load\n(copies/ml)") + 
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), legend.position = "none") + 
  scale_y_continuous(labels = scales::label_math(expr = 10^.x), limits = c(-1,6)) 

non_rebound_immune <- table_predictions %>% 
  filter(id %in% c(15, 19)) %>% 
  merge(distinct(dual_regime, run, regime)) %>% 
  select(regime, id, time, contains("mode")) %>% 
  ggplot(aes(time,log10( P_mode), color = regime, group = regime)) + 
  geom_line(linewidth = 1) +
  facet_wrap(~id) + 
  scale_color_brewer(palette = "Set1", labels = c("Immune control regime", "Reduced reservoir regime")) +
  labs(x = "Time after ART stop (days)", y = "Precursor\nImmune Cells") + 
  theme_minimal() +
  theme(panel.border = element_rect(color = "black", fill = NA),
        strip.background = element_blank(), panel.grid = element_blank(),
        plot.title = element_text(face = "bold"), legend.position = "none") +
  scale_y_continuous(labels = scales::label_math(expr = 10^.x), limits = c(2,9))

bicc_plot_df <- tabestimates %>% 
  filter(parameter != "a") %>% 
  mutate(regime = ifelse(run %in% regime1, "latent reservoir", "immune control"),
         parameter = factor(parameter, levels = c("lambda_pop", "betas_pop", "ta_pop","m_pop", "p_pop", "NPs_pop",
                                                  "omega_lambda", "omega_betas", "omega_ta", "omega_m", "omega_NPs", "omega_p"),
                            labels = c("lambda['pop']", "beta['pop']", "t['a, pop']", "m['pop']", "p['pop']", "N['p, pop']",
                                       "Omega[lambda]", "Omega[beta]", "Omega['ta']", "Omega['m']", "Omega['NP']", "Omega['p']"))) 

dual_regime_bicc <- bicc_plot_df %>% 
  filter(parameter %in% c("Omega[lambda]", "Omega['p']")) %>% 
  ggplot(aes(estimate, bicc, color = regime)) + 
  
  bicc_plot_df %>% 
  filter(grepl("beta", parameter)) %>% 
  ggplot(aes(estimate, bicc, color = regime)) + 
  
  bicc_plot_df %>% 
  filter(parameter %in% c("Omega['ta']", "t['a, pop']")) %>% 
  ggplot(aes(estimate, bicc, color = regime)) + 
  
  bicc_plot_df %>% 
  filter(parameter %in% c("Omega['m']", "m['pop']")) %>% 
  ggplot(aes(estimate, bicc, color = regime)) + 
  
  bicc_plot_df %>% 
  filter(parameter %in% c("Omega['NP']", "N['p, pop']")) %>% 
  ggplot(aes(estimate, bicc, color = regime)) +
  
  plot_layout(ncol = 1, guides = "collect") &
  
  
  geom_point(alpha = 0.75, size = 3) & 
  geom_point(data = . %>% group_by(regime) %>% filter(bicc == min(bicc)),
             size = 6, shape = 23, color = "black", aes(fill = regime)) &
  facet_wrap(~parameter, scales = "free_x", ncol = 2, labeller = "label_parsed", dir = "v", strip.position = "bottom") & 
  scale_color_brewer(palette = "Set1") & 
  scale_fill_brewer(palette = "Set1") & 
  labs(y = "BICc") & 
  theme_minimal() &
  theme(legend.position = "none", strip.placement = "outside", axis.title.x = element_blank(), strip.text = element_text(size = 15))
 
study1_regimes_plot <- ((dual_regime_bicc & theme(axis.text = element_text(size = 12), axis.title.y = element_text(size = 15)) | 
                           dual_regime_distribution_plot & theme(axis.text.x = element_text(size = 12)) | 
                           dual_regime_EBE_plot & theme(axis.text = element_text(size = 12), aspect.ratio = 1)) + plot_layout(width = c(2.5, 1.5,1))) / 
  (non_rebound_fit + non_rebound_immune + plot_layout(ncol = 2) & 
     theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), strip.text = element_text(size = 15)))+
  plot_layout(ncol = 1, heights = c(5, 0.5)) +
  plot_annotation(tag_levels = "a") &
  theme(plot.tag = element_text(size = 20))

ggsave(here("paper_figures", "study1_regimes.pdf"), study1_regimes_plot, width = 15, height = 20)


## Assess identifiability ####
correlate_estimates(tabestimates, remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name)
compare_estimates(tabestimates, remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name)
estimate_boxplots(tabestimates, out_folder = out_folder, project_name = project_name)

## Look at low m regime only:
correlate_estimates(tabestimates %>% filter(run %in% regime1), remove_bicc_outliers = T,
                    out_folder = out_folder, project_name = project_name, file_addendum = "_low_m")
compare_estimates(tabestimates %>% filter(run %in% regime1), remove_bicc_outliers = T,
                  out_folder = out_folder, project_name = project_name, file_addendum = "_low_m")
estimate_boxplots(tabestimates %>% filter(run %in% regime1), 
                  out_folder = out_folder, project_name = project_name , file_addendum = "_low_m")

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

pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "parameter_plots_low_m.pdf"), width = 7, height = 8)
for(parameter in getPopulationParameterInformation()$name){
  if(parameter %in% names(fixed)) next
  print(parameter_plots(
    parameter, remove_bicc_outliers = F, remove_estimate_outliers = F,
    log10_tabestimates = NULL, tabestimates = tabestimates %>% filter(run %in% regime1),
    tabiters = tabiters %>% filter(run %in% regime1), phase = phase[regime1,], initialvalues = initialvalues  %>% filter(run %in% regime1)
  ))
}
dev.off()

## Plot individual fits
indivDF <- indivFit %>% mutate(id = as.character(id)) %>% left_join(indivParams)
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "indiv_trajectories_by_param.pdf"))
for(param in names(indivParams)[2:7]){
  print(
    indivDF %>% 
      ggplot(aes(time, indivPred_mode, color = (!!as.symbol(param)), group = interaction(id, run))) + 
      geom_point(data = . %>% filter(run == min(run)), alpha = 0.5, inherit.aes = F, aes(time, vl)) +
      geom_line(alpha = 0.5) + 
      geom_line(data = . %>% merge(distinct(tabestimates, run, bicc) %>% filter(bicc == min(bicc))),
                color = "black", lty = "dashed") +
      facet_wrap(~id) + 
      scale_color_viridis_c(direction = -1)
  )
}
dev.off()

## Plot population fits
pdf(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "population_trajectories_by_param.pdf"))
indivFit %>% 
  left_join(tabestimates %>% distinct(run, bicc)) %>%
  ggplot(aes(time, popPred, color = bicc, group = interaction(run))) + 
  geom_line(alpha = 0.5) + 
  # facet_wrap(~id) +
  scale_color_viridis_c(direction = -1) 

for(param in names(indivParams)[2:7]){
  if(param %in% names(fixed)) next
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

## Second regime for avererimic animals has unrealistically high P
table_predictions %>% 
  left_join(distinct(tabestimates, bicc, run)) %>% 
  ggplot(aes(time, P_mode, group = factor(run), color = run %in% regime1)) + 
  geom_line(alpha = 0.5) + 
  facet_wrap(~id) + 
  scale_y_log10() 

## Big RSE for m and Np and several random effects
tabestimates %>% filter(bicc == min(bicc))

## Big RSE for m and some random effects
tabestimates %>% filter(run %in% regime1) %>% filter(bicc == min(bicc))

## Save best estimates from regime 1 to use as initial conditions for SCM
best_estimates <- tabestimates %>% filter(run %in% regime1) %>% filter(bicc == min(bicc))
write_csv(best_estimates, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "best_estimates.csv"))
