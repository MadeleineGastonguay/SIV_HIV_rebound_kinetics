# Compare model fitting for eclipse phase model with approximated early exponential growth rate and the 
# delayed differential equation that assumes virus is proportional to infected cells

# Focusing on difference in timing between all three models

library(tidyverse)
library(here)
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# set up out folder

out_folder <- here("compare_delay_models_timing_perturb_ic_study2")
dir.create(out_folder)

# Use Study 2 data for all

data_file <- here("data", "derived", "STUDY2_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", 3))

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

fixed <- c("lambda_pop" = 70, "p_pop" = 0.7)

treatment_effects <- list("NPs" = c("vac", "intTLR7VAC"), "ta" = "vac")

NE <- 10000

# use study 2 data, start with ODE model:

create_project(
  here("model", "2regime_model_final.txt"), 
  data_file,
  headerTypes, error_model,
  exclude_random_effects = NULL, treatment_effects = treatment_effects, 
  normal_dist = NULL, 
  fixed= fixed)

saveProject(here(out_folder, "no_delay.mlxtran"))




igs <- read_csv(here('results', "study2", "IGS_model_comparison", "initial_values.csv"))

# Fit the model

no_delay_param <- no_delay_predictions <- no_delay_EBE <- no_delay_BIC <- NULL
no_delay_times <- numeric(20)

for(i in 1:20){
  
  # Set initial conditions
  igs %>% filter(run == i) %>% 
    rename(name = parameter) %>% 
    mutate(initialValue = signif(initialValue, 1)) %>% 
    select(-run) %>% 
    as.data.frame %>% 
    setPopulationParameterInformation()
  
  start_time <- Sys.time()
  
  runScenario()
  
  end_time <- Sys.time()
  no_delay_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  saveProject(here(out_folder, "no_delay.mlxtran"))
  
  no_delay_param <- rbind(
    no_delay_param,
    data.frame(estimate = getEstimatedPopulationParameters()) %>% 
      rownames_to_column("parameter") %>% 
      left_join(getEstimatedStandardErrors()$stochasticApproximation) %>% 
      mutate(run = i)
  )
  no_delay_predictions <- rbind(
    no_delay_predictions,
    read_csv(here(out_folder, "no_delay", "predictions.txt")) %>% mutate(run = i)
  )
  no_delay_EBE <- rbind(
    no_delay_EBE,
    getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
      mutate(I_art=4*10^(-5)/(ta*0.4),
             V_art=I_art*50000/23,
             P_art=ifelse(0.001>p*(1-0.9)*V_art/(V_art+NPs*1e5), m/(0.001-p*(1-0.9)*V_art/(V_art+NPs*1e5)), 1000*NE),
             E_art = p*0.9*P_art*V_art/(1*(V_art+NPs*1e5)),
             R0 = lambda*betas*1e-7*50000/((1+E_art/NE)*0.05*0.4*23),
             run = i)
  )
  no_delay_BIC <- rbind(
    no_delay_BIC,
    as.data.frame(t(getEstimatedLogLikelihood()$importanceSampling )) %>% mutate(run = i)
  )
  
}

saveProject(here(out_folder, "no_delay.mlxtran"))

save(no_delay_param, no_delay_predictions, no_delay_EBE, no_delay_BIC, no_delay_times,
     file = here(out_folder, "no_delay", "fitting_summaries.RData"))

print(no_delay_times)


## Now use eclipse phase model

create_project(
  here("model", "2regime_model_eclipse.txt"), 
  data_file,
  headerTypes, error_model,
  exclude_random_effects = NULL, treatment_effects = treatment_effects, 
  normal_dist = NULL, 
  fixed = fixed)

saveProject(here(out_folder, "eclipse.mlxtran"))
setPopulationParameterInformation(igs)

# Fit the model

eclipse_param <- eclipse_predictions <- eclipse_EBE <- eclipse_BIC <- NULL
eclipse_times <- numeric(20)

for(i in 1:20){
  # Set initial conditions
  igs %>% filter(run == i) %>% 
    rename(name = parameter) %>% 
    mutate(initialValue = signif(initialValue, 1)) %>% 
    select(-run) %>% 
    as.data.frame %>% 
    setPopulationParameterInformation()
  
  start_time <- Sys.time()
  
  runScenario()
  
  end_time <- Sys.time()
  eclipse_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  saveProject(here(out_folder, "eclipse.mlxtran"))
  
  eclipse_param <- rbind(
    eclipse_param,
    data.frame(estimate = getEstimatedPopulationParameters()) %>% 
      rownames_to_column("parameter") %>% 
      left_join(getEstimatedStandardErrors()$stochasticApproximation) %>% 
      mutate(run = i)
  )
  eclipse_predictions <- rbind(
    eclipse_predictions,
    read_csv(here(out_folder, "eclipse", "predictions.txt")) %>% mutate(run = i)
  )
  eclipse_EBE <- rbind(
    eclipse_EBE,
    getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
      mutate(I_art=4*10^(-5)/(ta*0.4),
             V_art=I_art*50000/23,
             P_art=ifelse(0.001>p*(1-0.9)*V_art/(V_art+NPs*1e5), m/(0.001-p*(1-0.9)*V_art/(V_art+NPs*1e5)), 1000*NE),
             E_art = p*0.9*P_art*V_art/(1*(V_art+NPs*1e5)),
             R0 = lambda*betas*1e-7*50000/((1+E_art/NE)*0.05*0.4*23)*(1/(1+0.05)),
             run = i)
  )
  eclipse_BIC <- rbind(
    eclipse_BIC,
    as.data.frame(t(getEstimatedLogLikelihood()$importanceSampling )) %>% mutate(run = i)
  )
  
}

print(eclipse_times)

save(eclipse_param, eclipse_predictions, eclipse_EBE, eclipse_BIC, eclipse_times,
     file = here(out_folder, "eclipse", "fitting_summaries.RData"))

## Now use DDE model

create_project(
  here("model", "2regime_model_delay.txt"), 
  data_file,
  headerTypes, error_model,
  exclude_random_effects = NULL, treatment_effects = treatment_effects, 
  normal_dist = NULL, 
  fixed = fixed)

saveProject(here(out_folder, "DDE.mlxtran"))
setPopulationParameterInformation(igs)

# Fit the model

DDE_param <- DDE_predictions <- DDE_EBE <- DDE_BIC <- NULL
DDE_times <- numeric(20)

for(i in 1:20){
  
  # Set initial conditions
  igs %>% filter(run == i) %>% 
    rename(name = parameter) %>% 
    mutate(initialValue = signif(initialValue, 1)) %>% 
    select(-run) %>% 
    as.data.frame %>% 
    setPopulationParameterInformation()
  
  start_time <- Sys.time()
  
  runScenario()
  
  end_time <- Sys.time()
  DDE_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
  saveProject(here(out_folder, "DDE.mlxtran"))
  
  DDE_param <- rbind(
    DDE_param,
    data.frame(estimate = getEstimatedPopulationParameters()) %>% 
      rownames_to_column("parameter") %>% 
      left_join(getEstimatedStandardErrors()$stochasticApproximation) %>% 
      mutate(run = i)
  )
  DDE_predictions <- rbind(
    DDE_predictions,
    read_csv(here(out_folder, "DDE", "predictions.txt")) %>% mutate(run = i)
  )
  DDE_EBE <- rbind(
    DDE_EBE,
    getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
      mutate(I_art=4*10^(-5)/(ta*0.4),
             V_art=I_art*50000/23,
             P_art=ifelse(0.001>p*(1-0.9)*V_art/(V_art+NPs*1e5), m/(0.001-p*(1-0.9)*V_art/(V_art+NPs*1e5)), 1000*NE),
             E_art = p*0.9*P_art*V_art/(1*(V_art+NPs*1e5)),
             R0 = lambda*betas*1e-7*50000/((1+E_art/NE)*0.05*0.4*23)*exp(-0.05*1),
             run = i)
  )
  DDE_BIC <- rbind(
    DDE_BIC,
    as.data.frame(t(getEstimatedLogLikelihood()$importanceSampling )) %>% mutate(run = i)
  )
  
}


saveProject(here(out_folder, "DDE.mlxtran"))

save(DDE_param, DDE_predictions, DDE_EBE, DDE_BIC, DDE_times,
     file = here(out_folder, "DDE", "fitting_summaries.RData"))

print(no_delay_times)
print(eclipse_times)
print(DDE_times)

mean(no_delay_times)
mean(eclipse_times)
mean(DDE_times)

sd(no_delay_times)
sd(eclipse_times)
sd(DDE_times)

plot(no_delay_times, DDE_times)

# On Average, DDE takes 1.8x longer, but it can be up to 2x
summary(DDE_times/no_delay_times)
summary(DDE_times/eclipse_times)

# On Average, there is no difference in fitting times when using the eclipse phase model compared to no delay
summary(eclipse_times/no_delay_times)


data.frame(
  no_delay = no_delay_times,
  eclipse = eclipse_times,
  DDE = DDE_times
) %>% 
  pivot_longer(everything()) %>% 
  ggplot() + 
  geom_boxplot(aes(value, name)) +
  labs(x = "Fitting times (seconds)", y = "Model")  + 
  
  data.frame(
    no_delay = no_delay_times,
    eclipse = eclipse_times,
    DDE = DDE_times
  ) %>% 
  ggplot() + 
  geom_boxplot(aes(DDE/no_delay, "DDE vs No Delay")) + 
  geom_boxplot(aes(eclipse/no_delay, "Eclipse vs No Delay")) +
  geom_boxplot(aes(DDE/eclipse, "DDE vs Eclipse")) +
  xlim(c(0.5, 2)) + 
  labs(x = "Ratio of fitting times", y = "Comparison") 





## Compare consistency of fits

rbind(
  no_delay_predictions %>% mutate(model = "No Delay"),
  eclipse_predictions %>% mutate(model = "Eclipse phase"),
  DDE_predictions %>% mutate(model = "DDE")
) %>% 
  left_join(distinct(data, id, tlr7, vac)) %>% 
  arrange(vac, tlr7) %>% 
  mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_mode, color = model)) + 
  geom_point(data = data, aes(y = vl, shape = interaction(tlr7, vac)), alpha = 0.25, color = "black") +
  geom_line(alpha = 0.5,aes(group = interaction(model, run))) + 
  facet_wrap(~id, nrow = 3) + 
  theme_minimal() + 
  labs(x = "Time (days)", y = "Viral load (log10 copies/ml)", 
       title = "Fitting to Study 2 Data", shape = "TLR7:VAC",
       subtitle = "using TLR7 effect on lambda and ta") + 
  coord_cartesian(xlim = c(0,20))


# Compare best fits from each model
plot_df <- merge(
  rbind(
    no_delay_predictions %>% mutate(model = "No Delay"),
    eclipse_predictions %>% mutate(model = "Eclipse phase"),
    DDE_predictions %>% mutate(model = "DDE")
  ) ,
  
  rbind(
    no_delay_BIC %>% mutate(model = "No Delay"),
    eclipse_BIC %>% mutate(model = "Eclipse phase"),
    DDE_BIC %>% mutate(model = "DDE")
  ) 
) %>% 
  group_by(model) %>% filter(BICc == min(BICc)) %>% 
  left_join(distinct(data, id, tlr7, vac)) 

plot_df %>% distinct(model, BICc, run)

plot_df %>% 
  arrange(vac, tlr7) %>% 
  mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_mode, color = model)) + 
  geom_point(data = data, aes(y = vl, shape = interaction(tlr7, vac)), alpha = 0.25, color = "black") +
  geom_line(alpha = 0.5, aes(group = interaction(model, run)), linewidth = 1) + 
  facet_wrap(~id, nrow = 3) + 
  theme_minimal() + 
  labs(x = "Time (days)", y = "Viral load (log10 copies/ml)", 
       title = "Fitting to Study 2 Data", shape = "TLR7:VAC",
       subtitle = "using TLR7 effect on lambda and ta") + 
  coord_cartesian(xlim = c(0,20))


## Compare population parameter values/ treatment effects:

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(grepl("pop", parameter), !parameter %in% names(fixed)) %>%  #%>% pull(rse) %>% range(na.rm = T)
  ggplot(aes(run, estimate, color = model)) + 
  geom_errorbar(data = . %>% mutate(rse = ifelse(rse > 100, NA, rse)), 
                aes(ymin = estimate*(100-rse)/100, ymax = estimate*(100+rse)/100), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(x = "Run", y = "Population parameter estimate") 

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(grepl("beta_", parameter)) %>%
  ggplot(aes(run, exp(estimate), color = model)) + 
  geom_errorbar(#data = . %>% mutate(rse = ifelse(rse > 100, NA, rse)),
                aes(ymin = exp(estimate*(100-rse)/100), ymax = exp(estimate*(100+rse)/100)), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  facet_wrap(~parameter, scales = "free_y") +
  labs(x = "Run", y = "Treatment effect Fold-change") +
  ggh4x::facetted_pos_scales(
    y = list(
      scale_y_continuous(limits = c(0, 0.25)),
      scale_y_continuous(limits = c(0,0.13)),
      scale_y_continuous(limits = c(0,25))
    )
  ) 




# Compare EBE parameters from all three
plot_df <- rbind(
  no_delay_EBE %>% mutate(model = "No Delay"),
  eclipse_EBE %>% mutate(model = "Eclipse phase"),
  DDE_EBE %>% mutate(model = "DDE")
) %>% 
  pivot_longer(!c(id, model, run)) %>% 
  pivot_wider(names_from = model, values_from = value) 


limits <- plot_df %>% 
  # filter(name %in% c("lambda", "betas", "p", "NPs", "ta", "m")) %>% 
  select(-id, -run) %>% group_by(name) %>% 
  mutate(min = min(across(where(is.numeric))),
         max = max(across(where(is.numeric)))) %>% 
  distinct(name, min, max) %>% 
  column_to_rownames("name") %>% 
  as.matrix()

plot_df %>% 
  filter(name %in% c("lambda", "betas", "p", "NPs", "ta", "m")) %>% 
  left_join(distinct(data, id, tlr7, vac) %>% mutate(id = as.character(id))) %>% 
  ggplot(aes(`DDE`, `No Delay`, color = "No Delay", shape = interaction(tlr7, vac))) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_point(aes(DDE, `Eclipse phase`, color = "Eclipse Phase"), alpha = 0.5, size = 2) +
  geom_abline() + 
  facet_wrap(~name, scales = "free") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = limits["betas",]),
      scale_x_continuous(limits = limits["lambda",]),
      scale_x_continuous(limits = limits["m",], transform = "log10"),
      scale_x_continuous(limits = limits["NPs",], transform = "log10"),
      scale_x_continuous(limits = limits["p",]),
      scale_x_continuous(limits = limits["ta",], transform = "log10")
    ),
    y = list(
      scale_y_continuous(limits = limits["betas",]),
      scale_y_continuous(limits = limits["lambda",]),
      scale_y_continuous(limits = limits["m",], transform = "log10"),
      scale_y_continuous(limits = limits["NPs",], transform = "log10"),
      scale_y_continuous(limits = limits["p",]),
      scale_y_continuous(limits = limits["ta",], transform = "log10")
    )
  ) +
  labs(x = "EBE from fitting DDE to data", y = "EBE from fitting other model to data",
       color = "Other model", shape = "TLR7:VAC", 
       title = "Compare individual parameter estimates for Study 2 animals")

plot_df %>% 
  filter(name == "R0", !run %in% exclude) %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  ggplot(aes(`DDE`, `No Delay`, color = "No Delay", shape = factor(tlr7))) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_point(aes(DDE, `Eclipse phase`, color = "Eclipse Phase"), alpha = 0.5, size = 2) +
  geom_abline() +
  facet_wrap(~name) +
  labs(x = "R0 from fitting DDE to data", y = "R0 from fitting other model to data",
       color = "Other model", shape = "TLR7") +
  xlim(limits["R0",]) + ylim(limits["R0",])

# Compare ta EBEs from two delay methods versus R0:
plot_df %>% 
  filter(name == "ta") %>% 
  left_join(distinct(data, id, tlr7, vac) %>% mutate(id = as.character(id))) %>% 
  left_join(plot_df %>% filter(name == "R0") %>% group_by(id, run) %>% summarise(R0 = DDE) ) %>% 
  ggplot(aes(`DDE`, `Eclipse phase`,  shape = interaction(tlr7, vac), color = R0)) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_abline() + 
  facet_wrap(vac~tlr7, labeller = "label_both") + 
  scale_x_log10(limits = c(0.05, 10)) + scale_y_log10(limits = c(0.05, 10)) + 
  labs(x = "ta EBE from DDE", y = "ta EBE from Eclipse Phase model",
       title = "Compare ta EBEs")

# Looks like ta is larger in the Eclipse phase model than the DDE model for animals not on vaccine

rbind(
  eclipse_predictions %>% mutate(model = "Eclipse phase"),
  DDE_predictions %>% mutate(model = "DDE")
) %>% 
  left_join(distinct(data, id, tlr7, vac)) %>% 
  filter(vac == 0) %>% 
  arrange(tlr7) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_mode, group = interaction(run, model), 
             shape = factor(tlr7))) + 
  geom_line(alpha = 0.25, aes(color = model)) + 
  geom_point(data = . %>% filter(run == 1, model == "DDE"), aes(y = vl)) +
  facet_wrap(~id) + 
  coord_cartesian(xlim = c(0,30)) + 
  labs(title = "Study 2 animals not on vaccine")


## Compare EBEs from best fit of all three models

# Best runs
best_runs <- rbind(
  no_delay_BIC %>% mutate(model = "No Delay"),
  eclipse_BIC %>% mutate(model = "Eclipse phase"),
  DDE_BIC %>% mutate(model = "DDE")
)  %>% filter(!run %in% exclude) %>% 
  group_by(model) %>% filter(BICc == min(BICc)) %>% summarise(run = run, BICc = BICc, OFV = OFV)

plot_df %>% 
  pivot_longer(!c(id, run, name), names_to = "model") %>% 
  filter(interaction(model, run) %in% interaction(best_runs$model, best_runs$run)) %>% 
  select(-run) %>% 
  pivot_wider(names_from = model) %>% 
  filter(name %in% c("lambda", "betas", "p", "NPs", "ta", "m")) %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  ggplot(aes(`DDE`, `No Delay`, color = "No Delay", shape = factor(tlr7))) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_point(aes(DDE, `Eclipse phase`, color = "Eclipse Phase"), alpha = 0.5, size = 2) +
  geom_abline() + 
  facet_wrap(~name, scales = "free") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(1, 8.3)),
      scale_x_continuous(limits = c(54, 783)),
      scale_x_continuous(limits = c(0.4, 3038), transform = "log10"),
      scale_x_continuous(limits = c(0.6, 108), transform = "log10"),
      scale_x_continuous(limits = c(1.5, 6.5)),
      scale_x_continuous(limits = c(0.01, 84.5), transform = "log10")
    ),
    y = list(
      scale_y_continuous(limits = c(1, 8.3)),
      scale_y_continuous(limits = c(54, 783)),
      scale_y_continuous(limits = c(0.4, 3038), transform = "log10"),
      scale_y_continuous(limits = c(0.6, 108), transform = "log10"),
      scale_y_continuous(limits = c(1.5, 6.5)),
      scale_y_continuous(limits = c(0.01, 84.5), transform = "log10")
    )
  ) +
  labs(x = "EBE from fitting DDE to data", y = "EBE from fitting other model to data",
       color = "Other model", shape = "TLR7", 
       title = "Compare individual parameter estimates for Study 2 animals")  


rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(interaction(model, run) %in% interaction(best_runs$model, best_runs$run)) %>% 
  select(-run) %>% 
  filter(grepl("pop", parameter), !parameter %in% names(fixed)) %>%  #%>% pull(rse) %>% range(na.rm = T)
  ggplot(aes(model, estimate, color = model)) + 
  geom_errorbar(data = . %>% mutate(rse = ifelse(rse > 500, NA, rse)), 
                aes(ymin = estimate*(100-rse)/100, ymax = estimate*(100+rse)/100), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(x = "Run", y = "Population parameter estimate") 

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(interaction(model, run) %in% interaction(best_runs$model, best_runs$run)) %>% 
  select(-run) %>% 
  filter(grepl("beta_", parameter)) %>%  #%>% pull(rse) %>% range(na.rm = T)
  ggplot(aes(model, exp(estimate), color = model)) + 
  geom_errorbar(aes(ymin = exp(estimate*(100-rse)/100), ymax = exp(estimate*(100+rse)/100)), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(x = "Run", y = "Population parameter estimate") 

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(interaction(model, run) %in% interaction(best_runs$model, best_runs$run)) %>% 
  select(-run) %>% 
  filter(grepl("omega_", parameter)) %>%  #%>% pull(rse) %>% range(na.rm = T)
  ggplot(aes(model, exp(estimate), color = model)) + 
  geom_errorbar(aes(ymin = exp(estimate*(100-rse)/100), ymax = exp(estimate*(100+rse)/100)), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(x = "Run", y = "Population parameter estimate") 
