# Compare model fitting for eclipse phase model with approximated early exponential growth rate and the 
# delayed differential equation that assumes virus is proportional to infected cells

# Focusing on difference in timing between all three models

library(tidyverse)
library(here)
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# set up out folder

out_folder <- here("compare_delay_models_timing_perturb_ic")
dir.create(out_folder)

# Use Study 1 data for all

data_file <- here("data", "derived", "STUDY1_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", "catcov")

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

fixed <- c("lambda_pop" = 200, "p_pop" = 3)

treatment_effects <- list("lambda" = "tlr7", "ta" = "tlr7")

NE <- 10000

# use study 1 data, start with ODE model:

create_project(
  here("model", "2regime_model_final.txt"), 
  data_file,
  headerTypes, error_model,
  exclude_random_effects = NULL, treatment_effects = treatment_effects, 
  normal_dist = NULL, 
  fixed= fixed)

saveProject(here(out_folder, "no_delay.mlxtran"))




igs <- read_csv(here('results', "study1", "IGS_model_comparison", "initial_values.csv"))

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

save(eclipse_param, eclipse_predictions, eclipse_EBE, eclipse_BIC, eclipse_times,
     file = here(out_folder, "eclipse", "fitting_summaries.RData"))

print(eclipse_times)

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

mean(no_delay_times[-c(1:2)])
mean(eclipse_times[-c(1:2)])
mean(DDE_times[-c(1:2)])

sd(no_delay_times[-c(1:2)])
sd(eclipse_times[-c(1:2)])
sd(DDE_times[-c(1:2)])

plot(no_delay_times[-c(1:2)], DDE_times[-c(1:2)])



# On Average, DDE takes 1.4 - 1.5x longer, but it can be up to 2x
summary(DDE_times[-c(1:2)]/no_delay_times[-c(1:2)])
summary(DDE_times[-c(1:2)]/eclipse_times[-c(1:2)])

# On Average, there is no difference in fitting times when using the eclipse phase model compared to no delay
summary(eclipse_times[-c(1:2)]/no_delay_times[-c(1:2)])


data.frame(
  no_delay = no_delay_times[-c(1:2)],
  eclipse = eclipse_times[-c(1:2)],
  DDE = DDE_times[-c(1:2)]
) %>% 
  pivot_longer(everything()) %>% 
  ggplot() + 
  geom_boxplot(aes(value, name)) +
  labs(x = "Fitting times (seconds)", y = "Model")  + 
  
data.frame(
  no_delay = no_delay_times[-c(1:2)],
  eclipse = eclipse_times[-c(1:2)],
  DDE = DDE_times[-c(1:2)]
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
  left_join(distinct(data, id, tlr7)) %>% 
  arrange(tlr7) %>% 
  mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_mode, color = model)) + 
  geom_point(data = data, aes(y = vl, shape = as.factor(tlr7)), alpha = 0.25, color = "black") +
  geom_line(alpha = 0.5,aes(group = interaction(model, run))) + 
  facet_wrap(~id, nrow = 3) + 
  theme_minimal() + 
  labs(x = "Time (days)", y = "Viral load (log10 copies/ml)", 
       title = "Fitting to Study 1 Data", shape = "TLR7",
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
  left_join(distinct(data, id, tlr7)) 

plot_df %>% distinct(model, BICc, run)

plot_df %>% 
  arrange(tlr7) %>% 
  mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_mode, color = model)) + 
  geom_point(data = data, aes(y = vl, shape = as.factor(tlr7)), alpha = 0.25, color = "black") +
  geom_line(alpha = 0.5, aes(group = interaction(model, run)), linewidth = 1) + 
  facet_wrap(~id, nrow = 3) + 
  theme_minimal() + 
  labs(x = "Time (days)", y = "Viral load (log10 copies/ml)", 
       title = "Fitting to Study 1 Data", shape = "TLR7",
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
  filter(grepl("beta_", parameter)) %>%
  ggplot(aes(run, exp(estimate), color = model)) + 
  geom_errorbar(data = . %>% mutate(rse = ifelse(rse > 100, NA, rse)),
                aes(ymin = exp(estimate*(100-rse)/100), ymax = exp(estimate*(100+rse)/100)), width = 0.2, position = position_dodge(width = 0.5)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) + 
  facet_wrap(~parameter, scales = "free_y") +
  labs(x = "Run", y = "Treatment effect Fold-change") 




# Compare EBE parameters from all three
plot_df <- rbind(
  no_delay_EBE %>% mutate(model = "No Delay"),
  eclipse_EBE %>% mutate(model = "Eclipse phase"),
  DDE_EBE %>% mutate(model = "DDE")
) %>% 
  pivot_longer(!c(id, model, run)) %>% 
  pivot_wider(names_from = model, values_from = value) 

exclude <- plot_df %>% filter(name == "m", `No Delay` > 50000) %>% pull(run) %>% unique

plot_df %>% 
  filter(name %in% c("lambda", "betas", "p", "NPs", "ta", "m"), ! run %in% exclude) %>% 
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
       title = "Compare individual parameter estimates for Study 1 animals")

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
  xlim(c(0,43)) + ylim(c(0,43))

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
       title = "Compare individual parameter estimates for Study 1 animals")  


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
