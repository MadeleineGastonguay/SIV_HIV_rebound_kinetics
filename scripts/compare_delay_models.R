# Compare model fitting for eclipse phase model with approximated early exponential growth rate and the 
# delayed differential equation that assumes virus is proportional to infected cells

library(tidyverse)
library(here)
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# set up out folder

out_folder <- here("compare_delay_models")
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

# use study 1 data, start with ODE model:

create_project(
  here("model", "2regime_model_final.txt"), 
  data_file,
  headerTypes, error_model,
  exclude_random_effects = NULL, treatment_effects = treatment_effects, 
  normal_dist = NULL, 
  fixed= fixed)

saveProject(here(out_folder, "no_delay.mlxtran"))

igs <- getPopulationParameterInformation()
igs[3, 2] <- 2
igs[4, 2] <- 0.4
igs[6, 2] <- 8
igs[8, 2] <- 50


setPopulationParameterInformation(igs)

# Fit the model
runScenario()

NE <- 10000


no_delay_param <- data.frame(estimate = getEstimatedPopulationParameters()) %>% 
  rownames_to_column("parameter") %>% 
  left_join(getEstimatedStandardErrors()$stochasticApproximation)
no_delay_predictions <- read_csv(here(out_folder, "no_delay", "predictions.txt"))
no_delay_EBE <- getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
  mutate(I_art=4*10^(-5)/(ta*0.4),
         V_art=I_art*50000/23,
         P_art=ifelse(0.001>p*(1-0.9)*V_art/(V_art+NPs*1e5), m/(0.001-p*(1-0.9)*V_art/(V_art+NPs*1e5)), 1000*NE),
         E_art = p*0.9*P_art*V_art/(1*(V_art+NPs*1e5)),
         R0 = lambda*betas*1e-7*50000/((1+E_art/NE)*0.05*0.4*23))
no_delay_BIC <- getEstimatedLogLikelihood()$importanceSampling

saveProject(here(out_folder, "no_delay.mlxtran"))


## Now use eclipse phase model

create_project(
  here("model", "2regime_model_eclipse.txt"), 
  data_file,
  headerTypes, error_model,
  exclude_random_effects = NULL, treatment_effects = treatment_effects, 
  normal_dist = NULL, 
  fixed = fixed)

saveProject(here(out_folder, "eclipse_phase.mlxtran"))
setPopulationParameterInformation(igs)

# Fit the model
runScenario()

eclipse_param <- data.frame(estimate = getEstimatedPopulationParameters()) %>% 
  rownames_to_column("parameter") %>% 
  left_join(getEstimatedStandardErrors()$stochasticApproximation)
eclipse_predictions <- read_csv(here(out_folder, "eclipse_phase", "predictions.txt"))
eclipse_EBE <- getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
  mutate(I_art=4*10^(-5)/(ta*0.4),
         V_art=I_art*50000/23,
         P_art=ifelse(0.001>p*(1-0.9)*V_art/(V_art+NPs*1e5), m/(0.001-p*(1-0.9)*V_art/(V_art+NPs*1e5)), 1000*NE),
         E_art = p*0.9*P_art*V_art/(1*(V_art+NPs*1e5)),
         R0 = lambda*betas*1e-7*50000/((1+E_art/NE)*0.05*0.4*23)*(1/(1+0.05)))
eclipse_BIC <- getEstimatedLogLikelihood()$importanceSampling

saveProject(here(out_folder, "eclipse_phase.mlxtran"))

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
runScenario()

DDE_param <- data.frame(estimate = getEstimatedPopulationParameters()) %>% 
  rownames_to_column("parameter") %>% 
  left_join(getEstimatedStandardErrors()$stochasticApproximation)
DDE_predictions <- read_csv(here(out_folder, "DDE", "predictions.txt"))
DDE_EBE <- getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
  mutate(I_art=4*10^(-5)/(ta*0.4),
         V_art=I_art*50000/23,
         P_art=ifelse(0.001>p*(1-0.9)*V_art/(V_art+NPs*1e5), m/(0.001-p*(1-0.9)*V_art/(V_art+NPs*1e5)), 1000*NE),
         E_art = p*0.9*P_art*V_art/(1*(V_art+NPs*1e5)),
         R0 = lambda*betas*1e-7*50000/((1+E_art/NE)*0.05*0.4*23)*exp(-0.05*1))
DDE_BIC <- getEstimatedLogLikelihood()$importanceSampling

saveProject(here(out_folder, "DDE.mlxtran"))


## Compare model predictions:

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
  geom_line(alpha = 0.5, linewidth = 1) + 
  facet_wrap(~id, nrow = 3) + 
  theme_minimal() + 
  labs(x = "Time (days)", y = "Viral load (log10 copies/ml)", 
       title = "Fitting to Study 1 Data", shape = "TLR7",
       subtitle = "using TLR7 effect on lambda and ta")

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
  geom_line(alpha = 0.5, linewidth = 1) + 
  facet_wrap(~id, nrow = 3) + 
  theme_minimal() + 
  labs(x = "Time (days)", y = "Viral load (log10 copies/ml)", 
       title = "Fitting to Study 1 Data", shape = "TLR7",
       subtitle = "using TLR7 effect on lambda and ta") + 
  coord_cartesian(xlim = c(0,20)) 

rbind(
  no_delay_predictions %>% mutate(model = "No Delay"),
  eclipse_predictions %>% mutate(model = "Eclipse phase"),
  DDE_predictions %>% mutate(model = "DDE")
) %>% 
  select(model, id, time, indivPred_mode) %>% 
  pivot_wider(names_from = model, values_from = indivPred_mode) %>% 
  ggplot(aes(DDE, `Eclipse phase`, color = "Eclipse")) + 
  geom_point(alpha = 0.5) + 
  geom_point(aes(DDE, `No Delay`, color = "No Delay"), alpha = 0.5) +
  geom_abline() + 
  labs(x = "Viral load predicted from DDE", y = "Viral load predicted from other model",
       color = "Other model", title = "Comparing viral load predictions bewteen models")

# Compare EBE parameters from all three
plot_df <- rbind(
  no_delay_EBE %>% mutate(model = "No Delay"),
  eclipse_EBE %>% mutate(model = "Eclipse phase"),
  DDE_EBE %>% mutate(model = "DDE")
) %>% 
  pivot_longer(!c(id, model)) %>% 
  pivot_wider(names_from = model, values_from = value) 

plot_df %>% 
  filter(name %in% c("lambda", "betas", "p", "NPs", "ta", "m")) %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  ggplot(aes(`DDE`, `No Delay`, color = "No Delay", shape = factor(tlr7))) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_point(aes(DDE, `Eclipse phase`, color = "Eclipse Phase"), alpha = 0.5, size = 2) +
  geom_abline() + 
  facet_wrap(~name, scales = "free") +
  ggh4x::facetted_pos_scales(
    x = list(
      scale_x_continuous(limits = c(1.2, 6)),
      scale_x_continuous(limits = c(60, 700)),
      scale_x_continuous(limits = c(0.4, 1600), transform = "log10"),
      scale_x_continuous(limits = c(0.6, 103), transform = "log10"),
      scale_x_continuous(limits = c(2.5, 3.6)),
      scale_x_continuous(limits = c(0.01, 82), transform = "log10")
    ),
    y = list(
      scale_y_continuous(limits = c(1.2, 6)),
      scale_y_continuous(limits = c(60, 700)),
      scale_y_continuous(limits = c(0.4, 1600), transform = "log10"),
      scale_y_continuous(limits = c(0.6, 103), transform = "log10"),
      scale_y_continuous(limits = c(2.5, 3.6)),
      scale_y_continuous(limits = c(0.01, 82), transform = "log10")
    )
  ) +
  labs(x = "EBE from fitting DDE to data", y = "EBE from fitting other model to data",
       color = "Other model", shape = "TLR7", 
       title = "Compare individual parameter estimates for Study 1 animals")

plot_df %>% 
  filter(name == "R0") %>% 
  left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  ggplot(aes(`DDE`, `No Delay`, color = "No Delay", shape = factor(tlr7))) + 
  geom_point(alpha = 0.5, size = 2) + 
  geom_point(aes(DDE, `Eclipse phase`, color = "Eclipse Phase"), alpha = 0.5, size = 2) +
  geom_abline() +
  facet_wrap(~name) +
  labs(x = "R0 from fitting DDE to data", y = "R0 from fitting other model to data",
       color = "Other model", shape = "TLR7") + 
  xlim(c(0,40)) + ylim(c(0,40))

# Check that the differences in R0 is due to the difference in lambda*beta times the factor of 0.95
plot_df %>% 
  group_by(id, run) %>% 
  filter(name %in% c("lambda", "betas")) %>% 
  summarise(name = "lambda*beta", 
            `No Delay` = prod(`No Delay`),
            `Eclipse phase` = prod(`Eclipse phase`),
            DDE = prod(DDE)) %>%
  pivot_longer(!c(id, run, name), names_to = "model") %>% 
  rbind(plot_df %>% filter(name == "R0") %>% pivot_longer(!c(id, run, name), names_to = "model") ) %>% 
  pivot_wider(names_from = model) %>% 
  group_by(id, run, name) %>% 
  summarise(DDE_vs_no_delay = DDE/`No Delay`,
         eclipse_vs_no_delay = `Eclipse phase`/`No Delay`) %>% 
  pivot_longer(!c(id, run, name), names_to = "comparison") %>%
  pivot_wider() %>% 
  ggplot(aes(`lambda*beta`*0.95, R0, color = comparison)) + geom_point() + geom_abline()

## Compare population parameter values/ treatment effects:

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(grepl("pop", parameter) | parameter == "a") %>%
  ggplot(aes(model, estimate, color = model)) + 
  geom_errorbar(aes(ymin = estimate*(100-rse)/100, ymax = estimate*(100+rse)/100), width = 0.2) +
  geom_point(size = 2) + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) +
  labs(x = "Model", y = "Population parameter estimate") + 
  theme(legend.position = "none")

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(grepl("beta_", parameter)) %>%
  ggplot(aes(model, exp(estimate), color = model)) + 
  geom_errorbar(aes(ymin = exp(estimate*(100-rse)/100), ymax = exp(estimate*(100+rse)/100)), width = 0.2) +
  geom_point(size = 2) + 
  facet_wrap(~parameter, scales = "free_y") +
  labs(x = "Model", y = "Treatment effect Fold-change") + 
  theme(legend.position = "none")

rbind(
  no_delay_param %>% mutate(model = "No Delay"),
  eclipse_param %>% mutate(model = "Eclipse phase"),
  DDE_param %>% mutate(model = "DDE")
) %>% 
  filter(grepl("omega_", parameter)) %>%
  ggplot(aes(model, exp(estimate), color = model)) + 
  geom_errorbar(aes(ymin = exp(estimate*(100-rse)/100), ymax = exp(estimate*(100+rse)/100)), width = 0.2) +
  geom_point(size = 2) + 
  facet_wrap(~parameter, scales = "free_y") +
  labs(x = "Model", y = "Random effect Fold-change") + 
  theme(legend.position = "none")
