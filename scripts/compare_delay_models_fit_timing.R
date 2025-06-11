# Compare model fitting for eclipse phase model with approximated early exponential growth rate and the 
# delayed differential equation that assumes virus is proportional to infected cells

# Focusing on difference in timing between all three models

library(tidyverse)
library(here)
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# set up out folder

out_folder <- here("compare_delay_models_timing")
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

igs <- getPopulationParameterInformation()
igs[3, 2] <- 2
igs[4, 2] <- 0.4
igs[6, 2] <- 8
igs[8, 2] <- 50


setPopulationParameterInformation(igs)

# Fit the model

no_delay_param <- no_delay_predictions <- no_delay_EBE <- no_delay_BIC <- NULL
no_delay_times <- numeric(20)

for(i in 1:20){
  start_time <- Sys.time()
  
  runScenario()
  
  end_time <- Sys.time()
  no_delay_times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  
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

print(no_delay_times)
print(eclipse_times)
print(DDE_times)

mean(no_delay_times)
mean(eclipse_times)
mean(DDE_times)

sd(no_delay_times)
sd(eclipse_times)
sd(DDE_times)


# On Average, DDE takes 1.5x longer but it can be up to 4 times
summary(DDE_times/no_delay_times)
summary(DDE_times/eclipse_times)

# No difference in fitting time between eclipse phase and no delay model
summary(eclipse_times/no_delay_times)

