# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study1"
project_name <- "SCM.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "derived", str_interp("${toupper(study)}_DATA.txt"))
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

# Now run model building
## Create project
create_project(model_file, data_file, headerTypes, error_model,
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL, fixed = fixed)

# Use best estimates from IGS for initial conditions, to the 1st significant figure
# Use 1 for standard deviations (if it's less than this, it's difficult for SAEM to jump around, and if it's too big than it's difficult to converge)
igs <- read_csv(here(out_folder, "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
  rowwise %>% 
  mutate(estimate = ifelse(grepl("pop", parameter), estimate, 1),
         estimate = signif(estimate, 1)) %>% 
  select(parameter, estimate)

set_base(igs, fixed = fixed)

# saveProject(here(out_folder, project_name))

# Run SCM
relationships <- NULL # nothing locked in
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

# Best model: 421.1414, but it does not make sense to have the effect on beta and tlr7 on Np is insignificant


# Try alternative models ####

# tlr7 on lambda instead of beta

am1 <- list("lambda" = "tlr7",
            "ta" = "tlr7",
            "NPs" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1", 
  model_file, data_file, 
  igs, am1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 381.2489972   409.2489972   423.8723113   438.2976449     5.0000000     0.3246117 

# No treatment effects are signficant
# FIM condition number = 160

# Remove TLR7 on Np since it is not significant (use lambda instead of beta)

am2 <- list("lambda" = "tlr7",
            "ta" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 380.791757    406.791757    420.370549    434.795882      5.000000      0.285326 

# FIM condition number is 160 (RSE is big for m_pop, omega_beta, and omega_p)
# Both treatment effects are signficant

# Remove TLR7 on Np since it is not significant (keep tlr7 on beta)

am3 <- list("betas" = "tlr7",
            "ta" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3", 
  model_file, data_file, 
  igs, am3, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 380.4225462   406.4225462   420.0013379   434.4266714     5.0000000     0.3352546 

# Tlr7 effect on beta is not signficant
# Could not estimate eigenvalues of FIM


# Try alternative models without RE on beta and p  ####

exclude_random_effects <- c("betas", "p")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1_removeRE", 
  model_file, data_file, 
  igs, am1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 380.2541787   404.2541787   416.7884480   431.2137815     5.0000000     0.3480497 

# FIM condition number = 14
# TLR7 effect on Np is not significant

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2_removeRE", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 381.9511553   403.9511553   415.4409021   429.8662356     5.0000000     0.3266749 

# FIM condition number = 19
# Treatment effects remain significant
# RSE is reduced
# Very similar OFV to same model with random effects on beta and p

## Try RM models ####

# In SHIV study, we estimated an Ab effect on beta and Np and a tlr7 effect on m

rm1 <- list("m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1", 
  model_file, data_file, 
  igs, rm1, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 392.6126431   416.6126431   429.1469124   443.5722459     5.0000000     0.4632413 

# effect on m is not significant


rm2 <- list("m" = "tlr7",
            "ta" = "tlr7",
            "lambda" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2", 
  model_file, data_file, 
  igs, rm2, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 379.8120784   407.8120784   422.4353925   436.8607260     5.0000000     0.3280159 

# TLr7 effect on m is insignificant