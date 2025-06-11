# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "human"
project_name <- "SCM.mlxtran"
model_file <- here("model", "2regime_model_HIV_final_twash.txt")

## Results folder
out_folder <- here("results", study)

## Data file
data_file <- here("data", "derived", "HUMAN_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "ignore", "observation", "cens", "catcov")

## Random effects to exclude
exclude_random_effects <- "tw" # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 200, "p_pop" = 3, "tw_pop" = 3)

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
igs <- read_csv(here("results", "human", "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
  rowwise %>% 
  mutate(estimate = ifelse(grepl("pop", parameter), estimate, 1),
         estimate = signif(estimate, 1)) %>% 
  select(parameter, estimate)

set_base(igs, fixed = fixed)

saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run SCM
relationships <- NULL
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", linearization = FALSE,
                         relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 3)

# BICc: 1428.672	(< 2 from BICc of model without covariates)
# Effect on ta is insignificant 

## Test alternative models

# Previously, an effect of nnrti was inferred in twash
am1 <- list("tw" = "nnrti")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1", 
  model_file, data_file, 
  igs, am1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = "tw"
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 1366.2322585  1390.2322585  1417.0415366  1431.8265059     5.0000000     0.5312495 

# The covariate is signficant

# No covariates
am2 <- list()

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = "tw"
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 1370.7264071  1392.7264071  1417.3015787  1429.6223864     5.0000000     0.5011235 



