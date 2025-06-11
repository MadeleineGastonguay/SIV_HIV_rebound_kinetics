# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study1"
project_name <- "SCM_removed_RE.mlxtran"
model_file <- here("model", "2regime_model_final.txt")

## Results folder
out_folder <- here("results", study)

## Data file
data_file <- here("data", "derived",str_interp("${toupper(study)}_DATA.txt"))
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", "catcov")

## Random effects to exclude
exclude_random_effects <- c("betas", "p") # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 200, "p_pop" = 3)

## Error model to use
# Either "constant", "proportional", "combined1", or "combined2". See ?setErrorModel() for details
error_model <- list("constant" )
names(error_model) <- names(data)[headerTypes == "observation"]

# Now run model building
## Create project
create_project(model_file, data_file, headerTypes, error_model,
               exclude_random_effects, treatment_effects = NULL, normal_dist = NULL, fixed = fixed)

# Use best estimates from IGS for initial conditions, to the 1st significant figure
# Use 1 for standard deviations unless the best estimate for it is bigger
igs <- read_csv(here("results", "study1_control", "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
  rowwise %>% 
  mutate(estimate = ifelse(grepl("pop", parameter), estimate, max(estimate, 1)),
         estimate = signif(estimate, 1)) %>% 
  select(parameter, estimate)

set_base(igs, fixed = fixed)

saveProject(here(out_folder, project_name))

# Run SCM
relationships <- NULL
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", linearization = FALSE,
                         relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

# FIM condiiton number is 12
# TLR7 effect on Np is not significant
# BICc is 416
# There is no TLR7 effect on ta!

# Try locking ta ####
project_name <- "SCM_removed_RE_locked_ta.mlxtran"
saveProject(here(out_folder, project_name))
## Lock in tlr7 on ta:
relationships <- data.frame(parameters = "ta", covariates = "tlr7", locked = TRUE)
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", linearization = FALSE,
                         relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

# Tlr7 effect is small and insignificant
# BICc is 427

## AM1: TLR7 on lambda and ta

am1 <- list("lambda" = "tlr7",
            "ta" = "tlr7")

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1", 
        model_file, data_file, 
        igs, am1, 
        headerTypes, error_model,
        fixed = fixed,
        exclude_random_effects = exclude_random_effects)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 381.9511553   403.9511553   415.4409021   429.8662356     5.0000000     0.3266749 
