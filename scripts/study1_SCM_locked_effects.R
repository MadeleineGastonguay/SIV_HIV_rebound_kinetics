# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study1"
project_name <- "SCM_ta_locked.mlxtran"
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
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

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
## Lock in tlr7 on ta:
relationships <- data.frame(parameters = "ta", covariates = "tlr7", locked = TRUE)
# relationships <- data.frame(parameters = "lambda", covariates = "tlr7", locked = TRUE)
# relationships <- NULL
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", linearization = FALSE,
                         relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)


# tlr7 on Np is not significant (can't check significance for beta because theres no RE)
# Condition number = 100
# BICc = 415.42


# Try alternative models ####

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
# 381.2242979   405.2242979   417.7585672   432.1839007     5.0000000     0.3193926 

# Only tlr7 effect on lambda is significant
# Condition number = 110



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
# 381.2356319   403.2356319   414.7253788   429.1507123     5.0000000     0.2691841 


am3 <- list("NPs" = "tlr7",
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
# 375.3149388   397.3149388   408.8046856   423.2300191     5.0000000     0.2024884 

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

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2_removeRE", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3_removeRE", 
  model_file, data_file, 
  igs, am3, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)