# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM_tlr7_exclude_full_censored.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "derived", "VACABTLR7_DATA_exclude_full_censored.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", 3), rep("ignore", 3))

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 100, "p_pop" = 0.9)

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

if(!file.exists(here(out_folder, project_name))) saveProject(here(out_folder, project_name))

# Run SCM
relationships <- NULL # nothing locked in
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)

loadProject(here(out_folder, project_name))

# Best BICc: 859.826	

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 10, height = 7)

write_csv(out$summary, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_summary.csv"))


## Try Alternative Models:

# Model with Ab + Vac on Np instead of TLR7 + Ab
am1 <- list("lambda" = "ab",
            "betas" = "tlr7",
            "ta" = "vac",
            "NPs" = c("ab", "vac"),
            "m" = "ab"
            )

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
        model_file, data_file,
        igs, am1,
        headerTypes, error_model,
        fixed = fixed)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 787.0834335   821.0834335   850.2141586   863.0389054     5.0000000     0.3693764 


# Model with TLR7 on lambda instead of beta
am2 <- list("lambda" = c("ab", "tlr7"),
            "ta" = "vac",
            "NPs" = c("tlr7", "vac"),
            "m" = "ab"
)

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
        model_file, data_file,
        igs, am2,
        headerTypes, error_model,
        fixed = fixed)
# 
# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 810.0052999   844.0052999   873.1360250   885.9607718     5.0000000     0.3601058 

