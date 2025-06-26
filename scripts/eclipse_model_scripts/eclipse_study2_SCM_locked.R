# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study2"
project_name <- "SCM_locked.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "derived", "STUDY2_DATA.txt")
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

saveProject(here(out_folder, project_name))

# Run SCM
relationships <- data.frame(
  parameters = c("NPs", "ta", "lambda"), 
  covariates = c("vac"),
  locked = TRUE
)
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

# BICc: 513, FIM: 2500, ta effects insiginificant, vaccine reduces p

# If you remove TLR7*VAC effect on ta (AM3), the BICc is 516 but many effects are not significant and condition number is big

# If you remove vac effect on p (AM6), the BICc is 519 but many effects are not significant and condition number is big


### Compare model performance between model with and without p effect when p is fixed to a lower value:
fixed2 <- c("lambda_pop" = 100, "p_pop" = 0.7)

SCM_lowerp <- list("lambda" = "vac",
                   "ta" = c("vac", "intTLR7VAC"),
                   "NPs" = c("vac", "intTLR7VAC"),
                   "p" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_p",
  model_file, data_file,
  igs, SCM_lowerp,
  headerTypes, error_model,
  fixed = fixed2
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 442.3730372   476.3730372   502.3211661   515.1459129     5.0000000     0.2399918 

# The effect on p is now insignificant ... fixing it to 1 is slightly too big so we are inferring a treatment effect on it


SCM_lowerp2 <- list(
                   "ta" = c("vac", "intTLR7VAC"),
                   "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_p2",
  model_file, data_file,
  igs, SCM_lowerp2,
  headerTypes, error_model,
  fixed = fixed2
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 455.1095862   485.1095862   508.0049940   520.8297408     5.0000000     0.2375065 

SCM_lowerp3 <- list( "lambda" = "vac",
  "ta" = c("vac", "intTLR7VAC"),
  "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_p3",
  model_file, data_file,
  igs, SCM_lowerp3,
  headerTypes, error_model,
  fixed = fixed2
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 444.6910214   476.6910214   501.1127897   513.9375365     5.0000000     0.2283861 

# everything is significant, FIM condition number is 360


SCM_lowerp4 <- list( "lambda" = "vac",
                     "ta" = c("vac"),
                     "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_p4",
  model_file, data_file,
  igs, SCM_lowerp4,
  headerTypes, error_model,
  fixed = fixed2
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 448.7547422   478.7547422   501.6501501   514.4748969     5.0000000     0.3382489 

# FIM 230, lambda not significant

SCM_lowerp5 <- list( 
                     "ta" = c("vac"),
                     "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_p5",
  model_file, data_file,
  igs, SCM_lowerp5,
  headerTypes, error_model,
  fixed = fixed2
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 459.9355524   487.9355524   509.3045997   522.1293465     5.0000000     0.2970558 

SCM_lowerp_lowerl <- list( "ta" = c("vac"),
  "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_lowerl",
  model_file, data_file,
  igs, SCM_lowerp_lowerl,
  headerTypes, error_model,
  fixed = c("lambda_pop" = 70, "p_pop" = 0.7)
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 457.1537044   485.1537044   506.5227518   519.3474986     5.0000000     0.2395731 

SCM_lowerp_lowerl2 <- list( "ta" = c("vac"),
                           "NPs" = c("vac", "intTLR7VAC"),
                           "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "SCM_lower_lowerl2",
  model_file, data_file,
  igs, SCM_lowerp_lowerl2,
  headerTypes, error_model,
  fixed = c("lambda_pop" = 70, "p_pop" = 0.7)
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 448.4302230   478.4302230   501.3256309   514.1503777     5.0000000     0.2041456
