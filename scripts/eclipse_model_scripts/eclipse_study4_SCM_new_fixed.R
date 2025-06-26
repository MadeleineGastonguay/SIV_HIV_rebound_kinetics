# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM_new_fixed.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("ignore", 3), rep("catcov", 3))

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
       out$plot, width = 10, height = 7)

# BICc: 903.295, FIM 47, all significant

# The fourth covariate does not significantly improve BICc
out$summary %>% #group_by(n_covar) %>% filter(BICc == min(BICc)) %>% 
  mutate(d = threshold - BICc) %>% filter(d > 0, d < 2)

# AM1: up to model building step 4
am1 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
  model_file, data_file,
  igs, am1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 842.9299727   872.9299727   901.3072772   914.1320240     5.0000000     0.5482523 


# AM2: up to model building step 4 with the alternative treatment effect
am2 <- list("lambda" = "intTLR7AB",
            "p" = "intTLR7AB",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
  model_file, data_file,
  igs, am2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 844.7532917   874.7532917   903.1305962   915.9553430     5.0000000     0.4909552 

# All significant, but AB reduces p which does not make sense

# AM3: SCM before step 4
am3 <- list("lambda" = "intTLR7AB",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
  model_file, data_file,
  igs, am3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 848.6952373   876.6952373   903.1807215   916.0054683     5.0000000     0.5216833 

# AM4: final model with lambda and beta effects swapped
am4 <- list("betas" = "intTLR7AB",
            "lambda" = "intTLR7VAC",
            "ta" = "intTLR7VAC",
            "p" = c("intTLR7AB", "intTLR7VAC"),
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
  model_file, data_file,
  igs, am4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 850.9595749   886.9595749   921.0123403   933.8370871     5.0000000     0.4663922 


# AM5: SCM before step 4 with AB effect on beta instead of lambda
am5 <- list("betas" = "intTLR7AB",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5",
  model_file, data_file,
  igs, am5,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 865.5001906   893.5001906   919.9856748   932.8104216     5.0000000     0.5969083 

# AM6: SCM including step 4 with AB effect on beta instead of lambda
am6 <- list("betas" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6",
  model_file, data_file,
  igs, am6,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 860.7110928   890.7110928   919.0883973   931.9131441     5.0000000     0.5048768 

