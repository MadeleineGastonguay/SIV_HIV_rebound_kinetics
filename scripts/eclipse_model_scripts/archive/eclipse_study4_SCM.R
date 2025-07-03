# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM.mlxtran"
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
relationships <- NULL # nothing locked in
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 10, height = 7)

# BICc: 905 FIM 450

# AM1: 
am1 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "betas" = "intTLR7VAC",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
  model_file, data_file,
  igs, am1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 835.6213095   867.6213095   897.8904343   910.7151810     5.0000000     0.6063784 

am1.1 <- list("ta" = "intTLR7VAC",
            "betas" = c("intTLR7VAC", "intTLR7AB"),
            "NPs" = c("intTLR7AB", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1.1",
  model_file, data_file,
  igs, am1.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 844.1224755   876.1224755   906.3916002   919.2163470     5.0000000     0.4377525 


# AM2 
am2 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "betas" = "intTLR7VAC",
            "NPs" = c("intTLR7AB"),
            "p" = "intTLR7AB")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
  model_file, data_file,
  igs, am2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 843.3082511   875.3082511   905.5773758   918.4021226     5.0000000     0.4154383 

# AM3 
am3 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "p" = "intTLR7AB",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
  model_file, data_file,
  igs, am3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 839.2223070   871.2223070   901.4914318   914.3161786     5.0000000     0.3508755 

# AM3 without p:
am3.1 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3.1",
  model_file, data_file,
  igs, am3.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 845.2160090   875.2160090   903.5933135   916.4180603     5.0000000     0.5293002 

# AM3 without p, swap TLR7*Ab on beta instead of lambda:
am3.2 <- list("betas" = "intTLR7AB",
              "ta" = "intTLR7VAC",
              "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3.2",
  model_file, data_file,
  igs, am3.2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 858.9694768   888.9694768   917.3467813   930.1715281     5.0000000     0.6446101 

# AM4
am4 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "m" = "intTLR7AB",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
  model_file, data_file,
  igs, am4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 838.1677216   870.1677216   900.4368463   913.2615931     5.0000000     0.4079118

# AM5 - is "intTLR7AB" significant? what about instead of lambda? - based on prior models, do I want to have beta in here?
am5 <- list("lambda" = "intTLR7AB",
            "ta" = c("intTLR7VAC", "intTLR7AB"),
            "betas" = "intTLR7VAC",
            "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5",
  model_file, data_file,
  igs, am5,
  headerTypes, error_model,
  fixed = fixed
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 832.0699059   866.0699059   898.2308510   911.0555978     5.0000000     0.4165849 

# AM5.1 - is "intTLR7AB" significant? what about instead of lambda? - based on prior models, do I want to have beta in here?
am5.1 <- list("lambda" = "intTLR7AB",
            "ta" = c("intTLR7VAC", "intTLR7AB"),
            "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5.1",
  model_file, data_file,
  igs, am5.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 841.154475    873.154475    903.423600    916.248346      5.000000      0.716868 

# AM5.2 - is "intTLR7AB" significant? what about instead of lambda? - based on prior models, do I want to have beta in here?
am5.2 <- list("betas" = "intTLR7AB",
              "ta" = c("intTLR7VAC", "intTLR7AB"),
              "NPs" = c("intTLR7AB", "intTLR7VAC"))
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5.2",
  model_file, data_file,
  igs, am5.2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 847.4947802   879.4947802   909.7639049   922.5886517     5.0000000     0.4052333 


# AM6 - do I want to have beta in here?
am6 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "betas" = "intTLR7VAC",
            "NPs" = "intTLR7ABVAC")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6",
  model_file, data_file,
  igs, am6,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 844.7534175   874.7534175   903.1307220   915.9554688     5.0000000     0.5481257 


# AM6.1 no beta
am6.1 <- list("lambda" = "intTLR7AB",
            "ta" = "intTLR7VAC",
            "NPs" = "intTLR7ABVAC")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6.1",
  model_file, data_file,
  igs, am6.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 857.6352699   885.6352699   912.1207541   924.9455009     5.0000000     0.6263893 

# AM6.2 ab on beta instead of lambda
am6.2 <- list("betas" = "intTLR7AB",
              "ta" = "intTLR7VAC",
              "NPs" = "intTLR7ABVAC")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6.2",
  model_file, data_file,
  igs, am6.2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 872.7961920   900.7961920   927.2816762   940.1064230     5.0000000     0.5489666 
