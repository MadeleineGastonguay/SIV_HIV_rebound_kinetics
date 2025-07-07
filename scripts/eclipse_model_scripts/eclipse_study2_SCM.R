# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study2"
project_name <- "SCM.mlxtran"
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
fixed <- c("lambda_pop" = 50, "p_pop" = 0.7)

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


## AM1: Final model without tlr7*vac on ta because it doesn't significantly improve BICc

am1 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"),
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
  model_file, data_file,
  igs, am1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 450.5111415   480.5111415   503.4065493   516.2312961     5.0000000     0.2451323 


## AM2: final model without marginally significant effects

am2 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
  model_file, data_file,
  igs, am2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 458.030789    486.030789    507.399836    520.224583      5.000000      0.263612 

## AM3: tlr7*vac on ta instead of Np, keep vaccine effect on lamdba
am3 <- list("ta" = c("vac", "intTLR7VAC"),
            "NPs" = "vac",
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
  model_file, data_file,
  igs, am3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 454.3850019   484.3850019   507.2804098   520.1051566     5.0000000     0.2352656 

## AM3: tlr7*vac on ta instead of Np, no vaccine effect on lamdba
am3.1 <- list("ta" = c("vac", "intTLR7VAC"),
            "NPs" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3.1",
  model_file, data_file,
  igs, am3.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 459.8923310   487.8923310   509.2613784   522.0861252     5.0000000     0.2412635 

## AM4: tlr7 instead ot tlr7*vac on Np, keep vaccine effect on lambda
am4 <- list("ta" = "vac",
            "NPs" = c("tlr7", "vac"),
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
  model_file, data_file,
  igs, am4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 454.3960248   484.3960248   507.2914327   520.1161795     5.0000000     0.2478875 


## AM4.1: tlr7 instead of tlr7*vac on Np, no vac on lambda
am4.1 <- list("ta" = "vac",
            "NPs" = c("tlr7", "vac"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.1",
  model_file, data_file,
  igs, am4.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 458.9896446   486.9896446   508.3586920   521.1834387     5.0000000     0.2620964 