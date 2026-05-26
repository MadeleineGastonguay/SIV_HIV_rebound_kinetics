# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM_tlr7_exclude_tlr7_beta.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
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
# relationships <- NULL # nothing locked in
relationships <-  getModelBuildingSettings()$relationships
relationships[1,] <- c("betas", "tlr7", FALSE) # Exclude TLR7 effect from beta
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)

loadProject(here(out_folder, project_name))

# BEst model BICc: 916.4945	

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 10, height = 7)

write_csv(out$summary, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_summary.csv"))


## Try Alternative Models:

# Model with Ab + Vac on Np instead of TLR7 
am1 <- list("lambda" = "ab",
            "m" = "ab",
            "NPs" = c("ab", "vac")
            )

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
        model_file, data_file,
        igs, am1,
        headerTypes, error_model,
        fixed = fixed)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 845.8115551   875.8115551   904.1888596   917.0136064     5.0000000     0.5974002 

# Model with Vac on ta instead of ab on m
am2 <- list("lambda" = "ab",
            "ta" = "vac",
            "NPs" = c("ab", "vac")
)

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
        model_file, data_file,
        igs, am2,
        headerTypes, error_model,
        fixed = fixed)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 842.9299727   872.9299727   901.3072772   914.1320240     5.0000000     0.5482523 

# Model with Vac on ta in addition to ab on m
am2.1 <- list("lambda" = "ab",
              "ta" = "vac",
              "m" = "ab",
              "NPs" = c("ab", "vac")
)

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2.1",
        model_file, data_file,
        igs, am2.1,
        headerTypes, error_model,
        fixed = fixed)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 835.4198773   867.4198773   897.6890021   910.5137489     5.0000000     0.5285438 
# vac on ta insignificant

# Model with ab on beta instead of lambda
am3 <- list("betas" = "ab",
            "ta" = "vac",
            "NPs" = c("ab", "vac")
)

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
        model_file, data_file,
        igs, am3,
        headerTypes, error_model,
        fixed = fixed)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 860.7110928   890.7110928   919.0883973   931.9131441     5.0000000     0.5048768 

# Model with tlr7 on lambda
am4 <- list("lambda" = c("ab", "tlr7"),
            "ta" = "vac",
            "NPs" = c("ab", "vac")
)

alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
        model_file, data_file,
        igs, am4,
        headerTypes, error_model,
        fixed = fixed)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 843.8557769   875.8557769   906.1249017   918.9496485     5.0000000     0.4797469 
