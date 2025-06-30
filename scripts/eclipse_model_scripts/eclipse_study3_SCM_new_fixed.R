# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study3"
project_name <- "SCM_new_fixed.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "source", "ABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", ncol(data) - 4))

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 500, "p_pop" = 1)

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
igs <- read_csv(here(out_folder, "IGS_no_cov_fixed_param_new", "best_estimates.csv")) %>% 
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



# AM1 ab on Np instead of p

am1 <- list("betas" = "ab",
            "NPs" = "ab",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
  model_file, data_file,
  igs, am1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree         standardError 
# 687.8296195   715.8296195   740.8082743   754.0035610     5.0000000 0.2936267 


# AM2 (identifiability): ab on lambda instead of beta

am2 <- list("lambda" = "ab",
            "p" = "ab",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
  model_file, data_file,
  igs, am2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 687.6076929   715.6076929   740.5863477   753.7816344     5.0000000     0.3465158 


# AM3: tlr7 on lambda since we inferred that in study 1

am3 <- list("betas" = "ab",
            "NPs" = "ab",
            "lambda" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
  model_file, data_file,
  igs, am3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 693.0077239   721.0077239   745.9863788   759.1816655     5.0000000     0.2912992 

# AM3: tlr7 on lambda since we inferred that in study 1

am3.1 <- list("betas" = "ab",
            "NPs" = "ab",
            "lambda" = "tlr7",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3.1",
  model_file, data_file,
  igs, am3.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 691.754006    721.754006    748.516851    761.712137      5.000000      0.281861 

# BH1: ab reduces size of latent reservoir
bh1 <- list("ta" = "ab",
            "NPs" = "ab",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh1",
  model_file, data_file,
  igs, bh1,
  headerTypes, error_model,
  fixed = fixed
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 695.0289143   723.0289143   748.0075692   761.2028558     5.0000000     0.2747726 

# BH1.1: ab reduces size of latent reservoir
bh1.1 <- list("ta" = "ab",
            "betas" = "ab",
            "NPs" = "ab",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh1.1",
  model_file, data_file,
  igs, bh1.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 688.7128736   718.7128736   745.4757182   758.6710048     5.0000000     0.3149708 



# BH2: TLR7 reduces size of latent reservoir (shown in study 1)
bh2 <- list("betas" = "ab",
            "NPs" = "ab",
            "ta" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh2",
  model_file, data_file,
  igs, bh2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 696.326367    724.326367    749.305022    762.500308      5.000000      0.315193 


# BH2: TLR7 reduces size of latent reservoir (shown in study 1)
bh2.1 <- list("betas" = "ab",
            "NPs" = "ab",
            "ta" = "tlr7",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh2.1",
  model_file, data_file,
  igs, bh2.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 692.8304986   722.8304986   749.5933431   762.7886297     5.0000000     0.2710956 


# BH3: TLR7 might interact with Ab to boost immune system, which we see in study 2
bh3 <- list("betas" = "ab",
            "NPs" = c("inttlr7ab", "ab"),
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh3",
  model_file, data_file,
  igs, bh3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 689.90382     719.90382     746.66666     759.86195       5.00000       0.29476 


# BH3.1: TLR7 might interact with Ab to boost immune system, which we see in study 2
bh3.1 <- list("betas" = "ab",
            "NPs" = c("inttlr7ab"),
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh3.1",
  model_file, data_file,
  igs, bh3.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 699.0423987   727.0423987   752.0210536   765.2163402     5.0000000     0.2390052 

