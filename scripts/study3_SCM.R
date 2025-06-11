# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study3"
project_name <- "SCM.mlxtran"
model_file <- here("model", "2regime_model_final.txt")

## Results folder
out_folder <- here("results", study)

## Data file
data_file <- here("data", "source", "ABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov",3))

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
# Use 1 for standard deviations unless the best estimate for it is bigger
igs <- read_csv(here("results", "study3_control", "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
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
       out$plot, width = 5, height = 7)

# BICc: 753	
# OFV: 687
# All covariates are signficant
# Condition number 12


# Try alternative models ####

# ab on beta instead of lambda

am1 <- list("betas" = "ab",
            "m" = "tlr7",
            "NPs" = "ab")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1", 
  model_file, data_file, 
  igs, am1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 685.3767570   713.3767570   738.3554118   751.5506985     5.0000000     0.3193617 

# Swap ab effect on Np for p

am2 <- list("betas" = "ab",
            "m" = "tlr7",
            "p" = "ab")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 686.5674003   714.5674003   739.5460552   752.7413418     5.0000000     0.2978999 

# Equal to ab on Np

## Try immune effects all on m

am3 <- list("betas" = "ab",
            "m" = c("tlr7", "ab"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3", 
  model_file, data_file, 
  igs, am3, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 707.8922311   735.8922311   760.8708860   774.0661726     5.0000000     0.3699722 

# Significantly worse fit

## Try immune effects all on Np

am4 <- list("betas" = "ab",
            "NPs" = c("tlr7", "ab"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4", 
  model_file, data_file, 
  igs, am4, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 697.3984894   725.3984894   750.3771443   763.5724309     5.0000000     0.3601044 

# Significantly worse fit

## Try immune effects all on p

am5 <- list("betas" = "ab",
            "p" = c("tlr7", "ab"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5", 
  model_file, data_file, 
  igs, am5, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 700.9277227   728.9277227   753.9063776   767.1016642     5.0000000     0.3027202 

# Worse fit, cannot calculate FIM

## We inferred tlr7 on lambda in the other study, does it make sense here?

am6 <- list("betas" = "ab",
            "NPs" = "ab",
            "lambda" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6", 
  model_file, data_file, 
  igs, am6, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 694.6409014   722.6409014   747.6195562   760.8148429     5.0000000     0.3702332 


## BH: AB influences latent reservoir 
bh1 <- list("ta" = "ab",
           "NPs" = "ab",
           "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh1", 
  model_file, data_file, 
  igs, bh1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 701.3071391   729.3071391   754.2857940   767.4810806     5.0000000     0.4295346 

# effects on m and ta are insignificant
# condition number > 100

## BH: TLR7 influences latent reservoir
bh2 <- list("ta" = "tlr7",
            "NPs" = "ab",
            "betas" = "ab")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh2", 
  model_file, data_file, 
  igs, bh2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 700.2922864   728.2922864   753.2709413   766.4662279     5.0000000     0.3871229 

# TLR7 on ta insignificant

## BH: TLR7 and AB interact to influence Np
bh3 <- list("betas" = "ab",
            "NPs" = c("ab", "inttlr7ab"),
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh3", 
  model_file, data_file, 
  igs, bh3, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 687.9582819   717.9582819   744.7211264   757.9164130     5.0000000     0.3094995 

# Interaction effect is not signficant
# FIM condition number < 100


## RM: effects inferred in vaccine study
rm1 <- list("lambda" = "tlr7",
            "NPs" = c("ab"),
            "ta" = c("ab", "tlr7")
            )

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1", 
  model_file, data_file, 
  igs, rm1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 691.4471560   721.4471560   748.2100005   761.4052871     5.0000000     0.3637249 

# All coefficients are significant
# FIM # < 100

# TLR7 effects are negative .... this is the opposite of what we see in the other study


## RM: effects inferred in vaccine study plus interaction on ta
rm1b <- list("lambda" = "tlr7",
            "NPs" = c("ab"),
            "ta" = c("ab", "tlr7", "inttlr7ab")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1b", 
  model_file, data_file, 
  igs, rm1b, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 690.7938955   722.7938955   751.3409297   764.5362163     5.0000000     0.3222026 

# Interaction effect is insiginifcant
# TLR7 effects are negative....

## RM: effects inferred in vaccine study plus interaction on Np
rm2 <- list("lambda" = "tlr7",
             "NPs" = c("ab", "inttlr7ab"),
             "ta" = c("ab", "tlr7")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2", 
  model_file, data_file, 
  igs, rm2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 696.0827649   728.0827649   756.6297990   769.8250857     5.0000000     0.2692895 

# Several effects insignifiacnt, big condiiton number
# TLR7 effects are in the opposite direction

## RM: effects inferred in vaccine study plus interaction on Np and ta
rm2b <- list("lambda" = "tlr7",
            "NPs" = c("ab", "inttlr7ab"),
            "ta" = c("ab", "tlr7", "inttlr7ab")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2b", 
  model_file, data_file, 
  igs, rm2b, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 692.1910298   726.1910298   756.5222535   769.7175402     5.0000000     0.2918777 

# Interaction effects are insignifcant
# TLR7 effects are negative....