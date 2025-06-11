# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study2"
project_name <- "SCM.mlxtran"
model_file <- here("model", "2regime_model_eclipse.txt")

## Results folder
out_folder <- here("compare_delay_models_timing_perturb_ic_study2", "eclipse")

## Data file
data_file <- here("data", "derived",str_interp("${toupper(study)}_DATA.txt"))
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov",3))

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 70, "p_pop" = 0.7)

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
igs <- read_csv(here("results", "study2", "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
  rowwise %>% 
  mutate(estimate = ifelse(grepl("pop", parameter), estimate, 1),
         estimate = signif(estimate, 1)) %>% 
  select(parameter, estimate)

set_base(igs, fixed = fixed)

saveProject(here(out_folder, project_name))
# loadProject(here(out_folder, project_name))

# Run SCM

start_time <- Sys.time()

covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", linearization = FALSE)

end_time <- Sys.time()

print(as.numeric(difftime(end_time, start_time, units = "secs")))
# 19589.67 seconds = 5.4 hours

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb, ncol = 1)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

# BICc: 515
# The vaccine effect on beta is insignificant 

# The vaccine effect on beta is insignificant 
# Test alternate model ####

am1 <- list(
  "ta" = "vac",
  "NPs" = c("vac", "intTLR7VAC")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1", 
  model_file, data_file, 
  igs, am1, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 456.3426195   484.3426195   505.7116668   518.5364136     5.0000000     0.2880185 

# We inferred TLR7 effect on ta and lambda in the other study. Is there evidence for them here?

am2 <- list(
  "ta" = c("vac", "tlr7"),
  "NPs" = c("vac", "intTLR7VAC"),
  "lambda" = "tlr7"
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 454.9963740   486.9963740   511.4181424   524.2428892     5.0000000     0.2738282 

# TLR7 on lambda is insignificant

# What about just on ta?
am3 <- list(
  "ta" = c("vac", "tlr7"),
  "NPs" = c("vac", "intTLR7VAC")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3", 
  model_file, data_file, 
  igs, am3, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 454.2948768   484.2948768   507.1902847   520.0150315     5.0000000     0.2489209 

# What about just on lambda?
am3b <- list(
  "lambda" = "tlr7",
  "ta" = c("vac"),
  "NPs" = c("vac", "intTLR7VAC")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3b", 
  model_file, data_file, 
  igs, am3b, 
  headerTypes, error_model,
  fixed = fixed
)

# TLR7 on lambda still not significant
# BICc: 521.8501	

# Should interaction of vaccine and tlr7 on NP be additive?
am4 <- list(
  "ta" = c("vac"),
  "NPs" = c("vac", "tlr7")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4", 
  model_file, data_file, 
  igs, am4, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 459.1689085   487.1689085   508.5379558   521.3627026     5.0000000     0.2736671 

# TLR7 effect on Np is insignificant


# Should interaction of vaccine and tlr7 on NP be additive?
am5 <- list(
  "ta" = c("vac"),
  "NPs" = c("vac", "tlr7", "intTLR7VAC")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5", 
  model_file, data_file, 
  igs, am5, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 455.8222368   485.8222368   508.7176447   521.5423915     5.0000000     0.2903877 

# TLR7 on Np is insignificant, and the interaction remains significant 

# Would the vaccine effect on beta be significant on lambda?
am6 <- list(
  "ta" = c("vac"),
  "NPs" = c("vac", "intTLR7VAC"),
  "lambda" = "vac"
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6", 
  model_file, data_file, 
  igs, am6, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 449.8858153   479.8858153   502.7812232   515.6059700     5.0000000     0.3116782 

# The vaccine effect on lambda is significant

## Test BH models ###

# TLR7 and Vaccine interact to reduce latent reservoir
bh2 <- list("ta" = c("intTLR7VAC", "vac"),
            "NPs" = c("vac"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh2", 
  model_file, data_file, 
  igs, bh2, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 460.6873581   488.6873581   510.0564054   522.8811522     5.0000000     0.3067989 

# Vaccine effect on ta is not significant, but interaction effect is

bh2b <- list("ta" = c("intTLR7VAC", "vac"),
             "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh2b", 
  model_file, data_file, 
  igs, bh2b, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 452.3008694   482.3008694   505.1962772   518.0210240     5.0000000     0.2801614 

# All coefficients are singificant

## Test RM models ###

# In SHIV study, we estimated an Ab effect on beta and Np and a tlr7 effect on m

rm1 <- list("betas" = "vac",
            "NPs" = "vac",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1", 
  model_file, data_file, 
  igs, rm1, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 480.2753539   508.2753539   529.6444012   542.4691480     5.0000000     0.3135114 

# Worse fit, vaccine is estimated to reduce, not increase m

rm2 <- list("betas" = "vac",
            "NPs" = c("vac", "intTLR7VAC"),
            "m" = "tlr7",
            "ta" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2", 
  model_file, data_file, 
  igs, rm2, 
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 449.187136    481.187136    505.608905    518.433652      5.000000      0.251931 

# tlr7 on m and vac on betas are not signficant
# condition number > 500

