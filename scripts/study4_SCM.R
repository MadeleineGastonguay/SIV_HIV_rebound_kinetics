# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM.mlxtran"
model_file <- here("model", "2regime_model_final.txt")

## Results folder
out_folder <- here("results", study)

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("ignore", 3), rep("catcov", 3))

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 50, "p_pop" = 1)

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
igs <- read_csv(here("results", "study4_control", "IGS_no_cov_fixed_param", "best_estimates.csv")) %>% 
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

out <- plot_model_building(mb, ncol = 1)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

# BICc: 957.1811
# All effects are signficant

## Assess alternative models ####

final_model <- list(
  "lambda" = "inTLR7AB",
  "NPs" = c("intTLR7AB", "intTLR7VAC"),
  "m" = "intTLR7VAC"
)

# AM1 TLR7AB on betas instead of lambda

am1 <- final_model
am1$betas <- "intTLR7AB"
am1$lambda <- NULL

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1", 
  model_file, data_file, 
  igs, am1, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 886.1023608   916.1023608   944.4796653   957.3044121     5.0000000     0.4381983 

# Equivalent BICc, beta effect is signficant

# AM2 TLR7AB on ta

am2 <- am1
am2$ta <- "intTLR7AB"

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2", 
  model_file, data_file, 
  igs, am2, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 888.0028573   920.0028573   950.2719821   963.0967289     5.0000000     0.4421134 

# Effect on ta is insignificant


# AM3 TLR7ABVAC on Np instead of TLR7AB

am3 <- am1
am3$NPs <-  c("intTLR7ABVAC", "intTLR7VAC")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3", 
  model_file, data_file, 
  igs, am3, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 896.0891156   926.0891156   954.4664200   967.2911668     5.0000000     0.4928389 

# TLR7*VAC on Np is not significant


# AM4 TLR7VAC on ta instead of m

am4 <- am1
names(am4)[names(am4) == "m"] <- "ta"

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4", 
  model_file, data_file, 
  igs, am4, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 895.5570726   925.5570726   953.9343771   966.7591239     5.0000000     0.5880125 

# Ta effect is insignificant and it's estimated to be smaller ....

# AM5 Remove TLR7VAC on m
am5 <- am1
am5$m <- NULL

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5", 
  model_file, data_file, 
  igs, am5, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 890.7634941   918.7634941   945.2489782   958.0737250     5.0000000     0.4260288 

# Insignificant difference in BICc


# AM6 TLR7ABVAC on m instead of TLR7VAC
am6 <- am1
am6$m <- "intTLR7ABVAC"

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6", 
  model_file, data_file, 
  igs, am6, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 888.0796294   918.0796294   946.4569338   959.2816806     5.0000000     0.3216129 

# TLR7ABVAC effect is insigificant and estimates reduction in m

# AM7 TLR7AB on m
am7 <- am5
am7$m <- "intTLR7AB"

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am7", 
  model_file, data_file, 
  igs, am7, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 891.3752466   921.3752466   949.7525511   962.5772978     5.0000000     0.3168622 

# Insignificant, inferred to reduce m

# AM8 TLR7VAC on p instead of Np
am8 <- am5
am8$p <- "intTLR7VAC"
am8$NPs <- "intTLR7AB"

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am8", 
  model_file, data_file, 
  igs, am8, 
  headerTypes, error_model,
  fixed = fixed,
  exclude_random_effects = exclude_random_effects
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 894.6462560   922.6462560   949.1317402   961.9564870     5.0000000     0.4375549 

# Worse BICc, but significant effect on p

## Assess for effects on ta ####

# We previously inferred that the vaccine and TLR7 increases ta in SIV
# Going to run SCM with final model to see if adding any covariates on ta is significant
loadProject(here(out_folder, "SCM",  "am5.mlxtran"))
relationships <- data.frame(
  parameters = c("NPs", "NPs", "betas"),
  covariates = c("intTLR7AB", "intTLR7VAC", "intTLR7AB"),
  locked = TRUE
)
covariate_model_building(here(out_folder, "SCM_ta"), method = "SCM", criteria = "BIC", linearization = FALSE,
                         relationships = relationships, parameters = "ta")

# No evidence for any effects on ta!
