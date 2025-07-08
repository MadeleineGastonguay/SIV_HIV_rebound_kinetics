# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study3"
project_name <- "SCM.mlxtran"
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


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 688.4596544   716.4596544   741.4383093   754.6335959     5.0000000     0.2305609 

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
# 690.6165387   718.6165387   743.5951936   756.7904802     5.0000000     0.2393741 

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
# 691.645967    719.645967    744.624621    757.819908      5.000000      0.233143 



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
# 692.2264464   722.2264464   748.9892910   762.1845776     5.0000000     0.2239323 


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
# 695.2197544   723.2197544   748.1984093   761.3936959     5.0000000     0.3570143 




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
# 689.5344600   719.5344600   746.2973046   759.4925912     5.0000000     0.2967223 



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
# 697.0291410   725.0291410   750.0077959   763.2030825     5.0000000     0.3202687 

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
# 687.8372354   717.8372354   744.6000799   757.7953666     5.0000000     0.2833723 


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
# 688.9680914   718.9680914   745.7309359   758.9262226     5.0000000     0.2819871 

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
# 701.4744149   729.4744149   754.4530698   767.6483564     5.0000000     0.3039302 
