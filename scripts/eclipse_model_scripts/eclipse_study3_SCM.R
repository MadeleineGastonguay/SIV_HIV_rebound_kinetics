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
       out$plot, width = 5, height = 7)

# BICc: 756, all effects significant, condition number: 27

# AM1 (identifiability): ab on lambda instead of beta

am1 <- list("lambda" = "ab",
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
# 690.1720109   718.1720109   743.1506658   756.3459525     5.0000000     0.2510136 

# AM2 (identifiability): ab on p instead of Np

am2 <- list("betas" = "ab",
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
# 687.7602472   715.7602472   740.7389021   753.9341888     5.0000000     0.2497306 


# AM3 (SCM step similarity): ab*tlr7 instead of ab on beta

am3 <- list("betas" = "inttlr7ab",
            "p" = "ab",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
  model_file, data_file,
  igs, am3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 697.5241528   725.5241528   750.5028077   763.6980943     5.0000000     0.2744774 

# AM4: tlr7 on lambda instead of m since we inferred that in study 1
am4 <- list("betas" = "ab",
            "NPs" = "ab",
            "lambda" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
  model_file, data_file,
  igs, am4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 697.1583185   725.1583185   750.1369733   763.3322600     5.0000000     0.2305163

# AM4: tlr7 on lambda instead of m since we inferred that in study 1
am4.1 <- list("betas" = "ab",
            "NPs" = "ab",
            "lambda" = "tlr7",
            "m" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.1",
  model_file, data_file,
  igs, am4.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 689.181824    719.181824    745.944669    759.139955      5.000000      0.269306 

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
# 699.0973636   727.0973636   752.0760185   765.2713052     5.0000000     0.3167802 

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
# 688.2893644   718.2893644   745.0522089   758.2474956     5.0000000     0.3026112 

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
# 700.0644973   728.0644973   753.0431522   766.2384388     5.0000000     0.3363471

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
# 689.1118881   719.1118881   745.8747326   759.0700193     5.0000000     0.2549358 



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
# 686.7416003   716.7416003   743.5044448   756.6997315     5.0000000     0.1898599 