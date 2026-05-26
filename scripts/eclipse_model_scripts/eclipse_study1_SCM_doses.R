# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study1"
project_name <- "SCM_doses.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "derived", "STUDY1_DATA_DOSING.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", "ignore", 
                 "catcov", rep("contcov", 4), "ignore", "contcov")

## Random effects to exclude
exclude_random_effects <- NULL # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 200, "p_pop" = 3)

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
relationships <- NULL # nothing locked in
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)


loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_plot.png"),
       out$plot, width = 5, height = 7)

## Alternative model without m effect
am1 <- list("betas" = "max_dose",
            "ta" = "cum_dose")
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
  model_file, data_file,
  igs, am1,
  headerTypes, error_model,
  fixed = fixed
)


## Alternative model with effects on lambda instead of beta
am2 <- list("lambda" = "max_dose",
            "ta" = "cum_dose")
alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
  model_file, data_file,
  igs, am2,
  headerTypes, error_model,
  fixed = fixed
)


## Are random effects of m from model without treatment effect on m correlated with max dose?
indicator_model_RE <- read_csv(here("results_eclipse", "study1", "best_model", "IndividualParameters", "estimatedRandomEffects.txt"))
indicator_model_EBE <- read_csv(here("results_eclipse", "study1", "best_model", "IndividualParameters", "estimatedIndividualParameters.txt"))

left_join(distinct(data, id, max_dose),
          indicator_model_RE %>% select(id, eta_m_mode)) %>% 
  ggplot(aes(max_dose, eta_m_mode)) + 
  geom_point()

left_join(distinct(data, id, max_dose),
          indicator_model_EBE %>% select(id, m_mode)) %>% 
  ggplot(aes(max_dose, m_mode)) + 
  geom_point() + 
  scale_y_log10()

left_join(distinct(data, id, max_dose),
          indicator_model_RE %>% select(id, eta_m_mode)) %>% 
  summarise(cor(max_dose, eta_m_mode))
