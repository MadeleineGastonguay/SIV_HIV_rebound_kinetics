# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study2"
project_name <- "SCM_new_fixed.mlxtran"
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

# BICc: 517, condition number: 510, several insignificant effects

## AM1: SCM up to step 5 (no vaccine on beta)

am1 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"),
            "m" = "vac",
            "p" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
  model_file, data_file,
  igs, am1,
  headerTypes, error_model,
  fixed = fixed
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 449.6656374   481.6656374   506.0874058   518.9121526     5.0000000     0.2667942 


## AM2: SCM without insignificant effects

am2 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"),
            "p" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
  model_file, data_file,
  igs, am2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 458.054850    488.054850    510.950258    523.775005      5.000000      0.270129 

## AM3: vaccine effect on lambda instead of beta or m

am3 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"),
            "lambda" = "vac",
            "p" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
  model_file, data_file,
  igs, am3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 446.8947601   478.8947601   503.3165285   516.1412753     5.0000000     0.2832113


## AM4: AM3 without vaccine effect on p

am4 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"),
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
  model_file, data_file,
  igs, am4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 452.9871328   482.9871328   505.8825406   518.7072874     5.0000000     0.2281881 

## AM5: AM3 without vaccine effect on p or lambda

am5 <- list("ta" = "vac",
            "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5",
  model_file, data_file,
  igs, am5,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 465.725315    493.725315    515.094362    527.919109      5.000000      0.268229 


## AM6: test for tlr7*vac on ta 
am6 <- list("ta" = c("vac", "intTLR7VAC"),
            "NPs" = c("vac", "intTLR7VAC"),
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6",
  model_file, data_file,
  igs, am6,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 449.4745725   481.4745725   505.8963409   518.7210877     5.0000000     0.2510373 

## AM6.1: test for tlr7*vac on ta 
am6.1 <- list("ta" = c("vac", "intTLR7VAC"),
            "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6.1",
  model_file, data_file,
  igs, am6.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 458.7402880   488.7402880   511.6356958   524.4604426     5.0000000     0.2692892 

## AM6.2: test for tlr7*vac on ta 
am6.2 <- list("ta" = c("vac", "intTLR7VAC"),
            "NPs" = c("vac", "intTLR7VAC"),
            "m" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6.2",
  model_file, data_file,
  igs, am6.2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 450.2078235   482.2078235   506.6295919   519.4543387     5.0000000     0.2654896

## AM6.3: test for tlr7*vac on ta 
am6.3 <- list("ta" = c("vac", "intTLR7VAC"),
              "NPs" = c("vac", "intTLR7VAC"),
              "lambda" = "vac",
              "p" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6.3",
  model_file, data_file,
  igs, am6.3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 440.3472      474.3472      500.2953      513.1201        5.0000        0.1715 

## AM7: test for tlr7*vac on ta 
am7 <- list("ta" = c("intTLR7VAC"),
              "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am7",
  model_file, data_file,
  igs, am7,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 465.7983876   493.7983876   515.1674350   527.9921818     5.0000000     0.3502799


## AM7: test for tlr7*vac on ta 
am7.1 <- list("ta" = c("intTLR7VAC"),
            "NPs" = c("vac", "intTLR7VAC"),
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am7.1",
  model_file, data_file,
  igs, am7.1,
  headerTypes, error_model,
  fixed = fixed
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 465.3745535   495.3745535   518.2699614   531.0947082     5.0000000     0.3105186 


## AM7: test for tlr7*vac on ta 
am7.2 <- list("ta" = c("intTLR7VAC"),
              "NPs" = c("vac", "intTLR7VAC"),
              "m" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am7.2",
  model_file, data_file,
  igs, am7.2,
  headerTypes, error_model,
  fixed = fixed
)


# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 456.3280463   486.3280463   509.2234541   522.0482009     5.0000000     0.2205453 

## AM8: tlr7*vac on ta instead of Np
am8 <- list("ta" = c("intTLR7VAC", "vac"),
              "NPs" = c("vac"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am8",
  model_file, data_file,
  igs, am8,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 467.0504553   495.0504553   516.4195026   529.2442494     5.0000000     0.3468082 

## AM8: tlr7*vac on ta instead of Np
am8.1 <- list("ta" = c("intTLR7VAC", "vac"),
            "NPs" = c("vac"),
            "lambda" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am8.1",
  model_file, data_file,
  igs, am8.1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 458.1173874   488.1173874   511.0127952   523.8375420     5.0000000     0.2013012 
