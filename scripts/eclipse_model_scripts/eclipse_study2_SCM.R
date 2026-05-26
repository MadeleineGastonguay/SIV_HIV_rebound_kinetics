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

if(!file.exists(here(out_folder, project_name))) saveProject(here(out_folder, project_name))

# Run SCM
relationships <- NULL # nothing locked in
covariate_model_building(here(out_folder, project_name), method = "SCM", criteria = "BIC", 
                         linearization = FALSE, relationships = relationships)

loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

write_csv(out$summary, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_summary.csv"))

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

## AM5: immune effects on p instead of Np
am5 <- list("ta" = "vac",
              "p" = c("intTLR7VAC", "vac"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5",
  model_file, data_file,
  igs, am5,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 530.6418772   558.6418772   580.0109246   592.8356713     5.0000000     0.1783943 

## AM6: immune effects on m instead of Np
am6 <- list("ta" = "vac",
            "m" = c("intTLR7VAC", "vac"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am6",
  model_file, data_file,
  igs, am6,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 491.8627900   519.8627900   541.2318374   554.0565841     5.0000000     0.2809666 

## AM7: add TLR7 effects from study 1
am7 <- list("ta" = c("vac", "tlr7"),
            "NPs" = c("intTLR7VAC", "vac"),
            "lambda" = "tlr7")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am7",
  model_file, data_file,
  igs, am7,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 453.4134274   485.4134274   509.8351958   522.6599426     5.0000000     0.2660963


## rm1: vaccine effects estimated for SHIV
rm1 <- list("betas" = "vac", 
            "NPs" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1",
  model_file, data_file,
  igs, rm1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 478.9465817   504.9465817   524.7892685   537.6140153     5.0000000     0.2491585 

## rm2: vaccine effects estimated for SHIV plus estimates from final model
rm2 <- list("betas" = "vac", 
            "NPs" = c("vac", "intTLR7VAC"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2",
  model_file, data_file,
  igs, rm2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 470.5141533   498.5141533   519.8832007   532.7079474     5.0000000     0.2448972 

## Additional RMS:

bh2 <- list(
  "ta" = c("tlr7","vac"),
  "NPs" = c("vac", "intTLR7VAC")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "bh2",
  model_file, data_file,
  igs, bh2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 457.2281358   487.2281358   510.1235437   522.9482905     5.0000000     0.3018443

rm2b <- list(
  "lambda" = "vac",
  "NPs" = "vac")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2b",
  model_file, data_file,
  igs, rm2b,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 479.8830276   505.8830276   525.7257144   538.5504612     5.0000000     0.2747434 


rm3b <- list(
  "lambda" = "vac",
  "NPs" = c("vac", "intTLR7VAC")
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm3b",
  model_file, data_file,
  igs, rm3b,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 470.4516142   498.4516142   519.8206615   532.6454083     5.0000000     0.2752531 

rm4 <- list(
  "m" = "tlr7",
  "NPs" = c("vac", "intTLR7VAC"),
  "betas" = "vac"
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm4",
  model_file, data_file,
  igs, rm4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 468.138961    498.138961    521.034369    533.859116      5.000000      0.268917 

rm4b <- list(
  "betas" = "tlr7",
  "NPs" = c("vac", "intTLR7VAC"),
  "lambda" = "vac"
)

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm4b",
  model_file, data_file,
  igs, rm4b,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 469.2986109   499.2986109   522.1940188   535.0187655     5.0000000     0.3725877 