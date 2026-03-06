# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study3"
project_name <- "SCM_smaller_lambda.mlxtran"
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
fixed <- c("lambda_pop" = 200, "p_pop" = 1)

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
       out$plot, width = 10, height = 7)

write_csv(out$summary, here(out_folder, str_remove(project_name, "\\.mlxtran$"), "model_building_summary.csv"))

# BICc: 752.9313	


# AM1 ab on beta instead of lambda

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
# 688.6718082   716.6718082   741.6504631   754.8457498     5.0000000     0.3560193 

# AM2 ab on p instead of Np

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
# 687.4095179   715.4095179   740.3881727   753.5834594     5.0000000     0.2752094 

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
# 692.6669385   720.6669385   745.6455934   758.8408800     5.0000000     0.3290765 



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
# 687.5753661   717.5753661   744.3382106   757.5334972     5.0000000     0.3173582 

# AM4: AM1 with tlr7 on Np instead of m

am4 <- list("betas" = "ab",
            "NPs" = c("ab", "tlr7"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
  model_file, data_file,
  igs, am4,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 698.1615150   726.1615150   751.1401698   764.3354565     5.0000000     0.3417688 

# AM5: AM1 with ab effect on m instead of Np

am5 <- list("betas" = "ab",
            "m" = c("ab", "tlr7"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5",
  model_file, data_file,
  igs, am5,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 702.3129247   730.3129247   755.2915796   768.4868662     5.0000000     0.2588305 

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
# 697.0383263   725.0383263   750.0169811   763.2122678     5.0000000     0.3327012 

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
# 689.5745857   719.5745857   746.3374302   759.5327169     5.0000000     0.3235969 

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
# 698.7151849   726.7151849   751.6938398   764.8891264     5.0000000     0.3532935 

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
# 687.3965961   717.3965961   744.1594406   757.3547273     5.0000000     0.3291522 

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
# 688.4459367   718.4459367   745.2087813   758.4040679     5.0000000     0.3533916 

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
# 699.0507141   727.0507141   752.0293690   765.2246556     5.0000000     0.2401535 


# RM1: effects from SIV studies without interaction
rm1 <- list("lambda"= "tlr7",
            "ta" = c("tlr7", "ab"),
            "NPs" = "ab")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1",
  model_file, data_file,
  igs, rm1,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 689.0835191   719.0835191   745.8463637   759.0416503     5.0000000     0.2872809 

# RM1b: effects from SIV studies without interaction plus ta interaction
rm1b <- list("lambda"= "tlr7",
             "ta" = c("tlr7", "ab", "inttlr7ab"),
             "NPs" = "ab")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm1b",
  model_file, data_file,
  igs, rm1b,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 690.1809717   722.1809717   750.7280058   763.9232925     5.0000000     0.2637019 

# RM2: effects from SIV studies with interaction
rm2 <- list("lambda"= "tlr7",
            "ta" = c("tlr7", "ab"),
            "NPs" = c("ab", "inttlr7ab"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2",
  model_file, data_file,
  igs, rm2,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 690.2707296   722.2707296   750.8177637   764.0130504     5.0000000     0.2188556

# RM2b: Effects from SIV studies with ta and Np interaction
rm2b <- list("lambda"= "tlr7",
             "ta" = c("tlr7", "ab", "inttlr7ab"),
             "NPs" = c("ab", "inttlr7ab"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm2b",
  model_file, data_file,
  igs, rm2b,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 693.6445423   727.6445423   757.9757661   771.1710527     5.0000000     0.3220613 

# RM3: Vaccine effects from SIV + TLR7 effects from SHIV
rm3 <- list("m"= "tlr7",
             "ta" = "ab",
             "NPs" = "ab")

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm3",
  model_file, data_file,
  igs, rm3,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 697.0383263   725.0383263   750.0169811   763.2122678     5.0000000     0.3327012 

# RM3b: Vaccine effects from SIV + TLR7 effects from SHIV + interaction effect
rm3b <- list("m"= "tlr7",
            "ta" = "ab",
            "NPs" = c("ab", "inttlr7ab"))

alt_mod(
  here(out_folder, str_remove(project_name, "\\.mlxtran$")), "rm3b",
  model_file, data_file,
  igs, rm3b,
  headerTypes, error_model,
  fixed = fixed
)

# OFV           AIC           BIC          BICc  chosenDegree standardError 
# 699.366488    729.366488    756.129333    769.324619      5.000000      0.355372 