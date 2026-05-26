# Load libraries
library(here)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM_tlr7_exclude_beta_RE.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Define header types from the following options:
# See ?newProject() for details
headerTypes <- c("id", "time", "observation", "cens", rep("catcov", 3), rep("ignore", 3))

## Random effects to exclude
exclude_random_effects <- c('betas') # Either NULL or a character vector of parameters with no random effects

## Parameters to fix
fixed <- c("lambda_pop" = 100, "p_pop" = 0.9)

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


# ## Try Alternative Models:
# 
# # Model with Ab + Vac on Np instead of TLR7 + Ab
# am1 <- list("lambda" = "ab",
#             "betas" = "tlr7",
#             "ta" = "vac",
#             "p" = c("ab", "vac"),
#             "NPs" = c("ab", "vac")
#             )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am1",
#         model_file, data_file, 
#         igs, am1,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 814.6406807   850.6406807   884.6934461   897.5181929     5.0000000     0.3723626 
# 
# # If we cutoff SCM at steps that include BICc > 2
# am2 <- list("lambda" = "ab",
#             "betas" = "tlr7",
#             "ta" = "vac",
#             "NPs" = "tlr7",
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am2",
#         model_file, data_file, 
#         igs, am2,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 828.9889876   860.9889876   891.2581124   904.0828592     5.0000000     0.4752432 
# 
# 
# # If we cutoff SCM at steps that include BICc > 2 ... how do Np effects compare?
# am3 <- list("lambda" = "ab",
#             "betas" = "tlr7",
#             "ta" = "vac",
#             "NPs" = c("tlr7", "vac", "ab"),
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3",
#         model_file, data_file, 
#         igs, am3,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 826.2435154   862.2435154   896.2962808   909.1210276     5.0000000     0.4601308 
# 
# # m on ab insignificant
# # vac and tlr7 on Np insignificant
# # FIM cond number < 100
# 
# # If we cutoff SCM at steps that include BICc > 2 ... how do Np effects compare?
# am3.1 <- list("lambda" = "ab",
#             "betas" = "tlr7",
#             "ta" = "vac",
#             "NPs" = c("vac", "ab"),
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3.1",
#         model_file, data_file, 
#         igs, am3.1,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 823.5770471   857.5770471   889.7379922   902.5627390     5.0000000     0.3366905 
# 
# # Better BICc than TLR7 on Np alone ... ab on m is insigificant
# 
# # If we cutoff SCM at steps that include BICc > 2 ... how do Np effects compare?
# am3.2 <- list("lambda" = "ab",
#               "betas" = "tlr7",
#               "ta" = "vac",
#               "NPs" = c("vac", "ab")
#               # "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am3.2",
#         model_file, data_file, 
#         igs, am3.2,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 831.1170082   863.1170082   893.3861329   906.2108797     5.0000000     0.3316962 
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap effect on beta and lambda?
# am4 <- list("lambda" = "tlr7",
#               "betas" = "ab",
#               "ta" = "vac",
#               "NPs" =  c("vac", "ab"),
#               "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4",
#         model_file, data_file, 
#         igs, am4,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 855.4734440   889.4734440   921.6343891   934.4591359     5.0000000     0.5492714 
# 
# # the model is much worse!
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap effect on beta and lambda?
# am4.4 <- list("lambda" = "tlr7",
#             "betas" = "ab",
#             "ta" = "vac",
#             "NPs" =  c("vac", "ab")
#             # "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.4",
#         model_file, data_file, 
#         igs, am4.4,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 857.2791787   889.2791787   919.5483034   932.3730502     5.0000000     0.433185
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap ab effect from lambda to beta?
# am4.1 <- list(#"lambda" = "tlr7",
#             "betas" = c("ab", "tlr7"),
#             "ta" = "vac",
#             "NPs" = c("vac", "ab"),
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.1",
#         model_file, data_file, 
#         igs, am4.1,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 847.4629971   881.4629971   913.6239421   926.4486889     5.0000000     0.4308543 
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap ab effect from lambda to beta?
# am4.5 <- list(#"lambda" = "tlr7",
#   "betas" = c("ab", "tlr7"),
#   "ta" = "vac",
#   "NPs" = c("vac", "ab")
#   # "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.5",
#         model_file, data_file, 
#         igs, am4.5,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 850.1837768   882.1837768   912.4529016   925.2776484     5.0000000     0.5279838 
# 
# # FIM condition nimber big, ab effects insigifnicat
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap ab effect from lambda to beta?
# am4.6 <- list(#"lambda" = "tlr7",
#   "lambda" = c("ab", "tlr7"),
#   "ta" = "vac",
#   "NPs" = c("vac", "ab")
#   # "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.6",
#         model_file, data_file, 
#         igs, am4.6,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 843.8557769   875.8557769   906.1249017   918.9496485     5.0000000     0.4797469 
# 
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap ab effect from lambda to beta?
# am4.2 <- list(#"lambda" = "tlr7",
#   "betas" = c("ab"),
#   "ta" = "vac",
#   "NPs" = c("vac", "ab"),
#   "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.2",
#         model_file, data_file, 
#         igs, am4.2,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 868.0440767   900.0440767   930.3132014   943.1379482     5.0000000     0.5415846 
# 
# # If we cutoff SCM at steps that include BICc > 2 ... can we swap ab effect from lambda to beta?
# am4.3 <- list(#"lambda" = "tlr7",
#   "betas" = c("ab"),
#   "ta" = "vac",
#   "NPs" = c("vac", "ab")
#   # "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am4.3",
#         model_file, data_file, 
#         igs, am4.3,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # Vac on ta is not longer signifcant, worse BICc than model with lambda! 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 860.7110928   890.7110928   919.0883973   931.9131441     5.0000000     0.5048768 
# 
# # How does vaccine and antibody effect on beta compare to TLR7?
# 
# am5 <- list("lambda" = "ab",
#             "betas" = c("ab", "vac"),
#             "NPs" = c("ab", "vac"),
#             "ta" = "vac"
#             )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "am5",
#         model_file, data_file, 
#         igs, am5,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 836.1123033   870.1123033   902.2732484   915.0979952     5.0000000     0.5770401 
# 
# # Worse fit! 
# 
# ## Backwards elimination for thresholded model ####
# 
# # If we cutoff SCM at steps that include BICc > 2, would  tlr7 on Beta be removed during backwards selection?
# be1 <- list("lambda" = "ab",
#             "ta" = "vac",
#             "NPs" = "tlr7",
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "be1",
#         model_file, data_file, 
#         igs, be1,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 857.6260907   887.6260907   916.0033951   928.8281419     5.0000000     0.4956124 
# 
# # If we cutoff SCM at steps that include BICc > 2, would ab on m be removed during backwards selection?
# be2 <- list("lambda" = "ab",
#            "betas" = "tlr7",
#            "ta" = "vac",
#            "NPs" = "tlr7"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "be2",
#         model_file, data_file, 
#         igs, be2,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 841.4517168   871.4517168   899.8290213   912.6537681     5.0000000     0.3499422 
# 
# # If we cutoff SCM at steps that include BICc > 2, would tlr7 on Np be removed during backwards selection?
# be3 <- list("lambda" = "ab",
#             "betas" = "tlr7",
#             "ta" = "vac",
#             # "NPs" = "tlr7",
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "be3",
#         model_file, data_file, 
#         igs, be3,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 852.2210165   882.2210165   910.5983209   923.4230677     5.0000000     0.5276581 
# 
# # If we cutoff SCM at steps that include BICc > 2, would ab on lambda be removed during backwards selection?
# be4 <- list(#"lambda" = "ab",
#             "betas" = "tlr7",
#             "ta" = "vac",
#             "NPs" = "tlr7",
#             "m" = "ab"
# )
# 
# alt_mod(here(out_folder, str_remove(project_name, "\\.mlxtran$")), "be4",
#         model_file, data_file, 
#         igs, be4,
#         headerTypes, error_model, 
#         fixed = fixed)
# 
# # OFV           AIC           BIC          BICc  chosenDegree standardError 
# # 873.5817764   903.5817764   931.9590808   944.7838276     5.0000000     0.5728369 
# 
# 
# 
