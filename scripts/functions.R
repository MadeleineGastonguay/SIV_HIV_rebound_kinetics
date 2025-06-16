################################################################################
#### This script contains all functions required for NLME model building and evaluation
################################################################################

require(tidyverse)
require(here)
require(Rsmlx)
require(ggbeeswarm)
require(patchwork)
require(tune)
require(gridExtra)
require(GGally)
require(ggpubr)
require(ggpmisc)


# Define project settings and scenario tasks ####

set_project_settings <- function(){
  # Set scenario for Monolix
  scenario <- getScenario()
  scenario$tasks["standardErrorEstimation"] <- TRUE
  scenario$tasks["logLikelihoodEstimation"] <- TRUE
  scenario$linearization <- FALSE
  setScenario(scenario)
  
  # Set Parameter Estimation Settings
  setPopulationParameterEstimationSettings(
    nbBurningIterations=50,
    nbExploratoryIterations=1000,
    nbSmoothingIterations=500,
    smoothingAlpha=0.7,
    exploratoryAutoStop = TRUE,
    smoothingAutoStop = TRUE
  )
  
}

# Create a new project ####

create_project <- function(model_file, data_file, 
                           headerTypes, error_model, 
                           exclude_random_effects, treatment_effects = NULL, 
                           normal_dist = NULL, fixed = NULL){
  
  # Create a new project
  newProject(modelFile = model_file,
             data = list(dataFile = data_file,
                         headerTypes = headerTypes
             )
  )
  
  # Set Error model
  setErrorModel(error_model)
  
  # Remove random effects if needed
  if(!is.null(exclude_random_effects)){
    var_model <- getIndividualParameterModel()$variability$id
    var_model[exclude_random_effects] <- FALSE
    setIndividualParameterVariability(var_model)
  }
  
  if(!is.null(treatment_effects)){
    # treatment_effects should be a list of vectors of categorical variables to add to parameters
    covariateModel <- getIndividualParameterModel()$covariateModel
    for(parameter in names(treatment_effects)){
      covariateModel[[parameter]][treatment_effects[[parameter]]] <- TRUE
    }
    setCovariateModel(covariateModel)
  }
  
  if(!is.null(normal_dist)){
    distributions <- getIndividualParameterModel()$distribution
    distributions[normal_dist] <- "normal"
    setIndividualParameterDistribution(distributions)
  }
  
  if(!is.null(fixed)){
    popparams <- getPopulationParameterInformation()
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- fixed
    setPopulationParameterInformation(popparams)
  }
  
  
  print(getIndividualParameterModel())
  
  set_project_settings()
  
}

# Set base model with initial values ####
## Requires a project to be loaded before running

set_base <- function(initial_values, omega_under_1 = F, fixed = NULL){
  
  # Usually, we don't want an initial guess for omega parameters less than 1 to help with convergence of SAEM
  if(!omega_under_1){
    initial_values <- initial_values %>% 
      rowwise %>% 
      # don't use an initial condition less than 1 for random effects to help with convergence
      mutate(estimate = ifelse(grepl("omega", parameter), max(1, estimate), estimate)) 
  }
  
  initial_values <- initial_values %>% 
    column_to_rownames("parameter") 
  
  popparams <- getPopulationParameterInformation()
  popparams$initialValue <- initial_values[popparams$name, "estimate"]
  
  # Set values for fixed parameters
  if(!is.null(fixed)){
    popparams[popparams$name %in% names(fixed),"method"] <- "FIXED"
    popparams[popparams$name %in% names(fixed),"initialValue"] <- fixed
  }
  
  # If an initial condition for a treatment effect is not provided, set it to 0
  if(any(is.na(popparams[grepl("beta_", popparams$name), "initialValue"]))){
    popparams[grepl("beta_", popparams$name), "initialValue"] <- 0
  }
  
  setPopulationParameterInformation(popparams)
}


# Generate a report from a monolix file ####
# By default, the report will be saved as "report.html" in the project directory

generate_report <- function(project_file, output_file = "report.html"){
  output_file <- here(str_remove(project_file, "\\.mlxtran$"), output_file)
  rmarkdown::render(
    here('scripts', 'Monolix_template.Rmd'),
    params = list(project = project_file),
    output_file = output_file
  )
}

# initial guess sampling function ####
# Takes a tibble with ranges of values, sampling method, and number of iterations
## Requires a project to be loaded before running

initial_guess_sampling <- function(initialValue_ranges, log10 = NULL, 
                                   nits = 100, set_variance = NULL, set_treatment = NULL,
                                   record_shrinkage = FALSE,
                                   fixed = NULL, predefined_initialvalues = NULL,
                                   record_pvals = FALSE, record_indiv_params = FALSE,
                                   record_condition_number = FALSE){
  
  
  popparams <- getPopulationParameterInformation()
  
  estimates <- tabestimates <- tabiters <-  phase <- initialvalues <- shrinkage_df <- treatment_pvals <- indivParams <- predictions <- condition_numbers <- NULL;
  
  run_condition <- TRUE
  i <- 1
  while(run_condition){
    print(i)
    
    if(is.null(predefined_initialvalues)){
      
      # sample new initial estimates
      popini <- sapply(popparams$name, function(j){
        if(j %in% names(fixed)){
          fixed[j] 
        }else if(j %in% log10){
          # If sampling over orders of magnitude, sample from log10 space
          10^(runif(1, log10(initialValue_ranges[j,1]), log10(initialValue_ranges[j,2])))
        }else{
          if((grepl("omega", j) | j == "a") & !is.null(set_variance)){
            set_variance # If a value is provided for the standard deviations of random effects, do not randomly sample it
          }else if (grepl("beta_", j) & !is.null(set_treatment)){
            set_treatment # If a value is provided for the treatment effects, do not randomly sample it
          }else{
            runif(1, initialValue_ranges[j,1], initialValue_ranges[j,2])  
          }
        }
      })
      
      # set sampled values as new initial estimates
      newpopparams <- popparams
      newpopparams$initialValue <- popini
      initialvalues <- rbind(initialvalues, newpopparams %>% mutate(run = i) %>% rename(parameter = name))
      
      if(!is.null(fixed)){
        newpopparams[newpopparams$name %in% names(fixed),"method"] <- "FIXED"
        newpopparams[newpopparams$name %in% names(fixed),"initialValue"] <- fixed
      }
      
    }else{
      
      igs <- predefined_initialvalues %>% filter(run == i) %>% as.data.frame()
      newpopparams <- popparams %>% rowwise %>%  
        # Only set initial conditions for parameters that have saved initial conditions. Otherwise, use default values
        mutate(initialValue = ifelse(name %in% igs$parameter, igs[igs$parameter == name, "initialValue"], initialValue)) %>% 
        as.data.frame() 
      initialvalues <- rbind(initialvalues, newpopparams %>% mutate(run = i) %>% rename(parameter = name))
    }
    
    setPopulationParameterInformation(newpopparams)
    
    # run the estimation
    tryCatch({
      runScenario()
      
      # store the estimates and s.e. in a table
      estimates <- as.data.frame(getEstimatedPopulationParameters())
      names(estimates) <- "estimate"
      rses <- getEstimatedStandardErrors()$stochasticApproximation$rse
      names(rses) <- getEstimatedStandardErrors()$stochasticApproximation$parameter
      rses <- as.data.frame(rses)
      estimates <- merge(estimates, rses, by = "row.names")
      estimates$run <- i
      estimates$ofv <- getEstimatedLogLikelihood()$importanceSampling[1]
      estimates$aic <- getEstimatedLogLikelihood()$importanceSampling[2]
      estimates$bic <- getEstimatedLogLikelihood()$importanceSampling[3]
      estimates$bicc <- getEstimatedLogLikelihood()$importanceSampling[4]
      names(estimates)[names(estimates) == "Row.names"] <- "parameter"
      tabestimates <- rbind(tabestimates, estimates)
      
      # store the iterations
      iters <- getSAEMiterations()$estimates %>% mutate(iteration = 1:nrow(.))
      iters$run <- i
      tabiters <- rbind(tabiters, iters)
      
      # Get wald-test p-values for treatment effects;
      if(record_pvals){
        p <- getTests()$wald
        if(!is.null(p)){
          treatment_pvals <- rbind(
            treatment_pvals, 
            p %>% mutate(run = i)
          )
        }
      }
      
      # Store Phase Iteration
      phase <- rbind(phase, getSAEMiterations()$iterationNumbers)
      
      # Get Individual parameters if specified
      if(record_indiv_params){
        indivParams <- rbind(
          indivParams,
          getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode %>% 
            mutate(run = i)
        )
        
        predictions <- rbind(
          predictions, 
          read_delim(here(out_folder, gsub(".mlxtran", "", project_name), "predictions.txt")) %>% 
            mutate(run = i)
        )
      }
      
      # Get shrinkage for individual parameters if specified
      if(record_shrinkage){
        shrinkage <- getEtaShrinkage()
        if(any(!getIndividualParameterModel()$variability$id)){
          shrinkage_df <- rbind(
            shrinkage_df,
            data.frame(
              run = i,
              parameter = shrinkage$conditionalMean$parameters,
              conditionalMean = unlist(unname(shrinkage$conditionalMean[1,-1])),
              conditionalMode = unlist(unname(shrinkage$conditionalMode[1,-1])),
              conditionalDist = unlist(unname(shrinkage$conditionalDist[1,-1]))
            )
          )
        }else{
          shrinkage_df <- rbind(
            shrinkage_df,
            data.frame(
              run = i,
              parameter = shrinkage$conditionalMean$parameters,
              conditionalMean = shrinkage$conditionalMean$shrinkage,
              conditionalMode = shrinkage$conditionalMode$shrinkage,
              conditionalDist = shrinkage$conditionalDist$shrinkage
            )
          )
        }
      }
      
      # Get condition number if specified:
      if(record_condition_number){
        FIM <- getCorrelationOfEstimates()$stochasticApproximation
        eigenvalues <- eigen(FIM)$values
        condition_number <- max(eigenvalues)/min(eigenvalues)
        condition_numbers <- rbind(condition_numbers, data.frame(run = i, condition_number = condition_number))  
      }
      
      
    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
    i <- i + 1
    
    # Stop when we reach the number of successful iterations
    if(length(unique(tabestimates$run)) >= nits) run_condition <- FALSE
    
  }
  
  out <- list(tabestimates = tabestimates,
              tabiters = tabiters,
              phase = phase,
              initialvalues = initialvalues)
  
  if(record_shrinkage){
    out$shrinkage <- shrinkage_df
  }
  
  if(record_pvals){
    out$pvals <- treatment_pvals
  }
  
  if(record_indiv_params){
    out$indivParams <- indivParams
    out$indivFit <- predictions
  }
  
  if(record_condition_number){
    out$condition_number <- condition_numbers
  }
  
  # Rerun the model with the initial values that gave the best fit:
  best_run <- tabestimates %>% filter(bicc == min(bicc)) %>% pull(run) %>% unique
  igs <- initialvalues %>% filter(run == best_run) %>% rename(name = parameter) %>% select(-run)
  if(!is.null(fixed)){
    igs[igs$name %in% names(fixed),"method"] <- "FIXED"
    igs[igs$name %in% names(fixed),"initialValue"] <- fixed
  }
  setPopulationParameterInformation(igs)
  runScenario()
  
  return(out)
}

profile_likelihood <- function(parameter, range, igs, nits, use_estimates_as_ig = F, linearize = F, record_predictions = F){
  
  values <- seq(range[1], range[2], length.out = nits)
  
  initialValues <- igs
  if("estimate" %in% names(igs)){
    initialValues <- igs %>% select(name = parameter, initialValue = estimate) %>%
      mutate(method = "MLE", initialValue = ifelse(grepl("pop", name), initialValue, pmax(initialValue, 1))) 
  }
  if("parameter" %in% names(igs)){
    names(initialValues)[names(initialValues) == "parameter"] <- "name"
  }
  
  
  initialValues[initialValues$name == parameter, "method"] <- "FIXED"
  
  if(linearize){
    scenario <- getScenario()
    scenario$linearization <- TRUE
    setScenario(scenario)
  }
  
  tabestimates <- NULL
  predictions <- NULL
  for(i in 1:nits){
    print(i)
    
    # Set initial conditions
    value <- values[i]
    
    initialValues[initialValues$name == parameter, "initialValue"] <- value
    setPopulationParameterInformation(initialValues)
    
    # print(getPopulationParameterInformation())
    
    tryCatch({
      # Model fitting
      runScenario()
      
      # Get standard errors and likelihoods:
      if(linearize){
        SEs <- getEstimatedStandardErrors()$linearization
        likelihoods <- getEstimatedLogLikelihood()$linearization
      }else{
        SEs <- getEstimatedStandardErrors()$stochasticApproximation
        likelihoods <- getEstimatedLogLikelihood()$importanceSampling
      }
      
      
      # store the estimates and s.e. in a table
      estimates <- as.data.frame(getEstimatedPopulationParameters())
      names(estimates) <- "estimate"
      rses <- SEs$rse
      names(rses) <- SEs$parameter  
      rses <- as.data.frame(rses)
      estimates <- merge(estimates, rses, by = "row.names")
      estimates$param <- value
      estimates$ofv <- likelihoods[1]
      estimates$aic <- likelihoods[2]
      estimates$bic <- likelihoods[3]
      estimates$bicc <- likelihoods[4]
      names(estimates)[names(estimates) == "Row.names"] <- "parameter"
      
      print(paste0(parameter, ": ", value, " bicc: ", round(unique(estimates$bicc))))
      
      # Combine results
      tabestimates <- rbind(tabestimates, estimates)
      
      if(record_predictions){
        vl_pred <- getEstimatedPredictions()$tot  %>% mutate(id = as.character(id)) 
        EBE <- getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode  %>% mutate(id = as.character(id))
        table_pred <- read_delim(here(out_folder, gsub(".mlxtran", "", project_name), "tables.txt")) %>% 
          select(id, time, contains("mode"))  %>% mutate(id = as.character(id))
        predictions_temp <- right_join(
          EBE,
          left_join(
            vl_pred, table_pred
          )
        ) %>% mutate(parameter = value) %>% 
          select(parameter, everything())
        
        names(predictions_temp)[names(predictions_temp) == "parameter"] <- parameter
        
        # Combine resuts
        predictions <- rbind(predictions, predictions_temp)
      }
      
      
      if(use_estimates_as_ig){
        initialValues <- estimates %>% select(name = parameter, initialValue = estimate) %>% 
          mutate(method = "MLE", initialValue = ifelse(grepl("pop", name), initialValue, pmax(initialValue, 1))) 
        initialValues <- rbind(initialValues, data.frame(name = parameter, initialValue = value, method = "FIXED"))
      }
      
      
    }, error=function(e){cat("ERROR for ", i, ": ",conditionMessage(e), "\n")})
    
  }
  
  names(tabestimates)[names(tabestimates) == "param"] <- parameter
  
  if(record_predictions){
    return(list(tabestimates = tabestimates, predictions = predictions))
  }else{
    return(tabestimates)
  }
  
}

# Calculate BICw ####
calculate_BICw <- function(tabestimates){
  # Calculate BICc weights 
  tabestimates %>% 
    distinct(run, bicc) %>% 
    mutate(
      t = exp(-(bicc - min(bicc))/2),
      w = t/sum(t)
    ) %>% 
    select(-t) %>% 
    arrange(desc(w)) %>% 
    mutate(cumprob = cumsum(w))
}


# Find narrower parameter ranges from IGS ####
get_narrow_ranges <- function(tabestimates, p_threshold = NULL, bicc_cutoff = NULL){
  
  if(!is.null(p_threshold)){
    # First, get range of parameters that give "best" model
    bicw <- calculate_BICw(tabestimates)
    
    runs <- bicw %>% 
      filter(cumprob <= p_threshold) %>% 
      pull(run)
    
    bic_ranges <- tabestimates %>% 
      filter(run %in% runs) %>% 
      group_by(parameter) %>% 
      summarise(min = min(estimate), max = max(estimate), best = estimate[bicc == min(bicc)])
  }
  
  if(!is.null(bicc_cutoff)){
    bic_ranges <- tabestimates %>% 
      filter(bicc <= bicc_cutoff) %>% 
      group_by(parameter) %>% 
      summarise(min = min(estimate), max = max(estimate), best = estimate[bicc == min(bicc)])
  }
  
  # Now, calculate range as 1/2*best estimate to 2*best estimate, unless it is narrower than the above range
  search_ranges <- bic_ranges %>% 
    mutate(
      min = ifelse(best/2 > min, min, best/2),
      max = ifelse(best*2 < max, max, best*2)
    ) 
  
  return(search_ranges)
  
}

# Function to detect outliers ####
outlier_detector <- function(value, upper_only = F){
  if(upper_only){
    (value > 1.5*IQR(value, na.rm = T) +  quantile(value, probs = 0.75, na.rm = T)) 
  }else{
    (value > 1.5*IQR(value, na.rm = T) +  quantile(value, probs = 0.75, na.rm = T)) |
      (value < quantile(value, probs = 0.25, na.rm = T) - 1.5*IQR(value, na.rm = T) )  
  }
  
}

# Model building function ####
## Can run SCM or COSSAC
## By default, BICc is used for the slection criteria (criteria = "BIC")

covariate_model_building <- function(projectFile, method, criteria = "BIC", covariates = "all", 
                                     linearization = FALSE, relationships = NULL, parameters = NULL, lrt_thresold = c(0.05,0.01)){
  
  if(dir.exists(here(gsub(".mlxtran","",projectFile), "ModelBuilding"))){ 
    stop("Model building directory already exists. Aborting to avoid overwriting previous results.")
  }else if (!file.exists(projectFile)){
    saveProject(projectFile)
  }
  
  # get the default covariate search settings ####
  set <- getModelBuildingSettings()
  
  # modify the settings
  set$strategy <- method # "cossac" or "scm"
  set$useSambaBeforeCossac <- FALSE # if true, corresponds to "covSAMBA-COSSAC" in the GUI
  set$criterion <- criteria # "BIC" or "LRT"
  set$useLin <- linearization
  set$threshold$lrt <- lrt_thresold # first element is for forward selection, second for backward elimination
  
  if(!is.null(parameters)){
    set$parameters <- parameters
  }
  
  if(!is.null(relationships)){
    set$relationships <- relationships
  }
  
  if(!any(covariates == "all")){
    set$covariates <- covariates
  }
  
  #####
  # Run search
  runModelBuilding(settings=set)
  model_building_results <- getModelBuildingResults()
  
  #####
  # Generate report for best model
  # Get best model
  best_model <- which(unlist(model_building_results %>% lapply("[[", "bestModel")))
  best_model_file <- sprintf("ModelBuilding_%03d", best_model)
  best_cov <- model_building_results[[best_model]]$individualModels
  
  final_projectFile <- here(gsub(".mlxtran", "", projectFile), "ModelBuilding", sprintf("ModelBuilding_%03d", best_model), paste0(sprintf("ModelBuilding_%03d", best_model), ".mlxtran") )
  loadProject(final_projectFile)
  
  if(any(as.character(unlist(getLaunchedTasks())) == "FALSE")){
    scenario <- getScenario()
    scenario$tasks[] <- TRUE
    scenario$linearization <- FALSE
    setScenario(scenario)
    
    runScenario()
    computeChartsData()
    
    saveProject( here(gsub(".mlxtran", "", projectFile), "final_model.mlxtran"))
  }
  
  # Generate report with results
  generate_report(project_file = here(gsub(".mlxtran", "", projectFile), "final_model.mlxtran"), 
                  output_file = "SCM_final_model_report.html")
  
}

#####
# Generate report for an alternate model

alt_mod <- function(out_folder, project_name, 
                    model_file, data_file, 
                    initialvalues, treatment_effects, 
                    headerTypes, error_model, 
                    exclude_random_effects = NULL, normal_dist = NULL, fixed = NULL){
  
  
  create_project(model_file, data_file, headerTypes, error_model, exclude_random_effects, treatment_effects, normal_dist, fixed)
  
  set_base(initialvalues, fixed = fixed)
  
  saveProject(here(out_folder, paste0(project_name, ".mlxtran")))
  
  runScenario()
  
  print(getEstimatedLogLikelihood())
  
  # Generate report with results
  generate_report(project_file = here(out_folder, paste0(project_name, ".mlxtran")))
  
}




# Data formatting for visualization ####
filter_results <- function(tabestimates, remove_estimate_outliers, remove_bicc_outliers){
  
  ## Filter data
  caption <- ""
  # Runs with bicc outliers
  if(remove_bicc_outliers){
    remove_runs <- tabestimates %>% distinct(run, bicc) %>% filter(outlier_detector(bicc, upper_only = T)) %>% pull(run)
    if(length(remove_runs) >= 1){
      caption <- paste0(caption, "Runs with outlier BICcs removed\n")  
    }else{
      cat("No runs with outlier BICc to removed\n")
      caption <- caption
    }
    
  }else{
    remove_runs <- c()  
  }
  
  # Remove Estimate outliers 
  if(remove_estimate_outliers){
    tabestimates <- tabestimates %>% group_by(parameter) %>% filter(!outlier_detector(estimate)) %>% ungroup()
    caption <- paste0(caption, "Estimate outliers removed")
  }
  
  # Remove runs with bicc outliers
  tabestimates <- tabestimates %>% filter(!run %in% remove_runs)
  
  return(list(tabestimates = tabestimates, caption = caption))
}


################################################################################
# Plotting functions
################################################################################

# Plot BIC Versus estimates ####

plot_bicc_v_estimate <- function(tabestimates, out_folder = "", project_name = "", narrow = TRUE,
                                 remove_bicc_outliers = FALSE, remove_estimate_outliers = FALSE, log10_x = TRUE, 
                                 exclude_treatment_effects = FALSE, ylim = NULL,
                                 file_addendum = "",  parameter = NULL, BICc_threshold = NULL,
                                 color = "BICc", color_label = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  if(!is.null(parameter)){
    tabestimates <- tabestimates %>% filter(param == parameter)
  }else if(exclude_treatment_effects){
    tabestimates <- tabestimates %>% filter(!grepl("beta_", param))
  }
  
  if(!is.null(ylim)){
    tabestimates <- tabestimates %>% filter(bicc <= ylim[2])
  }
  
  ## Generate plot
  popparams <- getPopulationParameterInformation()
  
  if(any(color == "BICc")){
    tabestimates <- tabestimates %>% mutate(color = bicc)
  }else if(length(color) > 1){
    tabestimates <- tabestimates %>% mutate(color = run %in% color)
  }
  
  
  if(narrow){
    estimates_vs_bicc <- tabestimates %>%
      mutate(param = factor(param, levels = popparams$name)) %>%
      ggplot(aes(estimate, bicc, shape = bicc == min(bicc), fill = color)) +
      geom_point(size = 2, alpha = 0.75) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), size = 2)
  }else{
    estimates_vs_bicc <- tabestimates %>%
      mutate(param = factor(param, levels = popparams$name)) %>%
      ggplot(aes(estimate, bicc,  fill = color)) +
      geom_point(size = 2, alpha = 0.75, shape = 21) +
      geom_hline(yintercept = BICc_threshold, lty = "dashed", alpha = 0.5) + 
      labs(caption = paste0("Dashed line indicates BICc below which fall the most probable models\n", caption))
  }
  estimates_vs_bicc <- estimates_vs_bicc + 
    facet_wrap(~param, scales = "free_x") + 
    labs(x = "Final estimate", y = "BICc", shape = "Best fit", 
         title = "Comparing BIC and parameter estimates",
         caption = caption) + 
    scale_shape_manual(values = c(21, 24))
  
  if(any(color == "BICc")){
    estimates_vs_bicc <- estimates_vs_bicc + scale_fill_viridis_c(direction = -1) + labs(fill = "BICc")
  }else if(length(color) > 1){
    estimates_vs_bicc <- estimates_vs_bicc + 
      scale_fill_manual(values = c(`TRUE` = "darkred", `FALSE` = "gray")) + labs(fill = color_label) + 
      guides(fill = guide_legend(override.aes = list(color = c("gray", "darkred"))))
  }
  
  
  if(any(grepl("beta_", tabestimates$param))) estimates_vs_bicc <- estimates_vs_bicc + labs(subtitle = ifelse(!exclude_treatment_effects, "", "treatment effects excluded"))
  if(!is.null(ylim)) estimates_vs_bicc <- estimates_vs_bicc + ylim(ylim) + labs(caption = paste0(caption, "\nBICc greater than ", max(ylim), " removed"))
  if(!is.null(parameter)){
    if(!grepl("beta_", parameter) & log10_x) estimates_vs_bicc <- estimates_vs_bicc + scale_x_log10()
  }else if(log10_x & exclude_treatment_effects){ 
    estimates_vs_bicc <- estimates_vs_bicc + scale_x_log10()
  }
  
  
  if(out_folder != ""){
    file_name <- ifelse(narrow, 
                        ifelse(!exclude_treatment_effects, "BIC_vs_estimates_narrow_range_treatment_effects", "BIC_vs_estimates_narrow_range"),
                        ifelse(!exclude_treatment_effects, "BIC_vs_estimates_treatment_effects", "BIC_vs_estimates_log_range"))
    file_name <- paste0(file_name, file_addendum, ".pdf")
    
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), file_name), estimates_vs_bicc, width = 10, height = 7)
    
  }
  
  return(estimates_vs_bicc)
  
}

#####
# Plot initial value and final estimate distributions

plot_distributions <- function(initialvalues, tabestimates, out_folder = "", project_name = "", parameter_set = "pop",
                               remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                               save = T, parameter = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  subtitles <- c("pop" = "Population parameters",
                 "omega" = "Standard deviations of the random effect",
                 "beta" = "Treatment effects",
                 "error" = "Error model parameters")
  subtitle <- paste(subtitles[parameter_set], collapse = " and ")
  
  if(is.null(parameter)){
    params <- unlist(lapply(parameter_set, function(z){ 
      p <- c()
      if(z == "error"){
        p <- c(p, unique(tabestimates$param)[unique(tabestimates$param) %in% c("a", "b", "c")])
      }else{
        if(z == "beta") z <- paste0(z, "_")
        p <- c(p, unique(tabestimates$param)[grepl(z, unique(tabestimates$param))])
      }
      return(p)
    }))
    
    df <- merge(initialvalues %>% select(param = parameter, run, initialValue), tabestimates) %>% 
      filter(param %in% params)
  }else{
    df <-  merge(initialvalues %>% select(param = parameter, run, initialValue), tabestimates) %>% 
      filter(param == parameter) 
  }
  
  popparams <- getPopulationParameterInformation()
  
  plot <- df %>% 
    pivot_longer(c(initialValue, estimate)) %>% 
    mutate(name = factor(name, levels = c("initialValue", "estimate"), labels = c("Initial Value", "Final Estimate")),
           param = factor(param, levels = popparams$name)) %>% 
    ggplot(aes(name, value)) + geom_violin(scale = "width") + 
    facet_wrap(~param, scales = "free_y") +
    geom_boxplot(position = position_dodge(width=0.9), width = 0.1, fill = "white", outlier.shape = NA) +
    geom_beeswarm(aes( fill = bicc, color = bicc, group = name), size = 2, alpha = 0.75, shape = 21) + 
    scale_fill_viridis_c(direction = -1) +
    scale_color_viridis_c(direction = -1) +
    geom_point(data = . %>% filter(bicc == min(bicc)), aes(fill = bicc), shape = 24, size = 3, alpha = 0.75) +
    labs(x = "", y = "Parameter Value", fill = "BICc",
         title = "Comparing initial value and final estimate distributions",
         caption = caption,
         subtitle = subtitle) + 
    guides(color = "none")
  
  if(!is.null(parameter)){
    if(grepl("beta_", parameter) | grepl("omega_", parameter)){
      ylim <- layer_scales(plot)$y$range$range
      if(min(ylim) > 0) ylim[1] <- 0
      if(min(ylim) < 0 & max(ylim) < 0) ylim[2] <- 0
      plot <- plot + coord_cartesian(ylim = ylim)
    }
  }
  
  if(save){
    file_name <- paste0("iv_estimate_dist_", paste(parameter_set, collapse = "_"), ".pdf")
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), file_name), plot, width = 10, height = 7)  
  }
  
  return(plot)
}

#####
# Plot final estimates versus initial values

plot_iv_vs_estimates <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                                 remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                                 parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  df <- merge(initialvalues %>% select(param = parameter, run, initialValue), tabestimates) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    arrange(param)
  
  params <- as.list(unique(df$param))
  
  if(is.null(parameter)){
    plots <- lapply(params, function(p){
      df %>% filter(param == p) %>% 
        ggplot(aes(initialValue, estimate, fill = bicc, shape = bicc == min(bicc))) + 
        geom_point(size = 2, alpha = 0.75) + 
        geom_point(data = . %>% filter(bicc == min(bicc)), size = 2) + 
        coord_obs_pred() +
        labs(title = p) + 
        theme(plot.title = element_text(hjust = 0.5))
    })
    
    color_limit <- range(df$bicc)
    
    combined_plot <- wrap_plots(plots)
    
    plot <- combined_plot + 
      plot_layout(guides = "collect", axis_titles = "collect") &
      labs( color = "BICc", shape = "Best fit") &
      plot_annotation(title = "Comparing initial value and final estimates", caption = caption) &
      theme(axis.title = element_blank()) &
      scale_shape_manual(values = c(21, 24)) &
      scale_fill_viridis_c(direction = -1, limits = color_limit)
    
    plot <- grid.arrange(patchworkGrob(plot), left = "Final Estimate", bottom = "Initial Value")
  }else{
    color_limit <- range(df$bicc)
    
    plot <- df %>% 
      filter(param == parameter) %>% 
      ggplot(aes(initialValue, estimate, fill = bicc, shape = bicc == min(bicc))) + 
      geom_point(size = 2, alpha = 0.75) + 
      geom_point(data = . %>% filter(bicc == min(bicc)), size = 2) + 
      coord_obs_pred() +
      labs(color = "BICc", y = "Final Estimate", x = "Initial Value", shape = "Best fit") + 
      theme(plot.title = element_text(hjust = 0.5)) +
      scale_shape_manual(values = c(21, 24)) + 
      scale_fill_viridis_c(direction = -1, , limits = color_limit) + 
      ggpmisc::stat_correlation(method = "pearson",  label.x = 0, label.y = 1,
                                use_label("R", "P", sep = "\n"),
                                output.type = "text")
  }
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "estimates_vs_ivs.pdf"), plot, width = 12, height = 10)  
  }
  
  return(plot)
}

#####
# Correlate final and initial values

cor_iv_estimates <- function(tabestimates, initialvalues, out_folder, project_name,
                             remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  df <- merge(initialvalues %>% select(parameter, run, initialValue), tabestimates) %>% 
    mutate(parameter = factor(parameter, levels = popparams$name)) %>% 
    group_by(parameter) %>% 
    summarise(cor = cor.test(initialValue, estimate)$estimate, pval = cor.test(initialValue, estimate)$p.value)
  
  breaks <- round(range(log10(df$pval)))
  breaks <- seq(breaks[1], breaks[2], by = 2)
  
  plot <- df %>% 
    ggplot(aes(cor, fct_rev(parameter), fill = log10(pval), shape = pval <= 0.05)) + 
    geom_vline(xintercept = 0, lty = "dashed") + 
    geom_point(size = 3) + 
    xlim(c(-1,1)) + 
    labs(x = "Correlation", y = "Parameter", caption = caption, 
         title = "Correlation between initial values and final estimates") +
    scale_fill_viridis_c(option = "plasma", direction = -1, breaks = breaks) + 
    scale_shape_manual(values = c(21, 24))
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "cor_iv_estimates.pdf"), plot, width = 7, height = 8)  
  }
  
  return(plot)
}


#####
# Plot parameters against eachother

correlate_estimates <- function(tabestimates, method = "spearman", out_folder = "", project_name = "",
                                remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE,  file_addendum = ""){
  
  # Filter tabestimates to remove estimate outliers and/or BICc outliers
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  correlation_df <- tabestimates %>% 
    select(param, estimate, run)  %>% 
    pivot_wider(names_from = param, values_from = estimate) %>% 
    select(-run) %>% 
    cor(use = "complete.obs", method = method) %>% # only keep runs that have no estimate outliers (if estimate outliers were removed)
    as.data.frame %>% 
    rownames_to_column("parameter2") %>% 
    pivot_longer(!parameter2, names_to = "parameter1") 
  
  correlation_heatmap <- correlation_df %>% 
    mutate(parameter2 = factor(parameter2, levels = popparams$name),
           parameter1 = factor(parameter1, levels = popparams$name)) %>%
    ggplot(aes(parameter1, parameter2, fill = value)) + 
    geom_tile() + 
    geom_text(data = . %>% filter(value > 0.5 | value < -0.5),
              aes(label = round(value, 2))) +
    scale_fill_gradient2(limits = c(-1,1)) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title = element_blank()) + 
    labs(fill = "correlation", title = str_interp("${str_to_title(method)} Correlation of SAEM estimates"), caption = caption)
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), paste0("SAEM_correlations", file_addendum, ".pdf")), correlation_heatmap, width =10, height = 10)  
  }
  
  return(correlation_heatmap)
  
}

### Plot population estiamates against eachother
compare_estimates <- function(tabestimates, out_folder = "", project_name = "", pop_only = TRUE,
                              remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, file_addendum = ""){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  if(pop_only){
    df <- tabestimates %>% 
      select(param, estimate, run, bicc)  %>% 
      filter(grepl("pop", param)) 
  }else{
    df <- tabestimates %>% 
      select(param, estimate, run, bicc)  
  }
  
  plot <- df %>% 
    pivot_wider(names_from = param, values_from = estimate) %>% 
    GGally::ggpairs(columns = 3:ncol(.), 
                    lower = list(continuous = function(data, mapping){
                      ggplot(data, mapping) + 
                        geom_point(aes(fill = bicc, shape = bicc == min(bicc)), alpha = 0.75) + 
                        geom_point(data = . %>% filter(bicc == min(bicc)), aes(fill = bicc, shape = bicc == min(bicc)), size = 2) + 
                        scale_fill_viridis_c(direction = -1) + 
                        scale_shape_manual(values = c(21, 24)) + 
                        labs(fill = "BICc", shape = "Best fit") 
                    } 
                    )
    ) + 
    labs(caption = caption) 
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), paste0("paired_SAEM_population_estimates", file_addendum, ".pdf")), plot, width = 7, height = 7)  
  }
  
  return(plot)
  
}

#####
# Calculate coefficient of variance

calc_CV <- function(tabestimates, initialvalues, out_folder = "", project_name = "",
                    remove_estimate_outliers = FALSE, remove_bicc_outliers = FALSE, 
                    parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  cv <- merge(
    tabestimates, 
    initialvalues %>% rename(param = parameter)
  )  %>% 
    group_by(param) %>% 
    summarise(`Final Estimate` = sd(estimate)/mean(estimate),
              `Initial Value` = sd(initialValue)/mean(initialValue)) 
  
  if(!is.null(parameter)) cv <- cv %>% filter(param == parameter)
  
  cv <- cv %>% pivot_longer(!param) %>% 
    mutate(param = factor(param, levels = popparams$name)) %>% 
    ggplot(aes(factor(param), abs(value), fill = name)) + 
    geom_col(position = position_dodge(),alpha = 0.7) + 
    geom_text(aes(label = abs(round(value, 2))), position = position_dodge(width = 0.91), hjust = -0.1) + 
    labs(x = "Parameter", y = "Absolute Value of Coefficient of Variation", caption = caption, title = "Comparing CV of initial values and final estimates") + 
    coord_flip() + 
    theme(legend.position.inside = c(0.9,0.9), legend.title = element_blank(), 
          legend.position = "inside", legend.background = element_rect(color = "white")) + 
    scale_fill_brewer(palette = "Set1") 
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "coefficient_of_variation.pdf"), cv, width = 8.5, height = 10)  
  }
  
  return(cv)
}


###
# Plot SAEM iterations

plot_SAEM <- function(tabiters, phase, tabestimates, out_folder = "", project_name = "",
                      remove_bicc_outliers = FALSE, filter = NULL, parameter = NULL){
  
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  # caption <- filtered_results$caption
  
  popparams <- getPopulationParameterInformation()
  
  
  if(is.null(parameter)){
    plotList <- list()
    i <- 1
    for (param in c(popparams$name, "convergenceIndicator")) {
      if(param != "convergenceIndicator"){
        if (popparams[popparams$name == param, ]$method == "FIXED") next
      }
      p <- param
      
      temp <- tabiters %>% 
        left_join(data.frame(run = unique(tabiters$run), phase = phase[,1])) %>% 
        filter(run %in% unique(tabestimates$run)) %>% 
        left_join(tabestimates %>% distinct(run, bicc)) 
      
      if(!is.null(filter)) temp <- temp %>% filter(iteration-phase >= filter)
      
      plotList[[i]] <- temp %>% 
        ggplot( aes(x = iteration-phase, y = !!sym(param))) +
        geom_line(aes(group = run, color = bicc)) +
        theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
        # geom_vline(xintercept = changePhase, color = 1:length(changePhase)) +
        labs(title = param, x = NULL, y = NULL) + 
        # scale_y_log10() + 
        # coord_cartesian(xlim = c(-200, 314)) + 
        geom_vline(xintercept = 0, lty = "dashed", color = "black") + 
        geom_hline(yintercept = tabestimates %>% filter(bicc == min(bicc), param == p) %>% pull(estimate),
                   alpha = 0.5) + 
        scale_color_viridis_c(direction = -1) 
      i <- i + 1
    }
    
    plot <- grid.arrange(grobs = plotList, ncol = 5)
  }else{
    
    
    temp <- tabiters %>% 
      left_join(data.frame(run = unique(tabiters$run), phase = phase[,1])) %>% 
      filter(run %in% unique(tabestimates$run)) %>% 
      left_join(tabestimates %>% distinct(run, bicc)) 
    
    
    if(!is.null(filter)) temp <- temp %>% filter(iteration-phase >= filter)
    
    plot <- temp %>% 
      ggplot( aes(x = iteration-phase, y = !!sym(parameter))) +
      geom_line(aes(group = run, color = bicc)) +
      theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
      labs( x = "Iteration - phase switch", y = "Estimate") + 
      geom_vline(xintercept = 0, lty = "dashed", color = "black") + 
      geom_hline(yintercept = tabestimates %>% filter(bicc == min(bicc), param == parameter) %>% pull(estimate),
                 alpha = 0.5) + 
      scale_color_viridis_c(direction = -1) 
    
    ## adjust y axis if SAEM portion is much more variable than estimates
    us <- range(temp[parameter])
    ue <- range(tabestimates %>% filter(param == parameter) %>% pull(estimate))
    lower <- min(c(us, ue))
    if(log10(abs(us[2]/ue[2])) > 1) plot <- plot + coord_cartesian(ylim = c(lower, ue[2]*1.2))
  }
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "SAEM_trajectories.pdf"),
           plot,
           width = 15, height = 7)
  }
  
  return(plot)
}

plot_estimates <- function(tabestimates, out_folder = "", project_name = "",
                           remove_bicc_outliers = FALSE, parameter = NULL){
  
  filtered_results <- filter_results(tabestimates, remove_estimate_outliers = F, remove_bicc_outliers)
  
  tabestimates <- filtered_results$tabestimates %>% rename(param = parameter)
  
  popparams <- getPopulationParameterInformation()
  
  # Get plot y-limit, excluding unusually large rses:
  ylim_range <-  tabestimates %>% 
    group_by(param) %>%
    filter(!outlier_detector(rses)) %>% 
    mutate(upper = estimate * (1+as.numeric(rses)/100), lower = estimate * (1-as.numeric(rses)/100)) %>% 
    summarise(lower = min(lower), upper = max(upper)) %>% 
    column_to_rownames("param")
  
  if(is.null(parameter)){
    # plot population parameters
    plotList <- list()
    i <- 1
    for (param in popparams$name) {
      if (popparams[popparams$name == param, ]$method == "FIXED") next
      estimates <- tabestimates[tabestimates$param == param, ]
      plotList[[i]] <- ggplot(estimates %>% filter(run %in% unique(tabestimates$run)), 
                              aes(x = run, y = estimate, color = bicc)) +
        geom_point() +
        geom_errorbar(aes(ymax = estimate * (1+as.numeric(rses)/100), 
                          ymin = estimate * (1-as.numeric(rses)/100))) +
        theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
        labs(title = param, x = NULL, y = NULL) +
        geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate),
                   alpha = 0.5) + 
        geom_hline(yintercept = estimates %>% pull(estimate) %>% median(),
                   alpha = 0.5, lty = "dashed") + 
        scale_color_viridis_c(direction = -1) 
      i <- i + 1
    }
    
    plot <- grid.arrange(grobs = plotList, ncol = 5)
  }else{
    estimates <- tabestimates[tabestimates$param == parameter, ]
    plot <- ggplot(estimates %>% filter(run %in% unique(tabestimates$run)), 
                   aes(x = run, y = estimate, color = bicc)) +
      geom_point() +
      geom_errorbar(aes(ymax = estimate * (1+as.numeric(rses)/100), 
                        ymin = estimate * (1-as.numeric(rses)/100))) +
      theme(legend.position = "none", plot.title = element_text(hjust = .5)) +
      labs(x = NULL, y = "Final Estimate") +
      geom_hline(yintercept = estimates %>% filter(bicc == min(bicc)) %>% pull(estimate), alpha = 0.5) + 
      geom_hline(yintercept = estimates %>% pull(estimate) %>% median(),
                 alpha = 0.5, lty = "dashed") +
      scale_color_viridis_c(direction = -1) +
      coord_cartesian(ylim = unname(unlist(ylim_range[parameter,])))
    
    if(grepl("beta_", parameter) | grepl("omega_", parameter)){
      ylim <- layer_scales(plot)$y$range$range
      if(min(ylim) > 0) ylim[1] <- 0
      if(min(ylim) < 0 & max(ylim) < 0) ylim[2] <- 0
      plot <- plot + coord_cartesian(ylim = ylim)
    }
  }
  
  if(out_folder != ""){
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), "Convergence_estimates.pdf"),
           plot,
           width = 15, height = 7)
  }
  
  return(plot)
}

## 
# Function to plot normalized estimates 
estimate_boxplots <- function(tabestimates, out_folder = "", project_name = "", file_addendum = NULL){
  CV <- tabestimates %>% group_by(parameter) %>% 
    summarise(CV = signif(sd(estimate)/mean(estimate)*100, 2))
  
  df <- tabestimates %>% 
    group_by(parameter) %>% 
    mutate(estimate = estimate/median(estimate)) 
  
  limits <- range(df$estimate) 
  breaks <- c(1,2,5,10)
  bound <- ceiling(max(1/limits[1], limits[2]))
  breaks <- breaks[breaks <= bound]
  breaks <- sort(round(unique(c(1/breaks, breaks, bound, 1/bound)), 2))
  
  parameters <- unique(df$parameter)
  
  plot <- df %>% 
    mutate(parameter = factor(parameter, 
                              levels = c(parameters[grepl("pop", parameters)],
                                         parameters[grepl("omega", parameters)],
                                         parameters[grepl("beta_", parameters)],
                                         parameters[ parameters == "a"]))) %>% 
    ggplot(aes(factor(parameter), estimate)) + 
    geom_boxplot() + 
    geom_text(data = CV, aes(y = Inf, label = CV), vjust = 1) +
    annotate("text", -Inf, Inf, label = "CV (%)", vjust = 1, hjust = 0) +
    scale_y_log10(breaks = breaks, limits = range(breaks)) + 
    labs(x = "Parameter", y = "Normalized Estimate") 
  
  
  if(out_folder != ""){
    file <- ifelse(is.null(file_addendum), "normalized_estimate_boxplots.pdf", str_interp("normalized_estimate_boxplots${file_addendum}.pdf"))
    ggsave(here(out_folder, gsub(".mlxtran", "", project_name), file), plot, width = 12, height = 7)
    
  }
  
  return(plot)
}

###
# Function for combined plots

parameter_plots <- function(parameter, remove_bicc_outliers, remove_estimate_outliers,
                            log10_tabestimates, tabestimates, tabiters,
                            phase, initialvalues){
  
  
  tabestimates <- filter_results(tabestimates, remove_estimate_outliers, remove_bicc_outliers)$tabestimates
  tabestimates <- tabestimates %>% rename(param = parameter) %>% filter(param == parameter) %>% rename(parameter = param)
  
  if(!is.null(log10_tabestimates)){
    log10_tabestimates <- filter_results(log10_tabestimates, remove_estimate_outliers, remove_bicc_outliers)$tabestimates
    log10_tabestimates <- log10_tabestimates %>% rename(param = parameter) %>% filter(param == parameter) %>% rename(parameter = param)
  }
  
  SAEM_traj <- plot_SAEM(
    tabiters, phase, tabestimates, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    parameter = parameter#, filter = -200
  ) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  SAEM_convergence <- plot_estimates(
    tabestimates,
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    parameter = parameter
  ) +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    labs(x = "Run")
  
  narrow_range <- tabestimates %>% rename(param = parameter) %>% filter(param == parameter) %>% 
    pull(estimate) %>% range
  
  
  if(!is.null(log10_tabestimates)){
    wide_plot <- plot_bicc_v_estimate(
      log10_tabestimates, narrow = FALSE, 
      # remove_bicc_outliers =  remove_bicc_outliers, 
      # remove_estimate_outliers = remove_estimate_outliers,
      parameter = parameter
    ) + 
      scale_fill_viridis_c(direction = -1, option = "magma",
                           name = "BICc wide\nparameter\nsearch")  +
      theme(plot.title = element_blank(),
            plot.subtitle = element_blank(),
            strip.background = element_blank(),
            strip.text = element_blank()) + 
      annotate("rect", xmin = narrow_range[1], xmax = narrow_range[2], 
               ymin = -Inf, ymax = Inf,
               alpha = 0.3, fill = "lightgray")
  }
  
  
  narrow_plot <- plot_bicc_v_estimate(
    tabestimates,  narrow = TRUE, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    log10_x = F,  parameter = parameter
  ) + 
    labs(fill = "BICc narrow\nparameter\nsearch") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  iv_vs_final <- plot_iv_vs_estimates(
    tabestimates, initialvalues, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    # remove_estimate_outliers = remove_estimate_outliers,
    parameter = parameter
  ) +
    labs(fill = "BICc narrow\nparameter\nsearch") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank())
  
  # Get coefficients of variation
  cv <- merge(
    tabestimates, 
    initialvalues
  )  %>% 
    group_by(parameter) %>% summarise(`Final Estimate` = sd(estimate)/mean(estimate),
                                      `Initial Value` = sd(initialValue)/mean(initialValue)) %>% 
    pivot_longer(!parameter)
  
  distributions <- plot_distributions(
    initialvalues, tabestimates, 
    # remove_estimate_outliers = remove_estimate_outliers, 
    # remove_bicc_outliers =  remove_bicc_outliers, 
    parameter = parameter
  ) +
    labs(fill = "BICc narrow\nparameter\nsearch", 
         x = "Set") +
    theme(plot.title = element_blank(),
          plot.subtitle = element_blank(),
          strip.background = element_blank(),
          strip.text = element_blank()) + 
    geom_text(data = cv, aes(label = round(value, 2)*100, y = Inf), vjust = 1) + 
    annotate(geom = "text", label = "CV (%):", y = Inf, x = -Inf, vjust = 1, hjust = 0) + 
    scale_y_continuous(expand = expansion(mult = c(0, .1)))
  
  cv <- calc_CV(
    tabestimates, initialvalues, 
    # remove_estimate_outliers = remove_estimate_outliers, 
    # remove_bicc_outliers =  remove_bicc_outliers,
    parameter = parameter
  )
  
  caption <- ""
  if(remove_bicc_outliers) caption <- paste0(caption, "Runs with BICc outliers removed\n")
  if(remove_estimate_outliers) caption <- paste0(caption, "Estimate outliers removed")
  
  if(!is.null(log10_tabestimates)){
    final_plot <- wide_plot +  iv_vs_final +  SAEM_traj + narrow_plot + distributions + SAEM_convergence +
      plot_layout(nrow = 3, byrow = F, guides = "collect", widths = c(1, 1.1)) &
      plot_annotation(title = parameter, theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = caption) &
      labs(caption = "") 
  }else{
    final_plot <- plot_spacer() + iv_vs_final +  SAEM_traj + narrow_plot + distributions + SAEM_convergence +
      plot_layout(nrow = 3, byrow = F, guides = "collect", widths = c(1, 1.1)) &
      plot_annotation(title = parameter, theme = theme(plot.title = element_text(hjust = 0.5)),
                      caption = caption) &
      labs(caption = "") 
  }
  
  return(final_plot)
  
}

# Plot results of stepwise model building:

plot_model_building <- function(mb,  ...){
  mb_summary <- mb %>% 
    map_df(function(x) data.frame(
      param = x$tested[1], 
      treatment = x$tested[2],
      BICc = x$BICc, 
      LL = x$LL, 
      n_covar = sum(x$individualModels))
    )
  
  
  mb_summary_no_back_elim <- mb_summary %>% 
    mutate(id = 1:nrow(.), t = c(0,diff(n_covar))) %>% 
    filter(t >= 0) 
  
  best_mod_covar <- filter(mb_summary_no_back_elim, BICc == min(BICc))$n_covar
  
  mb_summary_no_back_elim <- mb_summary_no_back_elim %>% 
    filter(id <= max(id[n_covar == best_mod_covar])) %>% 
    select(-id, -t)
  
  thresholds <- mb_summary_no_back_elim %>% 
    group_by(n_covar) %>% 
    summarise(threshold = BICc[BICc == min(BICc)]) %>% 
    ungroup %>% 
    mutate(n_covar = n_covar + 1)
  
  df <- mb_summary_no_back_elim %>% 
    filter(n_covar != 0) %>% 
    left_join(thresholds)
  
  p <- df %>% 
    filter(param != "") %>% 
    group_by(n_covar) %>% mutate(best = BICc == min(BICc)) %>% 
    ungroup %>% group_by(param, treatment) %>% mutate(final_mod = any(best)) %>% 
    ggplot(aes(BICc, treatment, color = param)) + 
    geom_rect(data = . %>% filter(best), 
              aes(xmin = BICc, xmax = BICc + 2, ymin = -Inf, ymax = Inf), 
              fill = "grey", alpha = 0.25, color = NA) +
    geom_rect(data = . %>% filter(best), 
              aes(xmin = BICc, xmax = BICc + 5, ymin = -Inf, ymax = Inf), 
              fill = "grey", alpha = 0.15, color = NA) +
    geom_vline(aes(xintercept = threshold), alpha = 0.5) +
    geom_point(size = 2,  aes(shape = best), alpha = 0.25) +
    geom_point(data = . %>% filter(final_mod), size = 2,  aes(shape = best)) +
    ggrepel::geom_text_repel(
      data = . %>% filter(BICc <= threshold), aes(label = param), 
      hjust = 0, min.segment.length = 0, direction = "y",
      show.legend = F
    ) +
    facet_wrap(~n_covar, labeller = "label_both", dir = "v", ...) + 
    ggsci::scale_color_lancet() +
    guides(color = guide_legend(override.aes = list(alpha = 1)), 
           shape = guide_legend(override.aes = list(alpha = 1))) + 
    labs(color = "Parameter", shape = "Best covariate", caption = "Forward selection only") + 
    theme(panel.border = element_rect(fill = NA)) 
  
  # Add custom breaks while preserving the default range
  return(list(plot = p, summary = left_join(mb_summary, thresholds)))
  
}

# Plot ranges from each step ####
plot_parameter_ranges <- function(initialValue_ranges, log10_tabestimates, tabestimates = NULL, p_threshold = NULL, bicc_cutoff = NULL){
  iv_range <- initialValue_ranges %>% t() %>% as.data.frame %>% set_names(c("min", "max")) %>% rownames_to_column("parameter")
  log10_range <- log10_tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate)) %>% ungroup
  narrow_iv_range <- get_narrow_ranges(log10_tabestimates, p_threshold, bicc_cutoff) %>% select(-best)
  df <- rbind(
    iv_range %>% mutate(step = "1", range = "Initial Guess", parameter = paste(parameter, "pop", sep = "_")),
    log10_range %>% mutate(step = "2", range = "Estimate Range"),
    narrow_iv_range %>% mutate(step = "2", range = "Initial Guess")
  ) %>% 
    filter(grepl("pop", parameter))
  if(!is.null(tabestimates)){
    narrow_range <- tabestimates %>% group_by(parameter) %>% summarise(min = min(estimate), max = max(estimate)) %>% ungroup
    best <- tabestimates %>% filter(bicc == min(bicc), grepl("pop", parameter)) %>% select(parameter, estimate)
    df <- rbind(df, 
                narrow_range %>% mutate(step = "3", range = "Estimate Range")
    ) %>% 
      filter(grepl("pop", parameter))
  }
  
  plot <- ggplot(df, aes(x= step, ymin = min, ymax = max, color = range)) + 
    geom_linerange(linewidth = 2) + 
    # geom_linerange(data = filter(df, space == "Narrow Range"), linewidth = 2) + 
    facet_wrap(~parameter, scales = "free_y") + 
    scale_y_log10() + 
    scale_color_manual(values = c("gray", "darkred"))
  
  if(!is.null(tabestimates)){
    plot <- plot + 
      geom_point(data = best, aes(x = "3", y = estimate), 
                 color = "darkgreen", inherit.aes = F, size = 2) +
      geom_text(data = best, aes(x = "3", y = estimate, label = signif(estimate, 2)), 
                color = "darkgreen", inherit.aes = F, fontface = 2, vjust = 0, hjust = 1)
  }
  
  plot
  
}

# Plotting RSES ####
plot_vs_rses <- function(tabestimates, out_folder = "", project_name = "",  y = "bicc", file_addendum = NULL){
  names(tabestimates)[names(tabestimates) == y] <- "y"
  scales <- ifelse(y == "bicc", "free_x", "free")
  
  plot <- tabestimates %>% mutate(rses = abs(rses)) %>% 
    arrange(rses) %>% 
    mutate(color = fct_inorder(case_when(rses < 50~"<50", 
                                         rses >= 50 & rses < 100~ "50-100",
                                         rses > 100~ ">100"))) %>% 
    ggplot(aes(abs(rses), y, color = color)) +  
    geom_point() + 
    scale_x_log10() + 
    facet_wrap(~parameter, scales = scales) + 
    labs(x = "RSE (%)", y = y, color = "RSE (%)") +
    scale_color_brewer(palette = "Oranges")
  
  if(!out_folder == ""){
    file <- ifelse(is.null(file_addendum), str_interp("rses_vs_${y}.pdf"), str_interp("rses_vs_${y}${file_addendum}.pdf"))
    ggsave(here(out_folder, str_remove(project_name, "\\.mlxtran$"), file), 
           plot, width = 10, height = 7)
  }
  
  plot
}

# Shortcut functions ####
# functions to scale bicc color
color_bicc <- function() scale_color_viridis_c(direction = -1)
fill_bicc <- function() scale_fill_viridis_c(direction = -1)

