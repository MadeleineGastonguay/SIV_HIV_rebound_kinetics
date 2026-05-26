# Load libraries
library(here)
library(ggh4x)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))
source(here("scripts", "forward_simulation_functions.R"))
theme_set(theme_minimal())

# Project setup ####
study <- "study4"
project_name <- "SCM_tlr7.mlxtran"
model_file <- here("model", "2regime_model_delay.txt")

## Results folder
out_folder <- here("results_eclipse", study)

## Data file
data_file <- here("data", "source", "VACABTLR7_DATA.txt")
data <- read_delim(data_file, delim = "\t")

## Parameters to fix
fixed <- c("lambda_pop" = 100, "p_pop" = 0.9)


loadProject(here(out_folder, project_name))

## Get into the model building
mb <- getModelBuildingResults()

out <- plot_model_building(mb)

## Barplot figure with results:
SCM_summary <- out$summary
SCM_summary_figure_study4 <- SCM_summary %>%
  filter(n_covar != 0) %>%
  mutate(id = 1:nrow(.)) %>%
  filter(id < min(id[n_covar > 6])) %>%
  group_by(n_covar) %>%
  mutate(shaded = min(BICc) + 2,
         n_covar_label = factor(as.character(n_covar), levels = 1:7,
                                labels = paste(c("1st", "2nd", "3rd", "4th", "5th", "6th", "7th"), "treatment effect added")),
         label = ifelse(BICc < threshold, param, ""),
         param = factor(param, levels = c("lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         label = factor(label, levels = c("","lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("","lambda", "beta", "t['a']", "m", "p", "N['p']"))) %>%
  mutate(param = fct_rev(param),
         treatment = factor(treatment, levels = rev(c("tlr7", "vac", "ab", "intTLR7ABVAC")),
                            labels = rev(c("TLR7", "Vaccine", "Antibody", "Antibody*Vaccine")))
  ) %>%
  ggplot(aes(BICc - threshold, treatment, fill = param)) +
  geom_rect(data = ~filter(., BICc == min(BICc)),
            aes(xmin =BICc-threshold , xmax = shaded-threshold, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.5, inherit.aes = F) +
  geom_col(position = position_dodge(), color = NA) +
  geom_text(aes(BICc - threshold, label = label, color = param),
            parse = T, hjust = 1.1, position = position_dodge(width = 0.9), show.legend = F) +
  facet_wrap( ~ n_covar_label, ncol = 2, dir = "v") +
  # scale_y_discrete(labels = c("vac" = "Vaccine", "tlr7" = "TLR7", "inTLR7VAC" = "TLR7 * Vaccine")) +
  # scale_y_discrete(labels = rev(c( "Vaccine",  "TLR7", "TLR7 * Vaccine"))) +
  labs(x = "Change in BICc", y = "Treatment effect added to model", fill = "Parameter") +
  scale_fill_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0")), labels = function(x) parse(text = x)) +
  scale_color_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0"))) +
  guides(fill = guide_legend(reverse = T)) +
  # ggnewscale::new_scale_color() +
  # geom_col(data = . %>% filter(n_covar <= 3) %>% mutate(col = BICc == min(BICc)),
  #          aes(group = param, color = col), fill = NA, position = position_dodge(),  lty = "dashed", show.legend = F) +
  # scale_color_manual(values = c("#FFFFFF00", "black")) +
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        strip.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave(here(out_folder, str_remove(project_name, ".mlxtran"), "model_building_bar_plot.pdf"),
       SCM_summary_figure_study4 + coord_cartesian(xlim = c(-20,20)), width = 15, height = 15)

## Plot population trends at each step:
included_effects <- SCM_summary %>% 
  mutate(mod = 1:nrow(.)) %>% 
  filter(mod < min(mod[n_covar > 6])) %>% 
  group_by(n_covar) %>% 
  filter(BICc == min(BICc))

stepwise_preds <- stepwise_indivParam <- stepwise_popParam <-NULL
for(mod in included_effects$mod){
  num <- sprintf("%03d", mod)
  model_folder <- here(out_folder, "SCM_tlr7", "ModelBuilding", str_interp("ModelBuilding_${num}"), str_interp("ModelBuilding_${num}"))
  predictions <- read_csv(here(model_folder, "predictions.txt")) %>% select(-contains("mode"))
  indivParam <- read_csv(here(model_folder, "IndividualParameters", "estimatedIndividualParameters.txt")) %>% select(-contains("mode"))
  popParam <- read_csv(here(model_folder, "populationParameters.txt")) %>% select(-contains("mode")) %>% select(parameter, value)
  stepwise_preds <- rbind(stepwise_preds, predictions %>% mutate(model = num))
  stepwise_indivParam <- rbind(stepwise_indivParam, indivParam %>% mutate(model = num))
  stepwise_popParam <- rbind(stepwise_popParam, popParam %>% mutate(model = num))
}

write_csv(stepwise_preds, here(out_folder, "SCM_tlr7", "stepwise_predictions.csv"))
write_csv(stepwise_indivParam, here(out_folder, "SCM_tlr7", "stepwise_individualParameters.csv"))
write_csv(stepwise_popParam, here(out_folder, "SCM_tlr7", "stepwise_populationParameters.csv"))

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  ggplot(aes(time, popPred, group = interaction(id, model), color = treatment_effect)) + 
  geom_point(data = ~filter(.x, model == "001"), aes(time, vl), alpha = 0.5, color = "black") +
  geom_line(linewidth = 1) + 
  facet_wrap(ab ~ vac, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added")

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  filter(model %in% c("046", "060", "103")) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  ggplot(aes(time, popPred, group = interaction(id, model), color = factor(n_covar))) + 
  geom_point(data = ~filter(.x, model == "046"), aes(time, vl), alpha = 0.5, color = "black") +
  geom_line(linewidth = 1) + 
  facet_wrap(ab ~ vac, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
       title = "Predictions for candidate models")

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  filter(model != "001") %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(ab, vac) %>% 
  mutate(Treatment = fct_inorder(case_when(ab == 1 & vac == 1 ~ "TLR7 + Ab + Vac", 
                                           ab == 1 & vac == 0 ~ "TLR7 + Ab",
                                           ab == 0 & vac == 1 ~ "TLR7 + Vac",
                                           ab == 0 & vac == 0 ~ "Control"))) %>% 
  ggplot(aes(time, popPred, group = interaction(id, model), color = treatment_effect)) + 
  geom_point(aes(time, vl), alpha = 0.5, color = "black") +
  geom_line(linewidth = 1) + 
  facet_grid(treatment_effect ~ Treatment) + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added")

stepwise_preds %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, indivPred_SAEM, group = interaction(id, model), color = treatment_effect)) + 
  geom_hline(yintercept = 1.04, lty = "dashed") +
  geom_point(data = ~filter(.x, model == "001"), aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_wrap(~id, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6))

## Simulate with R: ####

fixed <- list(
  N_E = 1e4, f = 0.9, dTar = 0.05, dI_macaque = 0.4, dI_human = 1, dP = 1e-3, 
  dE = 1, k = 50000, c = 23, final_time=200, 
  I1_human = 4 * 10 ^ (-6), I1_monkey = 4 * 10 ^ (-5),  C =.7, rho=10,
  detection_threshold=log10(50), v_one_infected = 0.0724, reactivation_rate = 1e-3, 
  LR_decay = 2.7/365, alpha_factor = 1e-8, fano_factor = 10, whole_body_conversion=3*10^4, washout_time = 3,
  q = 1
)


N_E <- fixed$N_E #this is about what z approaches, in the dynamics between immune system and virus
#so we estimate at 1% of 10^6 CD4 cells per ml
f <- fixed$f
dTar <- fixed$dTar #target
# dI <- fixed$dI_macaque #infected
dP <- fixed$dP #long lived
dE <- fixed$dE #short memory
# k <- fixed$k_human #also k/c
c <- fixed$c
# I1 <- fixed$I1_human
C <- fixed$C
final_time <- fixed$final_time
rho <- fixed$rho
washout_time <- fixed$washout_time
k <- fixed$k
q <- fixed$q
I1 <- fixed$I1_monkey
dI <- fixed$dI_macaque


temp_df2 <- stepwise_indivParam %>% select(model, id, contains("SAEM"), vac, ab) %>% 
  mutate(treatment = case_when(ab == 0 & vac == 0 ~ "Control",
                               ab == 1 & vac == 0 ~ "TLR7 + Ab",
                               ab == 0 & vac == 1 ~ "TLR7 + Vac",
                               ab == 1 & vac == 1 ~ "TLR7 + Vac + Ab")) %>% 
  select(model, treatment, id, contains("SAEM")) %>% 
  rename_with(~gsub("_SAEM", "", .x), everything())

individual_trajectories <- data.frame()
for(i in 1:nrow(temp_df2)){
  print(i)
  params <- unlist(temp_df2[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T, delay = T)[[2]]
  sim <- cbind(temp_df2[i,1:3], sim)
  individual_trajectories <- rbind(individual_trajectories, sim )
}

save(individual_trajectories, file = here(out_folder, "SCM_tlr7", "stepwise_indiv_sims.RData"))


individual_trajectories %>% rename(Treatment = treatment) %>% 
  # filter(model %in% c("046", "060", "103")) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10), group = model, color = treatment_effect)) + 
  geom_hline(yintercept = 1.04, lty = "dashed") +
  geom_point(data = data, aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_wrap(~id, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6), xlim = c(0,168)) + 
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
       fill = "Treatment group") + 
  scale_fill_manual(values = colors, labels= c("Control", "TLR7+Antibody", "TLR7+Vaccine", "TLR7+Antibody+Vaccine"))


# individual_trajectories %>% rename(Treatment = treatment) %>% 
#   left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
#   left_join(distinct(data, id, ab, vac)) %>% 
#   mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
#   arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
#   ggplot(aes(time, V/log(10), group = model, color = treatment_effect)) + 
#   geom_hline(yintercept = 1.04, lty = "dashed") +
#   geom_point(data = data, aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
#   geom_line(linewidth = 1, alpha = 0.5) + 
#   facet_wrap(~id, labeller = "label_both") + 
#   scale_color_viridis_d( begin = 0.2, end = 1) +
#   coord_cartesian(ylim = c(-1, 6), xlim = c(0,168)) + 
#   labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
#        fill = "Treatment group") + 
#   scale_fill_manual(values = colors, labels= c("Control", "TLR7+Antibody", "TLR7+Vaccine", "TLR7+Antibody+Vaccine"))

# diff_id <- c(7132, 7133, 7137, 7138, 7142, 7145, 7147, 7149, 7152)
diff_id <- c(7132, 7133,  7138, 7142, 7145, 7147, 7149, 7151, 7152)

individual_trajectories %>% rename(Treatment = treatment) %>%
  filter(id %in% diff_id) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10), group = model, color = treatment_effect)) + 
  geom_hline(yintercept = 1.1, lty = "dashed") +
  geom_point(data = data %>% filter(id %in% diff_id), aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1) + 
  facet_wrap(~id, labeller = "label_both") + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6), xlim = c(0,168)) + 
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
       fill = "Treatment group") + 
  scale_fill_manual(values = colors, labels= c("TLR7+Antibody", "TLR7+Antibody+Vaccine"))

individual_trajectories %>% rename(Treatment = treatment) %>%
  filter(id %in% diff_id) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  ggplot(aes(time, V/log(10), group = model, color = treatment_effect)) + 
  geom_hline(yintercept = 1.1, lty = "dashed") +
  geom_point(data = data %>% filter(id %in% diff_id), aes(time, vl, fill = interaction(ab, vac)), shape = 21, inherit.aes = F) +
  geom_line(linewidth = 1) + 
  facet_grid(treatment_effect~id) + 
  scale_color_viridis_d( begin = 0.2, end = 1) +
  coord_cartesian(ylim = c(-1, 6), xlim = c(0,168)) + 
  labs(x = "Time since ART interruption (days)", y = "Viral load", color = "Treatment effect added",
       fill = "Treatment group") + 
  scale_fill_manual(values = colors, labels= c("TLR7+Antibody", "TLR7+Antibody+Vaccine"))



## Population trajectories
temp_df <- stepwise_popParam %>% 
  mutate(treatment = case_when(intTLR7AB == 0 & intTLR7VAC == 0 ~ "Control",
                               intTLR7AB == 1 & intTLR7VAC == 0 ~ "TLR7 + Ab",
                               intTLR7AB == 0 & intTLR7VAC == 1 ~ "TLR7 + Vac",
                               intTLR7AB == 1 & intTLR7VAC == 1 ~ "TLR7 + Vac + Ab")) %>% 
  select(model, treatment, id, contains("SAEM")) %>% 
  rename_with(~gsub("_SAEM", "", .x), everything())

individual_trajectories <- data.frame()
for(i in 1:nrow(temp_df2)){
  print(i)
  params <- unlist(temp_df2[i,c("lambda", "betas", "ta", "m", "p", "NPs")])
  params["betas"] <- params["betas"]*1e-7
  params["NPs"] <- params["NPs"]*1e5
  sim <- sim_traj_final(params, return_full = T, delay = T)[[2]]
  sim <- cbind(temp_df2[i,1:3], sim)
  individual_trajectories <- rbind(individual_trajectories, sim )
}


## Claculate difference between observation and prediction to see where we are over esitmating and estimating

stepwise_preds %>% 
  filter(id != 7122) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  mutate(error = ifelse(vl <= 1.1 & indivPred_SAEM <= 1.1, 0, indivPred_SAEM-vl)) %>% 
  ggplot(aes(time, error, color = treatment_effect)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(linewidth = 1, alpha = 0.5) + 
  facet_wrap(~id) +
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(y = "Difference between observed and predicted viral load", color = "Treatment effect added")

stepwise_preds %>% 
  filter(id %in% diff_id) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  arrange(vac, ab) %>% mutate(id = fct_inorder(as.character(id))) %>% 
  mutate(error = ifelse(vl <= 1.1 & indivPred_SAEM <= 1.1, 0, indivPred_SAEM-vl)) %>% 
  ggplot(aes(time, error, color = treatment_effect)) +
  geom_hline(aes(yintercept = 0)) +
  geom_line(linewidth = 1) + 
  facet_wrap(~id) +
  scale_color_viridis_d( begin = 0.2, end = 1) +
  labs(y = "Difference between observed and predicted viral load", color = "Treatment effect added")

## How do random effects change?

stepwise_popParam %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  filter(grepl("omega", parameter)) %>% 
  ggplot(aes(factor(n_covar), value, fill = treatment_effect)) + 
  geom_col() + 
  facet_wrap(~parameter) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added")

# how do treatment effects change?
stepwise_popParam %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  filter(grepl("beta_", parameter)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  ggplot(aes(factor(n_covar), value, fill = treatment_effect)) + 
  geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added")

# how do population parameters change?
stepwise_popParam %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, treatment, sep = ":")))) %>%
  filter(grepl("pop", parameter)) %>% 
  mutate(parameter = fct_inorder(parameter)) %>% 
  ggplot(aes(factor(n_covar), value, fill = treatment_effect)) + 
  geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~parameter, scales = "free_y", nrow = 2) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added")

## Compared to prior model, how does next model change?

compare_df <- merge(
  stepwise_preds %>% 
    select(id, model, time, popPred) %>% 
    filter(model != "001") %>% 
    left_join(mutate(included_effects, model = sprintf("%03d", mod))),
  stepwise_preds %>% 
    select(id, model, time, popPred_ref = popPred) %>% 
    filter(model != "122") %>% 
    left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>% 
    mutate(n_covar = n_covar + 1),
  by = c("n_covar", "id", "time")
) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, gsub("intTLR7", "int", treatment.x), sep = ":")))) 

compare_df %>%
  left_join(distinct(data, id, ab, vac)) %>% 
  # mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, treatment.x, sep = ":")))) %>%
  ggplot(aes(time, popPred-popPred_ref, color = treatment_effect, group = interaction(id, model.x))) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_line(linewidth = 1) + 
  facet_wrap(ab~vac, labeller = "label_both") + 
  scale_color_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Time", y = "Change from model at previous step", color = "Treatment effect added")

compare_df %>%
  left_join(distinct(data, id, ab, vac)) %>% 
  # mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, treatment.x, sep = ":")))) %>%
  arrange(ab, vac) %>% 
  mutate(Treatment = fct_inorder(case_when(ab == 1 & vac == 1 ~ "TLR7 + Ab + Vac", 
                                           ab == 1 & vac == 0 ~ "TLR7 + Ab",
                                           ab == 0 & vac == 1 ~ "TLR7 + Vac",
                                           ab == 0 & vac == 0 ~ "Control"))) %>% 
  ggplot(aes(time, popPred-popPred_ref, color = treatment_effect, group = interaction(id, model.x))) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_line(linewidth = 1) + 
  facet_grid(treatment_effect~Treatment) + 
  scale_color_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Time", y = "Change from model at previous step", color = "Treatment effect added")

compare_df %>%
  filter(n_covar <=4 ) %>% 
  left_join(distinct(data, id, ab, vac)) %>% 
  mutate(treatment_effect = fct_inorder(ifelse(param.x == "", "No treatment effects", paste(param.x, treatment.x, sep = ":")))) %>%
  ggplot(aes(time, popPred-popPred_ref, color = treatment_effect, group = interaction(id, model.x))) + 
  geom_hline(yintercept = 0, lty = "dashed") + 
  geom_line(linewidth = 1) + 
  facet_wrap(ab~vac, labeller = "label_both") + 
  scale_color_viridis_d(begin = 0.2, end = 1) + 
  labs(x = "Time", y = "Change from model at previous step", color = "Treatment effect added")


add_treatments <- function(df){
  df <- df %>% left_join(distinct(data, id, ab, vac)) %>%   
    arrange(ab, vac) %>% 
    mutate(Treatment = fct_inorder(case_when(ab == 1 & vac == 1 ~ "TLR7 + Ab + Vac", 
                                             ab == 1 & vac == 0 ~ "TLR7 + Ab",
                                             ab == 0 & vac == 1 ~ "TLR7 + Vac",
                                             ab == 0 & vac == 0 ~ "Control"))) 
  
  return(df)
}

stepwise_popParam_plotting <- stepwise_popParam %>% 
  rbind(data.frame(parameter = "beta_lambda_ab_1", value = 0, model = "001")) %>% 
  left_join(mutate(included_effects, model = sprintf("%03d", mod))) %>%  
  mutate(treatment_effect = fct_inorder(ifelse(param == "", "No treatment effects", paste(param, gsub("intTLR7", "int", treatment), sep = ":")))) 

study4_explore_mb <- compare_df %>% 
  add_treatments() %>%
  group_by(ab, vac) %>% 
  filter(id == min(id)) %>% 
  mutate(treatment_effect = factor(treatment_effect, levels = c("No treatment effects", levels(compare_df$treatment_effect)))) %>% 
  ggplot(aes(time, group = id)) + 
  geom_point(data= add_treatments(data), aes(time, vl), inherit.aes = FALSE, alpha = 0.25) + 
  geom_ribbon(aes(ymin = pmin(popPred_ref, popPred), ymax = pmax(popPred_ref, popPred)), alpha = 0.8, fill = "gray", color = "black") +
  # geom_ribbon(data = ~filter(.x, popPred_ref < popPred), aes(ymin = pmin(popPred_ref, popPred), ymax = pmax(popPred_ref, popPred), fill = "Increase"),  alpha = 0.75) + 
  # geom_ribbon(data = ~filter(.x, popPred_ref > popPred), aes(ymin = pmin(popPred_ref, popPred), ymax = pmax(popPred_ref, popPred), fill = "Decrease"), alpha = 0.75) + 
  # geom_line(aes(y = popPred_ref, color = "ref"), linewidth = 1) + 
  geom_line(data = stepwise_preds %>% 
              mutate(treatment_effect = factor("No treatment effects", levels = c("No treatment effects", levels(compare_df$treatment_effect)))) %>%  
              add_treatments() %>% filter(model == "001"), aes(time, popPred, color = treatment_effect), linewidth = 1.5) + 
  geom_line(aes(y = popPred, color = treatment_effect), linewidth = 1.5) + 
  facet_grid(treatment_effect~Treatment) +
  scale_color_viridis_d(begin = 0.2, end = 1)  +
  labs(x = "Time since ART interruption (days)", y = "Viral load") +
  # scale_fill_manual(values = c("Decrease" = "darkred", "Increase" = "darkgreen"))
  
  stepwise_popParam_plotting %>%
  filter(grepl("beta_", parameter)) %>% 
  arrange(treatment_effect, mod) %>%  
  mutate(parameter = fct_rev(fct_inorder(gsub("_intTLR7", ":", gsub("_1", "", gsub("beta_", "", parameter)))))) %>% 
  ggplot(aes(parameter, value, fill = treatment_effect)) + 
  # geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~treatment_effect, ncol = 1) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  # labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added") + 
  theme(legend.position = "none") + 
  plot_layout(guides = "collect", widths = c(3,1)) + coord_flip() +
  labs(x = "treatment effect", y = "Treatment effect") + 
  geom_hline(yintercept = 0) +
  
  stepwise_popParam_plotting %>% 
  filter(grepl("pop", parameter)) %>% 
  mutate(parameter = fct_rev(fct_inorder(gsub("_pop", "",parameter)))) %>% 
  group_by(parameter) %>% mutate(value = value/max(value)) %>% 
  ggplot(aes(parameter, value, fill = treatment_effect)) + 
  # geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~treatment_effect, scales = "free_y", ncol = 1) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  # labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added") + 
  theme(legend.position = "none") + 
  labs(x = "parameter", y = "Population parameter\nestimate normalized by maximum") +
  coord_flip() +
  
  stepwise_popParam_plotting %>% 
  filter(grepl("omega", parameter)) %>% 
  mutate(parameter = fct_rev(fct_inorder(gsub("omega_", "",parameter)))) %>% 
  ggplot(aes(parameter, value, fill = treatment_effect)) + 
  # geom_hline(yintercept = 0) +
  geom_col() + 
  facet_wrap(~treatment_effect, scales = "free_y", ncol = 1) + 
  scale_fill_viridis_d(begin = 0.2, end = 1) + 
  # labs(x = "Number of treatment effects", y = "Standard deviation of random effect", fill = "Treatment effect added") + 
  theme(legend.position = "none") + 
  labs(x = "parameter", y = "Random Effect") +
  coord_flip() +
  
  plot_layout(guides = "collect", widths = c(4,1,1,1), nrow = 1) &
  theme(legend.position = "none") &
  plot_annotation(tag_levels = "a")


ggsave(here(out_folder, str_remove(project_name, ".mlxtran"), "stepwise_model_building_plots.pdf"),
       study4_explore_mb, width = 15, height = 10)

stepwise_preds %>% 
  add_treatments() %>%
  group_by(ab, vac) %>% 
  filter( model == "039") %>% 
  # mutate(treatment_effect = factor(treatment_effect, levels = c("No treatment effects", levels(compare_df$treatment_effect)))) %>% 
  ggplot(aes(time, group = id)) + 
  # geom_ribbon(aes(ymin = pmin(popPred_ref, popPred), ymax = pmax(popPred_ref, popPred)), alpha = 0.75, fill = "gray") +
  # geom_ribbon(data = ~filter(.x, popPred_ref < popPred), aes(ymin = pmin(popPred_ref, popPred), ymax = pmax(popPred_ref, popPred), fill = "Increase"),  alpha = 0.75) + 
  # geom_ribbon(data = ~filter(.x, popPred_ref > popPred), aes(ymin = pmin(popPred_ref, popPred), ymax = pmax(popPred_ref, popPred), fill = "Decrease"), alpha = 0.75) + 
  geom_line(aes(y = popPred, color = Treatment, lty = Treatment), linewidth = 1, alpha = 0.5) + 
  # geom_line(aes(y = popPred_ref, color = "ref"), linewidth = 1) + 
  # geom_line(data = stepwise_preds %>% 
  #             mutate(treatment_effect = factor("No treatment effects", levels = c("No treatment effects", levels(compare_df$treatment_effect)))) %>%  
  #             add_treatments() %>% filter(model == "001"), aes(time, popPred, color = treatment_effect), linewidth = 1.5) + 
  facet_grid(Treatment~model) 
# scale_color_viridis_d(begin = 0.2, end = 1)  

stepwise_indivParam %>% 
  select(id, contains("SAEM"), starts_with("int"), "model") %>% 
  pivot_longer(contains("SAEM")) %>% 
  left_join(distinct(data,id, tlr7,ab,vac)) %>% 
  ggplot(aes(interaction(tlr7,ab,vac), value, color = model)) + 
  geom_violin() + 
  # geom_beeswarm() + 
  facet_wrap(~name, scales = "free_y") + 
  scale_y_log10()

## Get group parameters for final model: ####

popparams <- stepwise_popParam %>% 
  filter(model == 122, grepl("pop", parameter)) %>% 
  mutate(parameter = gsub("_pop", "", parameter)) %>% 
  select(parameter, control = value)

trt_effects <- stepwise_popParam %>% 
  filter(model == 122) %>% 
  filter(grepl("beta_", parameter)) %>% 
  separate(parameter, sep = "_", into = c("t", "parameter", "group", "t2")) %>% 
  mutate(group = group) %>% 
  select(-t, -t2, -model) %>% 
  pivot_wider(names_from = group, values_from = value, values_fill = 0) %>% 
  group_by(parameter) %>% 
  summarise(TLR7AB = ab + tlr7, TLR7VAC = tlr7+vac, TLR7ABVAC = ab + tlr7 + vac + intTLR7ABVAC) 

trt_effects <- rbind(trt_effects, 
                     data.frame(
                       parameter = popparams$parameter[!popparams$parameter %in% trt_effects$parameter], 
                       TLR7AB = 0, TLR7VAC = 0, TLR7ABVAC = 0
                       )
                     )

group_param <- left_join(popparams, pivot_longer(trt_effects, !parameter)) %>%
  mutate(value = control*exp(value)) %>% 
  pivot_wider()

write_csv(group_param, here(out_folder, str_remove(project_name, ".mlxtran"), "final_model_parameters_by_group.csv"))


