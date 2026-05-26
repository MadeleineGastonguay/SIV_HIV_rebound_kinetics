## Compare correlates of immune response and reservoir size to EBEs

library(tidyverse)
library(here)
library(ggsci)
library(survival)
library(patchwork)
library(ggbeeswarm)

source(here("scripts", "functions.R"))

theme_set(theme_bw() + 
            theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
                  panel.border = element_rect(colour = "black"), aspect.ratio=1) )
colors <- c("Control" = "black", "TLR7" = "#027dbb", "Vac" = "#fcba12", "Ab" ="#fa61c2", "TLR7 + Vac" = "#7ab148",  "TLR7 + Ab" = "#7f48b1", "TLR7 + Vac + Ab" =  "#a96e7c")


## Read in clinical covariates ####
# Total SIV DNA in CD4+ T memory cells for studies 1-3, intact provirus for study 4
reservoir_covar <- read_csv(here("data", "derived", "reservoir_correlates.csv")) %>% 
  mutate(Treatment = factor(Treatment, levels = names(colors)))
# Immune breadth against Gag, Pol,and Env (and total) available for Study 2 & 4
immune_covar <- read_csv(here("data", "derived", "immune_correlates.csv")) %>% 
  mutate(Treatment = factor(Treatment, levels = names(colors)))
# Size of immune response against Gag, Pol, and Env (and total) available for Study 4
immune_size_covar <- read_csv(here("data", "derived", "immune_size_correlates.csv"))


## Correlate Individual parameters with Reservoir size ####
## Start with EBEs: #####
EBEs <- read_csv(here("results_eclipse", "summary_single_fit_with_study4", "SIV_SHIV_individual_parameters.csv")) 
indivParams <- EBEs %>% select(data, id, contains("mode")) %>% mutate(Study = as.numeric(str_remove(data, "study")))

reservoir_combined_df <- reservoir_covar %>% 
  merge(indivParams) %>% 
  filter(!is.na(SIV_DNA)) %>% 
  mutate(log_SIV_DNA = log10(SIV_DNA), log_ta = log10(ta_mode))

# Plot SIV DNA versus ta in multiple celltypes
reservoir_plot <- reservoir_combined_df %>% 
  filter(Study != 4) %>% 
  ggplot(aes(ta_mode, SIV_DNA, color = factor(Study))) + 
  # geom_point(size = 2.5, aes(shape = SIV_DNA > 3)) +
  geom_point(size = 2.5, aes(shape = Censored == 0)) +
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  labs(x = "ta EBE", y = "SIV DNA (per million cells)", color = "Study") + 
  scale_color_bmj() +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1))  +
  guides(shape = "none")

# Plot intact provirus vs ta for study 4
reservoir_plot_s4 <- reservoir_combined_df %>% 
  filter(Study == 4) %>% 
  ggplot(aes(ta_mode, SIV_DNA, color = factor(Study))) + 
  # geom_point(size = 2.5, aes(shape = SIV_DNA > 3)) +
  geom_point(size = 2.5, aes(shape = Censored == 0)) +
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  labs(x = "ta EBE", y = "Intact provirus (per million cells)") + 
  scale_color_bmj() +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1))  +
  guides(shape = "none")

# Calculate correlation for each cell type with tobit regression accounting for censored values
# wrapper function for tobit model:
fit_tobit_ta <- function(df, celltype, censor_value = 3){
  df <- filter(df, cell_type == celltype)
  fit <- survreg(Surv(log_SIV_DNA, log_SIV_DNA > log10(censor_value), type = "left") ~ log_ta, 
                 data = df, dist = "gaussian")
  prediction <- tibble(
    Study = df$Study, id = df$id, cell_type = celltype, pval = summary(fit)$table[2,"p"],
    predicted = predict(fit, newdata = df)
  )
  return(list(fit = fit, prediction = prediction))
}


# Create a dataframe with predictions:
# wrapper function for predictions from tobit model:

# Studies 1-3
reservoir_predictions <- rbind(
  fit_tobit_ta(reservoir_combined_df, "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df, "GMMC")$prediction,
  fit_tobit_ta(reservoir_combined_df, "LNMC")$prediction
)


reservoir_plot <- reservoir_combined_df %>% 
  left_join(reservoir_predictions) %>% 
  filter(Study != 4) %>%
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
  ggplot(aes(ta_mode, SIV_DNA, color = factor(Study))) + 
  geom_point(size = 3, aes(shape = censored)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), color = "black", linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data =. %>% group_by(cell_type) %>% slice(1), 
            aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
            hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_log10() + scale_x_log10() + 
  # facet_wrap(~cell_type, ncol = 1) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_bmj() + 
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(x = "ta EBE", y = "SIV DNA (per million cells)", 
       color = "Study", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

reservoir_predictions_by_study <- rbind(
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 1), "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 1), "GMMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 1), "LNMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 2), "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 2), "LNMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 3), "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 3), "LNMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 4), "PBMC")$prediction
)

pvals <- reservoir_predictions_by_study %>% 
  filter(Study != 4) %>% 
  group_by(cell_type, Study) %>% slice(1) %>% 
  mutate(y = ifelse(Study == 1, 1e4, ifelse(Study == 2, 1e3, ifelse(Study == 3, 1e2, 1e1))),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC")))

reservoir_combined_df %>% 
  left_join(reservoir_predictions_by_study) %>% 
  filter(Study != 4) %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
  ggplot(aes(ta_mode, SIV_DNA, color = factor(Study))) + 
  geom_point(size = 3, aes(shape = censored)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data = pvals, 
            aes(Inf, y, label = paste0("p = ", signif(pval, 2))),
            hjust = 1, vjust = 1.5, size = 4, show.legend = F) +
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  scale_color_bmj() + 
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(x = "ta EBE", y = "SIV DNA (per million cells)", 
       color = "Study", title = "Correlates of Reservoir Size", shape = NULL) +
  theme(plot.title = element_text(hjust = 0.5))

study2_reservoir_correlates <- reservoir_combined_df %>% 
  left_join(reservoir_predictions_by_study) %>% 
  filter(Study == 2, cell_type == "PBMC") %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
  ggplot(aes(ta_mode, SIV_DNA)) + 
  geom_point(size = 3, aes(shape = censored, color = Treatment)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data = pvals %>% filter(Study == 2, cell_type == "PBMC"), 
            aes(Inf, y, label = paste0("p = ", signif(pval, 1))),
            hjust = 1.1, vjust = 1.5, size = 4, show.legend = F) +
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  scale_color_manual(values = (colors)) +
  labs(x = bquote("t"[a]~"EBE"), y = "SIV DNA (per million cells)", 
       color = NULL, title = "Correlates of Reservoir Size", shape = NULL) +
  theme(plot.title = element_text(hjust = 0.5)) 

study1_reservoir_correlates <- reservoir_combined_df %>% 
  left_join(reservoir_predictions_by_study) %>% 
  filter(Study == 1, cell_type == "PBMC") %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC")),
         Treatment = factor(Treatment, levels = c("Control", "TLR7"))) %>%
  ggplot(aes(ta_mode, SIV_DNA)) + 
  geom_point(size = 3, aes(shape = censored, color = Treatment)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data = pvals %>% filter(Study == 1, cell_type == "PBMC"), 
            aes(Inf, y, label = paste0("p = ", signif(pval, 1))),
            hjust = 1.1, vjust = 1.5, size = 4, show.legend = F) +
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  scale_color_manual(values = unname(colors)) +
  labs(x = bquote("t"[a]~"EBE"), y = "SIV DNA (per million cells)", 
       color = NULL, title = "Correlates of Reservoir Size", shape = NULL) +
  theme(plot.title = element_text(hjust = 0.5)) 

## Too many censored values to get a good fit
study3_reservoir_correlates <- reservoir_combined_df %>% 
  left_join(reservoir_predictions_by_study) %>% 
  filter(Study == 3, cell_type == "PBMC") %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC")),
         Treatment = factor(Treatment, levels = c("Control", "TLR7", "Ab", "Ab+TLR7"))) %>%
  ggplot(aes(ta_mode, SIV_DNA)) + 
  geom_point(size = 3, aes(shape = censored, color = Treatment)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data = pvals %>% filter(Study == 3, cell_type == "PBMC"), 
            aes(Inf, y, label = paste0("p = ", signif(pval, 1))),
            hjust = 1.1, vjust = 1.5, size = 4, show.legend = F) +
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  scale_color_manual(values = unname(colors)) +
  labs(x = bquote("t"[a]~"EBE"), y = "SIV DNA (per million cells)", 
       color = NULL, title = "Correlates of Reservoir Size", shape = NULL) +
  theme(plot.title = element_text(hjust = 0.5)) 


# Compare reservoir data between treatment groups:
treatment_v_reservoir <- reservoir_combined_df %>% 
  filter(cell_type == "PBMC", Study != 3) %>%
  mutate(Study = ifelse(Study == 1, "Study 1: SIV-TLR7 (Late ART)", "Study 2: SIV-TLR7+VAC (Early ART)")) %>% 
  ggplot(aes(Treatment, SIV_DNA)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 2.5, aes(color = factor(Study), shape = SIV_DNA == 3)) + 
  scale_y_log10() +
  scale_color_manual(values = ggsci::pal_bmj()(3)[1:2]) +
  labs(color = "Study", y = "SIV DNA (per million cells)") + 
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) + 
  guides(shape = 'none') +
  facet_wrap(~Study, scales = "free_x") + 
  theme(legend.position = "none") 

ggsave(here("results_eclipse", "summary_single_fit", "reservoir_correlates_treatment.pdf"),
       treatment_v_reservoir, width = 8, height = 5)

treatment_v_reservoir_study2 <- reservoir_combined_df %>% 
  filter(Study == 2, cell_type == "PBMC") %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "TLR7", "Vac", "Vac+TLR7"))) %>%
  ggplot(aes(Treatment, SIV_DNA, color = Treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 2.5, aes(shape = SIV_DNA == 3)) + 
  scale_y_log10() +
  scale_color_manual(values = ggsci::pal_bmj()(3)[1:2]) +
  labs(color = "Study", y = "SIV DNA (per million cells)") + 
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) + 
  scale_color_manual(values = unname(colors)) +
  guides(shape = 'none') +
  facet_wrap(~cell_type, scales = "free_x") + 
  theme(legend.position = "none") 

study2_res_plot <- study2_reservoir_correlates + theme(plot.title = element_blank()) + 
  treatment_v_reservoir_study2 + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(strip.text = element_blank(), axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

ggsave(here("presentation figures/study2_reservoir_correlates.pdf"), study2_res_plot, 
       width = 10, height = 5)

treatment_v_reservoir_study1 <- reservoir_combined_df %>% 
  filter(Study == 1, cell_type == "PBMC") %>%
  mutate(Treatment = factor(Treatment, levels = c("Control", "TLR7", "Vac", "Vac+TLR7"))) %>%
  ggplot(aes(Treatment, SIV_DNA, color = Treatment)) + 
  geom_boxplot(outlier.shape = NA) +
  geom_beeswarm(size = 2.5, aes(shape = SIV_DNA == 3)) + 
  scale_y_log10() +
  scale_color_manual(values = ggsci::pal_bmj()(3)[1:2]) +
  labs(color = "Study", y = "SIV DNA (per million cells)") + 
  scale_shape_manual(values = c("TRUE" = 1, "FALSE" = 16)) + 
  scale_color_manual(values = unname(colors)) +
  guides(shape = 'none') +
  facet_wrap(~cell_type, scales = "free_x") + 
  theme(legend.position = "none") 



study1_res_plot <- study1_reservoir_correlates + theme(plot.title = element_blank()) + 
  treatment_v_reservoir_study1 + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(strip.text = element_blank(), axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

ggsave(here("presentation figures/study1_reservoir_correlates.pdf"), study1_res_plot, 
       width = 10, height = 5)

# Compare to samples from the conditional distribution, not just EBE: ####

study1_conditional <- read_csv(here("results_eclipse", "study1", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))
study2_conditional <- read_csv(here("results_eclipse", "study2", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))
study3_conditional <- read_csv(here("results_eclipse", "study3", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))
study4_conditional <- read_csv(here("results_eclipse", "study4", "best_4.4", "IndividualParameters", "simulatedIndividualParameters.txt"))


all_conditional_param <- rbind(
  study1_conditional %>% mutate(Study = 1) %>% select(-tlr7),
  study2_conditional %>% mutate(Study = 2) %>% select(-c(intTLR7VAC, tlr7, vac)),
  study3_conditional %>% mutate(Study = 3) %>% select(-c(inttlr7ab, tlr7, ab)),
  study4_conditional %>% mutate(Study = 4) %>% select(-c(ab, tlr7, vac))
)

reservoir_conditional_df <- left_join(
  all_conditional_param %>% filter(Study != 4), 
  reservoir_covar %>% mutate(log_SIV_DNA = log10(SIV_DNA)),
  relationship = "many-to-many"
  ) %>% 
  filter(!is.na(SIV_DNA)) %>% 
  mutate(log_ta = log10(ta))

reservoir_conditional_df %>% 
  ggplot(aes(ta, SIV_DNA, color = factor(Study))) +
  geom_point() +
  facet_wrap(~cell_type) + 
  scale_x_log10() + scale_y_log10()

reservoir_conditional_df_S4 <- left_join(
  all_conditional_param %>% filter(Study == 4), 
  reservoir_covar %>% mutate(log_SIV_DNA = log10(SIV_DNA)),
  relationship = "many-to-many"
) %>% 
  filter(!is.na(SIV_DNA)) %>% 
  mutate(log_ta = log10(ta))

reservoir_conditional_df_S4 %>% 
  ggplot(aes(ta, SIV_DNA, color = factor(Study))) +
  geom_point() +
  facet_wrap(~cell_type) + 
  scale_x_log10() + scale_y_log10()

reservoir_conditional_predictions <- merge(
  reservoir_conditional_df,
  rbind(
    fit_tobit_ta(reservoir_conditional_df, "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df, "LNMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df, "GMMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.)))
  )
)

reservoir_plot_conditional <- reservoir_conditional_predictions %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
  ggplot(aes(log_ta, SIV_DNA, color = factor(Study))) + 
  geom_point(size = 3, aes(shape = censored)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), color = "black", linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data =. %>% group_by(cell_type) %>% slice(1), 
            aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
            hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_log10() + #scale_x_log10() + 
  scale_x_continuous(labels = scales::label_math(10^.x)) +
  coord_cartesian(ylim = c(1, 10000)) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_bmj() + 
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(x = bquote(t['a']), y = "SIV DNA (per million cells)", 
       color = "Study", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1)

reservoir_conditional_predictions_S4 <- merge(
  reservoir_conditional_df_S4,
  rbind(
    fit_tobit_ta(reservoir_conditional_df_S4, "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.)))
  )
)

reservoir_plot_conditional_S4 <- reservoir_conditional_predictions_S4 %>% 
  mutate(censored = ifelse(Censored == 1, "Censored", "Uncensored"),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
  ggplot(aes(ta, SIV_DNA, color = factor(Study))) + 
  geom_point(size = 3, aes(shape = censored)) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), color = "black", linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data =. %>% group_by(cell_type) %>% slice(1), 
            aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
            hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_log10() + #scale_x_log10() + 
  # scale_x_continuous(labels = scales::label_math(10^.x)) +
  coord_cartesian(ylim = c(1, 10000)) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_bmj() + 
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(x = bquote(t['a']), y = "SIV DNA (per million cells)", 
       color = "Study", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1)

# Break it up by study:

reservoir_conditional_predictions_by_study <- merge(
  rbind(reservoir_conditional_df, reservoir_conditional_df_S4),
  rbind(
    fit_tobit_ta(reservoir_conditional_df %>% filter(Study == 1), "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df %>% filter(Study == 2), "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df %>% filter(Study == 3), "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df_S4 %>% filter(Study == 4), "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df %>% filter(Study == 1), "LNMC" )$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df %>% filter(Study == 2), "LNMC" )$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df %>% filter(Study == 3), "LNMC" )$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.))),
    fit_tobit_ta(reservoir_conditional_df, "GMMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.)))
  )
)

make_reservoir_plot_conditional_by_study <- function(prediction_df, study, nrow = 1, breaks = NULL, scientific_x = FALSE){
  
  if(study == 3){
    treatment_levels <- c("Control", "TLR7", "Ab", "TLR7 + Ab")
  }else if(study == 4){
    treatment_levels <- c("Control", "TLR7 + Ab", "TLR7 + Vac", "TLR7 + Vac + Ab")
  }else{
    treatment_levels <- c("Control", "TLR7", "Vac", "TLR7 + Vac")
  }
  
  df <- prediction_df %>% 
    filter(Study == study) %>%
    mutate(#censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored"),
      cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC")),
      Treatment = factor(Treatment, levels = treatment_levels))
  
  if(scientific_x){
    plot <- ggplot(df, aes(log_ta , SIV_DNA, color = Treatment))  
  }else{
    plot <- ggplot(df, aes(ta , SIV_DNA, color = Treatment))
  }
  plot <- plot + geom_point(size = 3, aes(shape = factor(Censored))) + 
    # add tobit regression line
    geom_line(aes(y = 10^predicted), color = "black", linewidth = 1, alpha = 0.7) +
    # add pvalue
    geom_text(data =. %>% group_by(cell_type, Study) %>% slice(1), 
              aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
              hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
    scale_y_log10() + #scale_x_log10() + 
    facet_wrap(~cell_type, nrow = nrow) +
    scale_color_manual(values = (colors)) +
    scale_shape_manual(values = c("1" = 1, "0" = 16)) + 
    labs(x = bquote('Individual'~t['a']~'estimates'), y = "SIV DNA (per million cells)", 
         color = "Treatment", title = "Correlates of Reservoir Size", shape = NULL) +
    theme_bw() + 
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black"), aspect.ratio=1)
  
  if(!is.null(breaks) ){
    if(scientific_x){
      plot <- plot + scale_x_continuous(labels = scales::label_math(10^.x), breaks = breaks)   
    }else{
      plot <- plot + scale_x_continuous(breaks = breaks, trans = "log10") 
    }
  }else{
    if(scientific_x){
      plot <- plot + scale_x_continuous(labels = scales::label_math(10^.x))   
    }else{
      plot <- plot + scale_x_log10()
    }
  }
  plot
}


## Save figures with conditional individual parameters for each study:

make_treatment_v_reservoir <- function(reservoir_combined_df, study){
  
  # if(study == 3){
  #   treatment_levels <- c("Control", "TLR7", "Ab", "TLR7+Ab")
  # }else{
  #   treatment_levels <- c("Control", "TLR7", "Vac", "Vac+TLR7")
  # }
  # 
  
  treatment_levels <- c("Control", "TLR7", "Vac", "TLR7 + Vac", "Ab", "TLR7 + Ab", "TLR7 + Vac + Ab")
  
  reservoir_combined_df %>% 
    filter(Study %in% study) %>%
    mutate(Treatment = factor(Treatment, levels = treatment_levels),
           cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
    ggplot(aes(Treatment, SIV_DNA, color = Treatment)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(size = 2.5, aes(shape = factor(Censored))) + 
    scale_y_log10() +
    scale_color_manual(values = ggsci::pal_bmj()(3)[1:2]) +
    labs(color = "Study", y = "SIV DNA (per million cells)") + 
    scale_shape_manual(values = c("1" = 1, "0" = 16)) + 
    scale_color_manual(values = c((colors), (colors)[3:4])) +
    guides(shape = 'none') +
    facet_wrap(~cell_type, scales = "free_x", ncol = 1) + 
    theme(legend.position = "none") 
}


study1_res_plot_conditional <- make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 1, scientific_x = T) +
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 1) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "LNMC"), 1, scientific_x = T) +
  theme(axis.title.x = element_blank()) +
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="LNMC"), 1) + theme(axis.title.y = element_blank() , axis.title.x = element_blank()) + 
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "GMMC"), 1, scientific_x = T) +
  theme(axis.title.y = element_blank()) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="GMMC"), 1) + theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(#strip.text = element_blank(), 
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 10),
    plot.title = element_blank()) &
  coord_cartesian(ylim = c(1,10000)) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

# ggsave(here("paper_figures/study1_reservoir_correlates_conditional_indivParam.pdf"), study1_res_plot_conditional, 
#        width = 7, height = 8)

study2_res_plot_conditional <- make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 2) +
  theme(axis.title.x = element_blank()) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 2) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "LNMC"), 2 ) +
  # theme(axis.title.x = element_blank()) +
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="LNMC"), 2) + 
  theme(axis.title.y = element_blank()) + # , axis.title.x = element_blank()) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(#strip.text = element_blank(), 
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 10),
    plot.title = element_blank()) &
  coord_cartesian(ylim = c(1,10000)) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

# ggsave(here("paper_figures/study2_reservoir_correlates_conditional_indivParam.pdf"), study2_res_plot_conditional, 
#        width = 7, height = 7)

study3_res_plot_conditional <- make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 3) +
  theme(axis.title.x = element_blank()) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 3) + 
  theme(axis.title.y = element_blank(), axis.title.x = element_blank()) + 
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "LNMC"), 3) +
  # theme(axis.title.x = element_blank()) +
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="LNMC"), 3) + 
  theme(axis.title.y = element_blank()) + # , axis.title.x = element_blank()) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(#strip.text = element_blank(), 
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 10),
    plot.title = element_blank()) &
  coord_cartesian(ylim = c(0.005,300)) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

# ggsave(here("paper_figures/study3_reservoir_correlates_conditional_indivParam.pdf"), study3_res_plot_conditional, 
#        width = 7, height = 7)

study4_res_plot_conditional <- make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 4) +
  labs(y = "Intact provirus (per million cells)") +
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 4) + 
  theme(axis.title.y = element_blank()) + 
  plot_layout(ncol = 2, guides = "collect") & 
  theme(#strip.text = element_blank(), 
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 10),
    plot.title = element_blank()) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

## Make a simpler figure:

simple_reservoir_correlations <- reservoir_plot_conditional + labs(x = bquote("Individual"~t['a']~'estimate'))  +
  make_treatment_v_reservoir(reservoir_covar %>% mutate(study_lab = paste("Study",Study,"PBMC")) %>% filter(cell_type=="PBMC"), 1:3) + 
  facet_wrap(~study_lab, scales = "free_x") + 
  plot_layout(ncol = 1) &
  plot_annotation(tag_levels = "a") &
  theme(plot.title = element_blank())

ggsave(here("paper_figures", "reservoir_correlates_conditional_simple.pdf"), simple_reservoir_correlations, 
       width = 10, height = 7)

## Make a simpler figure:

simple_reservoir_correlations <- make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 1, scientific_x = T) +
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 2) +
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 3) +
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 4) 

# Correlate immune data with Np, p, and m EBEs ####

immune_combined_df <- immune_covar %>% filter(Study == 2) %>% 
  pivot_longer(6:11, names_to = "observation_type", values_to = "observation") %>% 
  left_join(EBEs %>% select(id, NPs_mode, m_mode, p_mode, tlr7:intTLR7VAC)) %>% 
  filter(!is.na(observation)) %>% 
  mutate(log_NPs = log10(NPs_mode))

# Look at observed data
immune_combined_df %>% 
  filter(grepl("breadth", observation_type)) %>%
  ggplot(aes(observation, Timing)) + 
  geom_boxplot() + geom_jitter() + 
  facet_wrap(~observation_type, scales = "free_x")

immune_combined_df %>% 
  filter(grepl("CD", observation_type)) %>%
  ggplot(aes(observation, Timing)) + 
  geom_boxplot() + geom_jitter() + 
  facet_wrap(~observation_type, scales = "free_x", ncol = 1) + 
  labs(x = "% IFN-gamma")

(immune_plot <- immune_combined_df %>% 
    filter(observation_type == "breadth_Total") %>%
    # filter(grepl("breadth", observation_type)) %>%
    mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
           observation_type = paste(gsub("breadth_", "", observation_type), "Breadth"),
           Study = factor(Study, levels = c(1,2,3,4))) %>%
    ggplot(aes(NPs_mode*1e5, observation, color = factor(Study))) + 
    geom_point(size = 2.5) + 
    scale_x_log10() + 
    # facet_wrap(observation_type~Timing, scales = "free", ncol = 2) +
    facet_wrap(observation_type~Timing,  ncol = 2) +
    labs(x = "NP EBE", y = "# peptide subpools giving IFN-g+ ELISPOT", color = "Study", title = "Immune Breadth vs NP") +
    scale_color_manual(values = ggsci::pal_bmj()(4)[c(2,4)]) +
    # scale_color_bmj(drop = F) + 
    theme_bw() + 
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black"), aspect.ratio=1) +
    stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black") +
    geom_smooth(method = "lm", se = F, color = "black", alpha = 0.7) 
)

ggsave(here("results_eclipse", "summary_single_fit", "immune_breadth_vs_NP.pdf"),
       immune_plot, width = 7, height = 5)

(immune_plot_pre_interrupt <- immune_combined_df %>% 
    filter(observation_type == "breadth_Total", Timing == "ART interruption") %>%
    # filter(grepl("breadth", observation_type)) %>%
    mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
           observation_type = paste(gsub("breadth_", "", observation_type), "Breadth"),
           Study = factor(Study, levels = c(1,2,3))) %>%
    ggplot(aes(NPs_mode*1e5, observation, color = factor(Study))) + 
    geom_point(size = 2.5) + 
    scale_x_log10() + 
    # facet_wrap(observation_type~Timing, scales = "free", ncol = 2) +
    facet_wrap(~Timing,  ncol = 2) +
    labs(x = "NP EBE", y = "# peptide subpools giving IFN-g+ ELISPOT", color = "Study", 
         title = "Total Immune Breadth vs NP") +
    scale_color_manual(values = ggsci::pal_bmj()(3)[2]) +
    # scale_color_bmj(drop = F) + 
    theme_bw() + 
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black"), aspect.ratio=1) +
    stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black") +
    geom_smooth(method = "lm", se = F, color = "black", alpha = 0.7) 
)

ggsave(here("results_eclipse", "summary_single_fit", "immune_breadth_vs_NP_ART_interrupt.pdf"),
       immune_plot_pre_interrupt, width = 4, height = 5)

immune_breadth_vs_treatment <- immune_combined_df %>% 
  filter(observation_type == "breadth_Total", Timing == "ART interruption") %>% 
  mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
         observation_type = paste(gsub("breadth_", "", observation_type), "Breadth"),
         Study = factor(Study, levels = c(1,2,3))) %>% 
  ggplot(aes(Treatment, observation)) + 
  geom_boxplot() + 
  geom_beeswarm(size = 2.5, aes(color = Study)) + 
  facet_wrap(~Timing) +
  scale_color_manual(values = ggsci::pal_bmj()(3)[2]) +
  labs(y = "# peptide subpools giving IFN-g+ ELISPOT")

ggsave(here("results_eclipse", "summary_single_fit", "immune_breadth_vs_treatment.pdf"),
       immune_breadth_vs_treatment, width = 4, height = 5)

immune_combined_df %>% 
  filter(observation_type == "breadth_Total") %>%
  mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
         observation_type = gsub("breadth_", "", observation_type)) %>%
  ggplot(aes(m_mode, observation, color = factor(Study))) + 
  geom_point(size = 2.5) + 
  scale_x_log10() + 
  facet_wrap(~Timing, scales = "free", nrow = 1) +
  labs(x = "m EBE", y = "# peptide subpools giving IFN-g+ ELISPOT", color = "Study", title = "Immune Breadth vs m") +
  scale_color_bmj() + 
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) +
  stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n") +
  geom_smooth(method = "lm", se = F)


immune_combined_df %>% 
  filter(observation_type == "breadth_Total") %>%
  mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
         observation_type = gsub("breadth_", "", observation_type)) %>%
  ggplot(aes(p_mode, observation, color = factor(Study))) + 
  geom_point(size = 2.5) + 
  # scale_x_log10() + 
  facet_wrap(~Timing, scales = "free", nrow = 1) +
  labs(x = "p EBE", y = "# peptide subpools giving IFN-g+ ELISPOT", color = "Study", title = "Immune Breadth vs m") +
  scale_color_bmj() + 
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) +
  stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n") +
  geom_smooth(method = "lm", se = F)




immune_plot_CD <- immune_combined_df %>% 
  filter(grepl("CD8", observation_type)) %>%
  mutate(Timing = paste(Timing, "(62 Weeks)")) %>% 
  ggplot(aes(NPs_mode*1e5, observation, color = factor(Study))) + 
  geom_point(size = 2.5) + 
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(Timing~., scales = "free", nrow = 2) +
  labs(x = "NP EBE", y = "CD8 %IFN-gamma+", color = "Study", title = "T-cell response vs NP") +
  # scale_color_bmj() +
  scale_color_manual(values = ggsci::pal_bmj()(3)[2]) +
  geom_smooth(method = "lm", se = F, color = "black", alpha = 0.5) +
  stat_cor(show.legend = F, label.x.npc = 0.7, label.y.npc = 0.9, label.sep = "\n", color = "black") +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

ggsave(here("results_eclipse", "summary_single_fit", "immune_CD8_vs_NP.pdf"),
       immune_plot_CD, width = 4, height = 5)

immune_combined_df %>% 
  filter(grepl("CD8", observation_type)) %>%
  mutate(Timing = paste(Timing, "(62 Weeks)")) %>% 
  ggplot(aes(p_mode, observation, color = factor(Study))) + 
  geom_point(size = 2.5) + 
  scale_y_log10() +
  facet_wrap(Timing~., scales = "free", nrow = 2) +
  labs(x = "p EBE", y = "CD8 %IFN-gamma+", color = "Study", title = "T-cell response vs p") +
  # scale_color_bmj() +
  scale_color_manual(values = ggsci::pal_bmj()(3)[2]) +
  geom_smooth(method = "lm", se = F, color = "black", alpha = 0.5) +
  stat_cor(show.legend = F, label.x.npc = 0.7, label.y.npc = 0.9, label.sep = "\n", color = "black") +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 



immune_combined_df %>% 
  filter(grepl("CD", observation_type)) %>%
  ggplot(aes(m_mode, observation, color = factor(Study))) + 
  geom_point(size = 3) + 
  scale_x_log10() + 
  scale_y_log10() +
  facet_wrap(Timing~observation_type, scales = "free", nrow = 2) +
  labs(x = "m EBE", y = "% Gag IFN-gamma+", color = "Study", title = "Immune Breadth vs m") +
  scale_color_bmj() + 
  
  # plot_layout(ncol = 2, guides = "collect") & 
  
  ## Check correlations ignoring censoring 
  
  # Spearman correlation
  immune_combined_df %>% 
  group_by(Timing, observation_type) %>% 
  summarise(
    cor_NPs = cor(NPs_mode, observation, method = "spearman"),
    # cor_logNPs = cor(log10(NPs_mode), observation, method = "spearman"),
    cor_m = cor(m_mode, observation, method = "spearman"),
    # cor_logm = cor(log10(m_mode), observation, method = "spearman"),
    cor_p = cor(p_mode, observation, method = "spearman")
  ) %>% 
  pivot_longer(contains("cor"), names_to = "parameter", values_to = "correlation") %>%
  mutate(parameter = fct_inorder(factor(gsub("cor_", "", parameter))),
         Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption", "Post-vaccination"))) %>%
  ggplot(aes(parameter, observation_type, fill = correlation)) + 
  geom_tile() + 
  geom_text(aes(label = round(correlation, 2))) +
  facet_wrap(~Timing, ncol = 3) +
  scale_fill_gradient2(limits = c(-1,1)) +
  labs(x = "Parameter", y = "Clinical Observation",
       title = "Spearman correlation of EBEs with clinical observations") 

## Figures of study 2 for presetnation

(study2_immune_plot_pre_interrupt <- immune_combined_df %>% 
    filter(observation_type == "breadth_Total", Timing == "ART interruption") %>%
    mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
           observation_type = paste(gsub("breadth_", "", observation_type), "Breadth"),
           Study = factor(Study, levels = c(1,2,3)),
           Treatment = factor(Treatment, levels = c("Control", "TLR7", "Vac", "Vac+TLR7"))) %>%
    ggplot(aes(NPs_mode*1e5, observation)) + 
    geom_point(size = 2.5, aes(color = Treatment)) + 
    scale_x_log10() + 
    # facet_wrap(observation_type~Timing, scales = "free", ncol = 2) +
    facet_wrap(~Timing,  ncol = 2) +
    labs(x = "NP EBE", y = "# peptide subpools giving IFN-g+ ELISPOT", color = NULL, 
         title = "Total Immune Breadth vs NP") +
    scale_color_manual(values = unname(colors)) +
    # scale_color_bmj(drop = F) + 
    theme_bw() + 
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black"), aspect.ratio=1) +
    stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black") +
    geom_smooth(method = "lm", se = F, color = "black", alpha = 0.7) 
)

study2_immune_breadth_vs_treatment <- immune_combined_df %>% 
  filter(observation_type == "breadth_Total", Timing == "ART interruption") %>% 
  mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
         observation_type = paste(gsub("breadth_", "", observation_type), "Breadth"),
         Study = factor(Study, levels = c(1,2,3)),
         Treatment = factor(Treatment, levels = c("Control", "TLR7", "Vac", "Vac+TLR7"))) %>%
  ggplot(aes(Treatment, observation, color = Treatment)) + 
  geom_boxplot(show.legend = F) + 
  geom_beeswarm(size = 2.5, show.legend = F) + 
  facet_wrap(~Timing) +
  scale_color_manual(values = unname(colors))  + 
  labs(y = "# peptide subpools giving IFN-g+ ELISPOT")

study2_immune_breadth_plot <- study2_immune_plot_pre_interrupt + 
  theme(plot.title = element_blank()) + 
  study2_immune_breadth_vs_treatment + 
  theme(axis.title.y = element_blank()) +
  plot_layout(ncol = 2, guides = "collect") & 
  theme(strip.text = element_blank(), axis.title = element_text(size = 15),
        axis.text = element_text(size = 10)) &
  plot_annotation(title = "Immune Breadth Pre-ART Interruption",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 

ggsave("presentation figures/study2_immune_breadth_correlates.pdf", study2_immune_breadth_plot, 
       width = 10, height = 5)

immune_combined_df %>% 
  filter(grepl("CD8", observation_type)) %>% 
  mutate(Study = factor(Study, levels = c(1,2,3))) %>%
  ggplot(aes(Treatment, observation)) + 
  geom_boxplot() + 
  geom_beeswarm(aes(color = Study), size = 2.5) + 
  # facet_wrap(~Timing) +
  scale_color_manual(values = ggsci::pal_bmj()(3)[2]) 
# scale_y_continuous(trans=scales::pseudo_log_trans(base = 10), limits = c(0,10), breaks = c(0.1,1,10)) 
# labs(y = "# peptide subpools giving IFN-g+ ELISOPT")

## Correlate conditional parameters with immune repsonse: ####

immune_combined_df_conditional <- left_join(
  rbind(study2_conditional %>% select(-intTLR7VAC, -tlr7, -vac) %>% mutate(Study = 2),
        study4_conditional %>% select(-ab, -tlr7, -vac) %>% mutate(Study = 4)
        ),
  immune_covar %>% #filter(Study == 2) %>% 
    pivot_longer(6:11, names_to = "observation_type", values_to = "observation"),
  relationship = "many-to-many"
) %>% 
  filter(!is.na(observation)) %>% 
  mutate(log_NPs = log10(NPs),
         log_m = log10(m),
         log_p = log10(p))

# Immune measures best correlate with NP:
immune_combined_df_conditional %>% 
  pivot_longer(c(p,log_NPs, log_m), values_to = "parameter_estimate", names_to = "parameter") %>% 
  filter(Timing != "Peak Immunity") %>% 
  ggplot(aes(parameter_estimate, observation, color = Treatment)) + 
  geom_point() + 
  facet_wrap(parameter ~interaction(observation_type, Timing), scales= "free", nrow = 3) + 
  stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black") +
  geom_smooth(method = "lm", se = FALSE, color = "black")

immune_combined_df_conditional_heatmaps <- immune_combined_df_conditional %>% 
  filter(grepl("breadth", observation_type)) %>% 
  pivot_longer(c(p, NPs, m, log_NPs, log_p, log_m), values_to = "parameter_estimate", names_to = "parameter") %>% 
  group_by(parameter, observation_type, Timing, Study) %>% 
  summarise(cor = cor(observation, parameter_estimate)) %>% 
  mutate(observation_type = case_when(grepl("breadth", observation_type) ~ (str_remove(observation_type, "breadth_")),
                                      grepl("CD4", observation_type) ~ "CD4 %IFN-gamma+",
                                      grepl("CD8", observation_type) ~ "CD8 %IFN-gamma+"
  ),
  # Timing = factor(Timing, labels = c("Peak Immunity (50 Weeks)",  "Pre-ART interruption (72 Weeks)"), 
  #                 levels = c("Peak Immunity", "ART interruption")),
  parameter = factor(parameter, levels = c("m", "log_m" ,"p", "log_p", "NPs", "log_NPs"), labels = c("m", "log10(m)", "p", "log10(p)",  "Np","log10(NPs)"))) %>%
  ggplot(aes(parameter, observation_type, fill = cor)) +
  geom_tile() + 
  geom_text(aes(label = round(cor,2))) +
  facet_wrap(paste("Study",Study)~Timing) +
  scale_fill_gradient2(limits = c(-1,1)) +
  labs(x = "Parameter", y = "Immune Breadth", fill = "Correlation")

ggsave(here("paper_figures", "immune_correlates_conditional_heatmap.pdf"), 
       immune_combined_df_conditional_heatmaps, width = 7, height = 5)

make_immune_combined_df_conditional <- function(immune_combined_df_conditional, obs_type, timing = "ART interruption", study = 2){
  immune_combined_df_conditional %>% 
    filter(observation_type == obs_type, Timing == timing, Study == study) %>% 
    mutate(observation_type = case_when(grepl("breadth", observation_type) ~ (str_remove(observation_type, "breadth_")),
                                        grepl("CD4", observation_type) ~ "CD4 %IFN-gamma+",
                                        grepl("CD8", observation_type) ~ "CD8 %IFN-gamma+"
    ),
    Timing = factor(Timing, labels = c("Peak Immunity (50 Weeks)",  "Pre-ART interruption (72 Weeks)"), 
                    levels = c("Peak Immunity", "ART interruption"))) %>% 
    ggplot(aes(log_NPs+5, observation, color = Treatment)) +
    geom_point(size = 3) + 
    facet_wrap(~observation_type, scales = "free_y", ncol = 1) +
    scale_x_continuous(labels = scales::label_math(10^.x)) +
    geom_smooth(method = "lm", color = "black", se = FALSE) +
    stat_cor(label.x.npc = 0.6, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black", size = 4) +
    # scale_x_log10() +
    labs(x = bquote("Individual"~N['p']~'estimate'), y = "# peptide subpools giving IFN-g+ ELISPOT", color = "Treatment") +
    scale_color_manual(values = (colors))
}

make_immune_boxplots <- function(immune_combined_df, obs_type, timing, study = 2){
  immune_combined_df_conditional %>% 
    filter(observation_type == obs_type, Timing == timing, Study == study) %>% 
    filter(rep == 1) %>% 
    # distinct(Study, id, Treatment, observation_type, observation ) %>% 
    mutate(observation_type = str_remove(observation_type, "breadth_")) %>% 
           # Treatment = factor(Treatment, levels = c("Control", "TLR7", "Vac", "Vac+TLR7"))) %>%
    ggplot(aes(Treatment, observation)) + 
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(aes(color = Treatment), size = 2.5, show.legend = F) + 
    facet_wrap(~observation_type, ncol = 1, scales = "free_y") + 
    scale_color_manual(values = (colors)) + 
    labs(y = "# peptide subpools giving IFN-g+ ELISPOT")
}

immune_correlations_conditional <- make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Total", "ART interruption") + 
  theme(axis.title.x = element_blank()) +
  make_immune_boxplots(immune_combined_df, "breadth_Total", "ART interruption") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Gag", "ART interruption") + 
  theme(axis.title.x = element_blank()) +
  make_immune_boxplots(immune_combined_df, "breadth_Gag", "ART interruption") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Env", "ART interruption") + 
  theme(axis.title.x = element_blank()) +
  make_immune_boxplots(immune_combined_df, "breadth_Env", "ART interruption") +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank()) +
  make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Pol", "ART interruption") + 
  make_immune_boxplots(immune_combined_df, "breadth_Pol", "ART interruption")  + 
  theme( axis.title.y = element_blank()) +
  plot_layout(ncol = 2, guides = "collect")

ggsave(here("paper_figures", "study2_immune_correlates_conditional_scatterplots.pdf"), 
       immune_correlations_conditional, width = 10, height = 15)

## Make a figure with study 2 and 4:

make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Total", "ART interruption") +
make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Total", "Post MVA vaccination") +
  theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = scales::label_math(10^.x), breaks = seq(-1,7, by = 2)) &
  coord_cartesian(ylim = c(0, 30))

make_immune_boxplots(immune_combined_df, "breadth_Total", "ART interruption") + 
  make_immune_boxplots(immune_combined_df, "breadth_Total", "Post MVA vaccination", study = 4) +
  theme(axis.title.y = element_blank()) &
  coord_cartesian(ylim = c(0, 30))

# simple_immune_correlations <-  immune_combined_df_conditional_heatmaps / 
#   (make_immune_combined_df_conditional(immune_combined_df_conditional %>% mutate(breadth_label = "Total Breadth"), "breadth_Total", "ART interruption") + 
#      theme(legend.position = "none") + coord_cartesian(ylim = c(0,30)) +
#      make_immune_boxplots(immune_combined_df %>% mutate(breadth_label = "Total Breadth"), "breadth_Total", "ART interruption") +
#      theme(axis.title.y = element_blank()) + coord_cartesian(ylim = c(0,30))  &
#      facet_wrap(~breadth_label)) + 
#   plot_layout(ncol = 1) &
#   plot_annotation(tag_levels = "a")
# 
# ggsave(here("paper_figures", "immune_correlates_conditional_simple.pdf"), simple_immune_correlations, 
#        width = 8, height = 7)


## Correlate initial CD4 counts with lambda EBE ####

cd4_correlates <- cd4_covar %>% 
  left_join(EBEs %>% select(id, lambda_mode)) %>% 
  filter(!is.na(Rebound_Initial_CD4)) %>% 
  ggplot(aes(lambda_mode, Rebound_Initial_CD4, color = factor(Study))) + 
  geom_point(size = 2.5) +
  labs(color = "Study", x = "Lambda EBE", y = "CD4 count (cells/ul)") + 
  coord_cartesian(ylim = c(0,1600)) + 
  scale_color_manual(values = ggsci::pal_bmj()(3)[1])  + 
  stat_cor(color = "black", label.x.npc = 0.75, label.y.npc = 0.9, label.sep = "\n") + 
  geom_smooth(method = "lm", color = "black", se = F) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1)  +
  
  cd4_covar %>% ggplot(aes(Treatment, Rebound_Initial_CD4)) + 
  geom_boxplot(outlier.shape = NA) + 
  ggbeeswarm::geom_beeswarm(size = 2.5, aes(color = factor(Study))) +
  scale_color_manual(values = ggsci::pal_bmj()(3)[1])  +
  labs(color = "Study",  y = "CD4 count (cells/ul)") + 
  coord_cartesian(ylim = c(0,1600)) + 
  
  plot_layout(guides = "collect", axes = "collect") &
  plot_annotation(title = "Initital CD4 Count at Rebound")

ggsave(here("results_eclipse", "summary_single_fit", "cd4_correlates.pdf"),
       cd4_correlates, width = 8, height = 5)

