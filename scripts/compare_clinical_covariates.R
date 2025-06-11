## Compare correlates of immune response and reservoir size to EBEs

library(tidyverse)
library(here)
library(ggsci)
library(survival)
library(patchwork)

theme_set(theme_bw() + 
            theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
                  panel.border = element_rect(colour = "black"), aspect.ratio=1) )

# Read in clinical covariates
reservoir_covar <- read_csv(here("data", "derived", "reservoir_correlates.csv"))
immune_covar <- read_csv(here("data", "derived", "immune_correlates.csv"))
cd4_covar <- read_csv(here("data", "derived", "cd4_correlates.csv"))

# Read in EBEs:
EBEs <- read_csv(here("results", "summary_single_fit", "SIV_SHIV_individual_parameters.csv")) 

indivParams <- EBEs %>% select(data, id, contains("mode")) %>% mutate(Study = as.numeric(str_remove(data, "study")))

# Correlate reservoir data with EBEs ####

reservoir_combined_df <- reservoir_covar %>% 
  pivot_longer(c(PBMC, LNMC, GMMC), names_to = "cell_type", values_to = "SIV_DNA") %>% 
  merge(indivParams) %>% 
  filter(!is.na(SIV_DNA)) %>% 
  mutate(log_SIV_DNA = log10(SIV_DNA), log_ta = log10(ta_mode))

reservoir_plot <- reservoir_combined_df %>% 
  ggplot(aes(ta_mode, SIV_DNA, color = factor(Study))) + 
  geom_point(size = 2.5, aes(shape = SIV_DNA > 3)) + 
  scale_y_log10() + scale_x_log10() + 
  facet_wrap(~cell_type, ncol = 1) +
  labs(x = "ta EBE", y = "SIV DNA (per million cells)", color = "Study") + 
  scale_color_bmj() +
  scale_shape_manual(values = c("TRUE" = 16, "FALSE" = 1))  +
  guides(shape = "none")

# Calculate correlation for each cell type with tobit regression
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

reservoir_predictions <- rbind(
  fit_tobit_ta(reservoir_combined_df, "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df, "GMMC")$prediction,
  fit_tobit_ta(reservoir_combined_df, "LNMC")$prediction
)


reservoir_plot <- reservoir_combined_df %>% 
  left_join(reservoir_predictions) %>% 
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
# theme(plot.title = element_text(hjust = 0.5)) 

# ggsave(here("results", "summary_single_fit", "reservoir_correlates.pdf"),
#        reservoir_plot, width = 6, height = 7.5)

ggsave(here("results", "summary_single_fit", "reservoir_correlates.pdf"),
       reservoir_plot, width = 10, height = 5)

reservoir_predictions_by_study <- rbind(
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 1), "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 1), "GMMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 1), "LNMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 2), "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 2), "LNMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 3), "PBMC")$prediction,
  fit_tobit_ta(reservoir_combined_df %>% filter(Study == 3), "LNMC")$prediction
)

pvals <- reservoir_predictions_by_study %>% 
  group_by(cell_type, Study) %>% slice(1) %>% 
  mutate(y = ifelse(Study == 1, 1e4, ifelse(Study == 2, 1e3, ifelse(Study == 3, 1e2))),
         cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC")))

reservoir_combined_df %>% 
  left_join(reservoir_predictions_by_study) %>% 
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

ggsave(here("results", "summary_single_fit", "reservoir_correlates_treatment.pdf"),
       treatment_v_reservoir, width = 8, height = 5)

# Correlate immune data with Np, p, and m EBEs ####

immune_combined_df <- immune_covar %>% 
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
           Study = factor(Study, levels = c(1,2,3))) %>%
    ggplot(aes(NPs_mode*1e5, observation, color = factor(Study))) + 
    geom_point(size = 2.5) + 
    scale_x_log10() + 
    # facet_wrap(observation_type~Timing, scales = "free", ncol = 2) +
    facet_wrap(observation_type~Timing,  ncol = 2) +
    labs(x = "NP EBE", y = "# peptide subpools giving IFN-g+ ELISOPT", color = "Study", title = "Immune Breadth vs NP") +
    scale_color_manual(values = ggsci::pal_bmj()(3)[2]) +
    # scale_color_bmj(drop = F) + 
    theme_bw() + 
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black"), aspect.ratio=1) +
    stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black") +
    geom_smooth(method = "lm", se = F, color = "black", alpha = 0.7) 
)

ggsave(here("results", "summary_single_fit", "immune_breadth_vs_NP.pdf"),
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
    labs(x = "NP EBE", y = "# peptide subpools giving IFN-g+ ELISOPT", color = "Study", 
         title = "Total Immune Breadth vs NP") +
    scale_color_manual(values = ggsci::pal_bmj()(3)[2]) +
    # scale_color_bmj(drop = F) + 
    theme_bw() + 
    theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
          panel.border = element_rect(colour = "black"), aspect.ratio=1) +
    stat_cor(label.x.npc = 0.7, label.y.npc = 0.9, show.legend = F, label.sep = "\n", color = "black") +
    geom_smooth(method = "lm", se = F, color = "black", alpha = 0.7) 
)

ggsave(here("results", "summary_single_fit", "immune_breadth_vs_NP_ART_interrupt.pdf"),
       immune_plot_pre_interrupt, width = 4, height = 5)

immune_combined_df %>% 
  filter(observation_type == "breadth_Total") %>%
  mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
         observation_type = gsub("breadth_", "", observation_type)) %>%
  ggplot(aes(m_mode, observation, color = factor(Study))) + 
  geom_point(size = 2.5) + 
  scale_x_log10() + 
  facet_wrap(~Timing, scales = "free", nrow = 1) +
  labs(x = "m EBE", y = "# peptide subpools giving IFN-g+ ELISOPT", color = "Study", title = "Immune Breadth vs m") +
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
  labs(x = "p EBE", y = "# peptide subpools giving IFN-g+ ELISOPT", color = "Study", title = "Immune Breadth vs m") +
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

ggsave(here("results", "summary_single_fit", "immune_CD8_vs_NP.pdf"),
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


# Compare immune breadth data between treatment groups

immune_breadth_vs_treatment <- immune_combined_df %>% 
  filter(observation_type == "breadth_Total", Timing == "ART interruption") %>% 
  mutate(Timing = factor(Timing, levels = c("Peak Immunity", "ART interruption"), labels = c("Peak Immunity (50 Weeks)", "Pre-ART interruption (72 Weeks)")),
         observation_type = paste(gsub("breadth_", "", observation_type), "Breadth"),
         Study = factor(Study, levels = c(1,2,3))) %>%
  ggplot(aes(Treatment, observation)) + 
  geom_boxplot() + 
  geom_beeswarm(aes(color = Study), size = 2.5) + 
  facet_wrap(~Timing) +
  scale_color_manual(values = ggsci::pal_bmj()(3)[2])  + 
  labs(y = "# peptide subpools giving IFN-g+ ELISOPT")

ggsave(here("results", "summary_single_fit", "immune_breadth_vs_treatment.pdf"),
       immune_breadth_vs_treatment, width = 4, height = 5)

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

ggsave(here("results", "summary_single_fit", "cd4_correlates.pdf"),
       cd4_correlates, width = 8, height = 5)

