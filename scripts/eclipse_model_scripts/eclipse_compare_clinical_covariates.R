## Compare correlates of immune response and reservoir size to EBEs

library(tidyverse)
library(here)
library(ggsci)
library(survival)
library(patchwork)
library(ggbeeswarm)
library(censReg)
library(NADA2)

source(here("scripts", "functions.R"))

theme_set(theme_bw() + 
            theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
                  panel.border = element_rect(colour = "black"), aspect.ratio=1) )
colors <- c("Control" = "black", "TLR7" = "#027dbb", "Vac" = "#fcba12", "Ab" ="#fa61c2", "TLR7 + Vac" = "#7ab148",  "TLR7 + Ab" = "#7f48b1", "TLR7 + Vac + Ab" =  "#a96e7c")
study_colors <- ggsci::pal_bmj()(4)
names(study_colors) <- paste("Study", 1:4)

## Read in clinical covariates ####
# Total SIV DNA in CD4+ T memory cells for studies 1-3, intact provirus for study 4
reservoir_covar <- read_csv(here("data", "derived", "reservoir_correlates.csv")) %>% 
  mutate(Treatment = factor(Treatment, levels = names(colors)))
# Immune breadth against Gag, Pol,and Env (and total) available for Study 2 & 4
immune_covar <- read_csv(here("data", "derived", "immune_correlates.csv")) %>% 
  mutate(Treatment = factor(Treatment, levels = names(colors)))
# Size of immune response against Gag, Pol, and Env (and total) available for Study 4
immune_size_covar <- read_csv(here("data", "derived", "immune_size_correlates.csv"))

## Read in sampled individual parameters from full conditional distribution ####
study1_conditional <- read_csv(here("results_eclipse", "study1", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))
study2_conditional <- read_csv(here("results_eclipse", "study2", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))
study3_conditional <- read_csv(here("results_eclipse", "study3", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))
study4_conditional <- read_csv(here("results_eclipse", "study4", "best_model", "IndividualParameters", "simulatedIndividualParameters.txt"))

all_conditional_param <- rbind(
  study1_conditional %>% mutate(Study = 1) %>% select(-tlr7),
  study2_conditional %>% mutate(Study = 2) %>% select(-c(intTLR7VAC, tlr7, vac)),
  study3_conditional %>% mutate(Study = 3) %>% select(-c(inttlr7ab, tlr7, ab)),
  study4_conditional %>% mutate(Study = 4) %>% select(-c(ab, tlr7, vac))
)

## Correlate with Reservoir size ####
# Because there is a different measure for study 4 than 1-3, I will look at it separately


reservoir_conditional_df <- left_join(
  all_conditional_param %>% filter(Study != 4), 
  reservoir_covar %>% mutate(log_SIV_DNA = log10(SIV_DNA)),
  relationship = "many-to-many"
) %>% 
  filter(!is.na(SIV_DNA)) %>% 
  mutate(log_ta = log10(ta))

reservoir_conditional_df %>% 
  mutate(cell_type = fct_rev(cell_type)) %>% 
  ggplot(aes(ta, SIV_DNA, color = paste("Study", Study), shape = ifelse(Censored, "Censored", "Uncensored"))) +
  geom_point(size = 3) +
  facet_wrap(~cell_type) + 
  scale_x_log10() + scale_y_log10() + 
  scale_color_manual("Study", values = study_colors) +
  scale_shape_manual("", values = c("Censored" = 1, "Uncensored" = 16)) +
  labs(x = "Individual ta estimate", y ="SIV/SHIV DNA (per million cells)")

reservoir_conditional_df_S4 <- left_join(
  all_conditional_param %>% filter(Study == 4), 
  reservoir_covar %>% mutate(log_SIV_DNA = log10(SIV_DNA)),
  relationship = "many-to-many"
) %>% 
  filter(!is.na(SIV_DNA)) %>% 
  mutate(log_ta = log10(ta))

# Lots of censored data and limited correlation between intact provirus and ta estimate
reservoir_conditional_df_S4 %>% 
  mutate(cell_type = fct_rev(cell_type)) %>% 
  ggplot(aes(ta, SIV_DNA, color = paste("Study", Study), shape = ifelse(Censored, "Censored", "Uncensored"))) +
  geom_point(size = 3) +
  facet_wrap(~cell_type) + 
  # scale_x_log10() + 
  scale_y_log10() + 
  scale_color_manual("Study", values = study_colors) +
  scale_shape_manual("", values = c("Censored" = 1, "Uncensored" = 16)) +
  labs(x = "Individual ta estimate", y ="Intact provirus (per million cells)")


# Calculate correlation for each cell type with tobit regression accounting for censored values ####
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
  labs(x = bquote(t['a']), y = "SIV DNA\n(per million cells)", 
       color = "Study", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1)

reservoir_conditional_predictions_S4 <- merge(
  reservoir_conditional_df_S4,
  fit_tobit_ta(reservoir_conditional_df_S4, "PBMC")$prediction %>% mutate(rep = rep(1:10, length.out = nrow(.)))
)

reservoir_plot_conditional_S4 <- reservoir_conditional_predictions_S4 %>% 
  mutate(censored = ifelse(Censored == 1, "Censored", "Uncensored"),
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
  labs(x = bquote(t['a']), y = "Intact provirus (per million cells)", 
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
  
  df <- prediction_df %>% 
    filter(Study == study) %>%
    mutate(cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC")),
           Treatment = factor(Treatment, levels = names(colors)))
  
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
    scale_color_manual(values = (colors), drop = F) +
    scale_shape_manual(values = c("1" = 1, "0" = 16)) + 
    labs(x = bquote('Individual'~t['a']~'estimates'), y = "SIV DNA\n(per million cells)", 
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
  
  reservoir_combined_df %>% 
    filter(Study %in% study) %>%
    mutate(Treatment = factor(Treatment, levels = names(colors)),
           cell_type = factor(cell_type, levels = c("PBMC", "LNMC", "GMMC"))) %>%
    ggplot(aes(Treatment, SIV_DNA, color = Treatment)) + 
    geom_hline(yintercept = 3, lty = "dashed", color = "gray") + 
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(size = 2, alpha = 0.7, aes(shape = factor(Censored)), cex=2, priority = "density") + 
    scale_y_log10() +
    scale_color_manual(values = ggsci::pal_bmj()(3)[1:2], drop = F) +
    labs(color = "Study", y = "SIV DNA\n(per million cells)") + 
    scale_shape_manual(values = c("1" = 1, "0" = 16)) + 
    scale_color_manual(values = c((colors), (colors)[3:4])) +
    guides(shape = 'none') +
    facet_wrap(~cell_type, scales = "free_x", ncol = 1) + 
    theme(legend.position = "none") 
}


study_res_plots_conditional <- make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 1, scientific_x = T) +
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 2) +
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 3) +
  make_reservoir_plot_conditional_by_study(reservoir_conditional_predictions_by_study %>% filter(cell_type == "PBMC"), 4) +
  labs(y = "Intact provirus (per million cells)") +
  
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 1) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 2) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 3) + 
  make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 4) + 
  labs(y = "Intact provirus (per million cells)") +
  plot_layout(nrow = 2, guides = "collect") & 
  theme(#strip.text = element_blank(), 
    axis.title = element_text(size = 15),
    axis.text = element_text(size = 10),
    plot.title = element_blank()) &
  # coord_cartesian(ylim = c(1,10000)) &
  plot_annotation(title = "Correlates of Reservoir Size",
                  theme = theme(plot.title = element_text(hjust = 0.5, size = 20))) 



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
    labs(x = bquote("Individual"~N['p']~'estimate'), y = "Number of peptide subpools giving IFN-g+ ELISPOT", color = "Treatment") +
    scale_color_manual(values = (colors), drop =F)
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
    geom_beeswarm(aes(color = Treatment), alpha = 0.7,
                  size = 2, show.legend = F, cex=2, priority = "density") +
    facet_wrap(~observation_type, ncol = 1, scales = "free_y") + 
    scale_color_manual(values = (colors), drop = F) + 
    labs(y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT")
}


## Make a figure with study 2 and 4:

make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Total", "ART interruption") +
  make_immune_combined_df_conditional(immune_combined_df_conditional, "breadth_Total", "Post MVA vaccination", study = 4) +
  theme(axis.title.y = element_blank()) +
  scale_x_continuous(labels = scales::label_math(10^.x), breaks = seq(-1,7, by = 2)) &
  coord_cartesian(ylim = c(0, 30))

make_immune_boxplots(immune_combined_df, "breadth_Total", "ART interruption") + 
  make_immune_boxplots(immune_combined_df, "breadth_Total", "Post MVA vaccination", study = 4) +
  theme(axis.title.y = element_blank()) &
  coord_cartesian(ylim = c(0, 30))

## Paper Figure ####

reservoir_plot_conditional <- reservoir_conditional_predictions %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored")) %>%
  filter(cell_type == "PBMC", Study != 4) %>% 
  ggplot(aes(log_ta, SIV_DNA, color = Treatment)) + 
  geom_hline(yintercept = 3, lty = "dashed", color = "gray") + 
  geom_jitter(data = . %>% filter(censored == "Censored"), size = 2, aes(shape = censored, y = SIV_DNA - 1.5), show.legend = T, alpha = 0.7, 
              width = 0, height = 0.3) + 
  geom_point(data = . %>% filter(censored == "Uncensored"), size = 2, aes(shape = censored), show.legend = T, alpha = 0.7) + 
  # add tobit regression line
  geom_line(aes(y = 10^predicted), color = "black", linewidth = 1, alpha = 1) +
  # add pvalue
  geom_text(data =. %>% group_by(cell_type) %>% slice(1), 
            aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
            hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_y_log10() + #scale_x_log10() + 
  scale_x_continuous(labels = scales::label_math(10^.x), breaks = c(-2,0,2)) +
  coord_cartesian(ylim = c(1, 10000)) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_manual(values =colors, drop = F) +
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(x = bquote('Individual'~t['a']~'estimate'), y = "SIV DNA\n(per million cells)", 
       color = "Treatment", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

immune_plot_conditional <- immune_combined_df_conditional %>% 
  filter(Timing %in% c("Post MVA vaccination", "ART interruption"),
         observation_type == "breadth_Total") %>% 
  ggplot(aes(log_NPs+5, observation, color = Treatment)) +
  geom_jitter(data = . %>% filter(observation == 0), size = 2, show.legend = T, alpha = 0.7,
             width = 0, height = 1, aes(y = observation - 1)) + 
  geom_point(data = . %>% filter(observation != 0), size = 2, show.legend = T, alpha = 0.7) + 
  scale_x_continuous(labels = scales::label_math(10^.x), breaks = seq(-2,8, by = 3), limits = c(-4,9)) +
  geom_smooth(method = "lm", color = "black", se = FALSE, alpha = 0.7, fullrange = TRUE) +
  stat_cor(label.x.npc = 0.5, label.y.npc = 1, show.legend = F, label.sep = "\n", color = "black", size = 4) +
  # scale_x_log10() +
  labs(x = bquote("Individual"~N['p']~'estimate'), y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT", color = "Treatment") +
  scale_color_manual(values = (colors), drop =F) +
  coord_cartesian(ylim = c(-2, 30), xlim = c(-3,8)) 

immune_plot_conditiona_jl <- immune_combined_df_conditional %>% 
  filter(Timing %in% c("Post MVA vaccination", "ART interruption"),
         observation_type == "breadth_Total") %>% 
  ggplot(aes(log_NPs+5, observation, color = Treatment)) +
  geom_jitter( size = 2, show.legend = T, alpha = 0.7) + 
  scale_x_continuous(labels = scales::label_math(10^.x), breaks = seq(-2,8, by = 3), limits = c(-4,9)) +
  geom_smooth(method = "lm", color = "black", se = FALSE, alpha = 0.7, fullrange = TRUE) +
  stat_cor(label.x.npc = 0.5, label.y.npc = 1, show.legend = F, label.sep = "\n", color = "black", size = 4) +
  # scale_x_log10() +
  labs(x = bquote("Individual"~N['p']~'estimate'), y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT", color = "Treatment") +
  scale_color_manual(values = (colors), drop =F) +
  coord_cartesian(ylim = c(-2, 30), xlim = c(-3,8))

reservoir_trt_s1 <- make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 1) + coord_cartesian(ylim = c(3, 10000)) 
reservoir_trt_s2 <- make_treatment_v_reservoir(reservoir_covar %>% filter(cell_type=="PBMC"), 2)  + coord_cartesian(ylim = c(3, 10000)) 

immune_trt_s2 <- make_immune_boxplots(immune_combined_df, "breadth_Total", "ART interruption", study = 2) + coord_cartesian(ylim = c(0,30))
immune_trt_s4 <- make_immune_boxplots(immune_combined_df, "breadth_Total", "Post MVA vaccination", study = 4) + coord_cartesian(ylim = c(0,30))

r1 <-  reservoir_plot_conditional + 
  labs(title = "Reservoir Correlations", tag = "a") +
  theme(axis.title.x = element_text(vjust = 10)) +
  scale_y_log10(labels = label_log()) +
  coord_cartesian(ylim = c(3, 10000)) +
  reservoir_trt_s1 + theme(axis.title.y = element_blank()) + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) +
  labs(title = "SIV-TLR7", tag = "b") + 
  scale_y_log10(labels = label_log()) +
  reservoir_trt_s2 + theme(axis.title.y = element_blank()) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) +
  labs(title = "SIV-TLR7-Vac", tag = "c") + 
  scale_y_log10(labels = label_log()) + 
  plot_layout(ncol = 3) & 
  guides(color = "none", shape = "none")

r2 <- immune_plot_conditional +
  labs(title = "Immune Correlations", tag = "d", 
       y = "Number of peptide subpools\ngiving IFN-g+ ELISPOT") + 
  theme(axis.title.x = element_text(vjust = 10)) +
  immune_trt_s2 + theme(axis.title.y = element_blank()) +
  labs(title = "SIV-TLR7-Vac", tag = "e") + 
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) +
  immune_trt_s4 + theme(axis.title.y = element_blank()) +
  scale_x_discrete(labels = c("Control", "TLR7 + Vac", "TLR7 + Ab", "TLR7 +\nVac + Ab")) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1), axis.title.x = element_blank()) +
  labs(title = "SHIV-TLR7-Vac-Ab", tag = "f") + 
  plot_layout(ncol = 3) & 
  guides(color = "none", shape = "none")

paper_figure <- (r1 / get_legend(reservoir_plot_conditional + 
                                   theme(legend.position = "bottom",
                                         legend.box = "vertical",
                                         legend.spacing.y = unit(-0.1, "cm"),
                                         legend.background = element_blank()) +
                                   guides(color = guide_legend(nrow = 1, 
                                                               override.aes = list(alpha = 1),
                                                               order = 2),
                                          shape = guide_legend(override.aes = list(alpha = 1, order = 1))) ) / 
                   plot_spacer() /
                   r2) + 
  plot_layout(heights = c(1,0.1,0.05,1)) & 
  theme(strip.text = element_blank(),
         plot.title = element_text(size = 12),
        panel.grid.minor = element_blank()) &
  # plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(face = 2),
        plot.title = element_text(size = 12),
        aspect.ratio = NULL,
        plot.background = element_blank(),
        plot.margin = margin(5,5,0,5)) 

ggsave(here("paper_figures", "combined_immune_correlate_figure.pdf"), paper_figure,
       width = 8, height = 7)

reservoir_plot_conditional_flip <- reservoir_conditional_predictions %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored")) %>%
  filter(cell_type == "PBMC", Study != 4) %>% 
  ggplot(aes(SIV_DNA, log_ta, color = Treatment)) + 
  geom_vline(xintercept = 3, lty = "dashed", color = "gray") + 
  geom_point(size = 2, aes(shape = censored), show.legend = T) + 
  # add tobit regression line
  geom_line(aes(x = 10^predicted), color = "black", linewidth = 1, alpha = 0.7) +
  # add pvalue
  geom_text(data =. %>% group_by(cell_type) %>% slice(1), 
            aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
            hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  scale_x_log10() + #scale_x_log10() + 
  scale_y_continuous(labels = scales::label_math(10^.x)) +
  coord_cartesian(xlim = c(3, 10000)) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_manual(values =colors, drop = F) +
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(y = bquote('Individual'~t['a']~'estimate'), x = "SIV DNA\n(per million cells)", 
       color = "Treatment", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 

immune_plot_conditional_flip <- immune_combined_df_conditional %>% 
  filter(Timing %in% c("Post MVA vaccination", "ART interruption"),
         observation_type == "breadth_Total") %>% 
  ggplot(aes(observation, log_NPs+5, color = Treatment)) +
  geom_point(size = 2, show.legend = T) + 
  scale_y_continuous(labels = scales::label_math(10^.x)) +
  geom_smooth(method = "lm", color = "black", se = FALSE, alpha = 0.7, fullrange = TRUE) +
  stat_cor(label.x.npc = 0.6, label.y.npc = 1, show.legend = F, label.sep = "\n", color = "black", size = 4) +
  # scale_x_log10() +
  labs(y = bquote("Individual"~N['p']~'estimate'), x = "Number of peptide subpools\ngiving IFN-g+ ELISPOT", color = "Treatment") +
  scale_color_manual(values = (colors), drop =F) +
  scale_x_continuous(limits = c(0, 35)) +
  coord_cartesian(xlim = c(0,30))

c1 <- reservoir_plot_conditional_flip + 
         labs(title = "Reservoir Correlations", tag = "a") +
         reservoir_trt_s1 + theme(axis.title.y = element_blank()) + 
         labs(title = "SIV-TLR7", tag = "b") + 
         scale_x_discrete(limits = rev) +
         scale_y_log10(labels = label_log()) +
         coord_flip() +
         reservoir_trt_s2 + theme(axis.title.y = element_blank()) +
         labs(title = "SIV-TLR7-Vac", tag = "c") + 
         scale_x_discrete(limits = rev) +
         scale_y_log10(labels = label_log()) +
         coord_flip() + 
         plot_layout(ncol = 1, guides = "collect") &
  theme(legend.position = "none")

c2 <- immune_plot_conditional_flip + 
  labs(title = "Immune Correlations", tag = "d") + 
  immune_trt_s2 + theme(axis.title.y = element_blank()) +
  labs(title = "SIV-TLR7-Vac", tag = "e") + 
  scale_x_discrete(limits = rev) + coord_flip() +
  immune_trt_s4 + 
  labs(title = "SHIV-TLR7-Vac-Ab", tag = "f") + 
  # scale_x_discrete(limits = rev) + 
  theme(axis.title.y = element_blank()) + coord_flip() +
  scale_x_discrete(labels = rev(c("Control", "TLR7 + Vac", "TLR7 + Ab", "TLR7 +\nVac + Ab")), 
                   limits = rev) +
  plot_layout(ncol = 1) &
  theme(legend.position = "none")

paper_figure_vert <- (c1 | get_legend(reservoir_plot_conditional_flip + theme(legend.title = element_text(hjust = 0.5)))| c2) +   
  plot_layout(widths = c(1,0.5,1)) +
  theme(strip.text = element_blank(),
        plot.title = element_text(size = 12),
        panel.grid.minor = element_blank()) &
  # plot_annotation(tag_levels = "a") & 
  theme(plot.tag = element_text(face = 2))

ggsave(here("paper_figures", "combined_immune_correlate_figure_vertical.pdf"), paper_figure_vert,
       width = 8, height = 10)

##### 
# Try Tobit regresion with censReg ... provides same results as my method


fit_censReg_ta <- function(df, celltype, censor_value = 3){
  # browser()
  df <- filter(df, cell_type == celltype)
  fit <- survreg(Surv(log_SIV_DNA, log_SIV_DNA > log10(censor_value), type = "left") ~ log_ta, 
                 data = df, dist = "gaussian")
  
  fit <- censReg(log_SIV_DNA ~ log_ta, left = log10(censor_value), data = df)
  prediction <- tibble(
    Study = df$Study, id = df$id, cell_type = celltype, pval = summary(fit)$table[2,"p"],
    predicted = predict(fit, newdata = df)
  )
  return(list(fit = fit, prediction = prediction))
}

fit_censReg_ta(reservoir_conditional_df, "PBMC")

## Calculate censored correlation with Kendall's tau
df <- reservoir_conditional_df %>% filter(cell_type == "PBMC")

kendall_tau <- NADA2::cenken(
  y       = df$log_SIV_DNA,
  ycen    = as.logical(df$Censored),
  x       = df$log_ta
)

print(kendall_tau)
# Returns: tau, slope, intercept, p-value

# Optional: plot with the censored regression line
with(df, cenxyplot(x = log_ta, y = log_SIV_DNA, ycen = as.logical(Censored), xcen = F, lty = "solid"))    # x vs. y
lines(kendall_tau, col = "red")

# Visualize using ggplot2
reservoir_conditional_df %>% filter(cell_type == "PBMC") %>% 
  ggplot(aes(x = log_ta, y = log_SIV_DNA, shape = Censored)) +
  # Plot detected points (normal dots)
  geom_point(data = ~filter(.x, Censored == 0), aes(shape = "Uncensored"), size = 3) +
  # Plot censored bounds (downward arrows)
  geom_point(data = ~filter(.x, Censored == 1), aes(shape = "Censored"), size = 3) +
  scale_shape_manual(values = c("Uncensored" = 19, "Censored" = 6)) +
  # Add the Akritas-Theil-Sen (ATS) regression line 
  geom_abline(intercept = kendall_tau$intercept, slope = kendall_tau$slope, 
              color = "blue", linetype = "dashed", linewidth = 1) +
  labs(title = "Censored Data Kendall's Tau Regression",
       subtitle = "Blue line represents the Akritas-Theil-Sen line",
       shape = "Observation Type",
       y = "Dependent Variable (y)", 
       x = "Independent Variable (x)") +
  theme_minimal()


reservoir_plot_conditional_kendall <- reservoir_conditional_predictions %>% 
  mutate(censored = ifelse(SIV_DNA == 3, "Censored", "Uncensored")) %>%
  filter(cell_type == "PBMC", Study != 4) %>% 
  ggplot(aes(log_ta, log_SIV_DNA, color = Treatment)) + 
  geom_hline(yintercept = log10(3), lty = "dashed", color = "gray") + 
  geom_jitter(data = . %>% filter(censored == "Censored"), size = 2, aes(shape = censored, y = log_SIV_DNA - 0.4), show.legend = T, alpha = 0.7, 
              width = 0, height = 0.3) + 
  geom_point(data = . %>% filter(censored == "Uncensored"), size = 2, aes(shape = censored), show.legend = T, alpha = 0.7) + 
  # add tobit regression line
  # geom_line(aes(y = 10^predicted), color = "black", linewidth = 1, alpha = 1) +
  geom_abline(intercept = kendall_tau$intercept, slope = kendall_tau$slope, linewidth = 1) +
  # add pvalue
  # geom_text(data =. %>% group_by(cell_type) %>% slice(1), 
  #           aes(Inf, Inf, label = paste0("p = ", signif(pval, 2))),
  #           hjust = 1.1, vjust = 1.5, size = 4, color = "black") +
  # scale_y_log10() + #scale_x_log10() + 
  scale_x_continuous(labels = scales::label_math(10^.x), breaks = c(-2,0,2)) +
  scale_y_continuous(labels = function(breaks) 10^breaks, breaks = 0:4, limits = c(-0.5,4)) +
  coord_cartesian(ylim = c(-0.25, 4)) +
  facet_wrap(~cell_type, nrow = 1) +
  scale_color_manual(values =colors, drop = F) +
  scale_shape_manual(values = c("Censored" = 1, "Uncensored" = 16)) + 
  labs(x = bquote('Individual'~t['a']~'estimate'), y = "SIV DNA\n(per million cells)", 
       color = "Treatment", title = "Correlates of Reservoir Size", shape = NULL) +
  theme_bw() + 
  theme(strip.background = element_blank(), plot.title = element_text(hjust = 0.5),
        panel.border = element_rect(colour = "black"), aspect.ratio=1) 
