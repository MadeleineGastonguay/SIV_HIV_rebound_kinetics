library(tidyverse)
library(ggsci)
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())

## Example inference plots ####
loadProject("results/study1/best_model.mlxtran")

(popparam <- getEstimatedPopulationParameters())

indivParams <- getEstimatedIndividualParameters(method = "conditionalMode")$conditionalMode
dat <- read_delim("data/derived/STUDY1_DATA.txt", delim = "\t")

parameters <- read_csv("results/study1/best_model/IndividualParameters/estimatedIndividualParameters.txt")

indivParams %>% left_join(distinct(data, id, tlr7) %>% mutate(id = as.character(id))) %>% 
  pull(lambda) %>% range

lambda_density <- tibble(
  lambda = seq(100, 400, by = 1),
  density_tlr70 = dnorm(log(lambda), log(200), popparam[["omega_lambda"]]),
  density_tlr71 = dnorm(log(lambda), log(200) + popparam[["beta_lambda_tlr7_1"]] , popparam[["omega_lambda"]])
)

p1 <- lambda_density %>% 
  ggplot(aes(lambda)) + 
  geom_line(aes(y = density_tlr70, color = "Control"), linewidth = 1) +
  scale_color_manual("", values = unname(colors[1:2])) + 
  geom_vline(xintercept = 200, lty = "dashed") +
  theme_classic() +
  scale_y_continuous(limits = c(0,3.5), expand = c(0,0)) + 
  scale_x_continuous(limits = range(lambda_density$lambda), expand = c(0,0)) + 
  theme(axis.line.y.left = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 20)) + 
  labs(x = bquote(lambda[i]))  +
  theme(legend.title = element_blank())
# theme(legend.position.inside = c(1,1), legend.justification = c(1,1), legend.position = "inside",
#       legend.background = element_blank(), legend.title = element_blank())

p2 <- p1 + 
  geom_line(aes(y = density_tlr71, color = "TLR7"), linewidth = 1) +
  geom_vline(xintercept = 200*exp(popparam[["beta_lambda_tlr7_1"]]), lty = "dashed", color = colors[2])

p2.1 <- p2 + 
  geom_rug(data = parameters %>% filter(tlr7 == 0),
           aes(lambda_mode), linewidth = 0.7, length = unit(0.1, "npc"), alpha = 0.75) + 
  geom_rug(data = parameters %>% filter(tlr7 == 1),
           aes(lambda_mode), color = colors[2], linewidth = 0.7, length = unit(0.1, "npc"), alpha = 0.75)

p1.1 <-  p1 + geom_rug(data = parameters %>% filter(tlr7 == 0),
                       aes(lambda_mode), linewidth = 0.7, length = unit(0.1, "npc"), alpha = 0.75) 

ggsave("presentation figures/lambda_density.pdf", p1, width = 4, height = 2)
ggsave("presentation figures/lambda_density_data.pdf", p1.1, width = 4, height = 2)
ggsave("presentation figures/lambda_density_treatment.pdf", p2, width = 4, height = 2)
ggsave("presentation figures/lambda_density_treatment_data.pdf", p2.1, width = 4, height = 2)


## Example inference plots 2 ####
lambda_density <- tibble(
  lambda = seq(100, 500, by = 1),
  density_control = dnorm(log(lambda), log(200), 0.14),
  density_tlr7 = dnorm(log(lambda), log(200) , 0.14),
  density_vac = dnorm(log(lambda), log(200) +  0.30508546   , 0.14),
  density_tlr7_vac = dnorm(log(lambda), log(200) +  0.30508546  + 0.2 , 0.14)
)

sim_param <- tibble(id = 1:40, 
                    log_lambda = c(rnorm(20, log(200), 0.14), 
                                   rnorm(10, log(200) + 0.30508546, 0.14), 
                                   rnorm(10, log(200) + 0.30508546  + 0.2, 0.14)),
                    lambda = exp(log_lambda),
                    treatment = rep(c("Control", "TLR7", "Vaccine", "Vaccine + TLR7"), each = 10)
                    )

sim_param_control <- tibble(id = 1:40, 
                    log_lambda = rnorm(40, log(200), 0.14), 
                    lambda = exp(log_lambda),
)


p1 <- lambda_density %>% 
  ggplot(aes(lambda)) + 
  geom_line(aes(y = density_control, color = "Control"), linewidth = 1) +
  scale_color_manual("", values = unname(colors)) + 
  geom_vline(xintercept = 200, lty = "dashed") +
  theme_classic() +
  scale_y_continuous(limits = c(0,3.5), expand = c(0,0)) + 
  scale_x_continuous(limits = range(lambda_density$lambda), expand = c(0,0)) + 
  theme(axis.line.y.left = element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), 
        axis.title.y = element_blank(), axis.title.x = element_text(size = 20)) + 
  labs(x = bquote(lambda[i]))  +
  theme(legend.title = element_blank())
# theme(legend.position.inside = c(1,1), legend.justification = c(1,1), legend.position = "inside",
#       legend.background = element_blank(), legend.title = element_blank())

p2 <- p1 + 
  geom_line(aes(y = density_tlr7, color = "TLR7"), linewidth = 1, lty = "dashed") +
  geom_line(aes(y = density_vac, color = "Vaccine"), linewidth = 1) +
  geom_line(aes(y = density_tlr7_vac, color = "Vaccine + TLR7"), linewidth = 1) +
  geom_vline(xintercept = 200*exp(0.30508546), lty = "dashed", color = colors[3]) + 
  geom_vline(xintercept = 200*exp(0.30508546  + 0.2), lty = "dashed", color = colors[4]) + 
  guides(color = guide_legend(override.aes = list(lty = "solid")))

p2.1 <- p2 + 
  geom_rug(data = sim_param,
           aes(lambda, color = treatment), linewidth = 0.7, length = unit(0.1, "npc"), alpha = 0.75) 

p1.1 <-  p1 + geom_rug(data = sim_param_control,
                       aes(lambda), linewidth = 0.7, length = unit(0.1, "npc"), alpha = 0.75) 

ggsave("presentation figures/lambda_density.pdf", p1, width = 4, height = 2)
ggsave("presentation figures/lambda_density_data.pdf", p1.1, width = 4, height = 2)
ggsave("presentation figures/lambda_density_treatment.pdf", p2, width = 4, height = 2)
ggsave("presentation figures/lambda_density_treatment_data.pdf", p2.1, width = 5, height = 2)

## Example practical identifiability issues: ####

load("results_eclipse/study2/IGS_no_cov/IGS.RData")

estimate_comparison <- plot_estimates(tabestimates, parameter = "lambda_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(lambda[pop])) + 
  coord_cartesian(ylim = c(0, 100)) +
  plot_estimates(tabestimates, parameter = "betas_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(beta[pop])) + 
  coord_cartesian(ylim = c(0, 40)) +
  plot_estimates(tabestimates, parameter = "ta_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(ta[pop])) + 
  coord_cartesian(ylim = c(0, 2)) +
  plot_estimates(tabestimates, parameter = "NPs_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(NP[pop])) + 
  coord_cartesian(ylim = c(0, 10)) +
  plot_estimates(tabestimates, parameter = "p_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(p[pop])) + 
  coord_cartesian(ylim = c(0, 3)) +
  plot_estimates(tabestimates, parameter = "m_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(m[pop])) + 
  coord_cartesian(ylim = c(0, 100)) +
  # plot_estimates(tabestimates, parameter = "m_pop") + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) 


estimate_CV <- tabestimates %>% 
  # filter(parameter %in% c("lambda_pop", "betas_pop", "NPs_pop", "p_pop")) %>% 
  filter(grepl("pop", parameter)) %>% 
  group_by(parameter) %>% summarise(CV = sd(estimate)/mean(estimate)*100) %>% 
  mutate(parameter = factor(parameter, 
                            levels = rev(c("lambda_pop", "betas_pop", "ta_pop", "NPs_pop", "p_pop", "m_pop")))) %>% 
  # labels = rev(c(bquote(lambda[pop]), bquote(beta[pop]), bquote("t"[a][pop]), bquote(NP[pop]), bquote(p[pop]), bquote(m[pop]))))) %>%
  ggplot(aes(CV, factor(parameter))) + geom_col(fill = "darkblue", alpha = 0.75) + 
  xlim(c(0,100)) + 
  labs(x = "Coefficient of variation for estimates\nof population parameters (%)", y = "") + 
  scale_y_discrete(labels = rev(c(bquote(lambda[pop]), bquote(beta[pop]), bquote("ta"[pop]), bquote(NP[pop]), bquote(p[pop]), bquote(m[pop])))) +
  theme(panel.grid.major.y = element_blank(), 
        axis.text.y = element_text(size = 15), axis.text.x = element_text(size = 15), axis.title = element_text(size = 20))

estimate_correlations <- correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)), label_threshold = 0.7) + 
  theme(plot.title = element_blank()) + 
  scale_x_discrete(labels = c(bquote(lambda[pop]), bquote(beta[pop]), bquote("ta"[pop]), bquote(NP[pop]), bquote(p[pop]), bquote(m[pop]))) +
  scale_y_discrete(labels = c(bquote(lambda[pop]), bquote(beta[pop]), bquote("ta"[pop]), bquote(NP[pop]), bquote(p[pop]), bquote(m[pop]))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text = element_text(size = 15))

ggsave("presentation figures/estimate_comparison.pdf", estimate_comparison, width = 7, height = 8)
ggsave("presentation figures/estimate_CV.pdf", estimate_CV, width = 6, height = 5)
ggsave("presentation figures/estimate_correlation.pdf", estimate_correlations, width = 6, height = 5)

rbind(
  filter_results(tabestimates, remove_bicc_outliers = T, remove_estimate_outliers = F)$tabestimates %>% 
    filter(grepl("pop", parameter)) %>% mutate(method = "no fixed"),
  filter_results(tabestimates_fixed, remove_bicc_outliers = T, remove_estimate_outliers = F)$tabestimates %>% 
    filter(grepl("pop", parameter)) %>% mutate(method = "fixed")
) %>% 
  group_by(parameter) %>% mutate(norm_estimate = estimate/median(estimate)) %>% 
  ggplot(aes(norm_estimate, interaction(method, factor(parameter)))) + 
  geom_boxplot() + 
  # ggbeeswarm::geom_beeswarm(aes(color = bicc), alpha = 0.75) + 
  scale_x_log10(breaks = c(0.2, 0.5, 1, 2, 5), limits = c(0.2, 5)) + 
  color_bicc()


## Example corrected identifiability issues: ####

load("results_eclipse/study2/IGS_no_cov_fixed_param/IGS.RData")

estimate_comparison_fixed <- plot_estimates(distinct(tabestimates, bicc, run) %>% 
                                              mutate(estimate = 50, parameter = "lambda_pop", rses = NA), 
                                            parameter = "lambda_pop", remove_bicc_outliers = T) + 
  ggtitle(bquote(lambda[pop])) +
  coord_cartesian(ylim = c(0, 100)) +
  plot_estimates(tabestimates, parameter = "betas_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(beta[pop])) + 
  coord_cartesian(ylim = c(0, 40)) +
  plot_estimates(tabestimates, parameter = "ta_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(ta[pop])) + 
  coord_cartesian(ylim = c(0, 2)) +
  plot_estimates(tabestimates, parameter = "NPs_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(NP[pop])) + 
  coord_cartesian(ylim = c(0, 10)) +
  plot_estimates(distinct(tabestimates, bicc, run) %>% 
                   mutate(estimate = 0.7, parameter = "p_pop", rses = NA), 
                 parameter = "p_pop", remove_bicc_outliers = T) + 
  ggtitle(bquote("p"[pop])) +
  coord_cartesian(ylim = c(0, 3)) +
  plot_estimates(tabestimates, parameter = "m_pop", remove_bicc_outliers = T) +
  ggtitle(bquote(m[pop])) + 
  coord_cartesian(ylim = c(0, 100)) +
  # plot_estimates(tabestimates, parameter = "m_pop") + 
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right", plot.title = element_text(hjust = 0.5), 
        axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank()) 


estimate_CV_fixed <- estimate_CV + 
  geom_col(data = tabestimates %>% filter(grepl("pop", parameter)) %>% 
             group_by(parameter) %>% summarise(CV = sd(estimate)/mean(estimate)*100) %>% 
             mutate(parameter = factor(parameter, 
                                       levels = rev(c("lambda_pop", "betas_pop", "ta_pop", "NPs_pop", "p_pop", "m_pop")))),
           fill = "orange", alpha = 0.75
  )

estimate_correlations_fixed <- correlate_estimates(tabestimates %>% filter(grepl("pop", parameter)), label_threshold = 0.7) + 
  theme(plot.title = element_blank()) + 
  scale_x_discrete(labels = c(bquote(beta[pop]), bquote("ta"[pop]), bquote(NP[pop]), bquote(m[pop]))) +
  scale_y_discrete(labels = c(bquote(beta[pop]), bquote("ta"[pop]), bquote(NP[pop]), bquote(m[pop]))) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5), axis.text = element_text(size = 15))

compare_estimates(tabestimates %>% filter(grepl("pop", parameter))) + 
  theme(plot.title = element_blank())

ggsave("presentation figures/estimate_comparison_fixed.pdf", estimate_comparison_fixed, width = 7, height = 8)
ggsave("presentation figures/estimate_CV_fixed.pdf", estimate_CV_fixed, width = 6, height = 5)
ggsave("presentation figures/estimate_correlation_fixed.pdf", estimate_correlations_fixed, width = 6, height = 5)

## Example of SCM sensitivity: ####

summary <- read_csv("~/Downloads/model_building_summary.csv")

plot_df <- summary %>% 
  filter(n_covar == 3) %>% 
  mutate(treatment = factor(treatment, 
                            levels = c("intTLR7VAC", "vac", "tlr7"),
                            labels = c("TLR7*Vaccine", "Vaccine", "TLR7")))


plot_df <- plot_df %>% select(-LL) %>% pivot_wider(names_from = param, values_from = BICc) %>% 
  pivot_longer(!c(treatment, n_covar, threshold), names_to = "param", values_to = "BICc") %>% 
  mutate(BICc = ifelse(is.na(BICc), threshold, BICc),
         label = case_when(threshold - BICc <= 0 ~ "",
                           param == "NPs" ~ "'N'[p]",
                           param == "ta" ~ "'t'[a]",
                           .default = param),
         param = factor(param, levels = rev(c("lambda", "betas", "ta", "NPs", "p", "m")))) 

SCM_example_no_bar <- plot_df %>% 
  ggplot(aes(0, treatment, fill = param)) + 
  # annotate(geom = "rect", 
  #          xmax = max(plot_df$threshold - plot_df$BICc),
  #          xmin = max(plot_df$threshold - plot_df$BICc) - 2, 
  #          ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.25, color = NA) +
  geom_col(position = "dodge", width = 0.7) +
  # geom_text(aes(label = label, color = param), 
  #           position = position_dodge(width = 0.7), hjust = -0.2, 
  #           size = 3.5, show.legend = F, parse= T, fontface = 2) +
  geom_vline(xintercept = 0, color = "darkgray") +
  labs(x = "Improvement to model fit (BICc)", y = "Treatment effect added to model", fill = "Parameter") + 
  scale_fill_manual(values = rev(c("#e97131", "#f3b699", "#156082", "#1a6b25", "#7ca67d", "#bed3bf")),
                    guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values = rev(c("#e97131", "#f3b699", "#156082", "#1a6b25", "#7ca67d", "#bed3bf")),
                     guide = guide_legend(reverse = TRUE)) +
  theme(panel.grid.minor.x = element_blank(), axis.title = element_text(size = 15)) +
  coord_cartesian(xlim = c(-5, 5))

ggsave("presentation figures/SCM_example_no_bars.pdf", SCM_example_no_bar, width = 6, height = 5)

SCM_example <- plot_df %>% 
  ggplot(aes(threshold - BICc, treatment, fill = param)) + 
  annotate(geom = "rect", 
           xmax = max(plot_df$threshold - plot_df$BICc),
           xmin = max(plot_df$threshold - plot_df$BICc) - 2, 
           ymin = -Inf, ymax = Inf, fill = "gray", alpha = 0.25, color = NA) +
  geom_col(position = "dodge", width = 0.7) +
  geom_text(aes(label = label, color = param), 
            position = position_dodge(width = 0.7), hjust = -0.2, 
            size = 3.5, show.legend = F, parse= T, fontface = 2) +
  geom_vline(xintercept = 0, color = "darkgray") +
  labs(x = "Improvement to model fit (BICc)", y = "Treatment effect added to model", fill = "Parameter") + 
  scale_fill_manual(values = rev(c("#e97131", "#f3b699", "#156082", "#1a6b25", "#7ca67d", "#bed3bf")),
                    guide = guide_legend(reverse = TRUE)) + 
  scale_color_manual(values = rev(c("#e97131", "#f3b699", "#156082", "#1a6b25", "#7ca67d", "#bed3bf")),
                     guide = guide_legend(reverse = TRUE)) +
  theme(panel.grid.minor.x = element_blank(), axis.title = element_text(size = 15)) +
  coord_cartesian(xlim = c(-5, 5))


ggsave("presentation figures/SCM_example.pdf", SCM_example, width = 6, height = 5)
