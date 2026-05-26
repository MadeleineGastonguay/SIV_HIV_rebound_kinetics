library(tidyverse)
library(here)
library(patchwork)
library(ggh4x)
library(ggnewscale)

## Load practical identifiability data for study 2 (TLR7 + Vaccine in SIV)

# Source scripts for initializing monolix and functions to run the analysis
source(here("scripts", "initialize_monolix.R"))
source(here("scripts", "functions.R"))

theme_set(theme_minimal())


# Project setup ####
## Results folder
study <- "study2"
project_name <- "IGS_no_cov.mlxtran"
out_folder <- here("results_eclipse", study)
if(!dir.exists(out_folder)) dir.create(out_folder)

## Data file
data_file <- here("data", "derived", str_interp("${toupper(study)}_DATA.txt"))
data <- read_delim(data_file, delim = "\t")


# Load Project
loadProject(here(out_folder, project_name))  
# Load results from IGS over orders of magnitude
load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "log10_IGS.RData"))
log10_tabestimates <- tabestimates
log10_initialvalues <- initialvalues

## Somehow the tabestimates from narrow range got concatenated to the log10 tabestimates, so I will filter those out:
log10_tabestimates <- log10_tabestimates %>% mutate(run2 = rep(1:150, each = 13)) %>% filter(run2 <= 100) %>% select(-run2)

# Load results from IGS over narrow range
load(here(out_folder, str_remove(project_name, "\\.mlxtran$"), "IGS.RData"))

fixed_results <- new.env()
load(here(out_folder, "IGS_no_cov_fixed_param", "IGS.RData"), envir = fixed_results)

narrow_range <- read_csv(here(out_folder, "initialValue_ranges_narrow.csv"))

### For plotting:

labels <- c(lambda_pop = bquote(lambda['pop']),
            betas_pop = bquote(beta['pop']),
            ta_pop = bquote(ta['pop']),
            m_pop = bquote(m['pop']),
            p_pop = bquote(p['pop']),
            NPs_pop = bquote(NP['pop']),
            omega_lambda = bquote(Omega[lambda]),
            omega_betas = bquote(Omega[beta]),
            omega_m = bquote(Omega['m']),
            omega_p = bquote(Omega['p']),
            omega_ta = bquote(Omega['ta']),
            omega_NPs = bquote(Omega['NP']))

### Show selection of initial conditions with lambda as an example:
iv_range <- narrow_range %>% filter(parameter == "lambda_pop")
tab_range <- tabestimates %>% filter(parameter == "lambda_pop") %>% summarise(min = min(estimate), max = max(estimate))
# bicc_threshold <- log10_tabestimates %>% filter(parameter == "lambda_pop", between(estimate, iv_range$min, iv_range$max)) %>% 
#   pull(bicc) %>% max

plot_bicc_v_estimate(log10_tabestimates %>% filter(bicc < 640), parameter = "lambda_pop", log10_x = T) + 
  annotate(geom = "rect", xmin = iv_range$min, xmax = iv_range$max, ymin = -Inf, ymax = Inf, 
           fill = "lightgrey", alpha = 0.25) +
  scale_fill_viridis_c(option = "A", direction = -1) +
  labs(title = "Fitting results from initial guess sampling to find region of global minimum") +
  plot_bicc_v_estimate(tabestimates, parameter = "lambda_pop", log10_x = F) + 
  coord_cartesian(xlim = unlist(iv_range[1,c(2,3)])) +
  labs(title = "Fitting results from initial guess sampling around global minimum") +
  plot_layout(ncol = 1) &
  theme(strip.text = element_blank(), axis.title.x = element_text(size = 15)) &
  labs(x = bquote(lambda['pop']))



## Parameter identifiability:

# 1) Subset the above IGS to fits with similar BICc, according to BICw

BICw_log10 <- calculate_BICw(log10_tabestimates)
runs_log10 <- BICw_log10 %>% filter(cumprob <= 0.99) %>% pull(run) 


BICw <- calculate_BICw(tabestimates)
runs <- calculate_BICw(tabestimates) %>% filter(cumprob <= 0.95) %>% pull(run) 
BICw %>% ggplot(aes(bicc, cumprob, color = run %in% runs)) + geom_point() + 
  scale_color_manual(values =c( "gray","darkred")) + 
  labs(color = "Similar BICc to\nbest model")


log10_plot <-  plot_bicc_v_estimate(log10_tabestimates %>% mutate(estimate = log10(estimate)), 
                                    parameter = "lambda_pop", log10_x = F, color = runs_log10) + 
  scale_x_continuous(#transform = "log10", 
    labels = scales::label_math(expr = 10^.x),
    breaks = 1:5, limits = c(0.99,5.7)) +
  labs(title = "Initial guess sampling to find region of global minimum", x = bquote(lambda['pop']))  +
  guides(fill = guide_legend(position = "inside", override.aes = list(shape = 22, size = 10, alpha = 1)), shape = "none") +
  theme(strip.text = element_blank(), plot.title = element_blank(), legend.position.inside = c(1,1),legend.justification = c(1,1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.background = element_rect(color = "white"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  
  BICw_log10 %>% mutate(run = fct_inorder(as.character(run))) %>% 
  ggplot(aes(run, cumprob, fill = run %in% runs_log10, color = run %in% runs_log10)) + geom_col() + 
  scale_fill_manual(values =c( "gray","darkred")) + 
  scale_color_manual(values =c( "gray","darkred")) + 
  labs(fill = "Near global minimum", y = "Cumulative Probability") + 
  guides(color = "none") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none"
  )  + 
  
  distinct(log10_tabestimates , run, bicc) %>% 
  ggplot(aes("", bicc)) + 
  geom_quasirandom(aes(fill = run %in% runs_log10),
                   show.legend = F, color = "black", alpha = 0.75, shape = 21, size = 2) + 
  scale_fill_manual(values =c( "gray","darkred")) +
  labs(x = "", y = "BICc") +
  theme(axis.text.y  = element_text(size = 12), axis.title.y = element_text(size = 15)) +
  
  # plot_layout(nrow = 1, guides = "collect", widths = c(1.5,0.5,0.5)) &
  plot_layout(nrow = 1, widths = c(1.5,0.5,0.5)) &
  labs(fill = "Near global mininum") 



narrow_plot <-  plot_bicc_v_estimate(tabestimates, 
                                     parameter = "lambda_pop", log10_x = F, color = runs) + 
  # scale_x_continuous(#transform = "log10", 
  #   labels = scales::label_math(expr = 10^.x),
  #   breaks = 1:5, limits = c(0.99,5.7)) +
  coord_cartesian(xlim = unlist(iv_range[1,c(2,3)])) +
  scale_fill_manual(values =c( "gray","darkgreen")) + 
  labs(title = "Initial guess sampling to find region of global minimum", x = bquote(lambda['pop']))  +
  guides(fill = guide_legend(position = "inside", override.aes = list(shape = 22, size = 10, alpha = 1)), shape = "none") +
  theme(strip.text = element_blank(), plot.title = element_blank(), legend.position.inside = c(1,1),legend.justification = c(1,1),
        axis.text = element_text(size = 12), axis.title = element_text(size = 15), legend.background = element_rect(color = "white"),
        legend.title = element_text(size = 15), legend.text = element_text(size = 12)) +
  
  BICw %>% mutate(run = fct_inorder(as.character(run))) %>% 
  ggplot(aes(run, cumprob, fill = run %in% runs, color = run %in% runs)) + geom_col() + 
  scale_fill_manual(values =c( "gray","darkgreen")) + 
  scale_color_manual(values =c( "gray","darkgreen")) + 
  labs(fill = "Near global minimum", y = "Cumulative Probability") + 
  guides(color = "none") +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.position = "none"
  )  + 
  
  distinct(tabestimates, run, bicc) %>% 
  ggplot(aes("", bicc)) + 
  geom_quasirandom(aes(fill = run %in% runs),
                   show.legend = F, color = "black", alpha = 0.75, shape = 21, size = 2) + 
  scale_fill_manual(values =c( "gray","darkgreen")) +
  labs(x = "", y = "BICc")  +
  theme(axis.text.y  = element_text(size = 12), axis.title.y = element_text(size = 15)) +
  
  # plot_layout(nrow = 1, guides = "collect", widths = c(1.5,0.5,0.5)) &
  plot_layout(nrow = 1, widths = c(1.5,0.5,0.5)) &
  labs(fill = "Set of best models") 

## Define initial value ranges
initialValue_ranges <- tibble(
  lambda = c(10,1000), 
  betas = c(0.001,100), # Will be multiplied by 10^-7 in the model 
  m = c(0.01, 100),
  p = c(0.01,10), 
  NPs = c(10^2, 10^7)*10^-5, # Will be multiplied by 10^5 in the model
  ta = c(10^-2,10)
) 

parameter_range_plot <- plot_parameter_ranges_v2(initialValue_ranges, log10_tabestimates, tabestimates, p_threshold = 0.95) + 
  theme(axis.text = element_text(size = 12),  legend.text = element_text(size = 12), legend.title = element_blank(), strip.text = element_text(size = 15))


IC_plot <- log10_plot / narrow_plot / free(parameter_range_plot + theme(legend.position = "none"))  + plot_annotation(tag_levels = "a") 

ggsave(here("paper_figures", "setting_ICs_supplemental_figure_mockup.pdf"), IC_plot, width = 15, height = 10)

narrow_range <- read_csv(here("results_eclipse", "study2", "initialValue_ranges_narrow.csv"))

labels2 <- labels
labels2$betas_pop <- bquote(beta['pop']~'(x'~10^-7~')')
labels2$NPs_pop <- bquote(NP['pop']~'(x'~10^5~')')

IC_plot2 <- free(BICw_log10 %>% mutate(run = fct_inorder(as.character(run))) %>% 
  ggplot(aes(run, cumprob, fill = run %in% runs_log10, color = run %in% runs_log10)) + geom_col() + 
  scale_fill_manual(values =c( "gray","darkred")) + 
  scale_color_manual(values =c( "gray","darkred")) + 
  labs(fill = "Near global minimum", y = "Cumulative Probability", x = "Model fit") + 
    guides(color = "none", fill = guide_legend(position = "bottom")) +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_blank(),
          panel.grid.minor.y = element_blank(),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 15),
          legend.title.position = "top"
          # legend.position.inside = c(0,1), legend.justification = c(0,1)
    ))  + 
  
  
  plot_bicc_v_estimate(log10_tabestimates %>% 
                       mutate(estimate = log10(estimate),
                              # estimate = ifelse(parameter == "NPs_pop",log10(estimate*1e5), log10(estimate)),
                              parameter = factor(parameter, levels = names(labels), labels = unlist(labels2))) %>% 
                       filter(grepl("pop",parameter), bicc < 640), 
                     color = runs_log10, point_size = 3, add_best = FALSE, log10_x = FALSE) +
  facet_wrap(~param, nrow = 1, scales = "free_x", strip.position = "bottom", labeller = "label_parsed") +
  ggh4x::facetted_pos_scales( x = list(
    scale_x_continuous(breaks = 1:5, labels = scales::label_math(expr = 10^.x)),
    scale_x_continuous(breaks = -2:2, labels = scales::label_math(expr = 10^.x)),
    scale_x_continuous(breaks = log10(c(0.5, 1, 2)), labels = function(breaks) 10^breaks, limits = log10(c(0.5,2.2))),
    scale_x_continuous(breaks = seq(-8,4,by=4), labels = scales::label_math(expr = 10^.x), limits = c(-8.2,4)),
    scale_x_continuous(breaks = -1:2, labels = function(breaks) 10^breaks),
    scale_x_continuous(breaks = 0:2, labels = function(breaks) 10^breaks)
    )) +
  theme( axis.title.x = element_blank(), legend.position = "none") +
 

  free(BICw %>% mutate(run = fct_inorder(as.character(run))) %>% 
  ggplot(aes(run, cumprob, fill = run %in% runs, color = run %in% runs)) + geom_col() + 
  scale_fill_manual(values =c( "gray","darkgreen")) + 
  scale_color_manual(values =c( "gray","darkgreen")) + 
  labs(fill = "Set of best models", y = "Cumulative Probability", x = "Model fit") + 
  guides(color = "none", fill = guide_legend(position = "bottom")) +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.title = element_text(size = 15),
        legend.title.position = "top"
        # legend.position.inside = c(0,1), legend.justification = c(0,1)
  ))  + 
  
  
  plot_bicc_v_estimate(tabestimates %>% mutate(#estimate = ifelse(parameter == "NPs_pop",(estimate*1e5), (estimate)), 
                                               parameter = factor(parameter, levels = names(labels), labels = unlist(labels2))) %>% 
                         filter(grepl("pop",parameter)), 
                       color = runs, point_size = 3, log10_x = FALSE, add_best = FALSE) + 
  scale_fill_manual(values = c( "gray","darkgreen")) +
  geom_segment(data = mutate(filter(narrow_range, grepl("pop", parameter)),
                             # min = ifelse(parameter == "NPs_pop", min*1e5, min),
                             # max = ifelse(parameter == "NPs_pop", max*1e5, max),
                             param = factor(parameter, levels = names(labels), labels = unlist(labels2))),
               aes(y = Inf, yend = Inf, x = min, xend = max), linewidth = 3, inherit.aes = F, color = "darkred") +
  facet_wrap(~param, nrow = 1, scales = "free_x", strip.position = "bottom", labeller = "label_parsed") +
  theme( axis.title.x = element_blank(), legend.position = "none") +
  coord_cartesian(ylim = c(568,586)) +
  
  plot_layout(ncol = 2, widths = c(1,5)) &
  plot_annotation(tag_levels = "a") &
  theme(strip.placement = "outside", axis.text = element_text(size = 12), axis.title = element_text(size = 15),
        plot.title = element_blank(), strip.text = element_text(size = 15), plot.tag = element_text(size = 20))

ggsave(here("paper_figures", "setting_ICs_supplemental_figure_mockup2.pdf"), IC_plot2, width = 12, height = 7)

# Use 50% of runs from wide range
length(runs_log10)
# Use 70% of runs for narrow range
length(runs)

# get_ICs <- plot_bicc_v_estimate(log10_tabestimates %>% mutate(run = rep(1:150, each = 13)) %>% filter(bicc < 640), parameter = "lambda_pop", log10_x = T, color = runs_log10) + 
#   labs(title = "Initial guess sampling to find region of global minimum") +
#   plot_bicc_v_estimate(tabestimates, parameter = "lambda_pop", log10_x = F, color = runs) + 
#   coord_cartesian(xlim = unlist(iv_range[1,c(2,3)])) +
#   labs(title = "Initial guess sampling around global minimum") +
#   plot_layout(ncol = 1, guides = "collect") &
#   theme(strip.text = element_blank(), axis.title.x = element_text(size = 15)) &
#   labs(x = bquote(lambda['pop'])) &
#   guides(fill = "none")
# 
# ggsave(here("paper_figures", "setting_ICs_supplemental_figure_mockup.pdf"), get_ICs, width = 7, height = 7)


tabestimates_95 <- tabestimates %>% filter(run %in% runs)

# 2) Look at distribution of parameters for these estimates

df <- merge(
  tabestimates_95 %>% 
    filter(grepl("pop", parameter)) %>% 
    mutate(parameter = gsub("_pop", "", parameter)) %>% 
    select(run, parameter, mean = estimate, bicc),
  
  tabestimates_95 %>% 
    filter(grepl("omega", parameter)) %>% 
    mutate(parameter = gsub("omega_", "", parameter)) %>% 
    select(run, parameter, sd = estimate, bicc)
) 

density_data <- df %>%
  group_by(parameter) %>%
  group_modify(~ {
    x_vals <- exp(seq((min(log(.x$mean) - 4*.x$sd)),
                      (max(log(.x$mean) + 4*.x$sd)),
                      length.out = 500))
    do.call(rbind, lapply(1:nrow(.x), function(i) {
      data.frame(
        x = x_vals,
        y = dnorm(log(x_vals), mean = log(.x$mean[i]), sd = .x$sd[i]),
        run = as.factor(.x$run[i])
      )
    }))
  }) %>%
  ungroup()

distribution_plot <- density_data %>%
  left_join(tabestimates_95 %>% distinct(run, bicc) %>% mutate(run = factor(run))) %>%
  mutate(x = case_when(#parameter == "betas" ~ x*1e-7, 
    parameter == "NPs" ~ x*1e5,
    .default = x),
    parameter = factor(parameter, 
                       levels =  c("lambda","betas", "ta", "m", "p", "NPs"),
                       labels =  c("lambda","beta~'(x'~10^-7~')'", "t['a']", "m", "p", "N['p']"))) %>%
  ggplot(aes(x = log10(x), y = y, color = bicc, group = run)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free", nrow = 1, labeller = label_parsed, strip.position = "bottom") +
  labs(title = "Normal Distributions by Parameter",
       x = "Value", y = "Density") +
  theme_minimal() + color_bicc() +
  # scale_x_log10() +
  scale_x_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) +
  # facetted_pos_scales(
  #   x = list(
  #     "lambda" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(10, 30, 100, 300))),
  #     "beta" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(5, 10, 30, 50)))
  #   )
  # ) +
  facetted_pos_scales(
    x = list(
      "lambda" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(10, 30, 100, 300)), limits = log10(c(10,350))),
      "beta~'(x'~10^-7~')'" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(5, 10, 30, 50))),
      "t['a']" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,2)),
      "m" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,4,by=2), limits = c(-2.5,4.5)),
      "p" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-1,1))
    )
  ) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), panel.grid = element_blank(), axis.line.x = element_line(),
        axis.text.x = element_text(size = 15),
        strip.placement = "outside", strip.text = element_text(size = 20))

distribution_plot_v2 <- plot_estimates(tabestimates_95, parameter = "lambda_pop") + coord_cartesian(ylim = c(0, 100)) + 
  theme(axis.title.y = element_text(size = 15)) +
  labs(title = bquote(lambda["pop"])) + 
  plot_estimates(tabestimates_95, parameter = "betas_pop") + coord_cartesian(ylim = c(0, 40)) +
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(beta["pop"])) + 
  plot_estimates(tabestimates_95, parameter = "ta_pop") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(ta["pop"])) + 
  plot_estimates(tabestimates_95, parameter = "m_pop") + coord_cartesian(ylim = c(0, 100)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(m["pop"])) + 
  plot_estimates(tabestimates_95, parameter = "p_pop") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(p["pop"])) + 
  plot_estimates(tabestimates_95, parameter = "NPs_pop") + coord_cartesian(ylim = c(0, 7.5))  + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(NP["pop"])) + 
  plot_estimates(tabestimates_95, parameter = "omega_lambda") + coord_cartesian(ylim = c(0, 1)) + 
  theme(axis.title.y = element_text(size = 15)) +
  labs(title = bquote(Omega[lambda])) + 
  plot_estimates(tabestimates_95, parameter = "omega_betas") + coord_cartesian(ylim = c(0, 1)) +
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega[beta])) + 
  plot_estimates(tabestimates_95, parameter = "omega_ta") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega['ta'])) + 
  plot_estimates(tabestimates_95, parameter = "omega_m") + coord_cartesian(ylim = c(0, 5)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega['m'])) + 
  plot_estimates(tabestimates_95, parameter = "omega_p") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega['p'])) + 
  plot_estimates(tabestimates_95, parameter = "omega_NPs") + coord_cartesian(ylim = c(0, 4.5)) +
  labs(title = bquote(Omega['NP'])) + 
  theme(axis.title.y = element_blank()) +
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "right", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), axis.text.x = element_blank()) 



evaluation_plot <- free(correlate_estimates(tabestimates_95 %>% filter(parameter != "a"), label_threshold = 0.7, text_color = "white") +
  scale_x_discrete(labels = labels) +
  scale_y_discrete(labels = labels) +
  theme(legend.title = element_text(size = 15), legend.text = element_text(size = 12), aspect.ratio = 1,
        plot.margin = margin(0,0,0,0))) +
  # theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + 
  
  tabestimates_95 %>% select(run, bicc, parameter, estimate) %>% 
  pivot_wider(names_from = parameter, values_from = estimate) %>% 
  ggplot(aes(lambda_pop, betas_pop, color = bicc)) + 
  geom_point(size = 3) + 
  color_bicc() +
  labs(x = bquote(lambda['pop']), y = bquote(beta['pop']~'(x'~10^-7~')')) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  
  tabestimates_95 %>% select(run, bicc, parameter, estimate) %>% 
  pivot_wider(names_from = parameter, values_from = estimate) %>% 
  ggplot(aes(p_pop, NPs_pop, color = bicc)) + 
  geom_point(size = 3) + 
  color_bicc() + 
  labs(x = bquote(p['pop']), y = bquote(NP['pop']~'(x'~10^5~')')) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  
  tabestimates_95 %>% select(run, bicc, parameter, estimate) %>% 
  pivot_wider(names_from = parameter, values_from = estimate) %>% 
  ggplot(aes(lambda_pop, m_pop, color = bicc)) + 
  geom_point(size = 3) + 
  color_bicc() +
  labs(x = bquote(lambda['pop']), y = bquote(m['pop'])) +
  theme(axis.title = element_text(size = 20),
        axis.text = element_text(size = 15)) +
  
  plot_layout(nrow = 1, widths = c(1.5,1,1,1)) 
  

# 2b) look at individual parameter value consistency

indivParams %>% 
  pivot_longer(!c(id, run), names_to = "parameter", values_to = "EBE") %>% 
  group_by(id, parameter) %>% 
  summarise(CV = sd(EBE)/median(EBE)*100) %>% 
  ggplot(aes(id, CV)) + 
  geom_col() +facet_wrap(~parameter) + ylim(c(0,100))

# 3) Fixed results:

fixed_BICw <- calculate_BICw(fixed_results$tabestimates)
fixed_runs <- fixed_BICw %>% filter(cumprob <= 0.95) %>% pull(run) 
fixed_BICw %>% ggplot(aes(bicc, cumprob, color = run %in% fixed_runs)) + geom_point() + 
  scale_color_manual(values =c( "gray","darkred")) + 
  labs(color = "Similar BICc to\nbest model")

rbind(indivParams %>% filter(run %in% runs) %>% mutate(set = "no fixed"),
      fixed_results$indivParams %>% filter(run %in% fixed_runs) %>% mutate(set = "fixed")) %>% 
  pivot_longer(!c(id, run, set), names_to = "parameter", values_to = "EBE") %>% 
  group_by(id, parameter, set) %>% 
  summarise(CV = sd(EBE)/median(EBE)*100) %>% 
  ggplot(aes(id, CV, fill = set)) + 
  geom_col(position = position_dodge()) + facet_wrap(~parameter)

rbind(indivParams %>% filter(run %in% runs) %>% mutate(set = "no fixed"),
      fixed_results$indivParams %>% filter(run %in% fixed_runs) %>% mutate(set = "fixed")) %>% 
  pivot_longer(!c(id, run, set), names_to = "parameter", values_to = "EBE") %>% 
  filter(parameter %in% c("m")) %>%
  group_by(id, set, parameter) %>% mutate(norm_EBE = EBE/median(EBE)) %>% 
  ggplot(aes(factor(id), norm_EBE, fill = set)) +
  geom_violin(scale = "width", color = NA) +
  # geom_boxplot(outlier.shape = NA) + 
  # geom_beeswarm(position = position_dodge()) + 
  scale_y_log10(limits = c(1/5, 5)) + facet_wrap(~parameter)

plot_bicc_v_estimate(fixed_results$tabestimates, log10_x = F, color = fixed_runs)

fixed_tabestimates_95 <- fixed_results$tabestimates %>% filter(run %in% fixed_runs)

# lingering correlation between ta and beta, but they seem to converge okay
correlate_estimates(fixed_tabestimates_95, label_threshold = 0.7) +
  
  fixed_tabestimates_95 %>% select(run, bicc, parameter, estimate) %>%
  pivot_wider(names_from = parameter, values_from = estimate) %>%
  ggplot(aes(ta_pop, betas_pop, color = bicc)) +
  geom_point(size = 2) +
  color_bicc() +
  
  plot_layout(nrow = 1, guides = "collect")

fixed_df <- merge(
  fixed_tabestimates_95 %>% 
    filter(grepl("pop", parameter)) %>% 
    mutate(parameter = gsub("_pop", "", parameter)) %>% 
    select(run, parameter, mean = estimate, bicc) %>% 
    rbind(distinct(fixed_tabestimates_95, run, bicc) %>% mutate(parameter = "lambda", mean = 70)) %>% 
    rbind(distinct(fixed_tabestimates_95, run, bicc) %>% mutate(parameter = "p", mean = 0.7)),
  
  fixed_tabestimates_95 %>% 
    filter(grepl("omega", parameter)) %>% 
    mutate(parameter = gsub("omega_", "", parameter)) %>% 
    select(run, parameter, sd = estimate, bicc)
) 

fixed_density_data <- fixed_df %>%
  group_by(parameter) %>%
  group_modify(~ {
    x_vals <- exp(seq((min(log(.x$mean) - 4*.x$sd)),
                      (max(log(.x$mean) + 4*.x$sd)),
                      length.out = 500))
    do.call(rbind, lapply(1:nrow(.x), function(i) {
      data.frame(
        x = x_vals,
        y = dnorm(log(x_vals), mean = log(.x$mean[i]), sd = .x$sd[i]),
        run = as.factor(.x$run[i])
      )
    }))
  }) %>%
  ungroup()

fixed_distribution_plot <- fixed_density_data %>%
  left_join(fixed_tabestimates_95 %>% distinct(run, bicc) %>% mutate(run = factor(run))) %>%
  mutate(x = case_when(#parameter == "betas" ~ x*1e-7, 
    parameter == "NPs" ~ x*1e5,
    .default = x),
    parameter = factor(parameter, 
                       levels =  c("lambda","betas", "ta", "m", "p", "NPs"),
                       labels =  c("lambda","beta~'(x'~10^-7~')'", "t['a']", "m", "p", "N['p']"))) %>%
  ggplot(aes(x = log10(x), y = y, color = bicc, group = run)) +
  geom_line(alpha = 0.5) +
  facet_wrap(~parameter, scales = "free", nrow = 1, labeller = label_parsed, strip.position = "bottom") +
  labs(title = "Normal Distributions by Parameter",
       x = "Value", y = "Density") +
  theme_minimal() + scale_color_viridis_c("BICc", option = "A", direction = -1, breaks = c(556, 557, 558, 559)) +
  scale_x_continuous(
    labels = scales::label_math(expr = 10^.x)
  ) +
  facetted_pos_scales(
    x = list(
      "lambda" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(10, 30, 100, 300)), limits = log10(c(10,350))),
      "beta~'(x'~10^-7~')'" = scale_x_continuous(labels = function(breaks) 10^breaks, breaks = log10(c(5, 10, 30, 50))),
      "t['a']" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,2)),
      "m" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-2,4,by=2), limits = c(-2.5,4.5)),
      "p" = scale_x_continuous(labels = scales::label_math(expr = 10^.x), breaks = seq(-1,1))
    )
  ) +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), panel.grid = element_blank(), 
        axis.line.x = element_line(), strip.placement = "outside", strip.text = element_text(size = 20))

fixed_distribution_plot_v2 <- plot_estimates(fixed_tabestimates_95 %>% filter(parameter == "a") %>% 
                                               mutate(parameter = "lambda_pop", estimate = 70, rses = 0) , 
                                             parameter = "lambda_pop") + coord_cartesian(ylim = c(0, 100)) + 
  theme(axis.title.y = element_text(size = 15)) +
  labs(title = bquote(lambda["pop"])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "betas_pop") + coord_cartesian(ylim = c(0, 40)) +
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(beta["pop"])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "ta_pop") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(ta["pop"])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "m_pop") + coord_cartesian(ylim = c(0, 100)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(m["pop"])) + 
  plot_estimates(fixed_tabestimates_95 %>% filter(parameter == "a") %>% 
                   mutate(parameter = "p_pop", estimate = 0.7, rses = 0) , 
                 parameter = "p_pop") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(p["pop"])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "NPs_pop") + coord_cartesian(ylim = c(0, 7.5))  + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(NP["pop"])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "omega_lambda") + coord_cartesian(ylim = c(0, 1)) + 
  theme(axis.title.y = element_text(size = 15)) +
  labs(title = bquote(Omega[lambda])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "omega_betas") + coord_cartesian(ylim = c(0, 1)) +
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega[beta])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "omega_ta") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega['ta'])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "omega_m") + coord_cartesian(ylim = c(0, 5)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega['m'])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "omega_p") + coord_cartesian(ylim = c(0, 2)) + 
  theme(axis.title.y = element_blank()) +
  labs(title = bquote(Omega['p'])) + 
  plot_estimates(fixed_tabestimates_95, parameter = "omega_NPs") + coord_cartesian(ylim = c(0, 4.5)) +
  labs(title = bquote(Omega['NP'])) + 
  theme(axis.title.y = element_blank()) +
  plot_layout(nrow = 2, guides = "collect") &
  theme(legend.position = "right", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), 
        panel.grid.minor.y = element_blank(), axis.text.x = element_blank()) &
  scale_color_viridis_c("BICc", option = "A", direction = -1, breaks = seq(556,559))



CV_plot <- rbind(tabestimates_95 %>% mutate(set = "free"), 
                 rbind(fixed_tabestimates_95) %>% mutate(set = "fixed")) %>% 
  filter(parameter != "a") %>% 
  select(set, parameter, estimate) %>% 
  group_by(set, parameter) %>% 
  summarise(CV = sd(estimate)/mean(estimate)*100) %>% 
  rbind(data.frame(set = "fixed", parameter = c("lambda_pop", "p_pop"), CV = 0)) %>% 
  mutate(set = factor(set, levels = c("free", "fixed")),
         parameter = factor(parameter, 
                            levels = c("lambda_pop","omega_lambda","betas_pop", "omega_betas", "ta_pop","omega_ta",  
                                       "m_pop", "omega_m", "p_pop", "omega_p",  "NPs_pop","omega_NPs"))) %>% 
  ggplot(aes(factor(parameter), CV, fill = set)) +
  geom_col(position = position_dodge()) + 
  scale_fill_manual(values = c("#2B9289", "#B43776"), labels = c(expression("No fixed parameters"), expression(lambda['pop']~"and"~p['pop']~"fixed"))) + 
  labs(x = "Parameter", y = "Coefficient of\nVariation (%)") +
  guides(fill = guide_legend(position = "inside")) + 
  theme(legend.position.inside = c(1,0.95), legend.justification = c(1,1), 
        legend.background = element_rect(color = NA), legend.title = element_blank(),
        axis.text.x = element_text(size = 15), panel.grid.minor.y = element_blank(), panel.grid.major.x = element_blank()) + 
  scale_x_discrete(labels = labels) 
# axis.text.x = element_text(angle = 45, hjust = 1)) 

# Add in parameter estimate boxplots:
estimate_boxplots_plot <- rbind(tabestimates_95 %>% mutate(set = "free"), 
                           rbind(fixed_tabestimates_95) %>% mutate(set = "fixed")) %>% 
  filter(parameter != "a") %>% 
  rbind(fixed_tabestimates_95 %>% filter(parameter == "betas_pop") %>% mutate(parameter = "lambda_pop", estimate = 70, set= "fixed")) %>% 
  rbind(fixed_tabestimates_95 %>% filter(parameter == "betas_pop") %>% mutate(parameter = "p_pop", estimate = 0.7, set= "fixed")) %>% 
  select(set, parameter, estimate) %>% 
  group_by(set, parameter) %>% mutate(norm_estimate = estimate/median(estimate)) %>% 
  mutate(set = factor(set, levels = c("free", "fixed"), labels = c("No Fixed Parameters", "Lambda and p means fixed")),
         parameter = factor(parameter, 
                            levels = c("lambda_pop","omega_lambda","betas_pop", "omega_betas", "ta_pop","omega_ta",  
                                       "m_pop", "omega_m", "p_pop", "omega_p",  "NPs_pop","omega_NPs"))) %>% 
  ggplot(aes(factor(parameter), norm_estimate, color = set)) + 
  geom_jitter(position = position_jitterdodge(), alpha = 0.5) + 
  geom_boxplot(fill = NA, outlier.shape = NA) +
  scale_y_log10(limits = c(1/5, 5), breaks = c( 1/5, 1/2, 1,  2, 5)) + 
  guides(color = guide_legend(position = "inside")) + 
  theme(panel.grid.minor.y = element_blank(),
        legend.position.inside = c(1,1), legend.justification = c(1,1), 
        legend.background = element_rect(color = NA), legend.title = element_blank(),
        legend.text = element_text(size = 12), axis.text.x = element_text(size = 15)) + 
  scale_x_discrete(labels = labels) +
  scale_color_manual(values = c("#2B9289", "#B43776"), labels = c(expression("No fixed parameters"), expression(lambda['pop']~"and"~p['pop']~"fixed")))  + 
  labs(x = "Parameter", y = "Normalized\nEstimate") 

plot1 <- distribution_plot + theme(axis.title.x =element_blank(), axis.text.x = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12), plot.margin = margin(0,0,0,0)) +
 ( (evaluation_plot) & guides(color = "none") & theme(axis.text = element_text(size = 15))) +
  fixed_distribution_plot + theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 15), legend.title = element_text(size = 15), legend.text = element_text(size = 12), plot.margin = margin(0,0,0,0)) +
  estimate_boxplots_plot + theme(axis.text = element_text(size = 15), axis.title = element_text(size = 20), axis.title.x = element_blank()) +
  CV_plot + theme(axis.text = element_text(size = 15),axis.title = element_text(size = 20), legend.position = "none") +
  plot_layout(ncol = 1,  heights = c(1,1.5,1,1,0.75)) & 
  theme(plot.title = element_blank(), plot.tag = element_text(size = 25)) &
  plot_annotation(tag_levels = 'a')

# ggsave(here("paper_figures", "practical_identifiability_supplemental_mockup_study2.pdf"),  plot1, width = 20, height = 17)
ggsave(here("paper_figures", "practical_identifiability_supplemental_mockup_study2_2.pdf"),  plot1, width = 20, height = 18)
# ggsave(here("paper_figures", "practical_identifiability_supplemental_mockup_study2.pdf"),  plot1, width = 13.33, height = 10)

plot2 <- (distribution_plot_v2 & labs(color = "BICc") & theme(axis.text.y = element_text(size = 12)) ) / 
  (evaluation_plot & guides(color = "none") & theme(axis.text = element_text(size = 12), plot.title = element_blank())) /
  (fixed_distribution_plot_v2 & theme(axis.text.y = element_text(size = 12)))  / 
  CV_plot + theme(axis.text = element_text(size = 12),axis.title = element_text(size = 15)) +
  plot_layout(ncol = 1,  heights = c(1.5,1.5,1.5,0.75), axes = "collect") & 
  theme(plot.tag = element_text(size = 20)) 
# plot_annotation(tag_levels = 'a')

ggsave(here("paper_figures", "practical_identifiability_supplemental_mockup_study2_version2.pdf"),  plot2, width = 20, height = 18)

## SCM steps: ####
# Study 1 ####
SCM_summary <- read_csv(here("results_eclipse", "study1", "SCM", "model_building_summary.csv")) 

SCM_summary_figure_study1 <- SCM_summary %>% 
  filter(n_covar != 0) %>% 
  mutate(id = 1:nrow(.)) %>% 
  filter(id < min(id[n_covar > 3])) %>% 
  group_by(n_covar) %>% 
  mutate(shaded = min(BICc) + 2, 
         n_covar_label = factor(as.character(n_covar), levels = 1:5, 
                                labels = paste(c("1st", "2nd", "3rd", "4th", "5th"), "treatment effect added")),
         label = ifelse(BICc < threshold, param, ""),
         param = factor(param, levels = c("lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         label = factor(label, levels = c("","lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("","lambda", "beta", "t['a']", "m", "p", "N['p']"))) %>% 
  mutate(param = fct_rev(param)) %>% 
  ggplot(aes(BICc - threshold, treatment, fill = param)) +
  geom_rect(data = ~filter(., BICc == min(BICc)),
            aes(xmin =BICc-threshold , xmax = shaded-threshold, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.5, inherit.aes = F) +
  geom_col(position = position_dodge(), color = NA) + 
  geom_text(aes(BICc - threshold, label = label, color = param),
            parse = T, hjust = 1.1, position = position_dodge(width = 0.9), show.legend = F) +
  facet_wrap( ~ n_covar_label, ncol = 1) +
  # scale_y_discrete(labels = c("vac" = "Vaccine", "tlr7" = "TLR7", "inTLR7VAC" = "TLR7 * Vaccine")) +
  scale_y_discrete(labels = rev(c(   "TLR7"))) +
  labs(x = "Change in BICc", y = "Treatment effect added to model", fill = "Parameter") +
  scale_fill_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0")), labels = function(x) parse(text = x)) + 
  scale_color_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0"))) +
  guides(fill = guide_legend(reverse = T)) + 
  new_scale_color() +
  # Add line around treatment effects in final model
  geom_col(data = . %>% filter(n_covar > 1) %>% mutate(col = (n_covar == 2 & param == "lambda") | (n_covar == 3 & param == "t['a']")), 
           aes(group = param, color = col), fill = NA, position = position_dodge(),  lty = "dashed", show.legend = F) +
  scale_color_manual(values = c("#FFFFFF00", "black")) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave(here("paper_figures", "SCM_steps_supplemental_figure_study1_mockup.pdf"), SCM_summary_figure_study1, width = 6.5, height = 8)


# Study 2 ####
SCM_summary <- read_csv(here("results_eclipse", "study2", "SCM", "model_building_summary.csv")) 

SCM_summary_figure_study2 <- SCM_summary %>% 
  filter(n_covar != 0) %>% 
  mutate(id = 1:nrow(.)) %>% 
  filter(id < min(id[n_covar > 5])) %>% 
  group_by(n_covar) %>% 
  mutate(shaded = min(BICc) + 2, 
         n_covar_label = factor(as.character(n_covar), levels = 1:5, 
                                labels = paste(c("1st", "2nd", "3rd", "4th", "5th"), "treatment effect added")),
         label = ifelse(BICc < threshold, param, ""),
         param = factor(param, levels = c("lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         label = factor(label, levels = c("","lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("","lambda", "beta", "t['a']", "m", "p", "N['p']"))) %>% 
  mutate(param = fct_rev(param),
         treatment = factor(treatment, levels = rev(c("tlr7", "vac", "intTLR7VAC")),
                            labels = rev(c("TLR7", "Vaccine", "TLR7 * Vaccine")))) %>% 
  ggplot(aes(BICc - threshold, treatment, fill = param)) +
  geom_rect(data = ~filter(., BICc == min(BICc)),
            aes(xmin =BICc-threshold , xmax = shaded-threshold, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.5, inherit.aes = F) +
  geom_col(position = position_dodge(), color = NA) + 
  geom_text(aes(BICc - threshold, label = label, color = param),
            parse = T, hjust = 1.1, position = position_dodge(width = 0.9), show.legend = F) +
  facet_wrap( ~ n_covar_label, ncol = 1) +
  # scale_y_discrete(labels = c("vac" = "Vaccine", "tlr7" = "TLR7", "inTLR7VAC" = "TLR7 * Vaccine")) +
  # scale_y_discrete(labels = rev(c( "Vaccine",  "TLR7", "TLR7 * Vaccine"))) +
  labs(x = "Change in BICc", y = "Treatment effect added to model", fill = "Parameter") +
  scale_fill_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0")), labels = function(x) parse(text = x)) + 
  scale_color_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0"))) +
  guides(fill = guide_legend(reverse = T)) + 
  new_scale_color() +
  geom_col(data = . %>% filter(n_covar <= 3) %>% mutate(col = BICc == min(BICc)), 
           aes(group = param, color = col), fill = NA, position = position_dodge(),  lty = "dashed", show.legend = F) +
  scale_color_manual(values = c("#FFFFFF00", "black")) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave(here("paper_figures", "SCM_steps_supplemental_figure_study2_mockup.pdf"), SCM_summary_figure_study2, width = 6.5, height = 10)


# Study 3 ####
SCM_summary <- read_csv(here("results_eclipse", "study3", "SCM_smaller_lambda", "model_building_summary.csv")) 

SCM_summary_figure_study3 <- SCM_summary %>% 
  filter(n_covar != 0) %>% 
  mutate(id = 1:nrow(.)) %>%  
  filter(id < min(id[n_covar > 3])) %>% 
  group_by(n_covar) %>% 
  mutate(shaded = min(BICc) + 2, 
         n_covar_label = factor(as.character(n_covar), levels = 1:5, 
                                labels = paste(c("1st", "2nd", "3rd", "4th", "5th"), "treatment effect added")),
         label = ifelse(BICc < threshold, param, ""),
         param = factor(param, levels = c("lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("lambda", "beta", "t['a']", "m", "p", "N['p']")),
         label = factor(label, levels = c("","lambda", "betas", "ta", "m", "p", "NPs"),
                        labels = c("","lambda", "beta", "t['a']", "m", "p", "N['p']"))) %>% 
  mutate(param = fct_rev(param),
         treatment = factor(treatment, levels = rev(c("tlr7", "ab", "inttlr7ab")),
                            labels = rev(c("TLR7", "Antibody", "TLR7 * Antibody")))) %>% 
  ggplot(aes(BICc - threshold, treatment, fill = param)) +
  geom_rect(data = ~filter(., BICc == min(BICc)),
            aes(xmin =BICc-threshold , xmax = shaded-threshold, ymin = -Inf, ymax = Inf),
            fill = "lightgrey", alpha = 0.5, inherit.aes = F) +
  geom_col(position = position_dodge(), color = NA) + 
  geom_text(aes(BICc - threshold, label = label, color = param),
            parse = T, hjust = 1.1, position = position_dodge(width = 0.9), show.legend = F) +
  facet_wrap( ~ n_covar_label, ncol = 1) +
  labs(x = "Change in BICc", y = "Treatment effect added to model", fill = "Parameter") +
  scale_fill_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0")), labels = function(x) parse(text = x)) + 
  scale_color_manual(values = rev(c("#E97130", "#F3B798", "#145F82","#1A6B25", "#7CA67C", "#BED3C0"))) +
  guides(fill = guide_legend(reverse = T)) + 
  new_scale_color() +
  geom_col(data = . %>% mutate(col = (n_covar == 1 & param == "beta" & treatment == "Antibody") | (n_covar == 2 & param == "N['p']" & treatment == "Antibody") | (n_covar == 3 & param == "m" & treatment == "TLR7")), 
           aes(group = param, color = col), fill = NA, position = position_dodge(),  lty = "dashed", show.legend = F) +
  scale_color_manual(values = c("#FFFFFF00", "black")) + 
  theme(axis.text = element_text(size = 12), axis.title = element_text(size = 15), 
        strip.text = element_text(size = 15), legend.text = element_text(size = 12),
        legend.title = element_text(size = 15))

ggsave(here("paper_figures", "SCM_steps_supplemental_figure_study3_mockup.pdf"), SCM_summary_figure_study3, width = 6.5, height = 8)
