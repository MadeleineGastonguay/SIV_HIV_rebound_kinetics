# Compare profile likelihoods for lambda
library(tidyverse)
library(here)
library(ggthemes)

theme_set(theme_minimal() + theme(plot.title = element_text(hjust = 0.5)))

studies <- paste0("study", 1:4)
lambda_profile <- lambda_table <- NULL
p_profile <- p_table <- NULL
m_profile <- m_table <- NULL
for(study in studies){
  lambda_profile <- rbind(
    lambda_profile,
    read_csv(here("results_eclipse", study, "lambda_profile_likelihood.csv")) %>%
      mutate(study = study, control = F),
    read_csv(here("results_eclipse", paste0(study, "_control"), "lambda_profile_likelihood.csv")) %>%
      mutate(study = study, control = T)
  )

  lambda_table <- rbind(
    lambda_table,
    read_csv(here("results_eclipse", study, "lambda_profile_likelihood_table.csv")) %>%
      mutate(study = study, control = F),
    read_csv(here("results_eclipse", paste0(study, "_control"), "lambda_profile_likelihood_table.csv")) %>%
      mutate(study = study, control = T)
  )

  p_profile <- rbind(
    p_profile,
    read_csv(here("results_eclipse", study, "p_profile_likelihood.csv")) %>%
      mutate(study = study, control = F),
    read_csv(here("results_eclipse", paste0(study, "_control"), "p_profile_likelihood.csv")) %>%
      mutate(study = study, control = T)
  )

  p_table <- rbind(
    p_table,
    read_csv(here("results_eclipse", study, "p_profile_likelihood_table.csv")) %>%
      mutate(study = study, control = F),
    read_csv(here("results_eclipse", paste0(study, "_control"), "p_profile_likelihood_table.csv")) %>%
      mutate(study = study, control = T)
  )

  m_profile <- rbind(
    m_profile,
    read_csv(here("results_eclipse", study, "m_profile_likelihood.csv")) %>% 
      mutate(study = study, control = F),
    read_csv(here("results_eclipse", paste0(study, "_control"), "m_profile_likelihood.csv")) %>% 
      mutate(study = study, control = T)
  )
  
  m_table <- rbind(
    m_table,
    read_csv(here("results_eclipse", study, "m_profile_likelihood_table.csv")) %>%
      mutate(study = study, control = F),
    read_csv(here("results_eclipse", paste0(study, "_control"), "m_profile_likelihood_table.csv")) %>%
      mutate(study = study, control = T)
  )
}


## Lambda results ####

# Fixing lambda to 100 makes sense for studies 1,2 & 4
# The likelihood for study 3 is flatter, but 100 should be okay
# interestingly, using control animals only indicates a bigger value of lamdba is ideal
# I think this is just becaause there is not enough data to estimate it well

# Plot results from all studies 
lambda_profile %>% 
  group_by(study, control) %>% 
  mutate(bicc_norm = bicc - min(bicc)) %>% 
  ggplot(aes(lambda_pop, bicc_norm, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 100, lty = "dashed") +
  geom_vline(xintercept = 1000, lty = "dashed") + 
  scale_x_log10()

# For study 1, there is a subset of lower BICc driven by high estimate of m that we will remove:
lambda_profile %>% filter(study == "study1", !control, parameter == "m_pop")  %>% 
  filter(estimate >= 250)  %>% pull(lambda_pop) -> big_m

# Plot full studies only
(p <- lambda_profile %>%
  filter(!control, 
         !(study == "study1" & !control & lambda_pop %in% big_m)
         ) %>% 
  ggplot(aes(lambda_pop, bicc)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 100, lty = "dashed") +
  labs(y = "BICc", title = "Profile Likelihood for Lambda")
)

# Plot full studies only (log scale)
p + scale_x_log10()

# Plot control studies only
lambda_profile %>%
  filter(!(study == "study1" & !control & lambda_pop %in% big_m), control) %>% 
  ggplot(aes(lambda_pop, bicc, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 100, lty = "dashed") +
  scale_x_log10()


# Look at consistency of lambda*beta
lambda_profile %>%
  filter(!control, !(study == "study1" & !control & lambda_pop %in% big_m), 
         parameter == "betas_pop") %>% 
  ggplot(aes(lambda_pop*estimate, bicc)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free") +
  geom_smooth() + 
  scale_x_log10()


lambda_profile %>%
  group_by(study, control) %>% 
  filter(!(study == "study1" & !control & lambda_pop %in% big_m)) %>% 
  filter(bicc == min(bicc)) %>% 
  summarise(unique(lambda_pop))


## p results ####

# Fixing p to 1 makes sense for studies 2-4. 
# The likelihood for study 1 is more flat, so 1 should be fine too

# Plot results from all studies 
p_profile %>% 
  group_by(study, control) %>% 
  mutate(bicc_norm = bicc - min(bicc)) %>% 
  ggplot(aes(p_pop, bicc_norm, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 1) +
  labs(y = "change in BICc compared to best model", color = "Control animals only") + 
  scale_color_calc() + 
  scale_x_log10()

# Plot full studies only
p_profile %>%
  filter(!control) %>% 
  ggplot(aes(p_pop, bicc)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 1, lty = "dashed") + 
  labs(y = "BICc", title = "Profile Likelihood for p") 

# Plot control studies only
p_profile %>%
  filter( control) %>% 
  ggplot(aes(p_pop, bicc, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 1, lty = "dashed") 

# Absolute minimum BICc
p_profile %>%
  group_by(study, control) %>% 
  filter(bicc == min(bicc)) %>% 
  summarise(unique(p_pop))

## m results ####

# These are all over the place, there is no clear minimum

# Plot results from all studies 
m_profile %>% 
  group_by(study, control) %>% 
  mutate(bicc_norm = bicc - min(bicc)) %>% 
  ggplot(aes(m_pop, bicc_norm, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  labs(y = "change in BICc compared to best model", color = "Control animals only") + 
  scale_color_calc() + 
  scale_x_log10()

# Plot full studies only
m_profile %>%
  filter(!control) %>% 
  ggplot(aes(m_pop, bicc)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  labs(y = "BICc") 

# Plot control studies only
m_profile %>%
  filter( control) %>% 
  ggplot(aes(m_pop, bicc, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() 

# Absolute minimum BICc
m_profile %>%
  group_by(study, control) %>% 
  filter(bicc == min(bicc)) %>% 
  summarise(unique(m_pop))

