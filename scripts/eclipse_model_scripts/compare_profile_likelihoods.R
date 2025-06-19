# Compare profile likelihoods for lambda
library(tidyverse)
library(here)

theme_set(theme_minimal())

studies <- paste0("study", 1:4)
lambda_profile <- lambda_table <- NULL
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
}

lambda_profile %>% 
  group_by(study, control) %>% 
  mutate(bicc_norm = bicc - min(bicc)) %>% 
  ggplot(aes(lambda_pop, bicc_norm, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 100) +
  geom_vline(xintercept = 1000) + 
  scale_x_log10()

lambda_profile %>% filter(study == "study1", !control, parameter == "m_pop")  %>% 
  filter(estimate >= 250)  %>% pull(lambda_pop) -> big_m

lambda_profile %>%
  filter(!control, !(study == "study1" & !control & lambda_pop %in% big_m)) %>% 
  ggplot(aes(lambda_pop, bicc)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 100, lty = "dashed") 

lambda_profile %>%
  group_by(study, control) %>% 
  mutate(bicc_norm = bicc - min(bicc)) %>% 
  filter(!(study == "study1" & !control & lambda_pop %in% big_m), control) %>% 
  ggplot(aes(lambda_pop, bicc_norm, color = control)) + 
  geom_point() + 
  facet_wrap(~study, scales = "free_y") +
  geom_smooth() + 
  geom_vline(xintercept = 100, lty = "dashed") +
  scale_x_log10()


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
  mutate(bicc_norm = bicc - min(bicc)) %>% 
  filter(!(study == "study1" & !control & lambda_pop %in% big_m), control) %>% 
  group_by(study) %>% filter(bicc == min(bicc)) %>% summarise(unique(lambda_pop))

