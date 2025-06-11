################################################################################
#### This script is for exploring the viral dynamics data and generating subsets
################################################################################

##### Load Libraries
library(tidyverse)
library(here)
library(patchwork)
library(RColorBrewer)

# working directory
here()

# set theme
theme_set(theme_bw())

# Read in data ####

SIV_data <- read_table(here("data", "source", "VACTLR7_DATA.txt")) 
SHIV_data <- read_table(here("data", "source", "ABTLR7_DATA.txt")) 
VACAB_data <- read_table(here("data", "source", "VACABTLR7_DATA.txt"))

# In SIV data, study = 1 corresponds to study 1, study = 0 corresponse to study 2
glimpse(SIV_data)
glimpse(SHIV_data)
glimpse(VACAB_data)

# Count number of animals in each study
SIV_data %>% distinct(id, study) %>% count(study)
SHIV_data %>% distinct(id) %>% nrow()
VACAB_data %>% distinct(id) %>% nrow()
# Study 1: 21, Study 2: 34, Study 3: 44, Study 4: 49

# Count animals in each treatment group:
SIV_data %>% distinct(id, study, tlr7, vac) %>% count(study, tlr7, vac)
#   study  tlr7   vac     n
# 1 0         0     0     9
# 2 0         0     1     9
# 3 0         1     0     8
# 4 0         1     1     8
# 5 1         0     0     8
# 6 1         1     0    13

SHIV_data %>% distinct(id, tlr7, ab) %>% count(tlr7, ab)
#    tlr7    ab     n
# 1     0     0    11
# 2     0     1    11
# 3     1     0    11
# 4     1     1    11

VACAB_data %>% distinct(id, tlr7, ab, vac) %>% count(tlr7,vac,ab)
#    tlr7   vac    ab     n
# 1     0     0     0    15
# 2     1     0     1    12
# 3     1     1     0    12
# 4     1     1     1    10


# Compare data collection time points ####

rbind(
  SIV_data %>%  mutate(study = abs(study - 2)) %>% distinct(study, id, time),
  SHIV_data %>%  mutate(study = 3) %>% distinct(study, id, time),
  VACAB_data %>%  mutate(study = 4) %>% distinct(study, id, time)
) %>% 
  arrange(study, id) %>% 
  mutate(y = fct_inorder(interaction(study, id))) %>% 
  ggplot(aes(time, y, color = factor(study))) + geom_point(size = 3, alpha = 0.5) +
  # facet_wrap(~study, ncol = 1, scales ="free_y") + 
  theme_classic() + 
  labs(x = "Days since treatment interruption", y = "Animal ID", color = "Study") 


# Visualize rebound trajectories ####
rbind(
  SIV_data %>% select(id, time, vl, tlr7, second_therapy = vac, study)  %>% 
    mutate(treatment = case_when(
      tlr7 == 0 & second_therapy == 0 ~ "Control",
      tlr7 == 1 & second_therapy == 0 ~ "TLR7",
      tlr7 == 0 & second_therapy == 1 ~ "VAC",
      tlr7 == 1 & second_therapy == 1 ~ "TLR7+VAC"
    ), data = "SIV", study = abs(study-2)),
  SHIV_data %>% select(id, time, vl, tlr7, second_therapy = ab)  %>% 
    mutate(treatment = case_when(
      tlr7 == 0 & second_therapy == 0 ~ "Control",
      tlr7 == 1 & second_therapy == 0 ~ "TLR7",
      tlr7 == 0 & second_therapy == 1 ~ "AB",
      tlr7 == 1 & second_therapy == 1 ~ "TLR7+AB"
    ), data = "SHIV", study = 3
    ),
  VACAB_data %>% filter(vac != 1) %>% select(id, time, vl, tlr7, second_therapy = ab) %>% 
    mutate(treatment = case_when(
      tlr7 == 0 & second_therapy == 0 ~ "Control",
      tlr7 == 1 & second_therapy == 0 ~ "TLR7",
      tlr7 == 0 & second_therapy == 1 ~ "AB",
      tlr7 == 1 & second_therapy == 1 ~ "TLR7+AB"
    ), data = "SHIV", study = 4),
  
  VACAB_data %>% filter(ab != 1) %>% select(id, time, vl, tlr7, second_therapy = vac) %>% 
    mutate(treatment = case_when(
      tlr7 == 0 & second_therapy == 0 ~ "Control",
      tlr7 == 1 & second_therapy == 0 ~ "TLR7",
      tlr7 == 0 & second_therapy == 1 ~ "VAC",
      tlr7 == 1 & second_therapy == 1 ~ "TLR7+VAC"
    ), data = "SHIV", study = 4)
  
  ) %>% mutate(combination = ifelse(second_therapy == 0, "Single Therapy", "Dual Therapy")) %>% 
  select(-second_therapy) %>% 
  rbind(VACAB_data %>% filter(ab == 1 & vac == 1) %>% select(id, time, vl, tlr7) %>% 
          mutate(treatment = "TLR7+VAC+AB", data = "SHIV", study = 4, combination = "Triple Therapy")) %>% 
  mutate(treatment = factor(treatment, levels = c("Control", "TLR7", "VAC", "AB", "TLR7+VAC", "TLR7+AB", 'TLR7+VAC+AB')),
         group = paste0(treatment, " (Study", study, ")"),
         combination = factor(combination, levels = c("Single Therapy", "Dual Therapy", "Triple Therapy"))) %>%
  arrange(treatment) %>% 
  mutate(group = fct_inorder(group)) %>% 
  ggplot(aes(time, vl, group = interaction(id, treatment, data), color = group)) + 
  geom_line( alpha = 0.7) + 
  facet_grid(combination~tlr7, label = "label_both") + 
  labs(x = "time (days)", y = "Log10 Viral Load", color = "Group") +
  scale_color_manual(
    values = c(brewer.pal(5, "Purples")[-1], brewer.pal(4, "Blues")[-1],
               brewer.pal(4, "Greens")[2], brewer.pal(4, "Oranges")[2],
               brewer.pal(4, "Greens")[3], brewer.pal(4, "Greens")[4], 
               brewer.pal(4, "Oranges")[3], brewer.pal(4, "Oranges")[4],
               "goldenrod4")
  )

SIV_data %>% 
  mutate(study = factor(ifelse(study == 1, "1", "2"), levels = as.character(1:3))) %>% 
  ggplot(aes(time, vl, group = id, color = study)) + 
  geom_line(show.legend = T, alpha = 0.7) + 
  facet_grid(vac~tlr7, label = "label_both") + 
  labs(x = "time (days)", y = "Log10 Viral Load",
       title = "SIV data", color = "Study") + 
  
  SHIV_data %>% 
  mutate(study = factor("3", levels = as.character(1:3))) %>% 
  ggplot(aes(time, vl, group = id, color = study)) + 
  geom_line(show.legend = T, alpha = 0.7) + 
  facet_grid(ab~tlr7, label = "label_both") + 
  labs(x = "time (days)", y = "Log10 Viral Load",
       title = "SHIV data", color = "Study") + 
  
  plot_layout(ncol = 1, guides = "collect") & 
  scale_color_brewer(palette = "Dark2", drop = F)



#####
# Create a data set with only STUDY1 data
study1 <- SIV_data %>% 
  filter(study == 1) %>% 
  select(id, time, vl, censor, tlr7) 

write_delim(study1, here("data", "derived", "STUDY1_DATA.txt"), delim = "\t")

#####
# Create a data set with only STUDY1 control data
study1_control <- SIV_data %>% 
  filter(study == 1, vac == 0, tlr7 == 0) %>% 
  select(id, time, vl, censor) 

write_delim(study1_control, here("data", "derived", "STUDY1_CONTROL_DATA.txt"), delim = "\t")

#####
# Create a data set with only STUDY1 data excluding non-rebounding animals
exclude <- SIV_data %>% filter(vl < 2, time > 50) %>% pull(id) %>% unique

SIV_data_rebound_s1 <- SIV_data %>% 
  filter(!(id %in% exclude), study == 1)

write_delim(SIV_data_rebound_s1, here("data", "derived", "STUDY1_REBOUND_DATA.txt"), delim = "\t")

#####
# Create a data set with only STUDY2 data
study2 <- SIV_data %>% 
  filter(study == 0) %>% 
  select(id, time, vl, censor, tlr7, vac, intTLR7VAC) 

write_delim(study2, here("data", "derived", "STUDY2_DATA.txt"), delim = "\t")

#####
# Create a data set with only STUDY2 control data
study2_control <- SIV_data %>% 
  filter(study == 0, tlr7 == 0 , vac == 0) %>% 
  select(id, time, vl, censor) 

write_delim(study2_control, here("data", "derived", "STUDY2_CONTROL_DATA.txt"), delim = "\t")

#####
# Create data set with only controls for SIV
SIV_data_control <- SIV_data %>% 
  filter(tlr7 == 0, vac == 0) %>% 
  select(id, time, vl, censor, study)

write_delim(SIV_data_control, here("data", "derived", "SIV_CONTROL_DATA.txt"), delim = "\t")


#####
# Create a data set excluding non-rebounding animals
exclude <- SIV_data %>% filter(vl < 2, time > 50) %>% pull(id) %>% unique

SIV_data_rebound <- SIV_data %>% 
  filter(!(id %in% exclude))

write_delim(SIV_data_rebound, here("data", "derived", "SIV_REBOUND_DATA.txt"), delim = "\t")


#####
# Create a dataset with only control + TLR7 from either SIV study
TLR7 <- SIV_data %>% 
  filter(vac != 1) %>% 
  select(id, time, vl, censor, tlr7, study, intTLR7STUD) 

write_delim(TLR7, here("data", "derived", "SIV_TLR7_DATA.txt"), delim = "\t")

#####
# Create a dataset with only control + VAC
VAC <- SIV_data %>% 
  filter(tlr7 != 1) %>% 
  select(id, time, vl, censor, vac) 

write_delim(VAC, here("data", "derived", "VAC_DATA.txt"), delim = "\t")

#####
# create table with dosing information

study1_dosing <- rbind(
  data.frame(id = 1:4, GS9620 = 0, num = 7, cum = 1.8, avg = 0.26, max = 0.3, ART_dur_wk = 76),
  data.frame(id = 5:10, GS9620 = 0, num = 0, cum = 0, avg = 0, max = 0, ART_dur_wk = 76),
  data.frame(id = 11:12, GS9620 = 0, num = 0, cum = 0, avg = 0, max = 0, ART_dur_wk = 115),
  data.frame(id = 13:15, GS9620 = 0, num = 19, cum = 1.9, avg = 0.1, max = 0.1, ART_dur_wk = 115),
  data.frame(id = 16:18, GS9620 = 1, num = 19, cum = 0.95, avg = 0.05, max = 0.05, ART_dur_wk = 115),
  data.frame(id = 19:21, GS9620 = 1, num = 10, cum = 1.5, avg = 0.15, max = 0.15, ART_dur_wk = 115)
)

study2_dosing <- rbind(
  data.frame(id = unique(filter(study2, tlr7 == 1)$id), GS9620 = 1, num = 11, cum = 1.65, avg = 0.15, max = 0.15, ART_dur_wk = 71),
  data.frame(id = unique(filter(study2, tlr7 != 1)$id), GS9620 = 0, num = 0, cum = 0, avg = 0, max = 0, ART_dur_wk = 71)
)

dosing <- rbind(study1_dosing , study2_dosing)

# Combine vl data with dosing

SIV_data_dosing <- left_join(SIV_data, dosing, by = "id") 

write_delim(SIV_data_dosing, here("data", "derived", "SIV_DATA_DOSING.txt"), delim = "\t")

#####
# Create a subset of only control animals for study 3 (SHIV)

study3_control <- SHIV_data %>% 
  filter(tlr7 == 0, ab == 0) %>% 
  select(id, time, vl, censored) 

write_delim(study3_control, here("data", "derived", "STUDY3_CONTROL_DATA.txt"), delim = "\t")

#####
# Create a subset of only control animals for study 4 (SHIV2)

study4_control <- VACAB_data %>% 
  filter(tlr7 == 0, ab == 0, tlr7 == 0) %>% 
  select(id, time, vl, censor) 

write_delim(study4_control, here("data", "derived", "STUDY4_CONTROL_DATA.txt"), delim = "\t")
