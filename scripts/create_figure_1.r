library(tidyverse)
library(here)
library(patchwork)
library(ggbreak)

## Read in data for each study: #####

# ART start and stop times for each study:
t_ART_start <- c(9*7, 7, 7, 9)
t_ART_stop <- c((76+9)*7, 72*7, 130*7, 85*7+9)

# Study 1
S1_ID <- "156-08"
STUDY1_DATA <- read_delim(here("data", "derived", "STUDY1_DATA.txt"), delim = "\t")
STUDY1_full_traj <- read_csv(here("data", "source", "study1_TLR7_before_rebound.csv")) %>% 
  filter(ID == S1_ID)

# Study 2
S2_ID <- "5886"
STUDY2_DATA <- read_delim(here("data", "derived", "STUDY2_DATA.txt"), delim = "\t")
STUDY2_b4_ATI <- read_csv(here("data", "source", "study2_TLR7VAC_average_before_rebound.csv")) %>% 
  mutate(ID = "Average")
STUDY2_full_traj <- rbind(STUDY2_b4_ATI, 
                          STUDY2_DATA %>% filter(id == S2_ID) %>% 
                            select(ID=id, time, vl) %>% mutate(time = time + t_ART_stop[2]))

# Study 3
S3_ID <- "6399"
STUDY3_DATA <- read_delim(here("data", "source", "ABTLR7_DATA.txt"), delim = "\t")
STUDY3_full_traj <- read_csv(here("data", "source", "study3_full_trajectories.csv")) %>% 
  filter(ID == S3_ID)

# Study 4
S4_ID <- 7148
STUDY4_DATA <- read_delim(here("data", "source", "VACABTLR7_DATA.txt"), delim = "\t") 
STUDY4_b4_ATI <- read_csv(here("data", "source", "study4_VACABTLR7_before_rebound.csv")) %>% mutate(time = time*7)
STUDY4_full_traj <-  rbind(
  STUDY4_b4_ATI,
  STUDY4_DATA %>% select(time, vl, ID= id) %>%
    mutate(time = time + t_ART_stop[4])) %>% 
  filter(ID == S4_ID)



## Set up plotting settings and functions ####
theme_set(theme_bw(base_size = 15) + theme(panel.grid = element_blank(), axis.ticks.length = unit(-2.75, "pt"), plot.title = element_text(hjust = 0.5), plot.background = element_blank()) )

add_axis_ticks <- function(plot, axis = "both", end = 1150, seq_by = 50){
  plot <- plot + scale_y_continuous(sec.axis = dup_axis(name = NULL, labels = NULL), breaks = 0:8, limits = c(0,8), expand = c(0,0)) +
    scale_x_continuous(sec.axis = dup_axis(name = NULL, labels = NULL), breaks = seq(0,end,by = seq_by), expand = c(0,0), limits = c(0,end)) 
  
  if(axis == "y"){
    plot <- plot + theme(axis.text.x = element_blank(), axis.title.x = element_blank())
  }else if(axis == "x"){
    plot <- plot + theme(axis.text.y = element_blank(), axis.title.y =  element_blank())
  }else if(axis == "none"){
    plot <- plot + theme(axis.text.y = element_blank(), axis.title.y =  element_blank(), axis.text.x = element_blank(), axis.title.x = element_blank())
  }
  
  return(plot + coord_cartesian(clip = "off"))
} 

# names are tlr7:vac:ab
# colors <- c("0.0.0" = "black", "1.0.0" = "#027dbb", "0.1.0" = "#fcba12", "0.0.1" ="#fc12a6", "1.1.0" = "#7ab148",  "1.0.1" = "#7f48b1", "1.1.1" =  "#a96e7c")
colors <- c("0.0.0" = "black", "1.0.0" = "#027dbb", "0.1.0" = "#fcba12", "0.0.1" ="#fa61c2", "1.1.0" = "#7ab148",  "1.0.1" = "#7f48b1", "1.1.1" =  "#a96e7c")


## Plot example timecourses ####

plot_example_timecourse <- function(traj_df, ART_start, ART_stop, end_time, break_seq = 50, x_break = NULL){
  cens <- min(traj_df$vl, na.rm = TRUE)
  # browser()
  if(!is.null(x_break)){
    traj_df <- traj_df %>% filter(time > x_break[2] | time < x_break[1])
    traj_df <- traj_df %>% mutate(time2 = ifelse(time > x_break[2], time - x_break[2] + x_break[1], time))
    labels <- traj_df$time2
    names(labels) <- traj_df$time
    traj_df$time <- traj_df$time2
    
    ART_stop <- ART_stop - x_break[2] + x_break[1]
  }
  
  plot <- traj_df %>% 
    ggplot(aes(time, vl)) + 
    geom_hline(yintercept = cens, linetype = "dashed", color = "gray") +
    geom_line() +
    geom_point() + 
    geom_vline(xintercept = ART_start, linetype = "dashed", color = "#d32a42") +
    geom_vline(xintercept = ART_stop, linetype = "dashed", color = "#d32a42") +
    labs(x = "Time since infection (days)", y = "Viral load\n(log10 copies/ml)") +
    theme(plot.margin = margin(0,10,10,10))
  
  plot <- add_axis_ticks(plot, end = end_time, seq_by = break_seq) 
  
  if(!is.null(x_break)){
    plot <- plot + 
      annotate(geom = "text", x =  x_break[1] + 25 , y = -Inf, label = "//", size = 5) + 
      scale_x_continuous(labels = function(lab){
      idx <- which(lab > x_break[1])
      lab[idx] <- lab[idx] + diff(x_break)
      return(lab)
    }, breaks = seq(0, 1150, by = 50),
    limits = c(0, end_time - x_break[2] + x_break[1]), expand = c(0,0)) 
  }
  return(plot)
}

S1_example <- plot_example_timecourse(STUDY1_full_traj, t_ART_start[1], t_ART_stop[1], 700, x_break = c(200, 400))
S2_example <- plot_example_timecourse(STUDY2_full_traj, t_ART_start[2], t_ART_stop[2], 700)
S3_example <- plot_example_timecourse(STUDY3_full_traj, t_ART_start[3], t_ART_stop[3], 1150, x_break = c(150,500))
S4_example <- plot_example_timecourse(STUDY4_full_traj, t_ART_start[4], t_ART_stop[4], 800)

## Make timing of interventions:

STUDY1_interventions <- tibble(
  intervention = factor(c("ART", "TLR7"), levels = c("ART", "TLR7", "d6", "d5", "d1", "d2", "d3", "d4")),
  label = c("ART", "TLR7-agonist"),
  start = c(9*7, 70*7),
  stop = c((76+9)*7, 82*7)
)

STUDY2_interventions <- tibble(
  intervention = factor(c("ART", "TLR7", "VAC","VAC"), levels = c("ART", "VAC", "TLR7", "d5", "d1", "d2", "d3", "d4")),
  label = c("ART", "TLR7-agonist", "Vaccine prime\n(Ad26)", "Vaccine boost\n(MVA)"),
  start = c(7, 50*7, 24*7, 48*7),
  stop = c(72*7, 70*7, 36*7, 60*7)
)

STUDY3_interventions <- tibble(
  intervention = factor(c("ART", "TLR7", "AB"), levels = c("ART", "AB", "TLR7", "d1", "d2", "d3", "d4", "d5")),
  label = c("ART", "TLR7-agonist", "Antibody (PGT121)"),
  start = c(7, 96*7, 106*7),
  stop = c(130*7, 114*7,114*7)
)

STUDY4_interventions <- tibble(
  intervention = factor(c("ART", "TLR7", "VAC","VAC", "AB"), levels = c("ART", "AB", "VAC", "TLR7", "d1", "d2", "d3", "d4")),
  label = c("ART", "TLR7-agonist", "Vaccine prime\n(Ad26)", "Vaccine boost\n(MVA)", "Antibody (PGT121)"),
  start = c(9, 50*7, 24*7, 48*7, 64*7),
  stop = c(85*7+9, 72*7, 36*7, 60*7, 72*7)
)

plot_timing <- function(timing_df, xend, x_break = NULL){
  colors <- c("#d32a42", unname(colors[c("0.0.1", "0.1.0", "1.0.0")]))
  names(colors) <- c("ART", "AB", "VAC", "TLR7")
  
  if(!is.null(x_break)){
    timing_df <- timing_df %>% mutate(start = ifelse(intervention != "ART", start - x_break[2] + x_break[1], start),
                                      stop = stop - x_break[2] + x_break[1])
    xend <- xend - x_break[2] + x_break[1]
  }
  
  timing_df %>% 
    ggplot(aes(xmin = start, xmax = stop, y = intervention, color = intervention)) + 
    geom_point(aes(start, intervention), alpha = 0) +
    geom_linerange() +
    theme_void() +
    scale_x_continuous(breaks = seq(0,xend,by = 50), expand = c(0,0), limits = c(0,xend)) +
    scale_y_discrete(drop = F) +
    coord_cartesian( clip = "off")  + guides(color = "none") + 
    scale_color_manual(values = colors)
}

S1_timing_plot <- plot_timing(STUDY1_interventions, 700, x_break = c(200, 400))
S2_timing_plot <- plot_timing(STUDY2_interventions, 700)
S3_timing_plot <- plot_timing(STUDY3_interventions, 1150, x_break = c(150,500))
S4_timing_plot <- plot_timing(STUDY4_interventions, 800)

# 
# S4_timing_plot <- STUDY4_interventions %>% 
#   ggplot(aes(xmin = start, xmax = stop, y = intervention, color = intervention)) + 
#   geom_point(aes(start, intervention), alpha = 0) +
#   geom_linerange(data = ~filter(.x, intervention == "TLR7"), linewidth = 1, lty = "dashed") +
#   geom_linerange(data = ~filter(.x, intervention %in% c("ART", "AB")), linewidth = 2) +
#   geom_text(data = ~filter(.x, intervention %in% c("TLR7", "AB")), aes(label = label, x = start), hjust = 1.05, fontface = 2) +
#   geom_text(data = ~filter(.x, intervention == "VAC"), aes(label = label, x = start + (stop-start)/2), hjust = 0.5, vjust = 0.5, fontface = 2) +
#   geom_segment(data = ~filter(.x, intervention == "VAC"), aes(x = start, y = 3, xend = start, yend = 2.7),
#                arrow = arrow(length = unit(8, "points"), type = "closed") )  + 
#   geom_segment(data = ~filter(.x, intervention == "VAC"), aes(x = stop, y = 3, xend = stop, yend = 2.7),
#                arrow = arrow(length = unit(8, "points"), type = "closed")) + 
#   theme_void() +
#   scale_x_continuous(breaks = seq(0,800,by = 50), expand = c(0,0), limits = c(0,800)) +
#   scale_y_discrete(drop = F) +
#   coord_cartesian( clip = "off")  + guides(color = "none") + 
#   scale_color_manual(values = c("#d32a42", unname(colors[c("0.0.1", "0.1.0", "1.0.0")])))


## Make trajectories after ATI:

plot_rebound_data <- function(data, TLR7, VAC, AB, axis, xend, seq_by){
  if(!"ab" %in% names(data)) data <- mutate(data, ab = 0)
  if(!"vac" %in% names(data)) data <- mutate(data, vac = 0)
  if(!"censor" %in% names(data)) data <- rename(data, censor = censored)
  
  cens <- max(data %>% filter(censor == 1) %>% pull(vl))
  trt <- as.character(interaction(TLR7, VAC, AB))
  
  name <- switch(trt,
                 "0.0.0" = "Control",
                 "1.0.0" = "TLR7",
                 "0.1.0" = "VAC",
                 "0.0.1" = "Ab",
                 "1.1.0" = "TLR7 + VAC",
                 "1.0.1" = "TLR7 + Ab",
                 "1.1.1" = "TLR7 + VAC + Ab"
  )
  
  data <- data %>% 
    filter(ab == AB, vac == VAC, tlr7 == TLR7) 
  
  N <- length(unique(data$id))
  
  plot <- data %>% 
    ggplot(aes(time, vl, group = id, color = interaction(tlr7, vac, ab))) +
    geom_hline(yintercept = cens, color = "gray", lty = "dashed") + 
    geom_line(data = ~filter(.x, !id %in% no_rebound_id), show.legend = FALSE) + 
    scale_color_manual(values = colors) +
    labs(x = "Time after ART stop (days)", y = "Viral load\n(log10 copies/ml)",
         title = name)  +
    theme(plot.title = element_text(color = colors[trt], hjust = 0.5, face = 2, vjust = -7), plot.margin = margin(0, 10, 10, 10))
  
  no_rebound_id <- data %>% group_by(id) %>% summarise(n_cens = sum(censor), total_obs = n()) %>% ungroup %>% filter(n_cens == total_obs) %>% pull(id) %>% unique
  if(length(no_rebound_id) > 0){
    no_rebound_data <- tibble(id = no_rebound_id, vl = cens, idx = 0:(length(no_rebound_id)-1), vl_jitter = vl - idx*0.1) %>% left_join(select(data, id, tlr7, vac, ab))  
    no_rebound_data <- rbind(no_rebound_data %>% mutate(time = 0), no_rebound_data %>% mutate(time = max(data$time)))
    # plot <- plot + geom_hline(data = no_rebound_data, aes(yintercept = vl_jitter, color = interaction(tlr7, vac, ab)), show.legend = FALSE) 
    plot <- plot + geom_line(data = no_rebound_data, aes(time, vl_jitter, color = interaction(tlr7, vac, ab)), show.legend = FALSE) 
  }
  
  plot <- plot + annotate(geom = "text", x= Inf, y= -Inf, label = str_interp("N=${N}"), hjust = 1.5, vjust = -0.5, size = 5) 
  
  return(add_axis_ticks(plot, axis, xend, seq_by))
}


S1_control <- plot_rebound_data(STUDY1_DATA, 0,0,0, "both", 120, 20) 
S1_TLR7 <- plot_rebound_data(STUDY1_DATA, 1,0,0, "x", 120, 20)

S2_control <- plot_rebound_data(STUDY2_DATA, 0,0,0, "y", 180, 30) 
S2_TLR7 <- plot_rebound_data(STUDY2_DATA, 1,0,0, "none", 180, 30)  
S2_VAC <-  plot_rebound_data(STUDY2_DATA, 0,1,0, "both", 180, 30)  
S2_TLR7VAC <- plot_rebound_data(STUDY2_DATA, 1, 1, 0, "x", 180, 30) 

S3_control <- plot_rebound_data(STUDY3_DATA, 0,0,0, "y", 200, 50) 
S3_TLR7 <- plot_rebound_data(STUDY3_DATA, 1,0,0, "none", 200, 50)  
S3_AB <- plot_rebound_data(STUDY3_DATA, 0,0,1, "both", 200, 50)  
S3_TLR7AB <- plot_rebound_data(STUDY3_DATA, 1, 0, 1, "x", 200, 50) 

S4_control <- plot_rebound_data(STUDY4_DATA, 0,0,0, "y", 180, 30) 
S4_TLR7AB <- plot_rebound_data(STUDY4_DATA, 1,0,1, "none", 180, 30)  
S4_TLR7VAC <- plot_rebound_data(STUDY4_DATA, 1,1,0, "both", 180, 30)  
S4_TLR7VACAB <- plot_rebound_data(STUDY4_DATA, 1, 1, 1, "x", 180, 30) 

figure1 <- S1_timing_plot + theme(plot.margin = margin(0,0,0,0)) + plot_spacer() + plot_spacer() + 
  S1_example + S1_control + S1_TLR7 + 
  S2_timing_plot + S2_control + theme(plot.margin = margin(0, 10, 0, 10)) + S2_TLR7 + theme(plot.margin = margin(0, 10, 0, 10)) +
  S2_example + S2_VAC + S2_TLR7VAC +   
  S3_timing_plot + S3_control + theme(plot.margin = margin(0, 10, 0, 10)) + S3_TLR7 + theme(plot.margin = margin(0, 10, 0, 10)) + 
  S3_example + S3_AB + S3_TLR7AB +   
  S4_timing_plot + S4_control + theme(plot.margin = margin(0, 10, 0, 10)) + S4_TLR7AB + theme(plot.margin = margin(0, 10, 0, 10)) + 
  S4_example + S4_TLR7VAC + S4_TLR7VACAB +   
  
  plot_layout(ncol = 3, widths = c(2.5, 1, 1), heights = c(0.5,1,1,1,1,1,1,1))

ggsave(here("paper_figures", "figure1_draft.pdf"), figure1 & theme(panel.background = element_blank()),  width = 14, height = 21, bg = "transparent")


