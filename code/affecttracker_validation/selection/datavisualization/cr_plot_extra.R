# CR dynamics by rating methods and videos/quadrants
# 1. Plotting means across participants
# 2. Individual plots
# 3. Compute CR descriptive statistics
# 4. Arousal vs. Valence
#
# Author: Antonin Fourcade
# Last version: 23.02.2024

# import packages
# comment: not sure what packages are really used here
library(tidyverse)
library(ggpmisc)
library(ggdist)
library(gghalves)
library(hrbrthemes)
library(scales)
library(viridis)
library(car)
library(readxl)
library(e1071)
library(irr)
library(gridExtra)

# set paths
#data_path <- "/Volumes/Extreme SSD/AffectiveVR/Phase_1/Data/"
#save_path <- "/Volumes/Extreme SSD/AffectiveVR/Phase_1/MBB2024_poster/"

data_path <- "E:/AffectiveVR/Phase_1/Data/AVR/"
save_path <- "E:/AffectiveVR/Phase_1/affectivevr/cr_plots/"

if (!file.exists(save_path)) {dir.create(save_path)}

# read data
data_cr <- read.csv(paste(data_path, "cr_rs_clean.csv", sep = ""))
data_cri_sr <- read.csv(paste(data_path, "cri_sr_clean.csv", sep = ""))

# set names of test sites
test_list <- c("Torino", "Berlin")
# set names of rating methods and baseline
rm_list <- c("Grid", "Flubber", "Proprioceptive")
base_name <- c("Baseline")
# set names of quadrants/videos
#quad_list <- c("HP", "LP", "LN", "HN")
quad_list <- c("HN", "HP", "LN", "LP")
# set dimensions
dimensions <- c("v", "a", "dist", "angle")
dimensions_names <- str_replace_all(dimensions, c("angle" = "Angle", "a" = "Arousal", "v" = "Valence", "dist" = "Distance"))

# format data
data_cr <- data_cr %>% 
  mutate(
    test_site = factor(test_site, levels = test_list),
    rating_method = factor(rating_method, levels = c(rm_list, base_name)),
    quadrant = factor(quadrant, levels = quad_list),
    sj_id = factor(sj_id, levels = unique(data_cr$sj_id)),
    cr_time = factor(cr_time, levels = unique(data_cr$cr_time)))

data_cri_sr <- data_cri_sr  %>% 
  mutate(
    test_site = factor(test_site, levels = test_list),
    rating_method = factor(rating_method, levels = c(rm_list, base_name)),
    quadrant = factor(quadrant, levels = quad_list),
    sj_id = factor(sj_id, levels = unique(data_cr$sj_id)))

# Count non-NaN datapoints per metric
colSums(!is.na(data_cr))

# Create dataframe with mean and std for each quadrant and rating method
dataSummary_cr <- data_cr  %>% 
  group_by(cr_time, quadrant, rating_method) %>% 
  dplyr::summarise(n = n(),
                   cr_v_mean = mean(cr_v, na.rm = T),
                   cr_v_sd = sd(cr_v, na.rm = T),
                   cr_v_se = cr_v_sd/sqrt(n),
                   cr_a_mean = mean(cr_a, na.rm = T),
                   cr_a_sd = sd(cr_a, na.rm = T),
                   cr_a_se = cr_a_sd/sqrt(n),
                   cr_dist_mean = mean(cr_dist, na.rm = T),
                   cr_dist_sd = sd(cr_dist, na.rm = T),
                   cr_dist_se = cr_dist_sd/sqrt(n),
                   cr_angle_mean = mean(cr_angle, na.rm = T),
                   cr_angle_sd = sd(cr_angle, na.rm = T),
                   cr_angle_se = cr_angle_sd/sqrt(n))
# rename quadrants
dataSummary_cr$quadrant <- str_replace(dataSummary_cr$quadrant, 'HP', 'High Arousal - Positive Valence (HP)')
dataSummary_cr$quadrant <- str_replace(dataSummary_cr$quadrant, 'HN', 'High Arousal - Negative Valence (HN)')
dataSummary_cr$quadrant <- str_replace(dataSummary_cr$quadrant, 'LP', 'Low Arousal - Positive Valence (LP)')
dataSummary_cr$quadrant <- str_replace(dataSummary_cr$quadrant, 'LN', 'Low Arousal - Negative Valence (LN)')

# Mean across participants  ------------------------------------------------
# Plot the mean CR time-series across SJs, for each quadrant and each rating_method
# -> one figure for each dimensions (v, a , dist, angle)
  
# change names of mean measurements
mean_names <- c('mean Valence', 'mean Arousal', 'mean Distance', 'mean Angle')
names(dataSummary_cr) <- str_replace_all(names(dataSummary_cr), c("cr_v_mean" = mean_names[1], "cr_a_mean" = mean_names[2], "cr_dist_mean" = mean_names[3], "cr_angle_mean" = mean_names[4]))

# run a loop over all means of the four different dimensions (valence, arousal, distance, angle)
for (d in mean_names){
  # find index of dimension mean in dataSummary_cr
  p <- match(d, names(dataSummary_cr))
  # create a df with only one of the four continuous measurements
  mean_df <- dataSummary_cr %>%
    dplyr::select('cr_time','quadrant','rating_method',all_of(p)) %>%
    group_by(quadrant,cr_time, rating_method)
  
  # save the dimension mean as title for figure 
  title <- d
  # rename name of mean for easier plotting
  names(mean_df) <- str_replace(names(mean_df),d,'cr')
  names(data_cr) <- 
  
  # create a plot that shows the average rating for each rating method, condition, and measurement.
  # Define the y limits depending on dimension
  lim_y <- case_when(d == 'mean Valence' ~ c(-1,1),
                     d == 'mean Arousal' ~ c(-1,1),
                     d == 'mean Distance' ~ c(-sqrt(2),sqrt(2)),
                     d == 'mean Angle' ~ c(-pi,pi))
  
  plot <- ggplot(mean_df, aes(x = cr_time, y = cr)) +
    geom_line(aes(color = rating_method, group = rating_method), linewidth = 1) +
    #geom_line(cr_dim, aes(x='cr_time', y='cr_' + rat_dim[dim], group='sub', color='sub'), size=0.25, alpha=0.25)
    facet_wrap(.~quadrant) +
    scale_x_discrete(breaks = c(0,5,seq(10, 60, 10))) +
    theme_bw(base_size = 20) +
    theme(plot.background = element_blank(), legend.background = element_blank()) +
    #ggtitle(paste(title,'by Rating Method and Video/Quadrant')) +
    ylim(lim_y) +
    xlab('Time (s)') +
    #ylab(paste('Continuous Rating ', title)) +
    ylab(title) +
    labs(color = 'Feedback')
  
  plot
  # save the file as vector-based graphic with type of measurement in the name
   ggsave(file= file.path(save_path, paste(sub(' ', '_',title),'CR.png',sep='_')), plot = plot, height=170,  width= 297, units = 'mm', dpi = 600)
   ggsave(file= file.path(save_path, paste(sub(' ', '_',title),'CR.svg',sep='_')), plot = plot, height=170,  width= 297, units = 'mm', dpi = 600)
   
}


# Individual plots - CR time-series per rating methods and videos/ --------

# loop over every participant
sjs <- unique(data_cr$sj_id)
for (i in sjs){
  
  # create a new df for that participant
  df_sj <- data_cr %>%
    filter(sj_id == i)%>%
    group_by(cr_time, quadrant, rating_method) 
  
  # exclude participants with no data
  if (nrow(df_sj) == 0) next # skip iteration if df for one participant is empty (e.g. subject no 3)
  
  # create file path for every participant so save corresponding plots
  dDir <- file.path(save_path, 'individual_plots',i)
  if (!file.exists(dDir)) {dir.create(dDir)}
  
  # loop over the four dimensions (valence, arousal, distance, angle) for every participant
  for (d in 1:length(dimensions)){
    # find index of dimension cr in df_sj
    p <- match(paste0('cr_',dimensions[d]), names(df_sj))
    # create a new df for one measurement 
    df_sj_dim <- df_sj %>%
      dplyr::select('cr_time','quadrant','rating_method',p) %>%
      group_by(quadrant,cr_time, rating_method)
   
    # get title for figure
    title <- dimensions_names[d]
    # change variable name for easier plotting
    names(df_sj_dim) <- str_replace(names(df_sj_dim),paste0('cr_',dimensions[d]),'cr')
    
    # Define the y limits depending on dimension
    lim_y <- case_when(dimensions[d] == 'v' ~ c(-1,1),
                       dimensions[d] == 'a' ~ c(-1,1),
                       dimensions[d] == 'dist' ~ c(-sqrt(2),sqrt(2)),
                       dimensions[d] == 'angle' ~ c(-pi,pi))
    
    # plot rating method, for every condition across time
    plot <- ggplot(df_sj_dim, aes(x = cr_time, y = cr)) +
      geom_line(aes(color = rating_method, group = rating_method), size = 1) +
      facet_wrap(.~quadrant) +
      ylim(lim_y) +
      scale_x_discrete(breaks = seq(5, 60, 10)) +
      theme_bw(base_size = 14)+
      ggtitle(paste(title,'Continuous Ratings by Rating Method and Video/Quadrant'))+
      xlab('Time (s)')+
      ylab(paste(title, 'Continuous Ratings'))+
      labs(color = 'Rating Method')
    
    #save file
    ggsave(file= file.path(dDir, paste(title,".png", sep ='')), plot = plot, height=210,  width= 297, units = 'mm')
    
  }
  
}


# Summary descriptive statistics ------------------------------------------

# create a summary dataframe of mean, sd, skewness and kurtosis for every rating method and dimension and per quadrant (condition)
digit <- 2 # round to 'digit' digits

mean_cr_time <- data_cr  %>% 
  group_by(sj_id, test_site, quadrant, rating_method) %>% 
  pivot_longer(c(cr_v:cr_angle)) %>% 
  group_by(sj_id, test_site, rating_method , quadrant, name) %>% 
  dplyr::summarise(mean = mean(value, na.rm = T)) %>% 
  dplyr::rename(var = name) %>% 
  pivot_longer(mean) %>% 
  unite(name, c(var, name)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate_if(is.numeric, ~round(., digit))

descriptiveSummary <- mean_cr_time  %>% 
  group_by(quadrant, rating_method, test_site) %>% #!grouped by test_site for following plot
  pivot_longer(c(cr_a_mean:cr_v_mean)) %>% 
  group_by(rating_method , quadrant, test_site, name) %>% 
  dplyr::summarise(mean = mean(value, na.rm = T),
            sd = sd(value, na.rm = T)) %>% 
  dplyr::rename(var = name) %>% 
  pivot_longer(mean:sd) %>% 
  unite(name, c(var, name)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate_if(is.numeric, ~round(., digit))

descriptiveSummary2 <- mean_cr_time  %>% 
  group_by(quadrant, rating_method) %>%
  pivot_longer(c(cr_a_mean:cr_v_mean)) %>% 
  group_by(rating_method , quadrant, name) %>% 
  dplyr::summarise(mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T)) %>% 
  dplyr::rename(var = name) %>% 
  pivot_longer(mean:sd) %>% 
  unite(name, c(var, name)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate_if(is.numeric, ~round(., digit))

  
#save summary df
write_csv(descriptiveSummary, file.path(save_path,"cr_descriptive_summary.csv"))

# Relationship between valence and arousal --------------------------------

# CR Mean across subjects
# Separated by testing sites and RMs
plot <- ggplot(descriptiveSummary, aes(x = cr_v_mean, y = cr_a_mean, color = rating_method, shape = test_site)) +
  geom_point(size = 2) +
  # geom_errorbar(aes(
  #   ymin = cr_a_mean - cr_a_sd,
  #   ymax = cr_a_mean + cr_a_sd,
  #   x = cr_v_mean,
  #   width = 0.1
  # )) +
  # geom_errorbarh(aes(
  #   xmin = cr_v_mean - cr_v_sd,
  #   xmax = cr_v_mean + cr_v_sd,
  #   y = cr_a_mean,
  #   height = 0.1
  # )) +
  theme_bw(base_size = 14) +
  ggtitle('Mean across participants - Arousal vs. Valence') +
  xlab('Mean Valence') +
  ylab('Mean Arousal') +
  labs(color = 'Rating Method', shape = 'Testing Site') +
  ylim(c(-1,1)) +
  xlim(c(-1,1)) +
  geom_vline(xintercept = 0, color = "black") +
  geom_hline(yintercept = 0, color = "black")
plot
#save file
ggsave(file= file.path(save_path, "mean_valence_arousal_test_sites.png"), plot = plot, height=210,  width= 297, units = 'mm')

# Separated by RMs
# plot <- ggplot(descriptiveSummary2, aes(x = cr_v_mean_mean, y = cr_a_mean_mean, fill = rating_method, shape = rating_method)) +
#   geom_point(data = mean_cr_time, aes(x = cr_v_mean, y = cr_a_mean), alpha = .5, size = 1.5, stroke = 0) +
#   geom_point(size = 4) +
#   scale_shape_manual(values = c(21, 23, 24)) +
#   # geom_errorbar(aes(
#   #   ymin = cr_a_mean_mean - cr_a_mean_sd,
#   #   ymax = cr_a_mean_mean + cr_a_mean_sd,
#   #   x = cr_v_mean_mean,
#   #   width = 0.1
#   # )) +
#   # geom_errorbarh(aes(
#   #   xmin = cr_v_mean_mean - cr_v_mean_sd,
#   #   xmax = cr_v_mean_mean + cr_v_mean_sd,
#   #   y = cr_a_mean_mean,
#   #   height = 0.1
#   # )) +
#   theme_bw(base_size = 20) +
#   #ggtitle('Mean across participants - Arousal vs. Valence') +
#   xlab('Mean Valence') +
#   ylab('Mean Arousal') +
#   labs(color = 'Rating Method (RM)') +
#   ylim(c(-1,1)) +
#   xlim(c(-1,1)) +
#   geom_vline(xintercept = 0, color = "black") +
#   geom_hline(yintercept = 0, color = "black") +
#   theme(plot.background = element_blank(), legend.background = element_blank())
# plot

# Separated by RMs
ticks <- data.frame(x = c(-.5, -.5, .5, .5, -1, -1, 1, 1, .05, -.05, .05, -.05, .05, -.05, .05, -.05),
                    y = c(.05, -.05, .05, -.05, .05, -.05, .05, -.05, -.5, -.5, .5, .5, -1, -1, 1, 1),
                    tick = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8),
                    x.label = c(-.5, -.5, .5, .5, -1, -1, 1, 1, .08, .08, -.08, -.08, .08, .08, -.08, -.08),
                    y.label = c(-.12, -.12, .12, .12, -.12, -.12, .12, .12, -.5, -.5, .5, .5, -1, -1, 1, 1),
                    label = c("-0.5", "-0.5", "0.5", "0.5", rep(NA, 4), #"-1.0", "-1.0", "1.0", "1.0"
                              "-0.5", "-0.5", "0.5", "0.5", rep(NA, 4)), #"-1.0", "-1.0", "1.0", "1.0"
                    hjust = c(rep(NA, 8), 0, 0, 1, 1, 0, 0, 1, 1))

violins <- mean_cr_time %>% 
  # mutate(y = case_when(cr_v_mean > 0 & cr_a_mean > 0 ~ 1.4, 
  #                      cr_v_mean < 0 & cr_a_mean > 0 ~ 1.4, 
  #                      cr_v_mean > 0 & cr_a_mean < 0 ~ -1.4,
  #                      cr_v_mean < 0 & cr_a_mean < 0 ~ -1.4),
  #        x = case_when(cr_a_mean > 0 & cr_v_mean > 0 ~ 1.4, 
  #                      cr_a_mean < 0 & cr_v_mean > 0 ~ 1.4,
  #                      cr_a_mean > 0 & cr_v_mean < 0 ~ -1.4, 
  #                      cr_a_mean < 0 & cr_v_mean < 0 ~ -1.4))
  mutate(y = case_when(quadrant %in% c("HP", "HN") ~ 1.9,
                       quadrant %in% c("LP", "LN") ~ -1.4),
         x = case_when(quadrant %in% c("HP", "LP") ~ 1.4,
                       quadrant %in% c("HN", "LN") ~ -1.9))

p.grid <- ggplot(descriptiveSummary2, aes(x = cr_v_mean_mean, y = cr_a_mean_mean, fill = rating_method, shape = quadrant)) +
  geom_point(data = mean_cr_time, aes(x = cr_v_mean, y = cr_a_mean), alpha = .2, size = 2, stroke = 0, show.legend = F) +
  geom_point(size = 4, show.legend = F) +
  scale_shape_manual(values = c(22, 21, 24, 23)) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  
  geom_segment(aes(x = -1, xend = 1, y = 1, yend = 1), linewidth = .25) +
  geom_segment(aes(y = -1, yend = 1, x = 1, xend = 1), linewidth = .25) +
  geom_segment(aes(x = -1, xend = 1, y = -1, yend = -1), linewidth = .25) +
  geom_segment(aes(y = -1, yend = 1, x = -1, xend = -1), linewidth = .25) +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linewidth = .25) +
  geom_segment(aes(y = -1, yend = 1, x = 0, xend = 0), linewidth = .25) +
  
  geom_line(data = ticks, inherit.aes = NULL, aes(x = x, y = y, group = tick)) +
  # geom_text(data = ticks %>% distinct(x.label, y.label, label, hjust), inherit.aes = NULL, aes(x = x.label, y = y.label, label = label, hjust = hjust)) +
  
  xlab('Mean Valence') +
  ylab('Mean Arousal') +
  labs(fill = 'Rating Method (RM)',
       shape = '360Â° Video') +
  theme_void(base_size = 20) +
  # theme(plot.background = element_blank(), 
  #       legend.background = element_blank(),
  #       panel.border = element_blank(),
  #       axis.ticks = element_blank(),
  #       axis.text = element_blank(),
  #       axis.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))


p.top <- ggplot(violins %>% filter(quadrant %in% c("HP", "HN")), aes(x = cr_v_mean, group = interaction(quadrant, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method, linetype = quadrant), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1))

p.bottom <- ggplot(violins %>% filter(quadrant %in% c("LP", "LN")), aes(x = cr_v_mean, group = interaction(quadrant, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method, linetype = quadrant), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_y_reverse()

p.right <- ggplot(violins %>% filter(quadrant %in% c("HP", "LP")), aes(y = cr_a_mean, group = interaction(quadrant, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method, linetype = quadrant), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_reverse()

p.left <- ggplot(violins %>% filter(quadrant %in% c("HN", "LN")), aes(y = cr_a_mean, group = interaction(quadrant, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method, linetype = quadrant), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  scale_linetype_manual(values = c("solid", "longdash")) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1))

grid.arrange(p.top, p.right, p.grid, p.left, p.bottom,
             widths = c(1, 3, 1),
             heights = c(1, 3, 1),
             layout_matrix = rbind(c(NA, 1, NA),
                                   c(2:4),
                                   c(NA, 5, NA)))

p.top1 <- ggplot(violins %>% filter(quadrant == "HP"), aes(x = cr_v_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1))
p.top2 <- ggplot(violins %>% filter(quadrant == "HN"), aes(x = cr_v_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1))

p.bottom1 <- ggplot(violins %>% filter(quadrant == "LP"), aes(x = cr_v_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_y_reverse()
p.bottom2 <- ggplot(violins %>% filter(quadrant == "LN"), aes(x = cr_v_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_y_reverse()

p.left1 <- ggplot(violins %>% filter(quadrant == "HP"), aes(y = cr_a_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1))
p.left2 <- ggplot(violins %>% filter(quadrant == "LP"), aes(y = cr_a_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1))

p.right1 <- ggplot(violins %>% filter(quadrant == "HN"), aes(y = cr_a_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_reverse()
p.right2 <- ggplot(violins %>% filter(quadrant == "LN"), aes(y = cr_a_mean, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_reverse()

p.all <- grid.arrange(p.top1, p.top2, p.right1, p.right2, p.grid, p.left1, p.left2, p.bottom1, p.bottom2,
             widths = c(1, 1, 6, 1, 1),
             heights = c(1, 1, 6, 1, 1),
             layout_matrix = rbind(c(NA, NA, 1, NA, NA),
                                   c(NA, NA, 2, NA, NA),
                                   c(4, 3, 5, 7, 6),
                                   c(NA, NA, 8, NA, NA),
                                   c(NA, NA, 9, NA, NA)))
#save
ggsave(file.path(save_path, "CR_mean_valence_arousal_RMs_grid.svg"), p.all, height = 160, width = 160, unit = "mm", dpi = 600)

###

plot <- ggplot(descriptiveSummary2, aes(x = cr_v_mean_mean, y = cr_a_mean_mean, fill = rating_method, shape = quadrant)) +
  geom_point(data = mean_cr_time, aes(x = cr_v_mean, y = cr_a_mean), alpha = .2, size = 2, stroke = 0) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(22, 21, 24, 23)) +
  
  geom_rect(aes(xmin = -2.3, ymin = 1.1, xmax = 2.3, ymax = 2.3), fill = "white") +
  geom_rect(aes(xmin = -2.3, ymin = -1.1, xmax = 2.3, ymax = -2.3), fill = "white") +
  geom_rect(aes(ymin = -2.3, xmin = 1.1, ymax = 2.3, xmax = 2.3), fill = "white") +
  geom_rect(aes(ymin = -2.3, xmin = -1.1, ymax = 2.3, xmax = -2.3), fill = "white") +
  
  geom_segment(aes(x = -1.1, xend = 1.1, y = 1.1, yend = 1.1), linewidth = .25) +
  geom_segment(aes(y = -1.1, yend = 1.1, x = 1.1, xend = 1.1), linewidth = .25) +
  geom_segment(aes(x = -1.1, xend = 1.1, y = -1.1, yend = -1.1), linewidth = .25) +
  geom_segment(aes(y = -1.1, yend = 1.1, x = -1.1, xend = -1.1), linewidth = .25) +
  geom_segment(aes(x = -1.1, xend = 1.1, y = 0, yend = 0), linewidth = .25) +
  geom_segment(aes(y = -1.1, yend = 1.1, x = 0, xend = 0), linewidth = .25) +
  
  # geom_density(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = cr_v_mean, y = after_stat(scaled), group = interaction(y, rating_method), fill = rating_method),
  #              color = NA, alpha = .4)
  # 
  
  
  
  # geom_violin(data = violins %>% filter(quadrant %in% c("HP", "HN")), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, quadrant)), 
  #             position = position_dodge(.5), width = .5, alpha = 1, color = "black", show.legend = F) + 
  # geom_violin(data = violins %>% filter(quadrant %in% c("LP", "LN")), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, quadrant)), 
  #             position = position_dodge(.5), width = .5, alpha = 1, color = "black", show.legend = F) + 
  # geom_violin(data = violins %>% filter(quadrant %in% c("HP", "LP")), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(y, quadrant)), 
  #             position = position_dodge(.5), width = .5, alpha = 1, color = "black", show.legend = F) + 
  # geom_violin(data = violins %>% filter(quadrant %in% c("HN", "LN")), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(y, quadrant)), 
  #             position = position_dodge(.5), width = .5, alpha = 1, color = "black", show.legend = F) + 
  
  geom_violin(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = cr_v_mean, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = cr_v_mean, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  
  geom_violin(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = x + .5, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = x + .5, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  
  geom_boxplot(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = cr_v_mean, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = cr_v_mean, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  
  geom_boxplot(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = x + .5, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = x + .5, y = cr_a_mean, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  
  # geom_violin(data = violins %>% filter(cr_v_mean > 0), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  # geom_violin(data = violins %>% filter(cr_v_mean < 0), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  # geom_violin(data = violins %>% filter(cr_a_mean > 0), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(x, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  # geom_violin(data = violins %>% filter(cr_a_mean < 0), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(x, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +

  # geom_boxplot(data = violins %>% filter(cr_v_mean > 0), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  # geom_boxplot(data = violins %>% filter(cr_v_mean < 0), inherit.aes = NULL, aes(x = cr_v_mean, y = y, group = interaction(y, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  # geom_boxplot(data = violins %>% filter(cr_a_mean > 0), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(x, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  # geom_boxplot(data = violins %>% filter(cr_a_mean < 0), inherit.aes = NULL, aes(x = x, y = cr_a_mean, group = interaction(x, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
 
  geom_line(data = ticks, inherit.aes = NULL, aes(x = x, y = y, group = tick)) + 
  geom_text(data = ticks %>% distinct(x.label, y.label, label, hjust), inherit.aes = NULL, aes(x = x.label, y = y.label, label = label, hjust = hjust)) +
  
  theme_bw(base_size = 20) +
  #ggtitle('Mean across participants - Arousal vs. Valence') +
  xlab('Mean Valence') +
  ylab('Mean Arousal') +
  labs(fill = 'Rating Method (RM)',
       shape = '360Â° Video') +
  # ylim(c(-1,1)) +
  # xlim(c(-1,1)) +
  # geom_vline(xintercept = 0, color = "black") +
  # geom_hline(yintercept = 0, color = "black") +
  theme(plot.background = element_blank(), 
        legend.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) + 
  coord_cartesian(xlim = c(-2.1, 2.1), ylim = c(-2.1, 2.1))

plot
#save file
ggsave(file= file.path(save_path, "CR_mean_valence_arousal_RMs.png"), plot = plot, height=210,  width= 297, units = 'mm', dpi = 600)
ggsave(file= file.path(save_path, "CR_mean_valence_arousal_RMs.svg"), plot = plot, height=210,  width= 297, units = 'mm', dpi = 600)

# All data points
# create dataframe for 'additive' colors
toPlot <- data_cr %>% 
  select(rating_method, cr_v, cr_a) %>% 
  mutate_if(is.numeric, ~round(., digit)) %>% # round values
  mutate(coordinates = paste0(cr_v, ",", cr_a)) %>% # get coordinates from cr_v and cr_a
  select(-c(cr_v, cr_a)) %>% 
  group_by(rating_method, coordinates) %>% 
  mutate(n = n()) %>% # count number of cases (i.e., how many times a coordinate is hit)
  distinct() %>% # only unique rows
  pivot_wider(names_from = rating_method, values_from = n) %>% # create variables from each rating method containing the number of cases
  mutate(ratings = paste0(ifelse(!is.na(Grid), "Grid/", ""), 
                          ifelse(!is.na(Flubber), "Flubber/", ""), 
                          ifelse(!is.na(Proprioceptive), "Proprioceptive", "")), # paste together names of ratings that have non-zero number of ratings
         ratings = case_when(str_sub(ratings, -1, -1) == "/" ~ str_sub(ratings, 1, -2), # drop last slash if necessary
                             T ~ ratings), # otherwise keep rating name
         nRatings = sum(c(Grid, Flubber, Proprioceptive), na.rm = T)) %>% # tally total number of ratings across rating methods
  select(ratings, coordinates, nRatings) %>% 
  separate(coordinates, c("cr_v", "cr_a"), sep = ",") %>% # separate coordinates into cr_v and cr_a again
  mutate(across(c(cr_v, cr_a), as.numeric),
         ratings = factor(ratings, levels = c("Grid", "Flubber", "Proprioceptive",
                                              "Grid/Flubber", "Grid/Proprioceptive", "Flubber/Proprioceptive",
                                              "Grid/Flubber/Proprioceptive"))) %>% 
  uncount(nRatings) # create nRatings number of identical rows for each case

# plot dataframe for 'additive' colors
plot <- ggplot(toPlot, aes(x = cr_v, y = cr_a)) +
  geom_point(aes(color = ratings), size = 0.5, alpha = 0.1) +
  theme_bw(base_size = 14) +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1))) +
  scale_color_manual(values = c(hue_pal()(3), "#7C9853", "#AD89B6", "#31AB9C", "#5C2C26")) +
  ggtitle('Arousal vs. Valence') +
  xlab('Valence') +
  ylab('Arousal') +
  labs(color = 'Rating Method')
plot
#save file
ggsave(file= file.path(save_path, "valence_arousal.png"), plot = plot, height=210,  width= 297, units = 'mm')

# SR mean
digit <- 2 # round to 'digit' digits
# Separated by RMs
sr <- data_cri_sr %>% select(sj_id, test_site, rating_method, quadrant, sr_v, sr_a)
sr_mean <- sr  %>% 
  group_by(quadrant, rating_method) %>%
  pivot_longer(c(sr_v:sr_a)) %>% 
  group_by(rating_method , quadrant, name) %>% 
  dplyr::summarise(mean = mean(value, na.rm = T),
                   sd = sd(value, na.rm = T)) %>% 
  dplyr::rename(var = name) %>% 
  pivot_longer(mean:sd) %>% 
  unite(name, c(var, name)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate_if(is.numeric, ~round(., digit))


ticks <- data.frame(x = c(-.5, -.5, .5, .5, -1, -1, 1, 1, .05, -.05, .05, -.05, .05, -.05, .05, -.05),
                    y = c(.05, -.05, .05, -.05, .05, -.05, .05, -.05, -.5, -.5, .5, .5, -1, -1, 1, 1),
                    tick = c(1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8),
                    x.label = c(-.5, -.5, .5, .5, -1, -1, 1, 1, .08, .08, -.08, -.08, .08, .08, -.08, -.08),
                    y.label = c(-.12, -.12, .12, .12, -.12, -.12, .12, .12, -.5, -.5, .5, .5, -1, -1, 1, 1),
                    label = c("-0.5", "-0.5", "0.5", "0.5", rep(NA, 4), #"-1.0", "-1.0", "1.0", "1.0"
                              "-0.5", "-0.5", "0.5", "0.5", rep(NA, 4)), #"-1.0", "-1.0", "1.0", "1.0"
                    hjust = c(rep(NA, 8), 0, 0, 1, 1, 0, 0, 1, 1))

violins <- sr %>% 
  # mutate(y = case_when(cr_v_mean > 0 & cr_a_mean > 0 ~ 1.4, 
  #                      cr_v_mean < 0 & cr_a_mean > 0 ~ 1.4, 
  #                      cr_v_mean > 0 & cr_a_mean < 0 ~ -1.4,
  #                      cr_v_mean < 0 & cr_a_mean < 0 ~ -1.4),
  #        x = case_when(cr_a_mean > 0 & cr_v_mean > 0 ~ 1.4, 
  #                      cr_a_mean < 0 & cr_v_mean > 0 ~ 1.4,
  #                      cr_a_mean > 0 & cr_v_mean < 0 ~ -1.4, 
  #                      cr_a_mean < 0 & cr_v_mean < 0 ~ -1.4))
  mutate(y = case_when(quadrant %in% c("HP", "HN") ~ 1.9,
                       quadrant %in% c("LP", "LN") ~ -1.4),
         x = case_when(quadrant %in% c("HP", "LP") ~ 1.4,
                       quadrant %in% c("HN", "LN") ~ -1.9))

p.grid <- ggplot(sr_mean, aes(x = sr_v_mean, y = sr_a_mean, fill = rating_method, shape = quadrant)) +
  geom_point(data = sr, aes(x = sr_v, y = sr_a), alpha = .2, size = 2, stroke = 0, show.legend = F) +
  geom_point(size = 4, show.legend = F) +
  scale_shape_manual(values = c(22, 21, 24, 23)) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  
  geom_segment(aes(x = -1, xend = 1, y = 1, yend = 1), linewidth = .25) +
  geom_segment(aes(y = -1, yend = 1, x = 1, xend = 1), linewidth = .25) +
  geom_segment(aes(x = -1, xend = 1, y = -1, yend = -1), linewidth = .25) +
  geom_segment(aes(y = -1, yend = 1, x = -1, xend = -1), linewidth = .25) +
  geom_segment(aes(x = -1, xend = 1, y = 0, yend = 0), linewidth = .25) +
  geom_segment(aes(y = -1, yend = 1, x = 0, xend = 0), linewidth = .25) +
  
  geom_line(data = ticks, inherit.aes = NULL, aes(x = x, y = y, group = tick)) +
  # geom_text(data = ticks %>% distinct(x.label, y.label, label, hjust), inherit.aes = NULL, aes(x = x.label, y = y.label, label = label, hjust = hjust)) +
  
  xlab('Mean Valence') +
  ylab('Mean Arousal') +
  labs(fill = 'Feedback',
       shape = '360° Video') +
  theme_void(base_size = 20) +
  # theme(plot.background = element_blank(), 
  #       legend.background = element_blank(),
  #       panel.border = element_blank(),
  #       axis.ticks = element_blank(),
  #       axis.text = element_blank(),
  #       axis.title = element_blank()) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1))

p.top1 <- ggplot(violins %>% filter(quadrant == "HP"), aes(x = sr_v, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1))
p.top2 <- ggplot(violins %>% filter(quadrant == "HN"), aes(x = sr_v, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1))

p.bottom1 <- ggplot(violins %>% filter(quadrant == "LP"), aes(x = sr_v, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_y_reverse()
p.bottom2 <- ggplot(violins %>% filter(quadrant == "LN"), aes(x = sr_v, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-1, 1)) +
  scale_y_reverse()

p.left1 <- ggplot(violins %>% filter(quadrant == "HP"), aes(y = sr_a, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1))
p.left2 <- ggplot(violins %>% filter(quadrant == "LP"), aes(y = sr_a, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1))

p.right1 <- ggplot(violins %>% filter(quadrant == "HN"), aes(y = sr_a, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_reverse()
p.right2 <- ggplot(violins %>% filter(quadrant == "LN"), aes(y = sr_a, group = interaction(y, rating_method))) +
  geom_density(aes(color = rating_method, fill = rating_method), alpha = .1, show.legend = F) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  scale_color_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  theme_void(base_size = 20) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(ylim = c(-1, 1)) +
  scale_x_reverse()

p.all <- grid.arrange(p.top1, p.top2, p.right1, p.right2, p.grid, p.left1, p.left2, p.bottom1, p.bottom2,
                      widths = c(1, 1, 6, 1, 1),
                      heights = c(1, 1, 6, 1, 1),
                      layout_matrix = rbind(c(NA, NA, 1, NA, NA),
                                            c(NA, NA, 2, NA, NA),
                                            c(4, 3, 5, 7, 6),
                                            c(NA, NA, 8, NA, NA),
                                            c(NA, NA, 9, NA, NA)))

#save
ggsave(file.path(save_path, "SR_mean_valence_arousal_RMs_grid.svg"), p.all, height = 160, width = 160, unit = "mm", dpi = 600)

###
violins <- sr %>% 
  # mutate(y = case_when(sr_v > 0 & sr_a > 0 ~ 1.4, 
  #                      sr_v < 0 & sr_a > 0 ~ 1.4, 
  #                      sr_v > 0 & sr_a < 0 ~ -1.4,
  #                      sr_v < 0 & sr_a < 0 ~ -1.4),
  #        x = case_when(sr_a > 0 & sr_v > 0 ~ 1.4, 
  #                      sr_a < 0 & sr_v > 0 ~ 1.4,
  #                      sr_a > 0 & sr_v < 0 ~ -1.4, 
  #                      sr_a < 0 & sr_v < 0 ~ -1.4))
  mutate(y = case_when(quadrant %in% c("HP", "HN") ~ 1.9,
                       quadrant %in% c("LP", "LN") ~ -1.4),
         x = case_when(quadrant %in% c("HP", "LP") ~ 1.4,
                       quadrant %in% c("HN", "LN") ~ -1.9))

plot <- ggplot(sr_mean, aes(x = sr_v_mean, y = sr_a_mean, fill = rating_method, shape = quadrant)) +
  geom_point(data = sr, aes(x = sr_v, y = sr_a), alpha = .2, size = 2, stroke = 0) +
  geom_point(size = 4) +
  scale_shape_manual(values = c(22, 21, 24, 23)) +
  
  geom_rect(aes(xmin = -2.3, ymin = 1.1, xmax = 2.3, ymax = 2.3), fill = "white") +
  geom_rect(aes(xmin = -2.3, ymin = -1.1, xmax = 2.3, ymax = -2.3), fill = "white") +
  geom_rect(aes(ymin = -2.3, xmin = 1.1, ymax = 2.3, xmax = 2.3), fill = "white") +
  geom_rect(aes(ymin = -2.3, xmin = -1.1, ymax = 2.3, xmax = -2.3), fill = "white") +
  
  geom_segment(aes(x = -1.1, xend = 1.1, y = 1.1, yend = 1.1), linewidth = .25) +
  geom_segment(aes(y = -1.1, yend = 1.1, x = 1.1, xend = 1.1), linewidth = .25) +
  geom_segment(aes(x = -1.1, xend = 1.1, y = -1.1, yend = -1.1), linewidth = .25) +
  geom_segment(aes(y = -1.1, yend = 1.1, x = -1.1, xend = -1.1), linewidth = .25) +
  geom_segment(aes(x = -1.1, xend = 1.1, y = 0, yend = 0), linewidth = .25) +
  geom_segment(aes(y = -1.1, yend = 1.1, x = 0, xend = 0), linewidth = .25) +
  
  geom_violin(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = sr_v, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = sr_v, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  
  geom_violin(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = x + .5, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  geom_violin(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = x + .5, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
              position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  
  geom_boxplot(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = sr_v, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = sr_v, y = y - .5, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  
  geom_boxplot(data = violins %>% filter(quadrant == "HP"), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "HN"), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LP"), inherit.aes = NULL, aes(x = x + .5, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  geom_boxplot(data = violins %>% filter(quadrant == "LN"), inherit.aes = NULL, aes(x = x + .5, y = sr_a, group = interaction(y, rating_method), fill = rating_method), 
               outlier.shape = NA, position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  
  # geom_violin(data = violins %>% filter(sr_v > 0), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) + 
  # geom_violin(data = violins %>% filter(sr_v < 0), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) + 
  # geom_violin(data = violins %>% filter(sr_a > 0), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(x, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  # geom_violin(data = violins %>% filter(sr_a < 0), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(x, rating_method), fill = rating_method), position = position_dodge(.5), width = .5, alpha = 1, color = NA, show.legend = F) +
  # 
  # geom_boxplot(data = violins %>% filter(sr_v > 0), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) + 
  # geom_boxplot(data = violins %>% filter(sr_v < 0), inherit.aes = NULL, aes(x = sr_v, y = y, group = interaction(y, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) + 
  # geom_boxplot(data = violins %>% filter(sr_a > 0), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(x, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  # geom_boxplot(data = violins %>% filter(sr_a < 0), inherit.aes = NULL, aes(x = x, y = sr_a, group = interaction(x, rating_method)), position = position_dodge(.5), width = .1, fill = "white", show.legend = F) +
  
  geom_line(data = ticks, inherit.aes = NULL, aes(x = x, y = y, group = tick)) + 
  geom_text(data = ticks %>% distinct(x.label, y.label, label, hjust), inherit.aes = NULL, aes(x = x.label, y = y.label, label = label, hjust = hjust)) +
  
  theme_bw(base_size = 20) +
  #ggtitle('Mean across participants - Arousal vs. Valence') +
  xlab('Mean Valence') +
  ylab('Mean Arousal') +
  labs(fill = 'Rating Method (RM)',
       shape = '360° Video') +
  # ylim(c(-1,1)) +
  # xlim(c(-1,1)) +
  # geom_vline(xintercept = 0, color = "black") +
  # geom_hline(yintercept = 0, color = "black") +
  theme(plot.background = element_blank(), 
        legend.background = element_blank(),
        panel.border = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank()) +
  scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF', '#C77CFF')) +
  guides(fill = guide_legend(override.aes = list(shape = 21))) +
  coord_cartesian(xlim = c(-2.1, 2.1), ylim = c(-2.1, 2.1))

plot
#save file
ggsave(file= file.path(save_path, "SR_mean_valence_arousal_RMs.png"), plot = plot, height=210,  width= 297, units = 'mm', dpi = 600)
ggsave(file= file.path(save_path, "SR_mean_valence_arousal_RMs.svg"), plot = plot, height=210,  width= 297, units = 'mm', dpi = 600)



# CR-SR correlation plots -------------------------------------------------

# Keep only mean as metric of interest
metric_select <- c('cr_mean')
metric_select_name <- str_remove(metrics, "cr_")
# keep only valence and arousal as dimension of interest
dim_select <- dimensions[1:2]
dim_select_names <- dimensions_names[1:2]

cor_labels <- array(c(c("r = 0.950", "r = 0.922", "r = 0.920"),c("r = 0.938", "r = 0.873", "r = 0.852")), dim = c(3,2))


# Loop over valence and arousal
for (d in 1:length(dim_select)){
  # metric with dimension
  metric_d <- c(paste(metric_select, dimensions[d], sep="_"))
  # sr with dimension
  sr_d <- c(paste('sr', dimensions[d], sep="_"))

  # subset data
  cri_sr_select <- data_cri_sr %>% filter(rating_method != base_name) %>%
    dplyr::select('sj_id','test_site', 'quadrant','rating_method', all_of(metric_d), all_of(sr_d)) %>%
    group_by(sj_id, test_site, quadrant, rating_method)
  
  # rename metric_d and sr_d for easier plotting
  names(cri_sr_select) <- str_replace(names(cri_sr_select),metric_d,'cr')
  names(cri_sr_select) <- str_replace(names(cri_sr_select),sr_d,'sr')
  
  p <- ggplot(cri_sr_select, aes(x = sr, y = cr, fill = rating_method)) +
    geom_point(aes(shape=quadrant)) +
    geom_smooth(method='lm', aes(color=rating_method)) +
    #facet_wrap(.~quadrant) +
    annotate("text",x=-0.7,y=c(0.9, 0.75, 0.6),label=cor_labels[,d], color = c('#F8766D', '#00BA38', '#619CFF'), size = 8) +
    scale_fill_manual(values = c('#F8766D', '#00BA38', '#619CFF')) +
    scale_shape_manual(values = c(22, 21, 24, 23)) +
    labs(fill = 'Rating Method (RM)',
         shape = '360Â° Video',
         color = 'Rating Method (RM)') +
    ggtitle(dim_select_names[d],) +
    xlab('SR') +
    ylab('CR mean') +
    theme_bw(base_size = 20) +
    ylim(c(-1,1)) +
    xlim(c(-1,1)) +
    theme(plot.title = element_text(hjust = 0.5), plot.background = element_blank(), legend.background = element_blank())
  p
  #save file
  ggsave(file= file.path(save_path, paste0("corr_cr_mean_sr_",dim_select_names[d], ".png")), plot = p, height=210,  width= 297, units = 'mm', dpi = 600)
  ggsave(file= file.path(save_path, paste0("corr_cr_mean_sr_",dim_select_names[d], ".svg")), plot = p, height=210,  width= 297, units = 'mm', dpi = 600)
   
}
