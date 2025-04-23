# CR indices (CRi) and SR - Correlation and cocor analyses
# 1. Correlation including all data
# 2. Correlation by rating method for specific CRi (metrics) and comparisons
# 3. Intraclass correlation coefficient (ICC) between the two testing sites
# 4. Comparison CRi between RMs
# 5. Comparison SR between RMs including Baseline
# Author: Antonin Fourcade
# Last version: 07.11.2024

# import packages
library(tidyverse)
library(ggpmisc)
library(ggdist)
library(gghalves)
library(hrbrthemes)
library(viridis)
library(car)
library(readxl)
library(cocor)
library(irr)
library(rstatix)
library(TOSTER)

# set paths
data_path <- "E:/AffectiveVR/Phase_1/Data/AVR/"

# read data
data <- read.csv(paste(data_path, "cri_sr_clean.csv", sep = ""))

# set names of test sites
test_list <- c("Torino", "Berlin")
# set names of rating methods and baseline
rm_list <- c("Grid", "Flubber", "Proprioceptive")
base_name <- c("Baseline")
# set names of quadrants/videos
quad_list <- c("HP", "LP", "LN", "HN")

# format data
data <- data %>% 
  mutate(
    test_site = factor(test_site, levels = test_list),
    rating_method = factor(rating_method, levels = c(rm_list, base_name)),
    quadrant = factor(quadrant, levels = quad_list),
    sj_id = factor(sj_id, levels = unique(data$sj_id)))

#Count non-NaN datapoints per metric
colSums(!is.na(data))

# summary ratings (SR)
sr <- data %>% dplyr::select(!starts_with("cr"))
# Affective dimensions
dimensions <- unlist(str_split(names(sr), "_"))[seq(9, length(sr)*2, 2)]
dimensions_names <- str_replace_all(dimensions, c("angle" = "Angle", "a" = "Arousal", "v" = "Valence", "dist" = "Distance"))

# Correlation including all data and pairwise comparison ------------------

corResults_path <- "E:/AffectiveVR/Phase_1/affectivevr/corr_results/"
if (!file.exists(corResults_path)) {dir.create(corResults_path)}
setwd(corResults_path)

# chose specific metrics to run the analysis with
metrics <- c('cr_last', 'cr_mean', 'cr_median', 'cr_mode', 'cr_min', 'cr_max',  
             'cr_std', 'cr_cv', 'cr_range', 'cr_iqr', 'cr_skew', 'cr_kurtosis',
             'cr_auc')

# get number of participants
n <- length(unique(data$sj_id))

# create summary matrices to save cor and cocor results
cor_res <- setNames(vector("list", length(dimensions_names)), dimensions_names)
cocor_res <- setNames(vector("list", length(dimensions_names)), dimensions_names)

# loop over dimensions
for (d in 1:ncol(sr)){
  
  # initialize cor_res[[d]] as dict with metrics as keys
  cor_res[[d]] <- setNames(vector("list", length(metrics)), metrics)
  # initialize cocor_res[[d]] as dict with metrics as keys
  cocor_res[[d]] <- setNames(vector("list", length(metrics)), metrics)
  
  # vector of metric and dimension names
  metrics_d <- paste(metrics, dimensions[d], sep="_")
  curSr <-  names(sr[d])
  
  # Loop over pairs of metrics
  for (m1 in 1:length(metrics_d)){
    
    # create results folder
    curDir <- file.path(corResults_path, metrics[m1])
    if (!file.exists(curDir)) {dir.create(curDir)}
    
    # Initialize cor_res[[d]][[m1]] as dict with r and pval as keys
    cor_res[[d]][[m1]] <- list(r = NA, pval = NA)
    # Initialize cocor_res[[d]][[m1]] as dict with metrics as keys
    cocor_res[[d]][[m1]] <- setNames(vector("list", length(metrics)), metrics)
    # Get the names of the variables
    j.label <- curSr
    k.label <- metrics_d[m1]
    # Compute the Pearson correlations
    cor.jk <- cor.test(data[,metrics_d[m1]], data[,curSr], method = "pearson")
    # Get the correlation coefficients
    r.jk <- cor.jk$estimate
    #save cor.jk in txt file
    sink(file = file.path(curDir, paste0("results_",metrics_d[m1],".txt")))
    cat(paste("Results" , metrics_d[m1], "\n"))
    print(cor.jk)
    sink()
    # store results in summary matrices
    cor_res[[d]][[m1]]$r <- as.numeric(r.jk)
    cor_res[[d]][[m1]]$pval <- cor.jk$p.value
    
    for (m2 in 1:length(metrics_d)){
      # Initialize cocor_res[[d]][[m1]][[m2]] as dict with z and pval as keys
      cocor_res[[d]][[m1]][[m2]] <- list(z = NA, pval = NA)
      
      # Get the names of the variables
      h.label <- metrics_d[m2]
      # Compute the Pearson correlations
      cor.jh <- cor.test(data[,metrics_d[m2]], data[,curSr], method = "pearson")
      cor.kh <- cor.test(data[,metrics_d[m1]], data[,metrics_d[m2]], method = "pearson")
      # Get the correlation coefficients
      r.jh <- cor.jh$estimate
      r.kh <- cor.kh$estimate
      
      # Compare dependent overlapping correlations, using cocor package
      diff.jk.jh <- cocor.dep.groups.overlap(r.jk, r.jh, r.kh, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(j.label, k.label, h.label), return.htest = TRUE)
      
      # store results in summary matrices
      cocor_res[[d]][[m1]][[m2]]$z <- as.numeric(diff.jk.jh$hittner2003$statistic)
      cocor_res[[d]][[m1]][[m2]]$pval <- as.numeric(diff.jk.jh$hittner2003$p.value)
    
    } 
  }
}

# Reformat cor_res
# Initialize empty lists to store results
result_list <- list()
# Loop through each category and extract the correlation and p-value
for (category in names(cor_res)) {
  # Initialize a temporary list for storing current category results
  temp_list <- list()
  # Loop through each metric in the category
  for (metric in names(cor_res[[category]])) {
    # Extract r and pval
    r_value <- cor_res[[category]][[metric]]$r
    p_value <- cor_res[[category]][[metric]]$pval
    
    # Store in a temporary list
    temp_list[[metric]] <- c(r = r_value, pval = p_value)
  }
  # Combine into a data frame and assign to the result list
  result_list[[category]] <- as.data.frame(do.call(rbind, temp_list))
}
# Combine all categories into one data frame
final_df <- do.call(cbind, result_list)
# Clean up column names
colnames(final_df) <- paste(rep(names(result_list), each = 2), rep(c("r", "pval"), times = length(result_list)), sep = "_")
# Stack r and pval for each category
stacked_df <- data.frame(Metric = rownames(final_df))
# Loop through each category to stack r and pval
for (category in names(result_list)) {
  stacked_df[[paste0(category, "_r")]] <- final_df[[paste0(category, "_r")]]
  stacked_df[[paste0(category, "_pval")]] <- final_df[[paste0(category, "_pval")]]
}
# Display the final stacked data frame
print(stacked_df)

# save the cor_res
write.csv(stacked_df, file.path(corResults_path, 'summary_corr_cri_sr.csv'))

# Reformat cocor_res
# Function to flatten cocor_res into a data frame
flatten_cocor_res <- function(cocor_res) {
  # Create an empty data frame to store results
  results <- data.frame()
  
  # Loop through each dimension
  for (d in names(cocor_res)) {
    for (m1 in names(cocor_res[[d]])) {
      for (m2 in names(cocor_res[[d]][[m1]])) {
        # Extract the relevant data
        row <- data.frame(
          Dimension = d,
          Metric1 = m1,
          Metric2 = m2,
          Z = cocor_res[[d]][[m1]][[m2]]$z,
          Pval = cocor_res[[d]][[m1]][[m2]]$pval
        )
        # Append to results
        results <- rbind(results, row)
      }
    }
  }
  
  return(results)
}

# Flatten cocor_res
cocor_df <- flatten_cocor_res(cocor_res)
write.csv(cocor_df, file.path(corResults_path, 'summary_cocor_results.csv'))
  
# Correlation by rating method for metrics: mean, std, skewness and kurtosis --------------------------------------------

cocorResults_path <- "E:/AffectiveVR/affectivevr/cocor_results"
if (!file.exists(cocorResults_path)) {dir.create(cocorResults_path)}
setwd(cocorResults_path)

# chose specific metrics to run the analysis with, and format names
metrics <- c('cr_mean', 'cr_std', 'cr_skew', 'cr_kurtosis')
metrics_names <- str_remove(metrics, "cr_")
metrics_names <- str_replace(metrics_names, "skew", "skewness")

# label cocor parameters
j.label <- paste0("CR_", rm_list[1])
k.label <- paste0("SR_", rm_list[1])
h.label <- paste0("CR_", rm_list[2])
m.label <- paste0("SR_", rm_list[2])
n.label <- paste0("CR_", rm_list[3])
q.label <- paste0("SR_", rm_list[3])
# create summary matrices to save cor and cocor results
cor_jk <- vector(mode = "list", length = ncol(sr))
cor_hm <- vector(mode = "list", length = ncol(sr))
cor_nq <- vector(mode = "list", length = ncol(sr))
cocor_jk_hm <- vector(mode = "list", length = ncol(sr))
cocor_jk_nq <- vector(mode = "list", length = ncol(sr))
cocor_hm_nq <- vector(mode = "list", length = ncol(sr))

# get number of participants
n <- length(unique(data$sj_id))

# loop to run the entire analysis by metric and dimension
for (d in 1:ncol(sr)){
  # create a file path for every dimension
  dDir <- file.path(cocorResults_path, dimensions_names[d])
  if (!file.exists(dDir)) {dir.create(dDir)}
  
  # vector of metric and dimension names
  metrics_d <- c(paste(metrics[1], dimensions[d], sep="_"), paste(metrics[2], dimensions[d], sep="_"), paste(metrics[3], dimensions[d], sep="_"), paste(metrics[4], dimensions[d], sep="_"))
  curSr <-  names(sr[d])
  for (m in 1:length(metrics_d)){
    # create results folder for every metric 
    curDir <- file.path(dDir, metrics_d[m])
    if (!file.exists(curDir)) {dir.create(curDir)}
    
    # subset data 
    curDim <- data %>% 
      dplyr::select(all_of(curSr),metrics_d[m], rating_method) 
    # separate data per rating method
    rm1 <- curDim %>% 
      filter(rating_method == rm_list[1])
    rm2 <- curDim %>% 
      filter(rating_method == rm_list[2])
    rm3 <- curDim %>% 
      filter(rating_method == rm_list[3])
    
    # Compute CRi-SR correlations and all combinations of variables
    r.jk <- cor(rm1[,metrics_d[m]], rm1[,curSr], use = "complete.obs")[1]
    r.hm <- cor(rm2[,metrics_d[m]], rm2[,curSr], use = "complete.obs")[1]
    r.nq <- cor(rm3[,metrics_d[m]], rm3[,curSr], use = "complete.obs")[1]
    
    r.jh <- cor(rm1[,metrics_d[m]], rm2[,metrics_d[m]], use = "complete.obs")[1]
    r.jm <- cor(rm1[,metrics_d[m]], rm2[,curSr], use = "complete.obs")[1]
    r.kh <- cor(rm1[,curSr], rm2[,metrics_d[m]], use = "complete.obs")[1]
    r.km <- cor(rm1[,curSr], rm2[,curSr], use = "complete.obs")[1]
    
    r.jn <- cor(rm1[,metrics_d[m]], rm3[,metrics_d[m]], use = "complete.obs")[1]
    r.jq <- cor(rm1[,metrics_d[m]], rm3[,curSr], use = "complete.obs")[1]
    r.kn <- cor(rm1[,curSr], rm3[,metrics_d[m]], use = "complete.obs")[1]
    r.kq <- cor(rm1[,curSr], rm3[,curSr], use = "complete.obs")[1]
    
    r.hn <- cor(rm2[,metrics_d[m]], rm3[,metrics_d[m]], use = "complete.obs")[1]
    r.hq <- cor(rm2[,metrics_d[m]], rm3[,curSr], use = "complete.obs")[1]
    r.mn <- cor(rm2[,curSr], rm3[,metrics_d[m]], use = "complete.obs")[1]
    r.mq <- cor(rm2[,curSr], rm3[,curSr], use = "complete.obs")[1]
  
    # Compare dependent non-overlapping correlations, using cocor package
    diff.jk.hm <- cocor.dep.groups.nonoverlap(r.jk, r.hm, r.jh, r.jm, r.kh, r.km, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(j.label, k.label, h.label, m.label), return.htest = TRUE)
    diff.jk.nq <- cocor.dep.groups.nonoverlap(r.jk, r.nq, r.jn, r.jq, r.kn, r.kq, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(j.label, k.label, n.label, q.label), return.htest = TRUE)
    diff.hm.nq <- cocor.dep.groups.nonoverlap(r.hm, r.nq, r.hn, r.hq, r.mn, r.mq, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(h.label, m.label, n.label, q.label), return.htest = TRUE)
    
    # store results in summary matrices
    cor_jk[[d]][m] <- round(r.jk, digits = 3)
    cor_hm[[d]][m] <- round(r.hm, digits = 3)
    cor_nq[[d]][m] <- round(r.nq, digits = 3)
    cocor_jk_hm[[d]][m] <- round(diff.jk.hm$hittner2003$p.value, digits = 3)
    cocor_jk_nq[[d]][m] <- round(diff.jk.nq$hittner2003$p.value, digits = 3)
    cocor_hm_nq[[d]][m] <- round(diff.hm.nq$hittner2003$p.value, digits = 3)
    
    # save cocor results in txt file (nicely formatted -> return.htest = FALSE)
    sink(file = file.path(curDir, paste0("cocor_results_", metrics_d[m],".txt")))
    cat(c(rm_list[1], "vs.", rm_list[2], "\n"))
    print(cocor.dep.groups.nonoverlap(r.jk, r.hm, r.jh, r.jm, r.kh, r.km, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(j.label, k.label, h.label, m.label), return.htest = FALSE))
    cat('\n')
    cat(c(rm_list[1], "vs.", rm_list[3], "\n"))
    print(cocor.dep.groups.nonoverlap(r.jk, r.nq, r.jn, r.jq, r.kn, r.kq, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(j.label, k.label, n.label, q.label), return.htest = FALSE))
    cat('\n')
    cat(c(rm_list[2], "vs.", rm_list[3], "\n"))
    print(cocor.dep.groups.nonoverlap(r.hm, r.nq, r.hn, r.hq, r.mn, r.mq, n, alternative = "two.sided", test = "all", alpha = 0.05, conf.level = 0.95, null.value = 0, var.labels = c(h.label, m.label, n.label, q.label), return.htest = FALSE))
    cat('\n')
    sink()
 }
}

### Create dataframe summarizing all cor and cocor results
# create dfs of each summary matrix
cor_jk <- data.frame(cor_jk, row.names = metrics_names)
names(cor_jk) <- dimensions_names
cor_hm <- data.frame(cor_hm, row.names = metrics_names)
names(cor_hm) <- dimensions_names
cor_nq <- data.frame(cor_nq, row.names = metrics_names)
names(cor_nq) <- dimensions_names
cocor_jk_hm <- data.frame(cocor_jk_hm, row.names = metrics_names)
names(cocor_jk_hm) <- dimensions_names
cocor_jk_nq <- data.frame(cocor_jk_nq, row.names = metrics_names)
names(cocor_jk_nq) <- dimensions_names
cocor_hm_nq <- data.frame(cocor_hm_nq, row.names = metrics_names)
names(cocor_hm_nq) <- dimensions_names

# create df with all cocor results
cocor_res <- list(cocor_jk_hm,cocor_jk_nq,cocor_hm_nq)
cocor_list <- map_dfr(cocor_res,function(i){
  c(unlist(i[1]),unlist(i[2]),unlist(i[3]), unlist(i[4]))})
cocor_list <- mutate_all(cocor_list, ~ paste("p =",.)) 
row_names <- c(paste(rm_list[1],"vs.",rm_list[2]), paste(rm_list[1],"vs.",rm_list[3]), paste(rm_list[2],"vs.",rm_list[3]))
cocor_df <- data.frame(cocor_list, row.names = row_names)

# create df with all cor results
cor_res <- list(cor_jk, cor_hm, cor_nq)
cor_list <- map_dfr(cor_res,function(i){
  c(unlist(i[1]),unlist(i[2]),unlist(i[3]), unlist(i[4]))})
row_names <- c(paste0("r_", rm_list[1]), paste0("r_", rm_list[2]), paste0("r_", rm_list[3]))
cor_df <- data.frame(cor_list, row.names = row_names)

# merge cor and cocor df in one summary df
summary_df <- rbind(cor_df,cocor_df)
cri_row <- rep(metrics_names, length(dimensions))
summary_df <- rbind(cri_row, summary_df)
row.names(summary_df)[1] <- c('CRi')
header <- c(dimensions_names[1],'','','',dimensions_names[2],'','','',dimensions_names[3],'','','',dimensions_names[4], '','','')
names(summary_df) <- header

# save the summary df
write.csv(summary_df, file.path(cocorResults_path, 'summary_cocor_cor_results.csv'))

# ICC - Test-retest reliability for CR between the 2 testing sites --------

iccResults_path <- "E:/AffectiveVR/affectivevr/icc_results/"
if (!file.exists(iccResults_path)) {dir.create(iccResults_path)}
setwd(iccResults_path)
file_sr <- file.path(iccResults_path, "sr_icc_results.txt")
if (file.exists(file_sr)) {file.remove(file_sr)}
file_cr <- file.path(iccResults_path, "cr_icc_results.txt")
if (file.exists(file_cr)) {file.remove(file_cr)}

digit <- 2 
dataSummary <- data  %>% 
  group_by(quadrant, rating_method, test_site) %>% #!grouped by test_site for following plot
  pivot_longer(c(sr_v:cr_cp_angle)) %>% 
  group_by(rating_method , quadrant, test_site, name) %>% 
  dplyr::summarise(mean = mean(value, na.rm = T)) %>% 
  dplyr::rename(var = name) %>% 
  pivot_longer(mean) %>% 
  unite(name, c(var, name)) %>% 
  pivot_wider(names_from = name, values_from = value) %>% 
  mutate_if(is.numeric, ~round(., digit))

site1_sr <- dataSummary %>% filter(test_site==test_list[1])  %>%
  ungroup() %>% 
  dplyr::select(starts_with('sr')) 
site2_sr <- dataSummary %>% filter(test_site==test_list[2])  %>%
  ungroup() %>%
  dplyr::select(starts_with('sr'))
site1_cr <- dataSummary %>% filter(test_site==test_list[1])  %>%
  ungroup() %>% 
  dplyr::select(starts_with('cr_mean')) 
site2_cr <- dataSummary %>% filter(test_site==test_list[2])  %>%
  ungroup() %>%
  dplyr::select(starts_with('cr_mean'))

for (d in 1:length(dimensions)){
  curDim_site1_sr <- site1_sr %>% 
    dplyr::select(ends_with(paste0(dimensions[d], '_mean')))
  curDim_site2_sr <- site2_sr %>% 
    dplyr::select(ends_with(paste0(dimensions[d], '_mean')))
  curDim_sr_icc <-data.frame('site1' = curDim_site1_sr, 'site2' = curDim_site2_sr)
  # save icc results in txt file
  sink(file = file_sr, append = TRUE)
  cat(c(str_to_upper(dimensions_names[d]), "- Test-retest", test_list[1], "-", test_list[2], "\n"))
  print(icc(curDim_sr_icc, model = "twoway", type = "consistency", unit = "average"))
  cat('\n')
  sink()
  
  curDim_site1_cr <- site1_cr %>% 
    dplyr::select(ends_with(paste0(dimensions[d], '_mean')))
  curDim_site2_cr <- site2_cr %>% 
    dplyr::select(ends_with(paste0(dimensions[d], '_mean')))
  curDim_cr_icc <-data.frame('site1' = curDim_site1_cr, 'site2' = curDim_site2_cr)
  # save icc results in txt file
  sink(file = file_cr, append = TRUE)
  cat(c(str_to_upper(dimensions_names[d]), "- Test-retest", test_list[1], "-", test_list[2], "\n"))
  print(icc(curDim_cr_icc, model = "twoway", type = "consistency", unit = "average"))
  cat('\n')
  sink()
  
}

  



# CRi comparison between RMs ----------------------------------------------

# select only CRi data
cri <- data %>% dplyr::select(!starts_with("sr"))

# Remove Baseline data
cri <- cri %>% filter(rating_method != base_name)

# chose specific metrics to run the analysis with, and format names
metrics <- c('cr_mean', 'cr_std', 'cr_skew', 'cr_kurtosis')

# number of decimal places for rounding numerical values
digit <- 3

# save path
criResults_path <- "E:/AffectiveVR/Phase_1/affectivevr/cri_comp_results/"
if (!file.exists(criResults_path)) {dir.create(criResults_path)}
setwd(criResults_path)

# IN PROGRESS: ANOVA with factors testing site and RMs
cri.anova <- anova_test(cri, dv = response, wid = sj_id, between = test_site, effect.size = "pes")
get_anova_table(test_site.anova)

# Paired t-tests and equivalence testing between each RM (Grid, Flubber, Proprioceptive) 
# for each dimension (Arousal, Valence, Distance, Angle) - TOST approach

# create summary matrix for TOST results
summary_tost <- vector(mode = "list", length = length(dimensions))

# Set the smallest effect size of interest (SESOI)
# Here we choose 0.125 point (raw effect size) on a [-1,1] rating for arousal and valence
# For distance: [0, -sqrt(2)] ratings -> 0.125*sqrt(2) = 0.177
# For angle: [-pi, pi] ratings -> pi*0.125 = 0.393
# -> ~6.25% change
sesoi_a_v <- c(-0.125,0.125)
sesoi_dist <- c(0,0.177)
sesoi_ang <- c(-0.393,0.393)
sesoi_type <- c('raw')
mu <- 0 # difference in means tested
alpha <- 0.05 # alpha level

# loop to run the entire analysis by metric and dimension
for (m in 1:length(metrics)){
  # create a file path for every metric
  mDir <- file.path(criResults_path, metrics[m])
  if (!file.exists(mDir)) {dir.create(mDir)}
  
  for (d in 1:length(dimensions)){
    # create results folder for every dimension 
    curDir <- file.path(mDir, dimensions_names[d])
    if (!file.exists(curDir)) {dir.create(curDir)}
    
    # current metric and dimension
    curMetric_d_name <- paste(metrics[m], dimensions[d], sep="_")

    # subset data 
    curMetrics_d <- cri %>% 
      dplyr::select(curMetric_d_name, "sj_id", "rating_method", "quadrant") 
    
    # rename cri_dimension into cri for easier processing
    names(curMetrics_d)[1] <- "cri"
    
    # Define the y limits depending on dimension
    lim_y <- case_when(dimensions[d] == 'v' ~ c(-1,1),
                       dimensions[d] == 'a' ~ c(-1,1),
                       dimensions[d] == 'dist' ~ c(0,sqrt(2)),
                       dimensions[d] == 'angle' ~ c(-pi,pi))
    
    # plotting per rating method
    p_rm <- ggplot(curMetrics_d, aes(x = rating_method, y = cri, fill = rating_method)) +
      geom_half_violin(position = position_nudge(x = 0.2), alpha = 0.5, side = "r") +
      geom_point(aes(color = as.numeric(as.character(sj_id))), position = position_jitter(width = .1, height = 0), size = 1.5) +
      geom_boxplot(width = .1, alpha = 0.5) +
      ylim(lim_y) +
      ylab(curMetric_d_name) +
      xlab("Rating method") +
      scale_color_continuous(name = "SJ_id") +
      theme_bw(base_size = 14)
    p_rm
    ggsave(file.path(curDir, paste0(curMetric_d_name, "_rating_method.png")), p_rm, height = 6, width = 8)
    
    # plotting per rating method and video/quadrant
    p_rm_quad <- ggplot(curMetrics_d, aes(x = rating_method, y = cri, fill = rating_method)) +
      geom_half_violin(position = position_nudge(x = 0.2), alpha = 0.5, side = "r") +
      geom_point(aes(color = as.numeric(as.character(sj_id))), position = position_jitter(width = .1, height = 0), size = 1.5) +
      geom_boxplot(width = .1, alpha = 0.5) +
      facet_wrap(.~quadrant) +
      theme_bw(base_size = 14) +
      ggtitle(paste0(curMetric_d_name, " by Rating Method and Video/Quadrant")) +
      ylim(lim_y) +
      xlab("Rating method") +
      ylab(curMetric_d_name) +
      labs(color = 'Rating Method')
    p_rm_quad
    ggsave(file.path(curDir, paste0(curMetric_d_name, "_quadrant_rating_method.png")), p_rm_quad, height = 6, width = 8)
    
    # Paired t-tests & Equivalence testing - TOST approach
    
    rm1 <- (curMetrics_d %>% filter(rating_method == rm_list[1]))$cri
    rm2 <- (curMetrics_d %>% filter(rating_method == rm_list[2]))$cri
    rm3 <- (curMetrics_d %>% filter(rating_method == rm_list[3]))$cri
    
    sesoi <- case_when(dimensions_names[d] == 'Valence' ~ sesoi_a_v,
                       dimensions_names[d] == 'Arousal' ~ sesoi_a_v,
                       dimensions_names[d] == 'Distance' ~ sesoi_dist,
                       dimensions_names[d] == 'Angle' ~ sesoi_ang)
    
    rm1_rm2_test <- tsum_TOST(m1=mean(rm1), sd1=sd(rm1), n1=length(rm1),
                               m2=mean(rm2), sd2=sd(rm2),  n2=length(rm2), 
                               r12 = cor(rm1, rm2),
                               hypothesis = c("EQU"),
                               paired = T,
                               var.equal = T,
                               eqb = sesoi, 
                               mu = mu,
                               eqbound_type = sesoi_type,
                               alpha = alpha,
                               bias_correction = TRUE,
                               rm_correction = T,
                               glass = NULL,
                               smd_ci = c( "nct"))
    rm1_rm3_test <- tsum_TOST(m1=mean(rm1), sd1=sd(rm1), n1=length(rm1),
                               m2=mean(rm3), sd2=sd(rm3), n2=length(rm3), 
                               r12 = cor(rm1,rm3),
                               hypothesis = c("EQU"),
                               paired = T,
                               var.equal = T,
                               eqb = sesoi,
                               mu = mu,
                               eqbound_type = sesoi_type,
                               alpha = alpha,
                               bias_correction = TRUE,
                               rm_correction = T,
                               glass = NULL,
                               smd_ci = c( "nct"))
    rm2_rm3_test <- tsum_TOST(m1=mean(rm2), sd1=sd(rm2), n1=length(rm2),
                               m2=mean(rm3), sd2=sd(rm3), n2=length(rm3), 
                               r12 = cor(rm2, rm3),
                               hypothesis = c("EQU"),
                               paired = T,
                               var.equal = T,
                               eqb = sesoi,
                               mu = mu,
                               eqbound_type = sesoi_type,
                               alpha = alpha,
                               bias_correction = TRUE,
                               rm_correction = T,
                               glass = NULL,
                               smd_ci = c( "nct"))
    
    
    # extract NHST p-values and correct for multiple comparison
    p_values = c(rm1_rm2_test$TOST$p.value[1], rm1_rm3_test$TOST$p.value[1], rm2_rm3_test$TOST$p.value[1])
    adjusted_p_values <- round(p.adjust(p_values, 'fdr'), digit)
    
    # create dataframe of TOST results
    all_tosts <- data.frame(rm1_rm2_test$TOST, rm1_rm3_test$TOST, rm2_rm3_test$TOST) %>% 
      mutate_if(is.numeric, ~round(., digit))
    # add NHST FDR corrected p-values
    for (i in seq(3, 1, -1)) {
      new_column <- data.frame("FDR adjusted p-value" = c(adjusted_p_values[i],"",""))
      all_tosts <- add_column(all_tosts, new_column, .after = 4*i)
    }
    # add header with post-hoc comparisons
    all_tosts_names <- names(all_tosts)
    header <- c(paste(rm_list[1],'vs',rm_list[2]),'','','','',paste(rm_list[1],'vs',rm_list[3]),'','','','',paste(rm_list[2],'vs',rm_list[3]),'','','','')
    names(all_tosts) <- header
    all_tosts <- rbind(all_tosts_names, all_tosts)
    row.names(all_tosts)[1] <- dimensions_names[d]
    
    # save dataframe in summary matrix
    summary_tost[[d]] <- all_tosts
    
    # save dataframe in csv file
    write.csv(all_tosts, file = file.path(curDir, paste0("TOST_", curMetric_d_name, ".csv")))
    
    #save results in txt file
    sink(file = file.path(curDir, paste0("results_", curMetric_d_name,".txt")))
    
    cat(c("Comparison", rm_list[1], "vs.", rm_list[2], "\n"))
    print(rm1_rm2_test)
    cat(c('FDR adjusted NHST p-value: ',adjusted_p_values[1],'\n'))
    cat('\n')
    
    cat(c("Comparison", rm_list[1], "vs.", rm_list[3], "\n"))
    print(rm1_rm3_test)
    cat(c('FDR adjusted NHST p-value: ',adjusted_p_values[2],'\n'))
    cat('\n')
    
    cat(c("Comparison", rm_list[2], "vs.", rm_list[3], "\n"))
    print(rm2_rm3_test)
    cat(c('FDR adjusted NHST p-value: ',adjusted_p_values[3],'\n'))
    cat('\n')
    
    sink()
    
  }

  # Reformat summary dataframe
  names(summary_tost) <- dimensions_names
  
  # save summary dataframe in csv file
  file_path = file.path(mDir, paste0(metrics[m], '_summary_tost.csv'))
  if (file.exists(file_path)) {file.remove(file_path)} # delete file if already exist (because append otherwise)
  lapply(summary_tost, function(x) write.table(x, file = file_path, append=T, sep=',', row.names=T, col.names=T ))
}

# SR comparison between RMs including Baseline -----------------------------------------------

# number of decimal places for rounding numerical values
digit <- 3

# save path
criResults_path <- "E:/AffectiveVR/Phase_1/affectivevr/cri_comp_results/"
if (!file.exists(criResults_path)) {dir.create(criResults_path)}
setwd(criResults_path)

# Paired t-tests and equivalence testing between each RM (Grid, Flubber, Proprioceptive) 
# for each dimension (Arousal, Valence, Distance, Angle) - TOST approach

# create summary matrix for TOST results
summary_tost <- vector(mode = "list", length = length(dimensions))

# Set the smallest effect size of interest (SESOI)
# Here we choose 0.125 point (raw effect size) on a [-1,1] rating for arousal and valence
# For distance: [0, -sqrt(2)] ratings -> 0.125*sqrt(2) = 0.177
# For angle: [-pi, pi] ratings -> pi*0.125 = 0.393
# -> ~6.25% change
sesoi_a_v <- c(-0.125,0.125)
sesoi_dist <- c(0,0.177)
sesoi_ang <- c(-0.393,0.393)
sesoi_type <- c('raw')
mu <- 0 # difference in means tested
alpha <- 0.05 # alpha level

# create a file path for every metric
mDir <- file.path(criResults_path, 'sr')
if (!file.exists(mDir)) {dir.create(mDir)}

for (d in 1:length(dimensions)){
  # create results folder for every dimension 
  curDir <- file.path(mDir, dimensions_names[d])
  if (!file.exists(curDir)) {dir.create(curDir)}
  
  # current metric and dimension
  curMetric_d_name <- paste('sr', dimensions[d], sep="_")
  
  # subset data 
  curMetrics_d <- sr %>% 
    dplyr::select(curMetric_d_name, "sj_id", "rating_method", "quadrant") 
  
  # rename sr_dimension into sr for easier processing
  names(curMetrics_d)[1] <- "sr"
  
  # Define the y limits depending on dimension
  lim_y <- case_when(dimensions[d] == 'v' ~ c(-1,1),
                     dimensions[d] == 'a' ~ c(-1,1),
                     dimensions[d] == 'dist' ~ c(0,sqrt(2)),
                     dimensions[d] == 'angle' ~ c(-pi,pi))
  
  # plotting per rating method
  p_rm <- ggplot(curMetrics_d, aes(x = rating_method, y = sr, fill = rating_method)) +
    geom_half_violin(position = position_nudge(x = 0.2), alpha = 0.5, side = "r") +
    geom_point(aes(color = as.numeric(as.character(sj_id))), position = position_jitter(width = .1, height = 0), size = 1.5) +
    geom_boxplot(width = .1, alpha = 0.5) +
    ylim(lim_y) +
    ylab(curMetric_d_name) +
    xlab("Rating method") +
    scale_color_continuous(name = "SJ_id") +
    theme_bw(base_size = 14)
  p_rm
  ggsave(file.path(curDir, paste0(curMetric_d_name, "_rating_method.png")), p_rm, height = 6, width = 8)
  
  # plotting per rating method and video/quadrant
  p_rm_quad <- ggplot(curMetrics_d, aes(x = rating_method, y = sr, fill = rating_method)) +
    geom_half_violin(position = position_nudge(x = 0.2), alpha = 0.5, side = "r") +
    geom_point(aes(color = as.numeric(as.character(sj_id))), position = position_jitter(width = .1, height = 0), size = 1.5) +
    geom_boxplot(width = .1, alpha = 0.5) +
    facet_wrap(.~quadrant) +
    theme_bw(base_size = 14) +
    ggtitle(paste0(curMetric_d_name, " by Rating Method and Video/Quadrant")) +
    ylim(lim_y) +
    xlab("Rating method") +
    ylab(curMetric_d_name) +
    labs(color = 'Rating Method')
  p_rm_quad
  ggsave(file.path(curDir, paste0(curMetric_d_name, "_quadrant_rating_method.png")), p_rm_quad, height = 6, width = 8)
  
  # Paired t-tests & Equivalence testing - TOST approach
  
  rm1 <- (curMetrics_d %>% filter(rating_method == rm_list[1]))$sr
  rm2 <- (curMetrics_d %>% filter(rating_method == rm_list[2]))$sr
  rm3 <- (curMetrics_d %>% filter(rating_method == rm_list[3]))$sr
  baseline <- (curMetrics_d %>% filter(rating_method == base_name))$sr
  
  sesoi <- case_when(dimensions_names[d] == 'Valence' ~ sesoi_a_v,
                     dimensions_names[d] == 'Arousal' ~ sesoi_a_v,
                     dimensions_names[d] == 'Distance' ~ sesoi_dist,
                     dimensions_names[d] == 'Angle' ~ sesoi_ang)
  
  rm1_rm2_test <- tsum_TOST(m1=mean(rm1), sd1=sd(rm1), n1=length(rm1),
                            m2=mean(rm2), sd2=sd(rm2),  n2=length(rm2), 
                            r12 = cor(rm1, rm2),
                            hypothesis = c("EQU"),
                            paired = T,
                            var.equal = T,
                            eqb = sesoi, 
                            mu = mu,
                            eqbound_type = sesoi_type,
                            alpha = alpha,
                            bias_correction = TRUE,
                            rm_correction = T,
                            glass = NULL,
                            smd_ci = c( "nct"))
  rm1_rm3_test <- tsum_TOST(m1=mean(rm1), sd1=sd(rm1), n1=length(rm1),
                            m2=mean(rm3), sd2=sd(rm3), n2=length(rm3), 
                            r12 = cor(rm1,rm3),
                            hypothesis = c("EQU"),
                            paired = T,
                            var.equal = T,
                            eqb = sesoi,
                            mu = mu,
                            eqbound_type = sesoi_type,
                            alpha = alpha,
                            bias_correction = TRUE,
                            rm_correction = T,
                            glass = NULL,
                            smd_ci = c( "nct"))
  rm2_rm3_test <- tsum_TOST(m1=mean(rm2), sd1=sd(rm2), n1=length(rm2),
                            m2=mean(rm3), sd2=sd(rm3), n2=length(rm3), 
                            r12 = cor(rm2, rm3),
                            hypothesis = c("EQU"),
                            paired = T,
                            var.equal = T,
                            eqb = sesoi,
                            mu = mu,
                            eqbound_type = sesoi_type,
                            alpha = alpha,
                            bias_correction = TRUE,
                            rm_correction = T,
                            glass = NULL,
                            smd_ci = c( "nct"))
  rm1_base_test <- tsum_TOST(m1=mean(rm1), sd1=sd(rm1), n1=length(rm1),
                                 m2=mean(baseline), sd2=sd(baseline), n2=length(baseline), 
                                 r12 = cor(rm1, baseline),
                                 hypothesis = c("EQU"),
                                 paired = T,
                                 var.equal = T,
                                 eqb = sesoi,
                                 mu = mu,
                                 eqbound_type = sesoi_type,
                                 alpha = alpha,
                                 bias_correction = TRUE,
                                 rm_correction = T,
                                 glass = NULL,
                                 smd_ci = c( "nct"))
  rm2_base_test <- tsum_TOST(m1=mean(rm2), sd1=sd(rm2), n1=length(rm2),
                                 m2=mean(baseline), sd2=sd(baseline), n2=length(baseline), 
                                 r12 = cor(rm2, baseline),
                                 hypothesis = c("EQU"),
                                 paired = T,
                                 var.equal = T,
                                 eqb = sesoi,
                                 mu = mu,
                                 eqbound_type = sesoi_type,
                                 alpha = alpha,
                                 bias_correction = TRUE,
                                 rm_correction = T,
                                 glass = NULL,
                                 smd_ci = c( "nct"))
  rm3_base_test <- tsum_TOST(m1=mean(rm3), sd1=sd(rm3), n1=length(rm3),
                                 m2=mean(baseline), sd2=sd(baseline), n2=length(baseline), 
                                 r12 = cor(rm3, baseline),
                                 hypothesis = c("EQU"),
                                 paired = T,
                                 var.equal = T,
                                 eqb = sesoi,
                                 mu = mu,
                                 eqbound_type = sesoi_type,
                                 alpha = alpha,
                                 bias_correction = TRUE,
                                 rm_correction = T,
                                 glass = NULL,
                                 smd_ci = c( "nct"))
  
  
  # extract NHST p-values and correct for multiple comparison
  p_values = c(rm1_base_test$TOST$p.value[1], rm2_base_test$TOST$p.value[1], rm3_base_test$TOST$p.value[1], rm1_rm2_test$TOST$p.value[1], rm1_rm3_test$TOST$p.value[1], rm2_rm3_test$TOST$p.value[1])
  adjusted_p_values <- round(p.adjust(p_values, 'fdr'), digit)
  
  # create dataframe of TOST results
  all_tosts <- data.frame(rm1_base_test$TOST, rm2_base_test$TOST, rm3_base_test$TOST, rm1_rm2_test$TOST, rm1_rm3_test$TOST, rm2_rm3_test$TOST) %>% 
    mutate_if(is.numeric, ~round(., digit))
  # add NHST FDR corrected p-values
  for (i in seq(6, 1, -1)) {
    new_column <- data.frame("FDR adjusted p-value" = c(adjusted_p_values[i],"",""))
    all_tosts <- add_column(all_tosts, new_column, .after = 4*i)
  }
  # add header with post-hoc comparisons
  all_tosts_names <- names(all_tosts)
  header <- c(paste(rm_list[1],'vs',base_name),'','','','',paste(rm_list[2],'vs',base_name),'','','','',paste(rm_list[3],'vs',base_name),'','','','',paste(rm_list[1],'vs',rm_list[2]),'','','','',paste(rm_list[1],'vs',rm_list[3]),'','','','',paste(rm_list[2],'vs',rm_list[3]),'','','','')
  names(all_tosts) <- header
  all_tosts <- rbind(all_tosts_names, all_tosts)
  row.names(all_tosts)[1] <- dimensions_names[d]
  
  # save dataframe in summary matrix
  summary_tost[[d]] <- all_tosts
  
  # save dataframe in csv file
  write.csv(all_tosts, file = file.path(curDir, paste0("TOST_", curMetric_d_name, ".csv")))
  
  #save results in txt file
  sink(file = file.path(curDir, paste0("results_", curMetric_d_name,".txt")))
  
  cat(c("Comparison", rm_list[1], "vs.", rm_list[2], "\n"))
  print(rm1_rm2_test)
  cat(c('FDR adjusted NHST p-value: ',adjusted_p_values[1],'\n'))
  cat('\n')
  
  cat(c("Comparison", rm_list[1], "vs.", rm_list[3], "\n"))
  print(rm1_rm3_test)
  cat(c('FDR adjusted NHST p-value: ',adjusted_p_values[2],'\n'))
  cat('\n')
  
  cat(c("Comparison", rm_list[2], "vs.", rm_list[3], "\n"))
  print(rm2_rm3_test)
  cat(c('FDR adjusted NHST p-value: ',adjusted_p_values[3],'\n'))
  cat('\n')
  
  sink()
  
}

# Reformat summary dataframe
names(summary_tost) <- dimensions_names

# save summary dataframe in csv file
file_path = file.path(mDir, 'sr_summary_tost.csv')
if (file.exists(file_path)) {file.remove(file_path)} # delete file if already exist (because append otherwise)
lapply(summary_tost, function(x) write.table(x, file = file_path, append=T, sep=',', row.names=T, col.names=T ))


