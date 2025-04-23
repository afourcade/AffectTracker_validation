# Script for analysis of questionnaires in Assessment phase
# note: for invasiveness questionnaire -> separate script: invasiveness.R
#
# Author: Antonin Fourcade
# Last version: 15.08.2023

# import packages
library(tidyverse)
library(ggpmisc)
library(ggdist)
library(gghalves)
library(hrbrthemes)
library(viridis)
library(car)
library(readxl)
library(rstatix)
library(TOSTER)

source("E:/AffectiveVR/Phase_1/affectivevr/util.R")

# set paths
data_path <- "E:/AffectiveVR/Phase_1/Data/AVR/"

# read assessment data
data_assess <- read.csv(paste(data_path, "assessment.csv", sep = ""))

# set names of test sites
test_list <- c("Torino", "Berlin")
# set names of rating methods and baseline
rm_list <- c("Grid", "Flubber", "Proprioceptive")
base_name <- c("Baseline")
# set names of questionnaires
q_list <- c("SUS", "invasive_presence", "Kunin")
# number of decimal places for rounding numerical values
digit <- 2

# format data
data_assess <- data_assess %>% 
  mutate(
    test_site = factor(test_site, levels = test_list),
    rating_method = factor(rating_method, levels = c(rm_list, base_name)),
    questionnaire = factor(questionnaire, levels = q_list),
    sj_id = factor(sj_id, levels = unique(data_assess$sj_id)))

# add sj_idx column
data_assess$sj_idx <- as.integer(data_assess$sj_id)

# Satisfaction (Kunin scale) ----------------------------------------------

# select data
satisfactq <- filter(data_assess, questionnaire=="Kunin")

# path for saving results
Results_path_satisfact <- "E:/AffectiveVR/Phase_1/affectivevr/assessment_results/satisfaction/"

# plotting per rating method
p_satisfact <- plot_rm_q(df=satisfactq, x=rating_method, y=response, fill=rating_method, xlab="Rating Method", ylab="Satisfaction", breaks_y=seq(0,6,1), limits_y = c(0,6), savepath=file.path(Results_path_satisfact, "satisfaction_rating_method.png"))
p_satisfact

# Compare ratings methods

# Set the smallest effect size of interest (SESOI)
# Here we choose 0.5 point (raw effect size) on a 7-point Likert scale
# -> ~7% change
sesoi <- c(-0.5,0.5) #c(-0.4,0.4)
sesoi_type <- c('raw')
mu <- 0 # difference in means tested
alpha <- 0.05 # alpha level
mc_cor <- "fdr"

tost_satisfact <- compare_rm_q(df=satisfactq, dv='response', rm_list=rm_list, base_name=base_name, sesoi=sesoi, sesoi_type=sesoi_type, mu=mu, alpha=alpha, mc_cor=mc_cor, digit=digit, savepath=Results_path_satisfact, save_suffix="satisfaction")

### Additional post-hoc paired t-tests
# Paired t-tests between each Feedback (Grid, Flubber, Proprioceptive) and Baseline
add_p_ttests_rm(df=satisfactq, dv='response', rm_list=rm_list, base_name=base_name, mc_cor=mc_cor, savepath=Results_path_satisfact)

### Test site comparison
# ANOVA with factor test_site and 2 levels (Torino, Berlin)
test_site.anova <- anova_test(satisfactq, dv = response, wid = sj_id, between = test_site, effect.size = "pes")
get_anova_table(test_site.anova)
# Plotting per test site
plot_rm_q(df=satisfactq, x=test_site, y=response, fill=test_site, xlab="Test_site", ylab="Satisfaction", breaks_y=seq(0,6,1), savepath=file.path(Results_path_satisfact, "satisfaction_testing_site.png"))
# Tost approach
# Extract data for each test site
test_site1 <- (satisfactq %>% filter(test_site == test_list[1]))$response
test_site2 <- (satisfactq %>% filter(test_site == test_list[2]))$response
#Paired t-tests & Equivalence testing
test_site1_test_site2_test <- tsum_TOST(m1=mean(test_site1), sd1=sd(test_site1), n1=length(test_site1),
                                        m2=mean(test_site2), sd2=sd(test_site2),  n2=length(test_site2), 
                                        hypothesis = c("EQU"),
                                        paired = F,
                                        var.equal = T,
                                        eqb = sesoi, 
                                        mu = mu,
                                        eqbound_type = sesoi_type,
                                        alpha = alpha,
                                        bias_correction = TRUE,
                                        rm_correction = T,
                                        glass = NULL,
                                        smd_ci = c( "nct"),
)
# save dataframe in csv file
write.csv(data.frame(test_site1_test_site2_test$TOST), file = file.path(Results_path_satisfact, "TOST_test_site.csv"))
#save results in txt file
sink(file = file.path(Results_path_satisfact, "results_compare_test_site.txt"))
cat(c("ANOVA with factor test_site and two levels:", test_list, "\n"))
print(get_anova_table(test_site.anova))
cat("\n")
sink()

# System Usability Scale (SUS) --------------------------------------------

# select data
susq <- filter(data_assess, questionnaire=="SUS")

# path for saving results
Results_path_sus <- "E:/AffectiveVR/Phase_1/affectivevr/assessment_results/SUS/"

# compute SUS score
# For each of the even numbered questions (positive), subtract 0 from the score.
# For each of the odd numbered questions (negative), subtract their value from 6.
# Take these new values and add up the total score. Then multiply this by 100/(7*6)=2.38
odd_q <- susq %>% filter(index %% 2 == 1) %>% group_by(sj_id, sj_idx, rating_method, test_site) %>% 
  dplyr::summarize(sus_score = sum(6-response))
even_q <- susq %>% filter(index %% 2 == 0) %>% group_by(sj_id, sj_idx, rating_method, test_site)%>% 
  dplyr::summarize(sus_score = sum(response))
sus_score <- data.frame(sj_id = odd_q$sj_id,
                        sj_idx = odd_q$sj_idx,
                        rating_method = odd_q$rating_method,
                        test_site = odd_q$test_site,
                        sus_score = (odd_q$sus_score+even_q$sus_score)*2.38)

# plotting per rating method
p_sus <- plot_rm_q(df=sus_score, x=rating_method, y=sus_score, fill=rating_method, xlab="Rating Method", ylab="SUS Score", breaks_y=seq(0,100,10), limits_y = c(0,100), savepath=file.path(Results_path_sus, "sus_rating_method.png"))
p_sus

# Compare ratings methods

# Set the smallest effect size of interest (SESOI)
# Here we choose 7 point (raw effect size) on a 0-100 scale
# -> ~7% change
sesoi <- c(-7,7) 
sesoi_type <- c('raw')
mu <- 0 # difference in means tested
alpha <- 0.05 # alpha level
mc_cor <- "fdr"

tost_sus <- compare_rm_q(df=sus_score, dv='sus_score', rm_list=rm_list, base_name=base_name, sesoi=sesoi, sesoi_type=sesoi_type, mu=mu, alpha=alpha, mc_cor=mc_cor, digit=digit, savepath=Results_path_sus, save_suffix="sus")

### Test site comparison
# ANOVA with factor test_site and 2 levels (Torino, Berlin)
test_site.anova <- anova_test(sus_score, dv = sus_score, wid = sj_id, between = test_site, effect.size = "pes")
get_anova_table(test_site.anova)
# Plotting per test site
plot_rm_q(df=sus_score, x=test_site, y=sus_score, fill=test_site, xlab="Test_site", ylab="SUS Score", breaks_y=seq(0,100,10), limits_y = c(0,100), savepath=file.path(Results_path_sus, "sus_testing_site.png"))
# Tost approach
# Extract data for each test site
test_site1 <- (sus_score %>% filter(test_site == test_list[1]))$sus_score
test_site2 <- (sus_score %>% filter(test_site == test_list[2]))$sus_score
#Paired t-tests & Equivalence testing
test_site1_test_site2_test <- tsum_TOST(m1=mean(test_site1), sd1=sd(test_site1), n1=length(test_site1),
                                        m2=mean(test_site2), sd2=sd(test_site2),  n2=length(test_site2), 
                                        hypothesis = c("EQU"),
                                        paired = F,
                                        var.equal = T,
                                        eqb = sesoi, 
                                        mu = mu,
                                        eqbound_type = sesoi_type,
                                        alpha = alpha,
                                        bias_correction = TRUE,
                                        rm_correction = T,
                                        glass = NULL,
                                        smd_ci = c( "nct"),
)
# save dataframe in csv file
write.csv(data.frame(test_site1_test_site2_test$TOST), file = file.path(Results_path_sus, "TOST_test_site.csv"))
#save results in txt file
sink(file = file.path(Results_path_sus, "results_compare_test_site.txt"))
cat(c("ANOVA with factor test_site and two levels:", test_list, "\n"))
print(get_anova_table(test_site.anova))
cat("\n")
sink()

# Presence and Representation (of inner emotions) ---------------------------------------------

# select data
presq <- filter(data_assess, questionnaire=="invasive_presence" & (index == 3 | index == 4))
repq <- filter(data_assess, questionnaire=="invasive_presence" & index == 2)

# path for saving results
Results_path_presrep <- "E:/AffectiveVR/Phase_1/affectivevr/assessment_results/presence_representation/"

# Compute Presence score
# reformat 'outside world' (negative question)
presq['response'][presq['index']==4]  <- 6 - presq['response'][presq['index']==4]
# Mean of the 2 questions
pres_score <- presq %>% group_by(sj_id, sj_idx, rating_method, test_site) %>% 
  dplyr::summarize(pres_score = mean(response))
pres_score <- as.data.frame(pres_score)

# plotting presence per rating method
p_pres <- plot_rm_q(df=pres_score, x=rating_method, y=pres_score, fill=rating_method, xlab="Rating Method", ylab="Presence Score", breaks_y=seq(0,6,1), limits_y = c(0,6), savepath=file.path(Results_path_presrep, "presence_rating_method.png"))
p_pres

# plotting representation per rating method
p_rep <- plot_rm_q(df=repq, x=rating_method, y=response, fill=rating_method, xlab="Rating Method", ylab="Representation inner emotions", breaks_y=seq(0,6,1), limits_y = c(0,6), savepath=file.path(Results_path_presrep, "representation_rating_method.png"))
p_rep

# Compare ratings methods for presence
# Set the smallest effect size of interest (SESOI)
# Here we choose 0.5 point (raw effect size) on a 7-point Likert scale
# -> ~7% change
sesoi <- c(-0.5,0.5)
sesoi_type <- c('raw')
mu <- 0 # difference in means tested
alpha <- 0.05 # alpha level
mc_cor <- "fdr"

tost_pres <- compare_rm_q(df=pres_score, dv='pres_score', rm_list=rm_list, base_name=base_name, sesoi=sesoi, sesoi_type=sesoi_type, mu=mu, alpha=alpha, mc_cor=mc_cor, digit=digit, savepath=Results_path_presrep, save_suffix="presence")

# Compare ratings methods for representation
# Set the smallest effect size of interest (SESOI)
# Here we choose 0.5 point (raw effect size) on a 7-point Likert scale
# -> ~7% change
sesoi <- c(-0.5,0.5) #c(-0.4,0.4)
sesoi_type <- c('raw')
mu <- 0 # difference in means tested
alpha <- 0.05 # alpha level
mc_cor <- "fdr"

tost_rep <- compare_rm_q(df=repq, dv='response', rm_list=rm_list, base_name=base_name, sesoi=sesoi, sesoi_type=sesoi_type, mu=mu, alpha=alpha, mc_cor=mc_cor, digit=digit, savepath=Results_path_presrep, save_suffix="representation")

### Test site comparison presence
# ANOVA with factor test_site and 2 levels (Torino, Berlin)
test_site_pres.anova <- anova_test(pres_score, dv = pres_score, wid = sj_id, between = test_site, effect.size = "pes")
get_anova_table(test_site_pres.anova)
# Plotting per test site
plot_rm_q(df=pres_score, x=test_site, y=pres_score, fill=test_site, xlab="Test_site", ylab="Presence Score", breaks_y=seq(0,6,1), savepath=file.path(Results_path_presrep, "presence_testing_site.png"))
# Tost approach
# Extract data for each test site
test_site1 <- (pres_score %>% filter(test_site == test_list[1]))$pres_score
test_site2 <- (pres_score %>% filter(test_site == test_list[2]))$pres_score
#Paired t-tests & Equivalence testing
test_site1_test_site2_test_pres <- tsum_TOST(m1=mean(test_site1), sd1=sd(test_site1), n1=length(test_site1),
                                        m2=mean(test_site2), sd2=sd(test_site2),  n2=length(test_site2), 
                                        hypothesis = c("EQU"),
                                        paired = F,
                                        var.equal = T,
                                        eqb = sesoi, 
                                        mu = mu,
                                        eqbound_type = sesoi_type,
                                        alpha = alpha,
                                        bias_correction = TRUE,
                                        rm_correction = T,
                                        glass = NULL,
                                        smd_ci = c( "nct"),
)
# save dataframe in csv file
write.csv(data.frame(test_site1_test_site2_test_pres$TOST), file = file.path(Results_path_presrep, "TOST_test_site_presence.csv"))
#save results in txt file
sink(file = file.path(Results_path_presrep, "results_compare_test_site_presence.txt"))
cat(c("ANOVA with factor test_site and two levels:", test_list, "\n"))
print(get_anova_table(test_site_pres.anova))
cat("\n")
sink()

### Test site comparison representation
# ANOVA with factor test_site and 2 levels (Torino, Berlin)
test_site_repq.anova <- anova_test(repq, dv = response, wid = sj_id, between = test_site, effect.size = "pes")
get_anova_table(test_site_repq.anova)
# Plotting per test site
plot_rm_q(df=repq, x=test_site, y=response, fill=test_site, xlab="Test_site", ylab="Representation inner emotions", breaks_y=seq(0,6,1), savepath=file.path(Results_path_presrep, "representation_testing_site.png"))
# Tost approach
# Extract data for each test site
test_site1 <- (repq %>% filter(test_site == test_list[1]))$response
test_site2 <- (repq %>% filter(test_site == test_list[2]))$response
#Paired t-tests & Equivalence testing
test_site1_test_site2_test_repq <- tsum_TOST(m1=mean(test_site1), sd1=sd(test_site1), n1=length(test_site1),
                                             m2=mean(test_site2), sd2=sd(test_site2),  n2=length(test_site2), 
                                             hypothesis = c("EQU"),
                                             paired = F,
                                             var.equal = T,
                                             eqb = sesoi, 
                                             mu = mu,
                                             eqbound_type = sesoi_type,
                                             alpha = alpha,
                                             bias_correction = TRUE,
                                             rm_correction = T,
                                             glass = NULL,
                                             smd_ci = c( "nct"),
)
# save dataframe in csv file
write.csv(data.frame(test_site1_test_site2_test_repq$TOST), file = file.path(Results_path_presrep, "TOST_test_site_representation.csv"))
#save results in txt file
sink(file = file.path(Results_path_presrep, "results_compare_test_site_representation.txt"))
cat(c("ANOVA with factor test_site and two levels:", test_list, "\n"))
print(get_anova_table(test_site_repq.anova))
cat("\n")
sink()
