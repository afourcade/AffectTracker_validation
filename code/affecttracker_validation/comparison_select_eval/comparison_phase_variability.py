########################################################################################################################
# Script for comparison analysis of CR variability from AVR experiment Selection (phase 1) and Evaluation (phase 2)
# Need csv files cri_sr data from Selection and Evaluation phases
# Output: csv files (dataframe) containing ANOVA results (and post-hoc t-tests if significant effect of phase) 
# and figures of variability
# Author: Antonin Fourcade
# Last version: 12.06.2025
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
import statsmodels.api as sm
from statsmodels.formula.api import ols
from pingouin import ttest, multicomp
import plotnine as plt9
import matplotlib.pyplot as plt

# %%
# Set paths and experiment parameters
avr_path = 'E:/AffectiveVR/'
phase2_data_path = avr_path + 'Phase_2/Data/'
bids_folder = 'AVR' # folder with data in bids format
phase1_data_path = avr_path + 'Phase_1/Data/'

# experiment parameters
# blocks = ['Practice', 'Experiment']
# logging_freq = ['CR', 'SR']
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2
cri_sr_select = ['cr_mean', 'cr_std'] # cri_sr of interest
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'} # dictionnary for rating dimensions
#anova_factors = ['video'] # ['phase', 'video']; factors for ANOVA, note: these are categorical variables
# position of interest (phase 2)
pos_select = ['seated']  # ['seated', 'standing']
# quadrants/videos (phase 1)
quadrants = ["HP", "HN", "LP", "LN"]
# feedback of interest (phase 1)
feedback_select = 'Flubber'

# path to cri_sr csv files from phase 1 and 2 and load data
filename_phase2 = 'cri_sr_clean.csv'
phase2_cri_sr_path = phase2_data_path + bids_folder + '/' + filename_phase2
phase2_cri_sr = pd.read_csv(phase2_cri_sr_path)
filename_phase1 = 'cri_sr_clean.csv'
phase1_cri_sr_path = phase1_data_path + bids_folder + '/' + filename_phase1
phase1_cri_sr = pd.read_csv(phase1_cri_sr_path)

#excluded_subs = ['sub-13'] # excluded subjects from phase 2
# remove excluded subjects from phase 2 cri_sr
phase2_cri_sr = phase2_cri_sr[~phase2_cri_sr['sub'].isin(['sub-13'])]

# save path
if len(pos_select) == 2:
    save_path = 'E:/AffectiveVR/affecttracker_validation/results/comparison_selection_evaluation/variability/all_positions/'
else:
    save_path = 'E:/AffectiveVR/affecttracker_validation/results/comparison_selection_evaluation/variability/' + pos_select[0] + '/'
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# %%
# List of all cri_sr of interest
cri_sr_list = [x + '_' + rat_dim[y] for x in cri_sr_select for y in rat_dim.keys()]

# %%
# Phase 1 data formatting
# Select only data from rating_method of interest in phase 1
phase1_cri_sr = phase1_cri_sr[phase1_cri_sr['rating_method'] == feedback_select]
# keep only columns in cri_sr_list and covariates
phase1_cri_sr = phase1_cri_sr[['sj_id', 'test_site', 'quadrant' ] + cri_sr_list]
# rename sj_id to sub and quadrant to video
phase1_cri_sr = phase1_cri_sr.rename(columns={'sj_id': 'sub', 'quadrant': 'video'})
# in test_site column, replace Torino by TOR and Berlin by BER
phase1_cri_sr.loc[phase1_cri_sr['test_site'] == 'Torino', 'test_site'] = 'TOR'
phase1_cri_sr.loc[phase1_cri_sr['test_site'] == 'Berlin', 'test_site'] = 'BER'
# add phase column
phase1_cri_sr['phase'] = 1

# %%
# Phase 2 data formatting
# select only position of interest
phase2_cri_sr = phase2_cri_sr[phase2_cri_sr['position'].isin(pos_select)]
# keep only columns in cri_sr_list and covariates
phase2_cri_sr = phase2_cri_sr[['sub', 'test_site'] + cri_sr_list]
# add video column
phase2_cri_sr['video'] = 'sequence'
# add phase column
phase2_cri_sr['phase'] = 2

# %%
# Merge phase 1 and 2 data
cri_sr_data = pd.concat([phase1_cri_sr, phase2_cri_sr])

# %%
# ANOVA analysis with factor video
# loop over questionnaires of interest
for cri_sr in cri_sr_list:
    # run ANOVA with statsmodels
    model = ols(cri_sr + ' ~ C(video)', data=cri_sr_data).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    anova_table.index.name = cri_sr # rename index column to cri_sr
    anova_table.reset_index(inplace=True) # turn index into column
    # save results
    if len(pos_select) == 2:
        save_name = 'anova_' + cri_sr + '.csv'
    else:
        save_name = 'anova_' + cri_sr + '_' + pos_select[0] + '.csv'
    anova_table.round(3).to_csv(save_path + save_name, index=False)

    # post-hoc analysis with t-tests if significant effects
    # if significant effect of video, run t-tests between videos
    if anova_table['PR(>F)'][0] < 0.05:
        print('Significant effect of video for ' + cri_sr)
        # run t-tests between each pair of videos
        videos = cri_sr_data['video'].unique()
        ttest_results = pd.DataFrame()
        for i in range(len(videos)):
            for j in range(i+1, len(videos)):
                video1 = cri_sr_data[cri_sr_data['video'] == videos[i]][cri_sr]
                video2 = cri_sr_data[cri_sr_data['video'] == videos[j]][cri_sr]
                t_test = ttest(video1, video2, paired=False, alternative='two-sided') # two-sided t-test
                t_test.index.name = cri_sr # rename index column to cri_sr
                t_test.reset_index(inplace=True) # turn index into column
                t_test[cri_sr] = videos[i] + ' vs. ' + videos[j] # add t-test description
                ttest_results = pd.concat([ttest_results, t_test])
        # rename p-val column to p-unc
        ttest_results.rename(columns={'p-val': 'p-unc'}, inplace=True)
        # correct p-values for multiple comparisons
        ttest_multicomp = multicomp(ttest_results['p-unc'], alpha=0.05, method='bonferroni')
        ttest_results['p-bonf'] = ttest_multicomp[1]
        ttest_results['significance'] = ttest_multicomp[0]
        # save results
        if len(pos_select) == 2:
            save_name = 'posthoc_ttest_video_' + cri_sr + '.csv'
        else:
            save_name = 'posthoc_ttest_video_' + cri_sr + '_' + pos_select[0] + '.csv'
        ttest_results.round(3).to_csv(save_path + save_name, index=False)

    

# %%
# Plot each cri_sr by video
plt9.options.figure_size = (10, 6) # set default figure size in inches
plt9.options.dpi = 300

# Colors
colors = ["#56B4E9", "#0072B2", "#6C6C6C", "#7a86d1", "#009E73"]  # light blue, dark blue, grey, blue, green
# colors = ["#F0E442", "#E69F00", "#D55E00", "#CC79A7", "#009E73"]  # yellow, light orange, dark orange, pink, green

# dictionary for video names
video_names = {'sequence': 'Evaluation Sequence', 'HP': 'Selection HP', 'HN': 'Selection HN', 'LP': 'Selection LP', 'LN': 'Selection LN'}
# change video names in cri_sr_data
cri_sr_data.replace({'video': video_names}, inplace=True)

# loop through each cri_sr
for cri_sr in cri_sr_list:
    # for plotting, create sub-id column by converting 'sub' to category and map categories to integers (starting from 1)
    cri_sr_data.loc[:,'sub-id'] = cri_sr_data.loc[:,'sub'].astype('category').cat.codes + 1
    # plot cri_sr distribution
    p = (
        plt9.ggplot(cri_sr_data, plt9.aes(x='video', y=cri_sr)) 
        + plt9.geom_violin(plt9.aes(fill='video'), position=plt9.position_nudge(x=0.2), style='right', alpha=0.5)
        # Add the boxplot layer
        + plt9.geom_boxplot(plt9.aes(fill='video'), position=plt9.position_nudge(x=0.2), width=0.1, alpha=1, outlier_shape='')
        # Add the points layer
        + plt9.geom_point(plt9.aes(color= 'sub-id'), position=plt9.position_jitter(width=0.1, height=0), size=1.5, alpha=0.5)
        #+ plt9.geom_point(plt9.aes(color= 'sub-id'), position=plt9.position_nudge(x=-0.2), size=1.5, alpha=0.5)  
        # Set the y-axis labels and breaks
        #+ plt9.scale_y_continuous(name = cri_sr, breaks= np.arange(range_responses[score][0], range_responses[score][1]+1, break_steps[score]), limits = range_responses[score])
        # Set the x-axis label
        #+ plt9.xlab('video')
        # set the fill color scale
        + plt9.scale_fill_manual(values=colors)
        # Set the color scale legend
        + plt9.scale_color_continuous(name = 'sub-id')
        # Set the theme to a white background with a base font size of 14
        + plt9.theme_bw(base_size = 14)
        + plt9.theme(svg_usefonts=True)
        )
    # save plot
    if len(pos_select) == 2:
            save_name = cri_sr + '_video.png'
    else:
            save_name = cri_sr + '_video_' + pos_select[0] + '.svg'
    p.save(filename=save_path + save_name)


# %%
