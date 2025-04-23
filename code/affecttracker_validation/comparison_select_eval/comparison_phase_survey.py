########################################################################################################################
# Script for comparison analysis of survey data from AVR experiment phase 1 and 2
# Need csv files survey data from phase 1 and 2
# Output: csv files (dataframe) containing ANOVA results (and post-hoc t-tests if significant effect of phase) 
# and figures of questionnaires
# Author: Antonin Fourcade
# Last version: 23.09.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
import statsmodels.api as sm
from statsmodels.formula.api import ols
from pingouin import ttest
import plotnine as plt9

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
score_list = ['invasiveness', 'presence_score', 'sus_score', 'satisfaction', 'emo_rep'] # scores of interest
# position of interest (phase 2)
pos_select = ['seated'] # ['seated', 'standing']
# feedback of interest (phase 1)
feedback_select = 'Flubber'

# path to survey csv files from phase 1 and 2 and load data
filename_phase2 = 'post_survey_preprocessed.csv'
phase2_survey_path = phase2_data_path + bids_folder + '/' + filename_phase2
phase2_survey = pd.read_csv(phase2_survey_path)
filename_phase1 = 'assessment_preprocessed.csv'
phase1_survey_path = phase1_data_path + bids_folder + '/' + filename_phase1
phase1_survey = pd.read_csv(phase1_survey_path)

# save path
if len(pos_select) == 2:
    save_path = avr_path + 'Phase_2/affectivevr/comparison_phase_survey/all_positions/'
else:
    save_path = avr_path + 'Phase_2/affectivevr/comparison_phase_survey/' + pos_select[0] + '/'
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# %%
# Phase 1 data formatting
# Select only Flubber rating_method data in phase 1
phase1_survey = phase1_survey[phase1_survey['rating_method'] == feedback_select]
# for all questionnaires but sus_score, we need to reformat responses from 0-6 to 1-7 to be comparable with phase 2
phase1_survey.loc[phase1_survey['questionnaire'] != 'sus_score', 'response'] = phase1_survey.loc[phase1_survey['questionnaire'] != 'sus_score', 'response'] + 1
# remove rating_method column
phase1_survey = phase1_survey.drop(columns=['rating_method'])
# rename sj_id to sub
phase1_survey = phase1_survey.rename(columns={'sj_id': 'sub'})
# in test_site column, replace Torino by TOR and Berlin by BER
phase1_survey.loc[phase1_survey['test_site'] == 'Torino', 'test_site'] = 'TOR'
phase1_survey.loc[phase1_survey['test_site'] == 'Berlin', 'test_site'] = 'BER'
# add phase column
phase1_survey['phase'] = 1

# %%
# Phase 2 data formatting
# select only position of interest
phase2_survey = phase2_survey[phase2_survey['position'].isin(pos_select)]
# remove position and gender columns
phase2_survey = phase2_survey.drop(columns=['position', 'gender'])
# add phase column
phase2_survey['phase'] = 2

# %%
# Merge phase 1 and 2 data
survey_data = pd.concat([phase1_survey, phase2_survey])

# %%
# ANOVA analysis with factor: phase (1, 2) and questionnaires as dependent variables
# loop over questionnaires of interest
for score in score_list:
    # select data for the score
    score_df = survey_data[survey_data['questionnaire'] == score]
    # run ANOVA with statsmodels
    model = ols('response ~ C(phase)', data=score_df).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    anova_table.index.name = score # rename index column to score
    anova_table.reset_index(inplace=True) # turn index into column
    # save results
    if len(pos_select) == 2:
        save_name = 'anova_' + score + '.csv'
    else:
        save_name = 'anova_' + score + '_' + pos_select[0] + '.csv'
    anova_table.round(2).to_csv(save_path + save_name, index=False)

    # post-hoc analysis with t-tests if significant effect of phase
    if anova_table['PR(>F)'][0] < 0.05:
        print('Significant effect of phase for ' + score)
        # run t-tests between phase 1 and 2
        phase1_score = score_df[score_df['phase'] == 1]['response']
        phase2_score = score_df[score_df['phase'] == 2]['response']
        t_test = ttest(phase1_score, phase2_score)
        t_test.index.name = score # rename index column to score
        t_test.reset_index(inplace=True) # turn index into column
        # save results
        if len(pos_select) == 2:
            save_name = 'posthoc_ttest_' + score + '.csv'
        else:
            save_name = 'posthoc_ttest_' + score + '_' + pos_select[0] + '.csv'
        t_test.round(2).to_csv(save_path + save_name, index=False)
    

# %%
# Plot each score by phase
plt9.options.figure_size = (6, 6) # set default figure size in inches
plt9.options.dpi = 300

# plot parameters
range_responses = {'invasiveness': [1,7], 'emo_rep': [1,7], 'presence_score': [1,7], 'sus_score': [0,100], 'satisfaction': [1,7]} # range of responses for each POST questionnaire
break_steps = {'invasiveness': 1, 'emo_rep': 1, 'presence_score': 1, 'sus_score': 10, 'satisfaction': 1} # step for breaks in y-axis

# for plotting, format phase column to category
survey_data['phase'] = survey_data['phase'].astype('category')

# loop through each score
for score in score_list:
    # select score
    score_df = survey_data[survey_data['questionnaire'] == score]
    # for plotting, create sub-id column by converting 'sub' to category and map categories to integers (starting from 1)
    score_df.loc[:,'sub-id'] = score_df.loc[:,'sub'].astype('category').cat.codes + 1
    # plot score distribution
    p = (
        plt9.ggplot(score_df, plt9.aes(x='phase', y='response')) 
        + plt9.geom_violin(plt9.aes(fill='phase'), position=plt9.position_nudge(x=0.2), style='right', alpha=0.5)
        # Add the boxplot layer
        + plt9.geom_boxplot(plt9.aes(fill='phase'), width=0.1, alpha=0.5)
        # Add the points layer
        + plt9.geom_point(plt9.aes(color= 'sub-id'), position=plt9.position_jitter(width=0.1, height=0), size=1.5) 
        # Set the y-axis labels and breaks
        + plt9.scale_y_continuous(name = score, breaks= np.arange(range_responses[score][0], range_responses[score][1]+1, break_steps[score]), limits = range_responses[score])
        # Set the x-axis label
        + plt9.xlab(xlab='phase')
        # Set the color scale legend
        + plt9.scale_color_continuous(name = 'sub-id')
        # Set the theme to a white background with a base font size of 14
        + plt9.theme_bw(base_size = 14)
        + plt9.theme(svg_usefonts=True)
        )
    # save plot
    if len(pos_select) == 2:
            save_name = score + '_phase.png'
    else:
            save_name = score + '_phase_' + pos_select[0] + '.png'
    p.save(filename=save_path + save_name)
        
# %%
