########################################################################################################################
# Script to plot (pre- and post-experiment) survey from AVR experiment
# Need csv file containing pre_ and post_survey data (from preprocess_PRE_survey.py and preprocess_POST_survey.py)
# Output: figures of questionnaires
# Author: Antonin Fourcade
# Last version: 08.02.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
import plotnine as plt9

# %%
# Set paths and experiment parameters
phase_path = 'D:/AffectiveVR/Phase_2/'
data_path = phase_path + 'Data/'
bids_folder = 'AVR' # folder with data in bids format

# experiment parameters
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2
covariates = ['test_site', 'position', 'gender'] # covariates for ANOVA, note: these are categorical variables

# survey parameters
range_responses_pre = {'vr_exp_score': [1,4], 
                       'maia_attention_regulation': [0,5], 'maia_body_listening': [0,5], 'maia_emotional_awareness': [0,5], 'maia_not_distracting': [0,5], 'maia_not_worrying': [0,5], 'maia_noticing': [0,5], 'maia_self_regulation': [0,5], 'maia_trust': [0,5],
                       'ssq_disorientation': [1,7], 'ssq_nausea': [1,7], 'ssq_oculomotor': [1,7], 'ssq_total_score': [1,21],
                       'tas_difficulty_describing_feelings': [1,35], 'tas_difficulty_identifying_feelings': [1,25], 'tas_externally_oriented_thinking': [1,40], 'tas_total_score': [0,100]} # range of responses for each PRE questionnaire
break_step_pre = {'vr_exp_score': 1, 
                  'maia_attention_regulation': 1, 'maia_body_listening': 1, 'maia_emotional_awareness': 1, 'maia_not_distracting': 1, 'maia_not_worrying': 1, 'maia_noticing': 1, 'maia_self_regulation': 1, 'maia_trust': 1,
                  'ssq_disorientation': 1, 'ssq_nausea': 1, 'ssq_oculomotor': 1, 'ssq_total_score': 3,
                  'tas_difficulty_describing_feelings': 5, 'tas_difficulty_identifying_feelings': 5, 'tas_externally_oriented_thinking': 5, 'tas_total_score': 10} # step for breaks in y-axis
range_responses_post = {'invasiveness': [1,7], 'emo_rep': [1,7], 'presence_score': [1,7], 'sus_score': [0,100], 'satisfaction': [1,7]} # range of responses for each POST questionnaire
break_step_post = {'invasiveness': 1, 'emo_rep': 1, 'presence_score': 1, 'sus_score': 10, 'satisfaction': 1} # step for breaks in y-axis

# path to survey csv files and load data
filename_pre = 'pre_survey_preprocessed.csv'
filename_post = 'post_survey_preprocessed.csv'
pre_path = data_path + bids_folder + '/' + filename_pre
post_path = data_path + bids_folder + '/' + filename_post
pre_survey = pd.read_csv(pre_path)
post_survey = pd.read_csv(post_path)

# save paths
save_path = phase_path + 'affectivevr/survey_plots/'
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

plt9.options.figure_size = (6, 6) # set default figure size in inches
plt9.options.dpi = 300 # set default dpi

# %%
# Loop through pre and post survey data
for survey, range_responses, break_steps in zip([pre_survey, post_survey], [range_responses_pre, range_responses_post], [break_step_pre, break_step_post]):
    # Pre survey data
    score_list = survey['questionnaire'].unique()
    # factors to group by
    factors = ['test_site', 'gender', 'position']

    # loop through each score
    for score in score_list:
        # save path
        save_score_path = save_path + score + '/'
        if not os.path.exists(save_score_path): # create folder if it doesn't exist
            os.makedirs(save_score_path)
        # select score
        score_df = survey[survey['questionnaire'] == score]
        # for plotting, create sub-id column by converting 'sub' to category and map categories to integers (starting from 1)
        score_df.loc[:,'sub-id'] = score_df.loc[:,'sub'].astype('category').cat.codes + 1

        # loop through each factor
        for factor in factors:
            # plot score distribution
            p = (
            plt9.ggplot(score_df, plt9.aes(x=factor, y='response'))
            # Add the violin plot layer
            + plt9.geom_violin(plt9.aes(fill=factor), position=plt9.position_nudge(x=0.2), style='right', alpha=0.5)
            # Add the points layer
            + plt9.geom_point(plt9.aes(color= 'sub-id'), position=plt9.position_jitter(width=0.1, height=0), size=1.5) 
            # Add the boxplot layer
            + plt9.geom_boxplot(plt9.aes(fill=factor), width=0.1, alpha=0.5)
            # Set the y-axis labels and breaks
            + plt9.scale_y_continuous(name = score, breaks= np.arange(range_responses[score][0], range_responses[score][1]+1, break_steps[score]), limits = range_responses[score])
            # Set the x-axis label
            + plt9.xlab(xlab=factor)
            # Set the color scale legend
            + plt9.scale_color_continuous(name = 'sub-id')
            # Set the theme to a white background with a base font size of 14
            + plt9.theme_bw(base_size = 14)
            )
            p.save(filename=save_score_path + score + '_' + factor + '.png')
        
        # plot score distribution for position in BER only
        p = (
        plt9.ggplot(score_df[score_df['test_site'] == 'BER'], plt9.aes(x='position', y='response'))
        # Add the violin plot layer
        + plt9.geom_violin(plt9.aes(fill='position'), position=plt9.position_nudge(x=0.2), style='right', alpha=0.5)
        # Add the points layer
        + plt9.geom_point(plt9.aes(color= 'sub-id'), position=plt9.position_jitter(width=0.1, height=0), size=1.5)
        # Add the boxplot layer
        + plt9.geom_boxplot(plt9.aes(fill='position'), width=0.1, alpha=0.5)
        # Set the y-axis labels and breaks
        + plt9.scale_y_continuous(name = score, breaks= np.arange(range_responses[score][0], range_responses[score][1]+1, break_steps[score]), limits = range_responses[score])
        # Set the x-axis label
        + plt9.xlab(xlab='BER position')
        # Set the color scale legend
        + plt9.scale_color_continuous(name = 'sub-id')
        # Set the theme to a white background with a base font size of 14
        + plt9.theme_bw(base_size = 14)
        )
        p.save(filename=save_score_path + score + '_position_BER.png')
        
# %%
