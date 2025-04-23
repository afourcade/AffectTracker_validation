########################################################################################################################
# Script to preprocess Post Experiment survey data from AVR experiment
# Output: csv file (dataframe): post_survey[_corrected]_preprocessed.csv
# Author: Antonin Fourcade
# Last version: 07.02.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os

# %%
# Set paths and experiment parameters
phase_path = 'D:/AffectiveVR/Phase_2/'
data_path = phase_path + 'Data/'
bids_folder = 'AVR' # folder with data in bids format

# survey parameters
nb_q = {'invasiveness': 1, 'emo_rep': 1, 'presence': 2, 'sus': 7, 'satisfaction': 1} # nb of questions for each questionnaire
questionnaires = ['invasiveness']*nb_q['invasiveness'] + ['emo_rep']*nb_q['emo_rep'] + ['presence']*nb_q['presence'] + ['sus']*nb_q['sus'] + ['satisfaction']*nb_q['satisfaction'] # list of questionnaires responses types
idx_q = list(range(1, nb_q['invasiveness']+1)) + list(range(1, nb_q['emo_rep']+1)) + list(range(1, nb_q['presence']+1)) + list(range(1, nb_q['sus']+1)) + list(range(1, nb_q['satisfaction']+1)) # list of idx for each question in questionnaires  

debug = False  # debug mode

# path to post survey csv file and load data
filename_post = 'post_survey_corrected.csv'
post_path = data_path + bids_folder + '/' + filename_post
post_survey = pd.read_csv(post_path)

# %%
# Presence score
# Select presence questionnaire
presence = post_survey[(post_survey['questionnaire'] == 'presence')]
# Reverse code response (1 to 7) with index 2
presence.loc[presence['index'] == 2, 'response'] = 8 - presence.loc[presence['index'] == 2, 'response']
# Compute presence score: mean of the responses
presence_score = presence.groupby('sub')['response'].mean().reset_index()
# add columns for questionnaire, site, position, gender to presence_score
presence_score['questionnaire'] = 'presence_score'
presence_score['test_site'] = presence.groupby('sub')['test_site'].first().values
presence_score['position'] = presence.groupby('sub')['position'].first().values
presence_score['gender'] = presence.groupby('sub')['gender'].first().values

# %%
# SUS score
# Select SUS questionnaire
sus = post_survey[(post_survey['questionnaire'] == 'sus')]
# For each of the even numbered questions (positive), subtract 1 from the score.
sus.loc[sus['index'] % 2 == 0, 'response'] = sus.loc[sus['index'] % 2 == 0, 'response'] - 1
# For each of the odd numbered questions (negative), subtract their value from 7.
sus.loc[sus['index'] % 2 == 1, 'response'] = 7 - sus.loc[sus['index'] % 2 == 1, 'response']
# Take these new values and add up the total score. Then multiply this by 100/(7*6)=2.38
sus_score = sus.groupby('sub')['response'].sum().reset_index()
sus_score['response'] = sus_score['response']*2.38
# add columns for questionnaire, site, position, gender to sus_score
sus_score['questionnaire'] = 'sus_score'
sus_score['test_site'] = sus.groupby('sub')['test_site'].first().values
sus_score['position'] = sus.groupby('sub')['position'].first().values
sus_score['gender'] = sus.groupby('sub')['gender'].first().values

# %%
# create new dataframe with presence and sus scores, as well as invasiveness, emotion representation and satisfaction responses
post_survey_preprocessed = pd.concat([presence_score, sus_score], ignore_index=True)
# add invasiveness, emotion representation and satisfaction responses (remove index column)
post_survey_preprocessed = pd.concat([post_survey[(post_survey['questionnaire'] == 'invasiveness') | (post_survey['questionnaire'] == 'emo_rep') | (post_survey['questionnaire'] == 'satisfaction')].drop(columns=['index']), post_survey_preprocessed], ignore_index=True)
# sort by sub and questionnaire
post_survey_preprocessed.sort_values(by=['sub', 'questionnaire'], inplace=True)
# save dataframe
filename = data_path + 'AVR/' + 'post_survey_preprocessed.csv'
post_survey_preprocessed.to_csv(filename, na_rep='NaN', index=False)
print('Dataframe saved as ' + filename)

