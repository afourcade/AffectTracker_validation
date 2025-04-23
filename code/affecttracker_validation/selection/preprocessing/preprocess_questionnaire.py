########################################################################################################################
# Script to preprocess assessment data from AVR experiment
# Output: csv file (dataframe): assessment_preprocessed.csv
# Author: Antonin Fourcade
# Last version: 16.02.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os

# %%
# Set paths and experiment parameters
phase_path = 'D:/AffectiveVR/Phase_1/'
data_path = phase_path + 'Data/'

# survey parameters
nb_q = {'invasiveness': 1, 'emo_rep': 1, 'presence': 2, 'sus': 7, 'satisfaction': 1} # nb of questions for each questionnaire

debug = False  # debug mode

# path to post survey csv file and load data
filename_assess = 'assessment.csv'
assess_path = data_path + filename_assess
assessment = pd.read_csv(assess_path)

# %%
# Invasiveness
# Select invasiveness questionnaire: questionnaire = invasive_presence and index = 1
invasiveness = assessment[(assessment['questionnaire'] == 'invasive_presence') & (assessment['index'] == 1)]
# rename elements in questionnaire column to invasiveness
invasiveness['questionnaire'] = 'invasiveness'
# drop index column
invasiveness = invasiveness.drop(columns=['index'])

# %%
# Satisfaction
# Select satisfaction questionnaire: questionnaire = Kunin
satisfaction = assessment[(assessment['questionnaire'] == 'Kunin')]
# rename elements in questionnaire column to satisfaction
satisfaction['questionnaire'] = 'satisfaction'
# drop index column
satisfaction = satisfaction.drop(columns=['index'])

# %%
# Emotion representation
# Select emotion representation questionnaire: questionnaire = invasive_presence and index = 2
emo_rep = assessment[(assessment['questionnaire'] == 'invasive_presence') & (assessment['index'] == 2)]
# rename elements in questionnaire column to emo_rep
emo_rep['questionnaire'] = 'emo_rep'
# drop index column
emo_rep = emo_rep.drop(columns=['index'])

# %%
# Presence score
# Select presence questionnaire: questionnaire = invasive_presence and index = 3 and 4
presence = assessment[(assessment['questionnaire'] == 'invasive_presence') & (assessment['index'] >= 3)]
# change index 3 and 4 to 1 and 2
presence.loc[presence['index'] == 3, 'index'] = 1
presence.loc[presence['index'] == 4, 'index'] = 2
# Reverse code response (0 to 6) with index 2
presence.loc[presence['index'] == 2, 'response'] = 6 - presence.loc[presence['index'] == 2, 'response']
# Compute presence score: mean of the responses
presence_score = presence.groupby(['sj_id', 'test_site', 'rating_method'])['response'].mean().reset_index()
# add column for questionnaire
presence_score['questionnaire'] = 'presence_score'

# %%
# SUS score
# Select SUS questionnaire, range 0 to 6
sus = assessment[(assessment['questionnaire'] == 'SUS')]
# For each of the even numbered questions (positive), subtract 0 from the score.
# For each of the odd numbered questions (negative), subtract their value from 6.
sus.loc[sus['index'] % 2 == 1, 'response'] = 6 - sus.loc[sus['index'] % 2 == 1, 'response']
# Take these new values and add up the total score. Then multiply this by 100/(7*6)=2.38
sus_score = sus.groupby(['sj_id', 'test_site', 'rating_method'])['response'].sum().reset_index()
sus_score['response'] = sus_score['response']*2.38
# add columns for questionnaire, site, position, gender to sus_score
sus_score['questionnaire'] = 'sus_score'

# %%
# create new dataframe with presence and sus scores, as well as invasiveness, emotion representation and satisfaction responses
assessment_preprocessed = pd.concat([invasiveness, satisfaction, emo_rep, presence_score, sus_score], ignore_index=True)
# sort by sub and questionnaire
assessment_preprocessed.sort_values(by=['sj_id', 'test_site', 'rating_method', 'questionnaire'], inplace=True)
# save dataframe
filename = data_path + 'assessment_preprocessed.csv'
assessment_preprocessed.to_csv(filename, na_rep='NaN', index=False)
print('Dataframe saved as ' + filename)
