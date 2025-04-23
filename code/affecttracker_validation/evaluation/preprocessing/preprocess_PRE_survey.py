########################################################################################################################
# Script to preprocess Pre Experiment survey data from AVR experiment
# Output: csv file (dataframe): pre_survey_preprocessed.csv
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
#nb_q = {'vr_exp': 2, 'maia': 37, 'ssq': 6, 'tas': 20} # nb of questions for each questionnaire
#range_responses_q = {'vr_exp': [1,4], 'maia': [0,5], 'ssq': [1,7] - [0,3], 'tas': [1,5]} # range of responses for each questionnaire

debug = True  # debug mode

# path to pre survey csv file and load data
filename_pre = 'pre_survey.csv'
pre_path = data_path + bids_folder + '/' + filename_pre
pre_survey = pd.read_csv(pre_path)

# %%
# MAIA scores
# dictionary for MAIA scales and items
maia_scales = {'maia_noticing':np.arange(1,5), 'maia_not_distracting':np.arange(5,11), 
               'maia_not_worrying':np.arange(11,16), 'maia_attention_regulation':np.arange(16,23), 
               'maia_emotional_awareness':np.arange(23,28), 'maia_self_regulation':np.arange(28,32), 
               'maia_body_listening':np.arange(32,35), 'maia_trust':np.arange(35,38)}
# Select MAIA questionnaire
maia = pre_survey[(pre_survey['questionnaire'] == 'maia')]
# Compute MAIA scores for each scale (mean of the responses for each sub)
maia_scores = pd.DataFrame()
for scale, items in maia_scales.items():
    maia_scores[scale] = maia[maia['index'].isin(items)].groupby('sub')['response'].mean()
# add columns for test_site, position, gender to maia_scores
maia_scores['test_site'] = maia.groupby('sub')['test_site'].first()
maia_scores['position'] = maia.groupby('sub')['position'].first()
maia_scores['gender'] = maia.groupby('sub')['gender'].first()
# pivot dataframe to long format
maia_scores = maia_scores.reset_index().melt(id_vars=['sub', 'test_site', 'position', 'gender'], value_vars=list(maia_scales.keys()), var_name='questionnaire', value_name='response')
    

# %%
# TAS scores
# dictionary for TAS scales and items
tas_scales = {'tas_difficulty_identifying_feelings':[1,3,6,7,9,13,14], 
              'tas_difficulty_describing_feelings':[2,4,11,12,17],
              'tas_externally_oriented_thinking':[5,8,10,15,16,18,19,20]}
# list of TAS reverse coded items
tas_reverse = [4,5,10,18,19]
# Select TAS questionnaire
tas = pre_survey[(pre_survey['questionnaire'] == 'tas')]
# Reverse code response (1 to 5) for reverse coded items
tas.loc[tas['index'].isin(tas_reverse), 'response'] = 6 - tas.loc[tas['index'].isin(tas_reverse), 'response']
# Compute TAS scores for each scale
tas_scores = pd.DataFrame()
for scale, items in tas_scales.items():
    tas_scores[scale] = tas[tas['index'].isin(items)].groupby('sub')['response'].sum()
# add column with total score
tas_scores['tas_total_score'] = tas_scores['tas_difficulty_identifying_feelings'] + tas_scores['tas_difficulty_describing_feelings'] + tas_scores['tas_externally_oriented_thinking']
# add columns for test_site, position, gender to tas_scores
tas_scores['test_site'] = tas.groupby('sub')['test_site'].first()
tas_scores['position'] = tas.groupby('sub')['position'].first()
tas_scores['gender'] = tas.groupby('sub')['gender'].first()
# pivot dataframe to long format
tas_scores = tas_scores.reset_index().melt(id_vars=['sub', 'test_site', 'position', 'gender'], value_vars=list(tas_scales.keys()) + ['tas_total_score'], var_name='questionnaire', value_name='response')

# %%
# SSQ scores
# dictionary for SSQ questionnaire
ssq_q = {'general_discomfort':1, 'nausea':2, 
              'dizziness':3, 'headache':4, 
              'blurred_vision':5, 'difficulty_concentrating':6}
# dictionary for SSQ scales
ssq_scales = {'ssq_nausea':['general_discomfort', 'nausea', 'difficulty_concentrating'], 
              'ssq_oculomotor':['general_discomfort', 'headache', 'difficulty_concentrating'], 
              'ssq_disorientation':['nausea', 'dizziness', 'blurred_vision']}
# Select SSQ questionnaire
ssq = pre_survey[(pre_survey['questionnaire'] == 'ssq')]
# Compute SSQ scores for each scale
ssq_scores = pd.DataFrame()
for scale, items in ssq_scales.items():
    # get index of questions for each scale
    q_idx = []
    for item in items:
        q_idx.append(ssq_q[item])
    # compute sum of the responses for each sub
    ssq_scores[scale] = ssq[ssq['index'].isin(q_idx)].groupby('sub')['response'].sum()
# add column with total score
ssq_scores['ssq_total_score'] = ssq_scores['ssq_nausea'] + ssq_scores['ssq_oculomotor'] + ssq_scores['ssq_disorientation']
# add columns for test_site, position, gender to ssq_scores
ssq_scores['test_site'] = ssq.groupby('sub')['test_site'].first()
ssq_scores['position'] = ssq.groupby('sub')['position'].first()
ssq_scores['gender'] = ssq.groupby('sub')['gender'].first()
# pivot dataframe to long format
ssq_scores = ssq_scores.reset_index().melt(id_vars=['sub', 'test_site', 'position', 'gender'], value_vars=list(ssq_scales.keys()) + ['ssq_total_score'], var_name='questionnaire', value_name='response')

# %%
# VR experience scores
# Select VR experience questionnaire
vr_exp = pre_survey[(pre_survey['questionnaire'] == 'vr_exp')]
# Compute VR experience scores by averaging the responses of the two questions, for each sub
vr_exp_scores = vr_exp.groupby('sub')['response'].mean().reset_index()
# add columns for questionnaire, test_site, position, gender to vr_exp_scores
vr_exp_scores['questionnaire'] = 'vr_exp_score'
vr_exp_scores['test_site'] = vr_exp.groupby('sub')['test_site'].first().values
vr_exp_scores['position'] = vr_exp.groupby('sub')['position'].first().values
vr_exp_scores['gender'] = vr_exp.groupby('sub')['gender'].first().values

# %%
# create new dataframe with all pre survey scores
pre_survey_preprocessed = pd.concat([maia_scores, tas_scores, ssq_scores, vr_exp_scores], ignore_index=True)
# sort dataframe by sub and questionnaire
pre_survey_preprocessed = pre_survey_preprocessed.sort_values(by=['sub', 'questionnaire']).reset_index(drop=True)
# save pre_survey_preprocessed to csv
pre_survey_preprocessed.to_csv(data_path + bids_folder + '/pre_survey_preprocessed.csv', index=False)

# %%
