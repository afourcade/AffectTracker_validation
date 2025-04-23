########################################################################################################################
# Script to import Pre Experiment survey data from AVR experiment
# Output: 2 csv files (dataframes): participants.csv (demographics) & pre_survey.csv (questionnaires)
# Author: Antonin Fourcade
# Last version: 03.05.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os

# %%
## Set paths and experiment parameters
data_path = 'E:/AffectiveVR/Phase_1/Data/' # path to data files

# experiment parameters
file_survey = 'survey.csv' # files for each site

debug = True  # debug mode

# set questionnaires parameters
nb_q = {'vr_exp': 2, 'ssq_pre': 6, 'rm_pref': 1, 'ssq_post': 6} # nb of questions for each questionnaire
questionnaires = ['vr_exp']*nb_q['vr_exp'] + ['ssq_pre']*nb_q['ssq_pre'] + ['rm_pref']*nb_q['rm_pref'] + ['ssq_post']*nb_q['ssq_post'] # list of questionnaires responses types
idx_q = list(range(1, nb_q['vr_exp']+1)) + list(range(1, nb_q['ssq_pre']+1)) + list(range(1, nb_q['rm_pref']+1)) + list(range(1, nb_q['ssq_post']+1)) # list of idx for each question in questionnaires  

#range_responses_q = {'vr_exp': [0,3], 'ssq_pre': [0,6], 'rm_pref': [1,3], 'ssq_post': [0,6]} # range of responses for each questionnaire

# dictionaries for survey data
survey_code = {'sj_id': 'ID', 'vr_pre_exp': 'Pre_VR exp', 'game_prof': 'Games_VR_prof', 'rm_pref': 'Preferred'}
vr_pre_exp_code = {0: 'never', 1: '1-3h', 2: '6-15h', 3: '15h+'}
game_prof_code = {0: 'not_at_all', 1: 'beginner', 2: 'medium', 3: 'expert'}

# SSQ scores
# dictionary for SSQ questionnaire
ssq_q = {'general_discomfort':1, 'nausea':2, 
              'dizziness':3, 'headache':4, 
              'blurred_vision':5, 'difficulty_concentrating':6}
# dictionary for SSQ scales
ssq_scales = {'ssq_nausea':['general_discomfort', 'nausea', 'difficulty_concentrating'], 
              'ssq_oculomotor':['general_discomfort', 'headache', 'difficulty_concentrating'], 
              'ssq_disorientation':['nausea', 'dizziness', 'blurred_vision']}

# %%
# Initialize variables for dataframe creation
sub = []
q_type = []
q_index = []
q_response = []

# %%
# debug
if debug:
    sj = 2
    sj_id = 'Y45K2'
    t_s = 'TOR'

# %%

# import data
survey = pd.read_csv(data_path + 'AVR/' + file_survey)

# get sj list for this site, change into bids-format sub-xx and print nb
sj_list_site = survey[survey_code['sj_id']].values 
print("Number of SJ in this site: " + str(len(sj_list_site)))
# MISSING TWO SJs and ids not correct!!!

# create survey_preprocessed dataframe
survey_preprocessed = pd.DataFrame({'sj_id': [], 'vr_exp': [], 'pre_ssq_nausea': [], 'pre_ssq_oculomotor': [], 'pre_ssq_disorientation': [], 'rm_pref': [], 'post_ssq_nausea': [], 'post_ssq_oculomotor': [], 'post_ssq_disorientation': []})

# Questionnaires dataframe
# Loop over participants
for sj, sj_id in enumerate(sj_list_site):
    # select questionnaires data for sj (drop demographics columns)
    survey_sj = survey[survey[survey_code['sj_id']] == sj_id].drop(survey.columns[-3:], axis=1)
    # get labels of sj responses  
    colnames = survey_sj.columns.values.tolist()  # get column names
    responses_labels_sj = [colnames[x] for x in np.where(np.array(survey_sj) == checked_code['checked'])[1].tolist()] # checked answers
    # get sj responses values  
    responses_sj = [int(x[-1]) for x in responses_labels_sj]
    # check nb of responses
    if len(responses_sj) != len(questionnaires):
        print("Warning: nb of responses for SJ " + str(sj_id) + " is " + str(len(responses_sj)) + " instead of " + str(len(questionnaires)))


    sub = np.append(sub, survey_sj[survey_code['sj_id']].values[0])
    q_type = np.append(q_type, questionnaires[i_r])
    q_index = np.append(q_index, idx_q[i_r])
    q_response = np.append(q_response, resp)
    # show progress
    print(sub_dict[t_s][sj_id] + " done") 

# %%
# create pre_survey dataframe with Questionnaire data
d = {'sub': sub, 'test_site': site, 'position': position, 'gender': gender, 'questionnaire': q_type, 'index': q_index, 'response': q_response}
filename = data_path + 'AVR/' + 'pre_survey.csv'
pre_survey = pd.DataFrame(data=d)
pre_survey.sort_values(by=['sub', 'questionnaire', 'index'], inplace=True) # sort by sub

# %%
# save dataframe
pre_survey.to_csv(filename, na_rep='NaN', index=False)