########################################################################################################################
# Script to import Post Experiment survey data from AVR experiment
# Output: 3 csv files (dataframes): post_survey.csv, post_survey_corrected.csv (with reformat to 7-pt Likert scale) and
# open_feedback.csv
# Author: Antonin Fourcade
# Last version: 30.01.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os

# %%
## Set paths and experiment parameters
data_path = 'D:/AffectiveVR/Phase_2/Data/' # path to data files

nb_q = {'invasiveness': 1, 'emo_rep': 1, 'presence': 2, 'sus': 7, 'satisfaction': 1} # nb of questions for each questionnaire
questionnaires = ['invasiveness']*nb_q['invasiveness'] + ['emo_rep']*nb_q['emo_rep'] + ['presence']*nb_q['presence'] + ['sus']*nb_q['sus'] + ['satisfaction']*nb_q['satisfaction'] # list of questionnaires responses types
idx_q = list(range(1, nb_q['invasiveness']+1)) + list(range(1, nb_q['emo_rep']+1)) + list(range(1, nb_q['presence']+1)) + list(range(1, nb_q['sus']+1)) + list(range(1, nb_q['satisfaction']+1)) # list of idx for each question in questionnaires  

# experiment parameters
test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2
files_site = {'TOR': 'data_AVRpostTOR_2024-01-19_15-17.csv',  
             'BER': 'data_AVRpostBER_2024-01-19_15-06.csv'} # files for each site

debug = False  # debug mode

# dictionaries for survey data
survey_code = {'sj_id': 'A102s', 'satisfaction': 'D101', 'open_feedback': 'F101s'}
checked_code = {'not checked': 1, 'checked': 2}

# dictionnary for sj_ids in bids format
sub_dict = {'BER': {15: 'sub-01', 17: 'sub-02', 19: 'sub-03', 21: 'sub-04', 23: 'sub-05', 25: 'sub-06', 27: 'sub-07', 29: 'sub-08', 31: 'sub-09', 33: 'sub-10', 35: 'sub-11', 37: 'sub-12', 39: 'sub-13', 41: 'sub-14', 43: 'sub-15', 45: 'sub-16', 47: 'sub-17', 49: 'sub-18', 51: 'sub-19', 53: 'sub-20', 55: 'sub-21', 57: 'sub-22', 59: 'sub-23', 61: 'sub-24', 63: 'sub-25', 65: 'sub-26', 67: 'sub-27', 69: 'sub-28', 71: 'sub-29', 73: 'sub-30', 75: 'sub-31', 77: 'sub-32', 79: 'sub-33', 81: 'sub-34', 83: 'sub-35', 85: 'sub-36', 87: 'sub-37', 89: 'sub-38', 91: 'sub-39', 93: 'sub-40', 95: 'sub-41', 97: 'sub-42'}, 
            'TOR': {'ZVJGB': 'sub-43', 'Y45K2': 'sub-44', 'Y4ZSQ': 'sub-45', 'XS2L5': 'sub-46', 'VJ0O2': 'sub-47', 'UWON3': 'sub-48', 'TYECS': 'sub-49', 'TDT7C': 'sub-50', 'RKWE8': 'sub-51', 'QSVKO': 'sub-52', 'Q5WEM': 'sub-53', 'OR8FW': 'sub-54', 'NQDJS': 'sub-55', 'MW7PU': 'sub-56', 'MQNKL': 'sub-57', 'M9ZP2': 'sub-58', 'KDC5Y': 'sub-59', 'I3W3C': 'sub-60', 'H52P0': 'sub-61', 'H4K8O': 'sub-62', 'G93XP': 'sub-63', 'FRIZD': 'sub-64', 'D1TP2': 'sub-65', 'BXD21': 'sub-66', 'BQLT9': 'sub-67', '79985': 'sub-68', '6747F': 'sub-69', '96WU9': 'sub-70', '78FK0': 'sub-71', '18GBC': 'sub-72', '9P6N0': 'sub-73', '8INPH': 'sub-74', '8D8KN': 'sub-75', '8D01A': 'sub-76', '7HEIW': 'sub-77', '6MRCB': 'sub-78', '5KNP8': 'sub-79', '3CSH0': 'sub-80', '2QHD2': 'sub-81', '2PSC5': 'sub-82', '1UZ0Q': 'sub-83'}}

# participant file containing position data
participant_file_path = 'AVR/participants.csv'
participant_file = pd.read_csv(data_path + participant_file_path, sep=',', header=0)

# %%
# Initialize variables for dataframe creation
sub = []
q_type = []
q_index = []
q_response = []
site = []
position = []
gender = []
open_feedback = []
sub_list = []	

# %%
# debug
if debug:
    sj = 0
    sj_id = 15
    t_s = 'Berlin'

# %%
# Loop over sites
for t_s in test_site: 
    # import data
    assessment_site = pd.read_csv(data_path + 'AVR-' + t_s + '/' + files_site[t_s], sep=',', header=0, encoding='utf-16')

    # get sj list for this site, change into bids-format sub-xx and print nb
    sj_list_site = assessment_site[survey_code['sj_id']].values 
    sub_list_site = assessment_site[survey_code['sj_id']].replace(sub_dict[t_s]).values 
    print("Number of SJ in this site: " + str(len(sub_list_site)))

    # get open feedback data
    open_feedback_site = assessment_site[survey_code['open_feedback']].values
    # prepare open feedback and list of sub for dataframe
    open_feedback = np.append(open_feedback, open_feedback_site)
    sub_list = np.append(sub_list, sub_list_site)
    # get satisfaction data 
    satisfaction_site = assessment_site[survey_code['satisfaction']].values

    # a bit of cleaning - keep only relevant data (rest of the questionnaires)
    assessment_site = assessment_site.drop(assessment_site.columns[:7], axis=1) # remove first 7 columns
    assessment_site = assessment_site.drop(assessment_site.columns[-15:], axis=1) # remove last 15 columns
    assessment_site = assessment_site.drop(columns=[survey_code['satisfaction']]) # remove satisfaction columns

    # Loop over participants
    for sj, sj_id in enumerate(sj_list_site):
        # select data for sj
        assessment_site_sj = assessment_site[assessment_site[survey_code['sj_id']] == sj_id]
        # get labels of sj responses  
        colnames = assessment_site_sj.columns.values.tolist()  # get column names
        responses_labels_sj = [colnames[x] for x in np.where(np.array(assessment_site_sj) == checked_code['checked'])[1].tolist()] # checked answers
        # get sj responses values  
        responses_sj = [int(x[-1]) for x in responses_labels_sj]
        # add satisfaction response  
        responses_sj = np.append(responses_sj, satisfaction_site[sj])  
        # check nb of responses
        if len(responses_sj) != len(questionnaires):
            print("Warning: nb of responses for SJ " + str(sj_id) + " is " + str(len(responses_sj)) + " instead of " + str(len(questionnaires)))

        # Loop over responses to create dataframe
        for i_r, resp in enumerate(responses_sj):
            sub = np.append(sub, sub_dict[t_s][sj_id]) # add sub
            q_type = np.append(q_type, questionnaires[i_r]) # add questionnaire type
            q_index = np.append(q_index, idx_q[i_r]) # add questionnaire index
            q_response = np.append(q_response, resp) # add response
            site = np.append(site, t_s) # add test site
            position_sub = participant_file[participant_file['sub'] == sub_dict[t_s][sj_id]]['position'].values # get position of sub
            position = np.append(position, position_sub) # add position
            gender_sub = participant_file[participant_file['sub'] == sub_dict[t_s][sj_id]]['gender'].values # get gender of sub
            gender = np.append(gender, gender_sub) # add gender
        # show progress
        print(sub_dict[t_s][sj_id] + " done") 

# %%
# create and save dataframe with Questionnaire data
d = {'sub': sub, 'test_site': site, 'position': position, 'gender': gender, 'questionnaire': q_type, 'index': q_index, 'response': q_response}
filename = data_path + 'AVR/' + 'post_survey.csv'
post_survey = pd.DataFrame(data=d)
post_survey.sort_values(by=['sub', 'questionnaire', 'index'], inplace=True) # sort by sub, questionnaire and index
post_survey.to_csv(filename, na_rep='NaN', index=False)

# %%
# Reformat responses with 8-pt Likert scale into 7-pt Likert scale
# Emotion representation questionnaire
emo_rep_df = post_survey[post_survey['questionnaire'] == 'emo_rep']
# Linearly map 1-8 to 1-7
emo_rep_df['response'] = ((emo_rep_df['response'] - 1) / (8 - 1)) * (7 - 1) + 1
# Presence questionnaire
presence_df = post_survey[post_survey['questionnaire'] == 'presence']
# Linearly map 1-8 to 1-7
presence_df['response'] = ((presence_df['response'] - 1) / (8 - 1)) * (7 - 1) + 1
# SUS questionnaire but only for Torino
sus_df = post_survey[(post_survey['questionnaire'] == 'sus') & (post_survey['test_site'] == 'TOR')]
# Linearly map 1-8 to 1-7
sus_df['response'] = ((sus_df['response'] - 1) / (8 - 1)) * (7 - 1) + 1
# Merge dataframes
post_survey[post_survey['questionnaire'] == 'emo_rep'] = emo_rep_df
post_survey[post_survey['questionnaire'] == 'presence'] = presence_df
post_survey[(post_survey['questionnaire'] == 'sus') & (post_survey['test_site'] == 'TOR')] = sus_df
# Save dataframe
filename = data_path + 'AVR/' + 'post_survey_corrected.csv'
post_survey.to_csv(filename, na_rep='NaN', index=False)

# %%
# create and save dataframe with open feedback data
d_of = {'sub': sub_list, 'open_feedback': open_feedback}
filename_of = data_path + 'AVR/' + 'open_feedback.csv'
open_feedback_df = pd.DataFrame(data=d_of)
open_feedback_df.sort_values(by=['sub'], inplace=True) # sort by sub
open_feedback_df.to_csv(filename_of, na_rep='NaN', index=False)

# %%
