########################################################################################################################
# Script to import Pre Experiment survey data from AVR experiment
# Output: 2 csv files (dataframes): participants.csv (demographics) & pre_survey.csv (questionnaires)
# Author: Antonin Fourcade
# Last version: 07.02.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os

# %%
## Set paths and experiment parameters
data_path = 'D:/AffectiveVR/Phase_2/Data/' # path to data files

# experiment parameters
test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# demographics = ['age', 'gender', 'education'] # list of demographics
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2
# position_sj = ['seated', 'standing']  # Important -> TOR: all sjs seated, BER: first n_standing_BER standing, rest seated
n_standing_BER = 21  # nb of standing sjs in BER
files_site = {'TOR': 'data_AVRpreTOR_2024-01-19_15-10.csv',  
             'BER': 'data_AVRpreBER_2024-01-19_15-05.csv'} # files for each site

debug = True  # debug mode

# set questionnaires parameters
nb_q = {'vr_exp': 2, 'maia': 37, 'ssq': 6, 'tas': 20} # nb of questions for each questionnaire
questionnaires = ['vr_exp']*nb_q['vr_exp'] + ['maia']*nb_q['maia'] + ['ssq']*nb_q['ssq'] + ['tas']*nb_q['tas'] # list of questionnaires responses types
idx_q = list(range(1, nb_q['vr_exp']+1)) + list(range(1, nb_q['maia']+1)) + list(range(1, nb_q['ssq']+1)) + list(range(1, nb_q['tas']+1)) # list of idx for each question in questionnaires  

# dictionaries for survey data
survey_code = {'sj_id': 'A102s', 'age': 'A103s', 'gender': 'A104', 'education': 'A105'}
gender_code = {1: 'female', 2: 'male', 3: 'non_binary', 4: 'other', -9: 'no_response'}
education_code = {1: 'elementary', 2: 'middle_school', 3: 'high_school', 4: 'bachelor', 5: 'master', 6: 'doctorate', -9: 'no_response'}
checked_code = {'not checked': 1, 'checked': 2}

# dictionnary for sj_ids in bids format
sub_dict = {'BER': {15: 'sub-01', 17: 'sub-02', 19: 'sub-03', 21: 'sub-04', 23: 'sub-05', 25: 'sub-06', 27: 'sub-07', 29: 'sub-08', 31: 'sub-09', 33: 'sub-10', 35: 'sub-11', 37: 'sub-12', 39: 'sub-13', 41: 'sub-14', 43: 'sub-15', 45: 'sub-16', 47: 'sub-17', 49: 'sub-18', 51: 'sub-19', 53: 'sub-20', 55: 'sub-21', 57: 'sub-22', 59: 'sub-23', 61: 'sub-24', 63: 'sub-25', 65: 'sub-26', 67: 'sub-27', 69: 'sub-28', 71: 'sub-29', 73: 'sub-30', 75: 'sub-31', 77: 'sub-32', 79: 'sub-33', 81: 'sub-34', 83: 'sub-35', 85: 'sub-36', 87: 'sub-37', 89: 'sub-38', 91: 'sub-39', 93: 'sub-40', 95: 'sub-41', 97: 'sub-42'}, 
            'TOR': {'ZVJGB': 'sub-43', 'Y45K2': 'sub-44', 'Y4ZSQ': 'sub-45', 'XS2L5': 'sub-46', 'VJ0O2': 'sub-47', 'UWON3': 'sub-48', 'TYECS': 'sub-49', 'TDT7C': 'sub-50', 'RKWE8': 'sub-51', 'QSVKO': 'sub-52', 'Q5WEM': 'sub-53', 'OR8FW': 'sub-54', 'NQDJS': 'sub-55', 'MW7PU': 'sub-56', 'MQNKL': 'sub-57', 'M9ZP2': 'sub-58', 'KDC5Y': 'sub-59', 'I3W3C': 'sub-60', 'H52P0': 'sub-61', 'H4K8O': 'sub-62', 'G93XP': 'sub-63', 'FRIZD': 'sub-64', 'D1TP2': 'sub-65', 'BXD21': 'sub-66', 'BQLT9': 'sub-67', '79985': 'sub-68', '6747F': 'sub-69', '96WU9': 'sub-70', '78FK0': 'sub-71', '18GBC': 'sub-72', '9P6N0': 'sub-73', '8INPH': 'sub-74', '8D8KN': 'sub-75', '8D01A': 'sub-76', '7HEIW': 'sub-77', '6MRCB': 'sub-78', '5KNP8': 'sub-79', '3CSH0': 'sub-80', '2QHD2': 'sub-81', '2PSC5': 'sub-82', '1UZ0Q': 'sub-83'}}


# %%
# Initialize variables for dataframe creation
demographic = pd.DataFrame({'sub': [], 'test_site': [], 'position': [], 'age': [], 'gender': [], 'education': []})
sub = []
q_type = []
q_index = []
q_response = []
site = []
position = []
gender = []

# %%
# debug
if debug:
    sj = 2
    sj_id = 'Y45K2'
    t_s = 'TOR'

# %%
# Loop over sites
for t_s in test_site: 
    # import data
    assessment_site = pd.read_csv(data_path + 'AVR-' + t_s + '/' + files_site[t_s], sep=',', header=0, encoding='utf-16')
    
    # a bit of cleaning - keep only relevant data
    assessment_site = assessment_site.drop(assessment_site.columns[:7], axis=1) # remove first 7 columns
    assessment_site = assessment_site.drop(assessment_site.columns[-15:], axis=1) # remove last 15 columns
    
    # get sj list for this site, change into bids-format sub-xx and print nb
    sj_list_site = assessment_site[survey_code['sj_id']].values 
    sub_list_site = assessment_site[survey_code['sj_id']].replace(sub_dict[t_s]).values 
    print("Number of SJ in this site: " + str(len(sub_list_site)))
    
    # Demographics dataframe
    # get demographics data
    age_site = assessment_site[survey_code['age']].values # get age data
    gender_site = assessment_site[survey_code['gender']].replace(gender_code).values # get gender data
    education_site = assessment_site[survey_code['education']].replace(education_code).values # get education data
    dem_test_site = [t_s] * len(sub_list_site) # add site to list
    if t_s == 'TOR':
        # all sj are seated in Torino
        position_site = ['seated'] * len(sub_list_site)
    elif t_s == 'BER':
        # first n_standing_BER sj are seated in Berlin, rest standing
        position_site = np.append(['standing'] * n_standing_BER, ['seated'] * (len(sj_list_site) - n_standing_BER))
    # create demographic dataframe for this site
    demographic_site = pd.DataFrame({'sub': sub_list_site, 'test_site': dem_test_site, 'position': position_site, 'age': age_site, 'gender': gender_site, 'education': education_site})

    # append demographic_site dataframe to demographic dataframe
    #demographic = demographic.append(demographic_site, ignore_index=True) # add to demographics dataframe
    demographic = pd.concat([demographic, demographic_site])

    # Questionnaires dataframe
    # Loop over participants
    for sj, sj_id in enumerate(sj_list_site):
        # select questionnaires data for sj (drop demographics columns)
        assessment_site_sj = assessment_site[assessment_site[survey_code['sj_id']] == sj_id].drop(assessment_site.columns[:5], axis=1)
        # get labels of sj responses  
        colnames = assessment_site_sj.columns.values.tolist()  # get column names
        responses_labels_sj = [colnames[x] for x in np.where(np.array(assessment_site_sj) == checked_code['checked'])[1].tolist()] # checked answers
        # get sj responses values  
        responses_sj = [int(x[-1]) for x in responses_labels_sj]
        # check nb of responses
        if len(responses_sj) != len(questionnaires):
            print("Warning: nb of responses for SJ " + str(sj_id) + " is " + str(len(responses_sj)) + " instead of " + str(len(questionnaires)))

        # Loop over responses to create dataframe
        for i_r, resp in enumerate(responses_sj):
            sub = np.append(sub, sub_dict[t_s][sj_id])
            q_type = np.append(q_type, questionnaires[i_r])
            q_index = np.append(q_index, idx_q[i_r])
            q_response = np.append(q_response, resp)
            site = np.append(site, t_s)
            position = np.append(position, position_site[sj])
            gender = np.append(gender, gender_site[sj])
        # show progress
        print(sub_dict[t_s][sj_id] + " done") 


# %%
# save Demographics dataframe with demographics data
filename = data_path + 'AVR/' + 'participants.csv'
demographic.sort_values(by=['sub'], inplace=True) # sort by sub
demographic.to_csv(filename, na_rep='NaN', index=False)

# %%
# create pre_survey dataframe with Questionnaire data
d = {'sub': sub, 'test_site': site, 'position': position, 'gender': gender, 'questionnaire': q_type, 'index': q_index, 'response': q_response}
filename = data_path + 'AVR/' + 'pre_survey.csv'
pre_survey = pd.DataFrame(data=d)
pre_survey.sort_values(by=['sub', 'questionnaire', 'index'], inplace=True) # sort by sub

# %%
# translate all responses to MAIA questionnaire from 1-6 to 0-5
maia = pre_survey[(pre_survey['questionnaire'] == 'maia')]
maia.loc[:, 'response'] = maia.loc[:, 'response'] - 1
pre_survey.loc[pre_survey['questionnaire'] == 'maia', 'response'] = maia.loc[:, 'response']

# %%
# save dataframe
pre_survey.to_csv(filename, na_rep='NaN', index=False)