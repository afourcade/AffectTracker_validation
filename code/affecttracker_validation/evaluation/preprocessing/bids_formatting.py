########################################################################################################################
# Script to create BIDS-formatted AVR project folder from AVR experiment
# Output: AVR folder with sub-xx folders for each participant
# Author: Antonin Fourcade
# Last version: 22.01.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
import shutil

# %%
data_path = 'D:/AffectiveVR/Phase_2/Data/'
sj_list_TOR = os.listdir(data_path + 'AVR-TOR/') # Torino participants
sj_list_BER = os.listdir(data_path + 'AVR-BER/') # Berlin participants
sj_list_BER = [int(i) for i in sj_list_BER] # convert to int
sj_list_BER.sort() # sort in ascending order

# exclude BER participants because no soSci survey recorded
to_be_removed = {1, 3, 5, 7, 9, 11, 13}
sj_list_BER = [item for item in sj_list_BER if item not in to_be_removed]

# %%
# dictionnary for sj_ids in bids format
sub_dict = {'BER': {15: 'sub-01', 17: 'sub-02', 19: 'sub-03', 21: 'sub-04', 23: 'sub-05', 25: 'sub-06', 27: 'sub-07', 29: 'sub-08', 31: 'sub-09', 33: 'sub-10', 35: 'sub-11', 37: 'sub-12', 39: 'sub-13', 41: 'sub-14', 43: 'sub-15', 45: 'sub-16', 47: 'sub-17', 49: 'sub-18', 51: 'sub-19', 53: 'sub-20', 55: 'sub-21', 57: 'sub-22', 59: 'sub-23', 61: 'sub-24', 63: 'sub-25', 65: 'sub-26', 67: 'sub-27', 69: 'sub-28', 71: 'sub-29', 73: 'sub-30', 75: 'sub-31', 77: 'sub-32', 79: 'sub-33', 81: 'sub-34', 83: 'sub-35', 85: 'sub-36', 87: 'sub-37', 89: 'sub-38', 91: 'sub-39', 93: 'sub-40', 95: 'sub-41', 97: 'sub-42'}, 
            'TOR': {'ZVJGB': 'sub-43', 'Y45K2': 'sub-44', 'Y4ZSQ': 'sub-45', 'XS2L5': 'sub-46', 'VJ0O2': 'sub-47', 'UWON3': 'sub-48', 'TYECS': 'sub-49', 'TDT7C': 'sub-50', 'RKWE8': 'sub-51', 'QSVKO': 'sub-52', 'Q5WEM': 'sub-53', 'OR8FW': 'sub-54', 'NQDJS': 'sub-55', 'MW7PU': 'sub-56', 'MQNKL': 'sub-57', 'M9ZP2': 'sub-58', 'KDC5Y': 'sub-59', 'I3W3C': 'sub-60', 'H52P0': 'sub-61', 'H4K8O': 'sub-62', 'G93XP': 'sub-63', 'FRIZD': 'sub-64', 'D1TP2': 'sub-65', 'BXD21': 'sub-66', 'BQLT9': 'sub-67', '79985': 'sub-68', '6747F': 'sub-69', '96WU9': 'sub-70', '78FK0': 'sub-71', '18GBC': 'sub-72', '9P6N0': 'sub-73', '8INPH': 'sub-74', '8D8KN': 'sub-75', '8D01A': 'sub-76', '7HEIW': 'sub-77', '6MRCB': 'sub-78', '5KNP8': 'sub-79', '3CSH0': 'sub-80', '2QHD2': 'sub-81', '2PSC5': 'sub-82', '1UZ0Q': 'sub-83'}}

# %%
# create BIDS-formatted project folder (if not already existing)
if not os.path.exists(data_path + 'AVR'):
    os.mkdir(data_path + 'AVR')
# Loop over testing sites
for site in ['BER', 'TOR']:
    # Loop over participants
    for sj_id in sub_dict[site].keys():
        # get path to sj data folder
        sj_path = data_path + 'AVR-' + site + '/' + str(sj_id) + '/S000'
        other_path = sj_path + '/other'
        trackers_path = sj_path + '/trackers'
        session_info_path = sj_path + '/session_info'
        # create BIDS-formatted sub folder
        sub_path = data_path + 'AVR/' + sub_dict[site][sj_id]
        beh_path = sub_path + '/beh'
        motion_path = sub_path + '/motion'
        unity_path = sub_path + '/unity'
        if not os.path.exists(sub_path):
            os.mkdir(sub_path, )
        if not os.path.exists(beh_path):
            os.mkdir(beh_path)
        shutil.copytree(trackers_path, motion_path, dirs_exist_ok=True)
        shutil.copytree(session_info_path, unity_path, dirs_exist_ok=True)
        # copy data files from sj data folder to BIDS-formatted sub folder
        shutil.copy2(sj_path + '/trial_results.csv', sub_path + '/trial_results.csv')
        shutil.copy2(other_path + '/rating_CR_T001.csv', beh_path + '/rating_CR_T001.csv')
        shutil.copy2(other_path + '/rating_CR_T002.csv', beh_path + '/rating_CR_T002.csv')
        shutil.copy2(other_path + '/rating_SR_T001.csv', beh_path + '/rating_SR_T001.csv')
        shutil.copy2(other_path + '/rating_SR_T002.csv', beh_path + '/rating_SR_T002.csv')
        shutil.copy2(other_path + '/executionOrder.csv', unity_path + '/executionOrder.csv')
        shutil.copy2(other_path + '/markerLog.csv', unity_path + '/markerLog.csv')
        # change sj_id in trial_results.csv to sub-xx
        trial_results = pd.read_csv(sub_path + '/trial_results.csv', sep=',', header=0)
        trial_results.replace({'ppid': {sj_id: sub_dict[site][sj_id]}}, inplace=True)
        # change paths in trial_results.csv to BIDS-formatted folder
        trial_results.replace(other_path.replace(data_path,'').replace('-' + site, ''), beh_path.replace(data_path,''), regex=True, inplace=True)
        trial_results.replace(trackers_path.replace(data_path,'').replace('-' + site, ''), motion_path.replace(data_path,''), regex=True, inplace=True)
        # save updated trial_results.csv
        trial_results.to_csv(sub_path + '/trial_results.csv', index=False)
