########################################################################################################################
# Script to import CR data from AVR experiment
# Preprocessing steps (optional): resampling and cut first seconds of data
# Output: 2 csv files (dataframes): cr[_rs_clean].csv and cr[_rs_clean]_nb_samples.csv
# Author: Antonin Fourcade
# Last version: 29.01.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

#%matplotlib qt # plot in a separate window

# %%
# Set paths and experiment parameters
data_path = 'D:/AffectiveVR/Phase_2/Data/'
bids_folder = 'AVR' # folder with data in bids format

# experiment parameters
# blocks = ['Practice', 'Experiment']
# logging_freq = ['CR', 'SR']
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2

cr_fs = 1/0.05  # sampling frequency CR in Hz
vid_len = 90*4 + 253 + 390 + 381  # length of the sequence of videos in s: 4xScifi + Invasion + Asteroids + Underwood

# preprocessing parameters
resample = True  # resample CR to cr_fs (if samples are not even). Note: does not interpolate NA values
debug = False  # debug mode
clean = True  # clean CR: remove first clean_s seconds of CR
clean_s = 5  # seconds to remove from CR

# participant file
participant_file_path = 'participants.csv'
participant_file = pd.read_csv(data_path + bids_folder + '/' + participant_file_path, sep=',', header=0)

# get list of participants
sub_list = participant_file['sub'].values

# %%
# Initialize variables
participant = []
site = []
position = []
gender = []
cr_v = []
cr_a = []
cr_dist = []
cr_angle = []
cr_time = []
cr_time_diff = []
nb_samples = []

if debug:
    rh_time_diff = []
    sub = 'sub-01'

# %%
# Read and preprocess data
for sub_idx, sub in enumerate(sub_list):

    # Read results in csv file
    sub_path = data_path + bids_folder + '/' + sub + '/'
    trial_results_filename = sub_path + 'trial_results.csv'
    trials_results = pd.read_csv(trial_results_filename)

    # Select Experiment data
    trials_exp = trials_results[trials_results['block_name'] == 'Experiment']

    # Debug mode
    if debug:
        # check sampling of right hand tracking (90Hz)
        rh_loc = 'righthand_movement_location_0'
        rh_path = data_path + trials_exp[rh_loc].item()
        rh = pd.read_csv(rh_path)
        rh_t = rh['time'].values
        rh_t_diff = np.append(np.nan, np.diff(rh_t))

    # get CR from file
    # location of CR file
    cr_loc = 'rating_CR_location_0'
    cr_path = data_path + trials_exp[cr_loc].item()
    # load CR file
    cr = pd.read_csv(cr_path)  
    # rescale time from 0 to end of videos
    cr['time'] = cr['time'] - cr['time'][0]
    # a little bit of cleaning: name of column
    cr.rename(columns={'arousal ': 'arousal'}, inplace=True)
    # get valence and arousal ratings
    cr_val = cr['valence'].values  
    cr_aro = cr['arousal'].values  
    # compute distance and angle for ratings
    cr_d = np.hypot(cr_val, cr_aro)  # distance
    cr_ang = np.arctan2(cr_aro, cr_val)  # angle
    # for check, get number and times of samples
    nb_cr = cr.__len__()
    cr_t = cr['time'].values

    # Resampling CR data
    if resample:
        # resample cr to cr_fs (if samples are not even). Note: does not interpolate NA values
        # TODO: better way to compute length of video?
        # vid_len_from_stamp = trials_results['stop_neutral04'][1] - trials_results['start_neutral01'][1] # maybe better way to do it
        # vid_len_from_t = cr_t[-1] - cr_t[0]

        # new samples (time points) for CR
        cr_t = np.arange(1/cr_fs, vid_len + 1/cr_fs, 1/cr_fs)  # TODO: Right way to do it?
        # compute (linear) interpolation functions
        cr_v_interp = interp1d(cr['time'], cr_val, kind='linear', bounds_error=False,
                               fill_value=(cr_val[0], cr_val[-1]))
        cr_a_interp = interp1d(cr['time'], cr_aro, kind='linear', bounds_error=False,
                               fill_value=(cr_aro[0], cr_aro[-1]))
        cr_d_interp = interp1d(cr['time'], cr_d, kind='linear', bounds_error=False,
                               fill_value=(cr_d[0], cr_d[-1]))
        cr_ang_interp = interp1d(cr['time'], cr_ang, kind='linear', bounds_error=False,
                                 fill_value=(cr_ang[0], cr_ang[-1]))
        # get cr values at new samples
        cr_val = cr_v_interp(cr_t)
        cr_aro = cr_a_interp(cr_t)
        cr_d = cr_d_interp(cr_t)
        cr_ang = cr_ang_interp(cr_t)
        nb_cr = cr_t.__len__() # update number of samples
        # debug plot
        if debug:
            plt.figure
            plt.plot(cr['time'], cr['arousal'].values)
            plt.plot(cr_t, cr_aro, 'r--')
            plt.show()

    # compute time difference between samples
    cr_t_diff = np.append(np.nan, np.diff(cr_t))

    # save CR data of this trial to create dataframe
    participant = np.append(participant, np.tile(sub, nb_cr)) # add sub
    site_sub = participant_file[participant_file['sub'] == sub]['test_site'].values # get site of sub
    site = np.append(site, np.tile(site_sub, nb_cr)) # add site
    position_sub = participant_file[participant_file['sub'] == sub]['position'].values # get position of sub
    position = np.append(position, np.tile(position_sub, nb_cr)) # add position
    gender_sub = participant_file[participant_file['sub'] == sub]['gender'].values # get gender of sub
    gender = np.append(gender, np.tile(gender_sub, nb_cr)) # add gender
    cr_v = np.append(cr_v, cr_val) # add valence
    cr_a = np.append(cr_a, cr_aro) # add arousal
    cr_dist = np.append(cr_dist, cr_d) # add distance
    cr_angle = np.append(cr_angle, cr_ang) # add angle
    cr_time = np.append(cr_time, cr_t) # add time

    # to check number of samples per trial
    cr_time_diff = np.append(cr_time_diff, cr_t_diff)
    nb_samples = np.append(nb_samples, nb_cr)

    # show progress
    print(sub + " done")  

# %%
# create dataframe with CR data
d_cr = {'sub': participant, 'test_site': site, 'position': position, 'gender': gender, 'cr_v': cr_v, 'cr_a': cr_a, 'cr_dist': cr_dist, 'cr_angle': cr_angle, 'cr_time': cr_time}
df_cr = pd.DataFrame(data=d_cr)

# cut first seconds of data
if clean:
    df_cr = df_cr[df_cr['cr_time'] >= clean_s]

# save dataframe in csv file
if clean and resample:
    filename_cr = data_path + bids_folder + '/' + 'cr_rs_clean.csv'
    filename_nbs = data_path + bids_folder + '/' + 'cr_rs_clean_nb_samples.csv'
elif resample:
    filename_cr = data_path + bids_folder + '/' + 'cr_rs.csv'
    filename_nbs = data_path + bids_folder + '/' + 'cr_rs_nb_samples.csv'
elif clean:
    filename_cr = data_path + bids_folder + '/' + 'cr_clean.csv'
    filename_nbs = data_path + bids_folder + '/' + 'cr_clean_nb_samples.csv'
else:
    filename_cr = data_path + bids_folder + '/' + 'cr.csv'
    filename_nbs = data_path + bids_folder + '/' + 'cr_nb_samples.csv'
df_cr.to_csv(filename_cr, na_rep='NaN', index=False)

# %%
# create and save dataframe with number of samples per trial
d_nbs = {'sub': sub_list, 'nb_samples': nb_samples}
df_nbs = pd.DataFrame(data=d_nbs)
df_nbs.to_csv(filename_nbs, na_rep='NaN', index=False)

# %%
# check time difference between samples (framerate)
plt.figure()
plt.hist(cr_time_diff, bins=50)
plt.show()
max_t_diff = np.nanmax(cr_time_diff)
min_t_diff = np.nanmin(cr_time_diff)
median_t_diff = np.nanmedian(cr_time_diff)
median_fs = 1/median_t_diff

if debug:
    # check the framerate of right hand tracking
    plt.figure()
    plt.hist(rh_time_diff, bins=500)
    plt.show()
    max_t_diff = np.nanmax(rh_time_diff)
    min_t_diff = np.nanmin(rh_time_diff)
    median_t_diff = np.nanmedian(rh_time_diff)
    median_fs = 1 / median_t_diff
    
# %%
