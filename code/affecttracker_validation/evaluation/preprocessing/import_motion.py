########################################################################################################################
# Script to import motion data from AVR experiment
# Preprocessing steps (optional): resampling and cut first seconds of data
# Output: 2 csv files (dataframes): cr[_rs_clean].csv and cr[_rs_clean]_nb_samples.csv
# Author: Antonin Fourcade
# Last version: 02.02.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy import signal
from scipy import stats

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

motion_fs = 90  # sampling frequency in Hz
motion_params = {'pos_x': 'x', 'pos_y': 'y', 'pos_z': 'z', 'rot_x': 'roll', 'rot_y': 'pitch' , 'rot_z': 'yaw'} # motion parameters

vid_len = 90*4 + 253 + 390 + 381  # length of the sequence of videos in s: 4xScifi + Invasion + Asteroids + Underwood

# preprocessing parameters
resample = True  # resample CR to motion_fs (if samples are not even). Note: does not interpolate NA values
new_fs = motion_fs/3
clean = False  # clean CR: remove first clean_s seconds of CR
clean_s = 5  # seconds to remove from CR
debug = False  # debug mode

# participant file
participant_file_path = 'participants.csv'
participant_file = pd.read_csv(data_path + bids_folder + '/' + participant_file_path, sep=',', header=0)

# get list of participants
sub_list = participant_file['sub'].values

# %%
# Initialize variables
xrcam_all = pd.DataFrame()

if debug:
    xrcam_time_diff = []
    sub = 'sub-01'
    sub_idx = 0

# %%
# Read and preprocess data
for sub_idx, sub in enumerate(sub_list):

    # Read results in csv file
    sub_path = data_path + bids_folder + '/' + sub + '/'
    trial_results_filename = sub_path + 'trial_results.csv'
    trials_results = pd.read_csv(trial_results_filename)

    # Select Experiment data
    trials_exp = trials_results[trials_results['block_name'] == 'Experiment']

    # get motion data from file
    # location of head movement file
    xrcam_loc = 'xrcam_movement_location_0'
    xrcam_path = data_path + trials_exp[xrcam_loc].item()
    xrcam = pd.read_csv(xrcam_path)
    # get start time of the first video
    start_vid_time = trials_exp['start_neutral01'].values[0]
    # cut data for duration of video sequence
    xrcam = xrcam[xrcam['time'] >= start_vid_time]
    xrcam = xrcam[xrcam['time'] <= start_vid_time + vid_len]
    #reset index
    xrcam = xrcam.reset_index(drop=True)
    # rescale time from 0 to end of videos
    xrcam['time'] = xrcam['time'] - xrcam['time'][0]
    # for check, get number and times of samples
    nb_samples = xrcam.__len__()
    xrcam_t = xrcam['time'].values

    if resample:
        # resample cr to motion_fs (if samples are not even). Note: does not interpolate NA values
        # TODO: better way to compute length of video?
        # vid_len_from_stamp = trials_results['stop_neutral04'][1] - trials_results['start_neutral01'][1] # maybe better way to do it
        # vid_len_from_t = cr_t[-1] - cr_t[0]

        #TODO: downsample?

        # new samples (time points) for CR
        xrcam_t = np.arange(1/new_fs, vid_len + 1/new_fs, 1/new_fs)  # TODO: Right way to do it?
        # update number of samples
        nb_samples = xrcam_t.__len__()
        # Compute (linear) interpolation functions
        # Create a dictionary to store the interpolation functions
        interp_funcs = {}
        xrcam_rs = {}
        # Loop over the motion_params to compute the interpolation functions and resample the data
        for motion_param in motion_params:
            interp_funcs[motion_param] = interp1d(xrcam['time'], xrcam[motion_param], kind='linear', bounds_error=False,
                                        fill_value=(xrcam[motion_param][0], xrcam[motion_param][len(xrcam[motion_param]) - 1]))
            # Get xrcam values at new samples
            #TODO: use roll, pitch, yaw instead of rot_x, rot_y, rot_z? -> xrcam_rs[motion_params[motion_param]]
            xrcam_rs[motion_param] = interp_funcs[motion_param](xrcam_t)
        
        # debug plot
        if debug:
            plt.figure
            plt.plot(xrcam['time'], xrcam['rot_y'].values)
            plt.plot(xrcam_t, xrcam_rs['rot_y'], 'r--')
            plt.show()
        
        # update xrcam with resampled data
        xrcam = pd.DataFrame(data=xrcam_rs)
        # complete dataframe for this sub
        participant = np.tile(sub, nb_samples)
        site = np.tile(participant_file['test_site'][sub_idx], nb_samples)
        position = np.tile(participant_file['position'][sub_idx], nb_samples)
        gender = np.tile(participant_file['gender'][sub_idx], nb_samples)
        df_to_add = pd.DataFrame({'sub': participant, 'test_site': site, 'position': position, 'gender': gender, 'time': xrcam_t})
        xrcam = pd.concat([df_to_add, xrcam], axis=1)
        # concatenate xrcam to xrcam_all
        xrcam_all = pd.concat([xrcam_all, xrcam], axis=0)

    # show progress
    print(sub + " done")  
        
# %%
# cut first seconds of data
if clean:
    xrcam_all = xrcam_all[xrcam_all['time'] >= clean_s]

# save dataframe in csv file
if clean and resample:
    filename_xrcam = data_path + bids_folder + '/' + 'head_motion_rs_clean.csv'
elif resample:
    filename_xrcam = data_path + bids_folder + '/' + 'head_motion_rs.csv'
elif clean:
    filename_xrcam = data_path + bids_folder + '/' + 'head_motion_clean.csv'
else:
    filename_xrcam = data_path + bids_folder + '/' + 'head_motion.csv'
xrcam_all.to_csv(filename_xrcam, na_rep='NaN', index=False)
# %%
