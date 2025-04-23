########################################################################################################################
# Script to plot motion data from AVR experiment
# Need csv file containing motion data (from import_motion.py), and participant file (from import_PRE_questionnaire.py)
# Output: figures of motion parameters
# Author: Antonin Fourcade
# Last version: 02.02.2024
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
# blocks = ['Practice', 'Experiment']
# logging_freq = ['CR', 'SR']
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2

# videos parameters
scfi_len = 90  # length of the sequence of scifi videos in s
invasion_len = 253  # length of the invasion video in s
asteroids_len = 390  # length of the asteroids video in s
underwood_len = 381  # length of the underwood video in s
vid_len = scfi_len*4 + invasion_len + asteroids_len + underwood_len  # length of the sequence of videos in s: 4xScifi + Invasion + Asteroids + Underwood
# create list of video time boundaries
video_sequence = [scfi_len, invasion_len, scfi_len, asteroids_len, scfi_len, underwood_len, scfi_len]
video_bounds = [0]
for video_duration in video_sequence:
    video_bounds.append(video_bounds[-1] + video_duration)
video_labels = ['Scifi', 'Invasion', 'Scifi', 'Asteroids', 'Scifi', 'Underwood', 'Scifi']

motion_fs = 90  # sampling frequency in Hz
motion_start = 0 # start motion data, clean -> 5, not cleaned -> 0

debug = True # debug mode (True/False)

# dictionnaries for motion parameters and plotting
motion_params = {'pos_x': 'x', 'pos_y': 'y', 'pos_z': 'z', 'rot_x': 'roll', 'rot_y': 'pitch' , 'rot_z': 'yaw'} # motion parameters
#TODO: what are the limits for position?
lim_y = {'pos_x': [-1,1], 'pos_y': [0,2], 'pos_z': [-1,1], 'rot_x': [0, 360], 'rot_y': [0, 360], 'rot_z': [0, 360]} 

# path to motion csv file and load data
filename_hm = 'head_motion_rs.csv'
hm_path = data_path + bids_folder + '/' + filename_hm
hm_all = pd.read_csv(hm_path)

# save path
save_path = phase_path + 'affectivevr/motion_plots/'
save_ind_plot_path = save_path + 'individual_plots/'
# create folders if not there
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)
if not os.path.exists(save_ind_plot_path): # create folder if it doesn't exist
    os.makedirs(save_ind_plot_path)

# get list of participants
sub_list = hm_all['sub'].unique()

if debug:
    motion_param = 'rot_y'
    sub = 'sub-01'

# %%
## Plots mean CR time series for motion parameter
plt9.options.figure_size = (8, 24) # set default figure size in inches
# Loop through motion parameters
for motion_param in motion_params:
    # select colums of interest
    hm_param = hm_all[['sub', 'test_site', 'position', 'gender', 'time', motion_param]]
    # get mean and std across cr_time
    hm_param_stats = hm_param.groupby('time')[motion_param].agg(['mean', 'std'])

    # plot motion parameter means and stds for each time
    plot_sub = (
        plt9.ggplot()
        + plt9.geom_line(hm_param, plt9.aes(x='time', y=motion_param, group='sub', color='sub'), size=0.25, alpha=0.25)
        + plt9.geom_line(hm_param_stats, plt9.aes(x=hm_param_stats.index, y='mean'), color='black', size=0.5)
        + plt9.scale_x_continuous(breaks = [motion_start] + video_bounds + np.arange(0,vid_len+1,60).tolist(), 
                                  labels = [motion_start] + video_labels + [vid_len] + np.arange(0,vid_len+1,60).tolist(),
                                  )
        #+ plt9.geom_errorbar(cr_dim_stats, plt9.aes(ymin='mean-std', ymax='mean+std'), width=0.01, alpha=0.1)
        + plt9.geom_vline(xintercept = video_bounds, color = "grey", linetype = "dashed", size = 1)
        + plt9.ylim(lim_y[motion_param])
        + plt9.labs(x='Time (s)', y= 'mean head ' + motion_param) 
        + plt9.theme_bw(base_size=15)
        + plt9.theme(legend_position='none')
    )
    # save plot
    plot_sub.save(filename='mean_head_' + motion_param + '.png', path=save_path, dpi=300, height=20, width=60, units='cm') 

# %%
