########################################################################################################################
# Main script for the project
# Author: Antonin Fourcade
# Last version: 06.02.2024
########################################################################################################################

# %%
# import packages
import bids_formatting
import import_cr
import import_motion
import import_PRE_survey
import import_POST_survey
import cri_sr
import descriptive_stats
import cr_plot
import sr_plot
import motion_plot

# experiment parameters
# blocks = ['Practice', 'Experiment']
# logging_freq = ['CR', 'SR']
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2
# cri_list = ['last', 'mean', 'median', 'mode', 'max', 'min', 'std', 'cv', 'range', 'iqr', 'skew', 'kurtosis', 'auc', 'cp']
# rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'} # dictionnary for rating dimensions
# cr_fs = 1/0.05  # sampling frequency CR in Hz

# videos parameters
# scfi_len = 90  # length of the sequence of scifi videos in s
# invasion_len = 253  # length of the invasion video in s
# asteroids_len = 390  # length of the asteroids video in s
# underwood_len = 381  # length of the underwood video in s
# video_labels = ['Scifi', 'Invasion', 'Scifi', 'Asteroids', 'Scifi', 'Underwood', 'Scifi']

# surveys parameters
# PRE_nb_q = {'vr_exp': 2, 'maia': 37, 'ssq': 6, 'tas': 20} # nb of questions for each questionnaire
# POST_nb_q = {'invasiveness': 1, 'emo_rep': 1, 'presence': 2, 'sus': 7, 'satisfaction': 1} # nb of questions for each questionnaire


# %%
# Format data into BIDS
bids_formatting

# %%
# import continuous ratings (CR) data
import_cr

# %%
# import motion data
import_motion

# %%
# import survey data
# PRE-experiment
import_PRE_survey
# POST-experiment
import_POST_survey

# %%
# preprocess CR into CR indices (CRi) and import SR data
cri_sr

# %%
# descriptive statistics
descriptive_stats

# %%
# plotting
# CR
cr_plot
# SR
sr_plot
# Motion
motion_plot