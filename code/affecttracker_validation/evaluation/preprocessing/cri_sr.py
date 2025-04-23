########################################################################################################################
# Script to preprocess CR data from AVR experiment into CR indices, and also get SR
# Need csv file containing CR data, from import_cr.py
# Output: csv file (dataframe) containing all CR indices and SR
# Author: Antonin Fourcade
# Last version: 29.01.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
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
# cri_list = ['last', 'mean', 'median', 'mode', 'max', 'min', 'std', 'cv', 'range', 'iqr', 'skew', 'kurtosis', 'auc', 'cp']

cr_fs = 1/0.05  # sampling frequency CR in Hz

# script parameters
plot = False  # plot individual CRs
debug = True  # debug mode

# path to CR csv file and load data
filename_cr = 'cr_rs_clean.csv'
cr_path = data_path + bids_folder + '/' + filename_cr
cr_all = pd.read_csv(cr_path)

# get list of participants
sub_list = cr_all['sub'].unique()

# %%
# Initialize variables
participant = []
site = []
position = []
gender = []
sr_v = []
cr_last_v = []
cr_mean_v = []
cr_median_v = []
cr_mode_v = []
cr_max_v = []
cr_min_v = []
cr_std_v = []
cr_cv_v = []
cr_range_v = []
cr_iqr_v = []
cr_skew_v = []
cr_kurtosis_v = []
cr_auc_v = []
cr_cp_v = []
sr_a = []
cr_last_a = []
cr_mean_a = []
cr_median_a = []
cr_mode_a = []
cr_max_a = []
cr_min_a = []
cr_std_a = []
cr_cv_a = []
cr_range_a = []
cr_iqr_a = []
cr_skew_a = []
cr_kurtosis_a = []
cr_auc_a = []
cr_cp_a = []
sr_dist = []
cr_last_dist = []
cr_mean_dist = []
cr_median_dist = []
cr_mode_dist = []
cr_max_dist = []
cr_min_dist = []
cr_std_dist = []
cr_cv_dist = []
cr_range_dist = []
cr_iqr_dist = []
cr_skew_dist = []
cr_kurtosis_dist = []
cr_auc_dist = []
cr_cp_dist = []
sr_angle = []
cr_last_angle = []
cr_mean_angle = []
cr_median_angle = []
cr_mode_angle = []
cr_max_angle = []
cr_min_angle = []
cr_std_angle = []
cr_cv_angle = []
cr_range_angle = []
cr_iqr_angle = []
cr_skew_angle = []
cr_kurtosis_angle = []
cr_auc_angle = []
cr_cp_angle = []

# debug
if debug:
    sub = 'sub-01'

# %%
# Read and preprocess data
# Loop over participants
for sub_idx, sub in enumerate(sub_list):
    # Read results in csv file
    sub_path = data_path + bids_folder + '/' + sub + '/'
    trial_results_filename = sub_path + 'trial_results.csv'
    trials_results = pd.read_csv(trial_results_filename)

    # SR
    # Get SR from trial_results
    sr_sub = trials_results[trials_results['block_name'] == 'Experiment'].loc[:, ['SR_valence', 'SR_arousal']]
    # Compute SR distance and angle
    sr_sub_dist = np.hypot(sr_sub['SR_valence'].values[0], sr_sub['SR_arousal'].values[0])
    sr_sub_angle = np.arctan2(sr_sub['SR_arousal'].values[0], sr_sub['SR_valence'].values[0])
    # Append SR to list to create dataframe later
    sr_v = np.append(sr_v, sr_sub['SR_valence'].values[0])
    sr_a = np.append(sr_a, sr_sub['SR_arousal'].values[0])
    sr_dist = np.append(sr_dist, sr_sub_dist)
    sr_angle = np.append(sr_angle, sr_sub_angle)

    # CR
    # Get CR from cr_all
    cr_sub = cr_all[cr_all['sub'] == sub]
    # Compute CR indices (CRi) and append to list to create dataframe later
    #TODO: check how to deal with NaNs
    # for now -> omit NaNs
    # last rating
    cr_sub_last_v = cr_sub['cr_v'][~cr_sub['cr_v'].isna()].iloc[-1]
    cr_sub_last_a = cr_sub['cr_a'][~cr_sub['cr_a'].isna()].iloc[-1]
    cr_sub_last_dist = cr_sub['cr_dist'][~cr_sub['cr_dist'].isna()].iloc[-1]
    cr_sub_last_angle = cr_sub['cr_angle'][~cr_sub['cr_angle'].isna()].iloc[-1]
    cr_last_v = np.append(cr_last_v, cr_sub_last_v)
    cr_last_a = np.append(cr_last_a, cr_sub_last_a)
    cr_last_dist = np.append(cr_last_dist, cr_sub_last_dist)
    cr_last_angle = np.append(cr_last_angle, cr_sub_last_angle)
    # mean
    cr_sub_mean_v = cr_sub['cr_v'].mean(skipna=True)
    cr_sub_mean_a = cr_sub['cr_a'].mean(skipna=True)
    cr_sub_mean_dist = np.nanmean(cr_sub['cr_dist'])
    cr_sub_mean_angle = np.nanmean(cr_sub['cr_angle'])
    cr_mean_v = np.append(cr_mean_v, cr_sub_mean_v)
    cr_mean_a = np.append(cr_mean_a, cr_sub_mean_a)
    cr_mean_dist = np.append(cr_mean_dist, cr_sub_mean_dist)
    cr_mean_angle = np.append(cr_mean_angle, cr_sub_mean_angle)
    # median
    cr_sub_median_v = cr_sub['cr_v'].median(skipna=True)
    cr_sub_median_a = cr_sub['cr_a'].median(skipna=True)
    cr_sub_median_dist = cr_sub['cr_dist'].median(skipna=True)
    cr_sub_median_angle = cr_sub['cr_angle'].median(skipna=True)
    cr_median_v = np.append(cr_median_v, cr_sub_median_v)
    cr_median_a = np.append(cr_median_a, cr_sub_median_a)
    cr_median_dist = np.append(cr_median_dist, cr_sub_median_dist)
    cr_median_angle = np.append(cr_median_angle, cr_sub_median_angle)
    # mode
    cr_sub_mode_v = cr_sub['cr_v'].mode(dropna=True)[0]
    cr_sub_mode_a = cr_sub['cr_a'].mode(dropna=True)[0]
    cr_sub_mode_dist = cr_sub['cr_dist'].mode(dropna=True)[0]
    cr_sub_mode_angle = cr_sub['cr_angle'].mode(dropna=True)[0]
    cr_mode_v = np.append(cr_mode_v, cr_sub_mode_v)
    cr_mode_a = np.append(cr_mode_a, cr_sub_mode_a)
    cr_mode_dist = np.append(cr_mode_dist, cr_sub_mode_dist)
    cr_mode_angle = np.append(cr_mode_angle, cr_sub_mode_angle)
    # max
    cr_sub_max_v = cr_sub['cr_v'].max(skipna=True)
    cr_sub_max_a = cr_sub['cr_a'].max(skipna=True)
    cr_sub_max_dist = cr_sub['cr_dist'].max(skipna=True)
    cr_sub_max_angle = cr_sub['cr_angle'].max(skipna=True)
    cr_max_v = np.append(cr_max_v, cr_sub_max_v)
    cr_max_a = np.append(cr_max_a, cr_sub_max_a)
    cr_max_dist = np.append(cr_max_dist, cr_sub_max_dist)
    cr_max_angle = np.append(cr_max_angle, cr_sub_max_angle)
    # min
    cr_sub_min_v = cr_sub['cr_v'].min(skipna=True)
    cr_sub_min_a = cr_sub['cr_a'].min(skipna=True)
    cr_sub_min_dist = cr_sub['cr_dist'].min(skipna=True)
    cr_sub_min_angle = cr_sub['cr_angle'].min(skipna=True)
    cr_min_v = np.append(cr_min_v, cr_sub_min_v)
    cr_min_a = np.append(cr_min_a, cr_sub_min_a)
    cr_min_dist = np.append(cr_min_dist, cr_sub_min_dist)
    cr_min_angle = np.append(cr_min_angle, cr_sub_min_angle)
    # std
    cr_sub_std_v = cr_sub['cr_v'].std(skipna=True)
    cr_sub_std_a = cr_sub['cr_a'].std(skipna=True)
    cr_sub_std_dist = cr_sub['cr_dist'].std(skipna=True)
    cr_sub_std_angle = cr_sub['cr_angle'].std(skipna=True)
    cr_std_v = np.append(cr_std_v, cr_sub_std_v)
    cr_std_a = np.append(cr_std_a, cr_sub_std_a)
    cr_std_dist = np.append(cr_std_dist, cr_sub_std_dist)
    cr_std_angle = np.append(cr_std_angle, cr_sub_std_angle)
    # cv: std/|mean|
    cr_sub_cv_v = cr_sub_std_v/np.fabs(cr_sub_mean_v)
    cr_sub_cv_a = cr_sub_std_a/np.fabs(cr_sub_mean_a)
    cr_sub_cv_dist = cr_sub_std_dist/np.fabs(cr_sub_mean_dist)
    cr_sub_cv_angle = cr_sub_std_angle/np.fabs(cr_sub_mean_angle)
    cr_cv_v = np.append(cr_cv_v, cr_sub_cv_v)
    cr_cv_a = np.append(cr_cv_a, cr_sub_cv_a)
    cr_cv_dist = np.append(cr_cv_dist, cr_sub_cv_dist)
    cr_cv_angle = np.append(cr_cv_angle, cr_sub_cv_angle)
    # range
    cr_sub_range_v = cr_sub_max_v - cr_sub_min_v
    cr_sub_range_a = cr_sub_max_a - cr_sub_min_a
    cr_sub_range_dist = cr_sub_max_dist - cr_sub_min_dist
    cr_sub_range_angle = cr_sub_max_angle - cr_sub_min_angle
    cr_range_v = np.append(cr_range_v, cr_sub_range_v)
    cr_range_a = np.append(cr_range_a, cr_sub_range_a)
    cr_range_dist = np.append(cr_range_dist, cr_sub_range_dist)
    cr_range_angle = np.append(cr_range_angle, cr_sub_range_angle)
    # interquartile range (iqr): Q3 -Q1
    cr_sub_iqr_v = stats.iqr(cr_sub['cr_v'].values, nan_policy='omit')
    cr_sub_iqr_a = stats.iqr(cr_sub['cr_a'].values, nan_policy='omit')
    cr_sub_iqr_dist = stats.iqr(cr_sub['cr_dist'], nan_policy='omit')
    cr_sub_iqr_angle = stats.iqr(cr_sub['cr_angle'], nan_policy='omit')
    cr_iqr_v = np.append(cr_iqr_v, cr_sub_iqr_v)
    cr_iqr_a = np.append(cr_iqr_a, cr_sub_iqr_a)
    cr_iqr_dist = np.append(cr_iqr_dist, cr_sub_iqr_dist)
    cr_iqr_angle = np.append(cr_iqr_angle, cr_sub_iqr_angle)
    # skewness
    cr_sub_skew_v = cr_sub['cr_v'].skew(skipna=True)
    cr_sub_skew_a = cr_sub['cr_a'].skew(skipna=True)
    cr_sub_skew_dist = cr_sub['cr_dist'].skew(skipna=True)
    cr_sub_skew_angle = cr_sub['cr_angle'].skew(skipna=True)
    cr_skew_v = np.append(cr_skew_v, cr_sub_skew_v)
    cr_skew_a = np.append(cr_skew_a, cr_sub_skew_a)
    cr_skew_dist = np.append(cr_skew_dist, cr_sub_skew_dist)
    cr_skew_angle = np.append(cr_skew_angle, cr_sub_skew_angle)
    # kurtosis
    cr_sub_kurtosis_v = cr_sub['cr_v'].kurtosis(skipna=True)
    cr_sub_kurtosis_a = cr_sub['cr_a'].kurtosis(skipna=True)
    cr_sub_kurtosis_dist = cr_sub['cr_dist'].kurtosis(skipna=True)
    cr_sub_kurtosis_angle = cr_sub['cr_angle'].kurtosis(skipna=True)
    cr_kurtosis_v = np.append(cr_kurtosis_v, cr_sub_kurtosis_v)
    cr_kurtosis_a = np.append(cr_kurtosis_a, cr_sub_kurtosis_a)
    cr_kurtosis_dist = np.append(cr_kurtosis_dist, cr_sub_kurtosis_dist)
    cr_kurtosis_angle = np.append(cr_kurtosis_angle, cr_sub_kurtosis_angle)
    # area under curve (auc)
    cr_sub_auc_v = np.trapz(cr_sub['cr_v'].dropna(), dx=1/cr_fs)
    cr_sub_auc_a = np.trapz(cr_sub['cr_a'].dropna(), dx=1/cr_fs)
    cr_sub_auc_dist = np.trapz(cr_sub['cr_dist'][~np.isnan(cr_sub['cr_dist'])], dx=1/cr_fs)
    cr_sub_auc_angle = np.trapz(cr_sub['cr_angle'][~np.isnan(cr_sub['cr_angle'])], dx=1/cr_fs)
    cr_auc_v = np.append(cr_auc_v, cr_sub_auc_v)
    cr_auc_a = np.append(cr_auc_a, cr_sub_auc_a)
    cr_auc_dist = np.append(cr_auc_dist, cr_sub_auc_dist)
    cr_auc_angle = np.append(cr_auc_angle, cr_sub_auc_angle)
    # cumulative power (cp)
    # OK to drop nans? or would be better to interpolate?
    f_v, psd_v = signal.welch(cr_sub['cr_v'].dropna(), cr_fs)
    f_a, psd_a = signal.welch(cr_sub['cr_a'].dropna(), cr_fs)
    f_dist, psd_dist = signal.welch(cr_sub['cr_dist'][~np.isnan(cr_sub['cr_dist'])], cr_fs)
    f_angle, psd_angle = signal.welch(cr_sub['cr_angle'][~np.isnan(cr_sub['cr_angle'])], cr_fs)
    cr_sub_cp_v = np.trapz(psd_v, f_v, dx=f_v[1])
    cr_sub_cp_a = np.trapz(psd_a, f_a, dx=f_a[1])
    cr_sub_cp_dist = np.trapz(psd_dist, f_dist, dx=f_dist[1])
    cr_sub_cp_angle = np.trapz(psd_angle, f_angle, dx=f_angle[1])
    cr_cp_v = np.append(cr_cp_v, cr_sub_cp_v)
    cr_cp_a = np.append(cr_cp_a, cr_sub_cp_a)
    cr_cp_dist = np.append(cr_cp_dist, cr_sub_cp_dist)
    cr_cp_angle = np.append(cr_cp_angle, cr_sub_cp_angle)
    
    if plot:
        plt.figure()
        plt.semilogy(f_v, psd_v, label='valence')
        plt.semilogy(f_a, psd_a, label='arousal')
        plt.semilogy(f_dist, psd_dist, label='distance')
        plt.semilogy(f_angle, psd_angle, label='angle')
        plt.legend()
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [V**2/Hz]')
        plt.title(sub)
        fig_path = sub_path + 'figures/'
        if not os.path.exists(fig_path):
            os.makedirs(fig_path)
        fig_name = sub + '_VA_psd.png'
        plt.savefig(fig_path + fig_name)
        plt.close()
    
    # Append participant info to list to create dataframe later
    participant = np.append(participant, sub)  # add sub
    site = np.append(site, cr_sub['test_site'].unique())  # add site
    position = np.append(position, cr_sub['position'].unique())  # add position
    gender = np.append(gender, cr_sub['gender'].unique())

    # print progress
    print(sub + " done")

# create dataframe with all CR indices
d_cri_sr = {'sub': participant, 'test_site': site, 'position': position, 'gender': gender,
     'sr_v': sr_v, 'cr_last_v': cr_last_v, 'cr_mean_v': cr_mean_v, 'cr_median_v': cr_median_v, 'cr_mode_v': cr_mode_v, 'cr_max_v': cr_max_v, 'cr_min_v': cr_min_v,
     'cr_std_v': cr_std_v, 'cr_cv_v': cr_cv_v, 'cr_range_v': cr_range_v, 'cr_iqr_v': cr_iqr_v, 'cr_skew_v': cr_skew_v, 'cr_kurtosis_v': cr_kurtosis_v, 'cr_auc_v': cr_auc_v, 'cr_cp_v': cr_cp_v,
     'sr_a': sr_a, 'cr_last_a': cr_last_a, 'cr_mean_a': cr_mean_a, 'cr_median_a': cr_median_a, 'cr_mode_a': cr_mode_a, 'cr_max_a': cr_max_a, 'cr_min_a': cr_min_a,
     'cr_std_a': cr_std_a, 'cr_cv_a': cr_cv_a, 'cr_range_a': cr_range_a, 'cr_iqr_a': cr_iqr_a, 'cr_skew_a': cr_skew_a, 'cr_kurtosis_a': cr_kurtosis_a, 'cr_auc_a': cr_auc_a, 'cr_cp_a': cr_cp_a,
     'sr_dist': sr_dist, 'cr_last_dist': cr_last_dist, 'cr_mean_dist': cr_mean_dist, 'cr_median_dist': cr_median_dist, 'cr_mode_dist': cr_mode_dist, 'cr_max_dist': cr_max_dist, 'cr_min_dist': cr_min_dist,
     'cr_std_dist': cr_std_dist, 'cr_cv_dist': cr_cv_dist, 'cr_range_dist': cr_range_dist, 'cr_iqr_dist': cr_iqr_dist, 'cr_skew_dist': cr_skew_dist, 'cr_kurtosis_dist': cr_kurtosis_dist, 'cr_auc_dist': cr_auc_dist, 'cr_cp_dist': cr_cp_dist,
     'sr_angle': sr_angle, 'cr_last_angle': cr_last_angle, 'cr_mean_angle': cr_mean_angle, 'cr_median_angle': cr_median_angle, 'cr_mode_angle': cr_mode_angle, 'cr_max_angle': cr_max_angle, 'cr_min_angle': cr_min_angle,
     'cr_std_angle': cr_std_angle, 'cr_cv_angle': cr_cv_angle, 'cr_range_angle': cr_range_angle, 'cr_iqr_angle': cr_iqr_angle, 'cr_skew_angle': cr_skew_angle, 'cr_kurtosis_angle': cr_kurtosis_angle, 'cr_auc_angle': cr_auc_angle, 'cr_cp_angle': cr_cp_angle
     }
df_cri_sr = pd.DataFrame(data=d_cri_sr)
# save dataframe
filename = data_path + bids_folder + '/cri_sr_clean.csv'
df_cri_sr.to_csv(filename, na_rep='NaN', index=False)
