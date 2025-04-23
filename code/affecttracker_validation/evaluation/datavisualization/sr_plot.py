########################################################################################################################
# Script to plot SRs from AVR experiment
# Need csv file containing SR data (from cri_sr.py)
# Output: figures of SRs
# Author: Antonin Fourcade
# Last version: 05.02.2024
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

# dictionnary for rating dimensions and plotting
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'}
lim_y = {'valence': [-1,1], 'arousal': [-1,1], 'distance': [0, np.sqrt(2)], 'angle': [-np.pi, np.pi]} 

debug = True # debug mode (True/False)

# path to SR csv file and load data
filename_sr = 'cri_sr_clean.csv'
sr_path = data_path + bids_folder + '/' + filename_sr
sr_all = pd.read_csv(sr_path)
# keep only SR data
sr_all = sr_all[['sub', 'test_site', 'position', 'gender', 'sr_v', 'sr_a', 'sr_dist', 'sr_angle']]

# save path
save_path = phase_path + 'affectivevr/sr_plots/'
# create folders if not there
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# %%
## Plots Arousal-Valence space
plt9.options.figure_size = (6, 4) # set default figure size in inches
# select colums of interest
sr_av = sr_all[['sub', 'test_site', 'position', 'gender', 'sr_v', 'sr_a']]

# %%
# plot SR for each participant in arousal-valence space
sr_av_mean = sr_av.groupby('sub')[['sr_v', 'sr_a']].agg(['mean'])
sr_av_mean.columns = ['_'.join(col) for col in sr_av_mean.columns.values] # rename columns
sr_av_mean = sr_av_mean.reset_index()
# add test_site, position and gender columns with concatenation
sr_av_mean = sr_av_mean.merge(sr_av[['test_site', 'position', 'gender']], left_index=True, right_index=True)

# plot
plot_av_sub = (
    plt9.ggplot(sr_av_mean, plt9.aes(x='sr_v_mean', y='sr_a_mean', fill='gender', color='test_site', shape='position'))
    + plt9.geom_point(size=2)
    + plt9.labs(x='Valence', y='Arousal')
    + plt9.theme_bw(base_size=15)
    + plt9.ylim(lim_y['arousal'])
    + plt9.xlim(lim_y['valence'])
    + plt9.theme(legend_position='right', aspect_ratio=1)
)
# save plot
plot_av_sub.save(filename='av_sub.png', path=save_path, dpi=300, height=20, width=25, units='cm')


# %%
# factors to group by
factors = ['test_site', 'gender']
# Loop over factors
for factor in factors:
    sr_av_stats = sr_av.groupby(factor)[['sr_v', 'sr_a']].agg(['mean', 'std'])
    sr_av_stats.columns = ['_'.join(col) for col in sr_av_stats.columns.values] # rename columns
    # add errorbars columns: sr_v_mean-sr_v_std, sr_v_mean+sr_v_std, sr_a_mean-sr_a_std, sr_a_mean+sr_a_std
    sr_av_stats['sr_v_mean-sr_v_std'] = sr_av_stats['sr_v_mean'] - sr_av_stats['sr_v_std']
    sr_av_stats['sr_v_mean+sr_v_std'] = sr_av_stats['sr_v_mean'] + sr_av_stats['sr_v_std']
    sr_av_stats['sr_a_mean-sr_a_std'] = sr_av_stats['sr_a_mean'] - sr_av_stats['sr_a_std']
    sr_av_stats['sr_a_mean+sr_a_std'] = sr_av_stats['sr_a_mean'] + sr_av_stats['sr_a_std']

    # plot mean and std for each factor
    plot_factor = (
        plt9.ggplot(sr_av_stats, plt9.aes(x='sr_v_mean', y='sr_a_mean', color=sr_av_stats.index))
        + plt9.geom_point(size=2)
        + plt9.geom_errorbar(plt9.aes(ymin='sr_a_mean-sr_a_std', ymax='sr_a_mean+sr_a_std'), alpha=0.5)
        + plt9.geom_errorbarh(plt9.aes(xmin='sr_v_mean-sr_v_std', xmax='sr_v_mean+sr_v_std'), alpha=0.5)
        + plt9.labs(x='Valence', y='Arousal', color=factor)
        + plt9.theme_bw(base_size=15)
        + plt9.ylim(lim_y['arousal'])
        + plt9.xlim(lim_y['valence'])
        + plt9.theme(legend_position='right', aspect_ratio=1)
    )
    # save plot
    plot_factor.save(filename='av_mean_' + factor + '.png', path=save_path, dpi=300, height=20, width=25, units='cm')

# %%
# plot position in arousal-valence space
# get mean and std for each position
sr_av_stats = sr_av.groupby('position')[['sr_v', 'sr_a']].agg(['mean', 'std'])
sr_av_stats.columns = ['_'.join(col) for col in sr_av_stats.columns.values] # rename columns
sr_av_stats = sr_av_stats.drop('standing') # remove standing position (redundant with sr_stats_site)
sr_av_stats.index = pd.MultiIndex.from_tuples([('all', 'seated')], names=['test_site', 'position']) # reformat index
# get mean and std for each position and each test_site
sr_av_stats_site = sr_av.groupby(['test_site', 'position'])[['sr_v', 'sr_a']].agg(['mean', 'std'])
sr_av_stats_site.columns = ['_'.join(col) for col in sr_av_stats_site.columns.values] # rename columns
# concatenate sr_av_stats with sr_av_stats_site as a row with index 'position' = 'mean_seated'
sr_av_stats_site = pd.concat([sr_av_stats_site, sr_av_stats], axis=0)
# add errorbars columns: sr_v_mean-sr_v_std, sr_v_mean+sr_v_std, sr_a_mean-sr_a_std, sr_a_mean+sr_a_std
sr_av_stats_site['sr_v_mean-sr_v_std'] = sr_av_stats_site['sr_v_mean'] - sr_av_stats_site['sr_v_std']
sr_av_stats_site['sr_v_mean+sr_v_std'] = sr_av_stats_site['sr_v_mean'] + sr_av_stats_site['sr_v_std']
sr_av_stats_site['sr_a_mean-sr_a_std'] = sr_av_stats_site['sr_a_mean'] - sr_av_stats_site['sr_a_std']
sr_av_stats_site['sr_a_mean+sr_a_std'] = sr_av_stats_site['sr_a_mean'] + sr_av_stats_site['sr_a_std']
# plot mean and std
plot_pos = (
    plt9.ggplot(sr_av_stats_site, plt9.aes(x='sr_v_mean', y='sr_a_mean', color=sr_av_stats_site.index.get_level_values('test_site'), shape=sr_av_stats_site.index.get_level_values('position')))
    + plt9.geom_point(size=2)
    + plt9.geom_errorbar(sr_av_stats_site, plt9.aes(ymin='sr_a_mean-sr_a_std', ymax='sr_a_mean+sr_a_std'), alpha=0.5)
    + plt9.geom_errorbarh(sr_av_stats_site, plt9.aes(xmin='sr_v_mean-sr_v_std', xmax='sr_v_mean+sr_v_std'), alpha=0.5)
    + plt9.labs(x='Valence', y='Arousal', color='position', shape='test_site')
    + plt9.theme_bw(base_size=15)
    + plt9.ylim(lim_y['arousal'])
    + plt9.xlim(lim_y['valence'])
    + plt9.theme(legend_position='right', aspect_ratio=1)
)
# save plot
plot_pos.save(filename='av_mean_position.png', path=save_path, dpi=300, height=20, width=25, units='cm')
