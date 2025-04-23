########################################################################################################################
# Script to plot CRs from AVR experiment
# Need csv file containing CR data (from import_cr.py)
# Output: figures of CRs
# Author: Antonin Fourcade
# Last version: 30.10.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
import plotnine as plt9

# %%
# Set paths and experiment parameters
phase_path = 'E:/AffectiveVR/Phase_2/'
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

# dictionnary for rating dimensions and plotting
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'}
lim_y = {'valence': [-1,1], 'arousal': [-1,1], 'distance': [0, np.sqrt(2)], 'angle': [-np.pi, np.pi]} 

cr_fs = 1/0.05  # sampling frequency CR in Hz
cr_start = 5 # start cr, clean -> 5, not cleaned -> 0
debug = True # debug mode (True/False)

# position of interest
pos_select = ['seated'] # ['seated', 'standing']

# path to CR csv file and load data
filename_cr = 'cr_rs_clean.csv'
cr_path = data_path + bids_folder + '/' + filename_cr
cr_all = pd.read_csv(cr_path)

# save path
save_path = phase_path + 'affectivevr/cr_plots/'
save_ind_plot_path = save_path + 'individual_plots/'
# create folders if not there
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)
if not os.path.exists(save_ind_plot_path): # create folder if it doesn't exist
    os.makedirs(save_ind_plot_path)

# get list of participants
sub_list = cr_all['sub'].unique()

if debug:
    dim = 'valence'
    sub = 'sub-01'

# %%
## Plots mean CR time series for each dimension
plt9.options.figure_size = (8, 24) # set default figure size in inches
plt9.options.font_size = 20 # set default font size
margin = 0 # space between y-axis and plot
# Loop through rating dimensions
for dim in rat_dim:
    # select colums of interest
    selected_columns_info = cr_all[['sub', 'test_site', 'position', 'gender', 'cr_time']]
    selected_columns_dim = cr_all.filter(like='_' + rat_dim[dim], axis=1)
    # concatenate the two dataframes along the columns
    cr_dim = pd.concat([selected_columns_info, selected_columns_dim], axis=1)
    # select only position of interest
    cr_dim = cr_dim[cr_dim['position'].isin(pos_select)]
    # get mean and std across cr_time
    cr_dim_stats = cr_dim.groupby('cr_time')['cr_' + rat_dim[dim]].agg(['mean', 'std'])

    # plot cr means and stds for each cr_time
    plot_sub = (
        plt9.ggplot()
        + plt9.geom_line(cr_dim, plt9.aes(x='cr_time', y='cr_' + rat_dim[dim], group='sub', color='sub'), size=0.25, alpha=0.35)
        + plt9.geom_line(cr_dim_stats, plt9.aes(x=cr_dim_stats.index, y='mean'), color='#00BA38', size=1)
        + plt9.scale_x_continuous(breaks = [cr_start] + video_bounds + np.arange(0,vid_len+1,60).tolist(), 
                                  labels = [cr_start] + video_labels + [vid_len] + np.arange(0,vid_len+1,60).tolist(),
                                  expand = [margin, margin],
                                  #xlim = [cr_start, vid_len],
                                  )
        #+ plt9.geom_errorbar(cr_dim_stats, plt9.aes(ymin='mean-std', ymax='mean+std'), width=0.01, alpha=0.1)
        + plt9.geom_vline(xintercept = video_bounds, color = "grey", linetype = "dashed", size = 1)
        + plt9.ylim(lim_y[dim])
        + plt9.labs(x='Time (s)', y= 'mean ' + dim.capitalize()) 
        + plt9.theme_bw(base_size=30)
        + plt9.theme(legend_position='none', svg_usefonts=True, axis_text_x=plt9.element_text(rotation=45, hjust=1))
    )
    # save plot
    if len(pos_select) == 1:
        save_pos_path = save_path + pos_select[0] + '/'
    else:
        save_pos_path = save_path + 'all_positions/'
    if not os.path.exists(save_pos_path): # create folder if it doesn't exist
        os.makedirs(save_pos_path)   
    plot_sub.save(filename='mean_' + dim + '_cr.png', path=save_pos_path, dpi=300, height=20, width=60, units='cm')
    plot_sub.save(filename='mean_' + dim + '_cr.svg', path=save_pos_path, dpi=600, height=20, width=60, units='cm')  

# %%
# Individual plots
plt9.options.figure_size = (8, 24) # set default figure size in inches
# Loop over participants
for sub in sub_list:
    cr_sub = cr_all[cr_all['sub'] == sub]
    # Loop through rating dimensions
    for dim in rat_dim:
        # select colums of dimension of interest and include cr_time
        cr_dim = cr_sub.filter(like='_' + rat_dim[dim], axis=1)
        cr_dim['cr_time'] = cr_sub['cr_time']
        # plot cr against cr_time
        plot_sub = (
            plt9.ggplot(cr_dim, plt9.aes(x='cr_time', y='cr_' + rat_dim[dim]))
            + plt9.geom_line(size=0.5, color='blue')
            + plt9.scale_x_continuous(breaks = [cr_start] + video_bounds + np.arange(0,vid_len+1,60).tolist(), 
                                  labels = [cr_start] + video_labels + str(vid_len) + np.arange(0,vid_len+1,60).tolist(),
                                  )
            + plt9.geom_vline(xintercept = video_bounds, color = "grey", linetype = "dashed", size = 1)
            + plt9.ylim(lim_y[dim])
            + plt9.labs(x='Time (s)', y= dim) 
            + plt9.theme_bw(base_size=15)
            + plt9.theme(legend_position='none', svg_usefonts=True)
        )
        # save plot
        save_sub_path = save_ind_plot_path + sub + '/'
        if not os.path.exists(save_sub_path): # create folder if it doesn't exist
            os.makedirs(save_sub_path)
        plot_sub.save(filename=sub + '_' + dim + '_cr.png', path=save_sub_path, dpi=300, height=20, width=60, units='cm')


# %%
## Plots Arousal-Valence space
plt9.options.figure_size = (6, 4) # set default figure size in inches
# select colums of interest
cr_av = cr_all[['sub', 'test_site', 'position', 'gender', 'cr_time', 'cr_v', 'cr_a']]

# %%
# get mean for each participant
cr_av_mean = cr_av.groupby('sub')[['cr_v', 'cr_a']].agg(['mean'])
cr_av_mean.columns = ['_'.join(col) for col in cr_av_mean.columns.values] # rename columns
# add test_site, position and gender columns
cr_av_mean['test_site'] = cr_av.groupby('sub')['test_site'].agg(['first'])
cr_av_mean['position'] = cr_av.groupby('sub')['position'].agg(['first'])
cr_av_mean['gender'] = cr_av.groupby('sub')['gender'].agg(['first'])

# %%
# factors to group by
factors = ['test_site', 'gender']
# Loop over factors
for factor in factors:
    # get mean and std for each factor
    cr_av_stats = cr_av_mean.groupby(factor)[['cr_v_mean', 'cr_a_mean']].agg(['mean', 'std'])
    cr_av_stats.columns = ['_'.join(col) for col in cr_av_stats.columns.values] # rename columns
    # add errorbars columns: cr_v_mean_mean-cr_v_mean_std, cr_v_mean_mean+cr_v_mean_std, cr_a_mean_mean-cr_a_mean_std, cr_a_mean_mean+cr_a_mean_std
    cr_av_stats['cr_v_mean_mean-cr_v_mean_std'] = cr_av_stats['cr_v_mean_mean'] - cr_av_stats['cr_v_mean_std']
    cr_av_stats['cr_v_mean_mean+cr_v_mean_std'] = cr_av_stats['cr_v_mean_mean'] + cr_av_stats['cr_v_mean_std']
    cr_av_stats['cr_a_mean_mean-cr_a_mean_std'] = cr_av_stats['cr_a_mean_mean'] - cr_av_stats['cr_a_mean_std']
    cr_av_stats['cr_a_mean_mean+cr_a_mean_std'] = cr_av_stats['cr_a_mean_mean'] + cr_av_stats['cr_a_mean_std']

    # plot mean and std for each factor
    plot_factor = (
        plt9.ggplot(cr_av_stats, plt9.aes(x='cr_v_mean_mean', y='cr_a_mean_mean', color=cr_av_stats.index))
        + plt9.geom_point(size=2)
        + plt9.geom_errorbar(plt9.aes(ymin='cr_a_mean_mean-cr_a_mean_std', ymax='cr_a_mean_mean+cr_a_mean_std'), alpha=0.5)
        + plt9.geom_errorbarh(plt9.aes(xmin='cr_v_mean_mean-cr_v_mean_std', xmax='cr_v_mean_mean+cr_v_mean_std'), alpha=0.5)
        + plt9.labs(x='Valence', y='Arousal', color=factor)
        + plt9.theme_bw(base_size=15)
        + plt9.ylim(lim_y['arousal'])
        + plt9.xlim(lim_y['valence'])
        + plt9.theme(legend_position='right', aspect_ratio=1, svg_usefonts=True)
    )
    # save plot
    plot_factor.save(filename='av_mean_' + factor + '.png', path=save_path, dpi=300, height=20, width=25, units='cm')

# %%
# plot position in arousal-valence space
# get mean and std for each position
cr_av_stats = cr_av_mean.groupby('position')[['cr_v_mean', 'cr_a_mean']].agg(['mean', 'std'])
cr_av_stats.columns = ['_'.join(col) for col in cr_av_stats.columns.values] # rename columns
cr_av_stats = cr_av_stats.drop('standing') # remove standing position (redundant with cr_stats_site)
cr_av_stats.index = pd.MultiIndex.from_tuples([('all', 'seated')], names=['test_site', 'position']) # reformat index
# get mean and std for each position and each test_site
cr_av_stats_site = cr_av_mean.groupby(['test_site', 'position'])[['cr_v_mean', 'cr_a_mean']].agg(['mean', 'std'])
cr_av_stats_site.columns = ['_'.join(col) for col in cr_av_stats_site.columns.values] # rename columns
# concatenate cr_av_stats with cr_av_stats_site as a row with index 'position' = 'all_seated'
cr_av_stats_site = pd.concat([cr_av_stats_site, cr_av_stats], axis=0)
# add errorbars columns: cr_v_mean_mean-cr_v_mean_std, cr_v_mean_mean+cr_v_mean_std, cr_a_mean_mean-cr_a_mean_std, cr_a_mean_mean+cr_a_mean_std
cr_av_stats['cr_v_mean_mean-cr_v_mean_std'] = cr_av_stats['cr_v_mean_mean'] - cr_av_stats['cr_v_mean_std']
cr_av_stats['cr_v_mean_mean+cr_v_mean_std'] = cr_av_stats['cr_v_mean_mean'] + cr_av_stats['cr_v_mean_std']
cr_av_stats['cr_a_mean_mean-cr_a_mean_std'] = cr_av_stats['cr_a_mean_mean'] - cr_av_stats['cr_a_mean_std']
cr_av_stats['cr_a_mean_mean+cr_a_mean_std'] = cr_av_stats['cr_a_mean_mean'] + cr_av_stats['cr_a_mean_std']
# plot mean and std
plot_pos = (
    plt9.ggplot(cr_av_stats_site, plt9.aes(x='cr_v_mean_mean', y='cr_a_mean_mean', color=cr_av_stats_site.index.get_level_values('test_site'), shape=cr_av_stats_site.index.get_level_values('position')))
    + plt9.geom_point(size=2)
    + plt9.geom_errorbar(cr_av_stats_site, plt9.aes(ymin='cr_a_mean_mean-cr_a_mean_std', ymax='cr_a_mean_mean+cr_a_mean_std'), alpha=0.5)
    + plt9.geom_errorbarh(cr_av_stats_site, plt9.aes(xmin='cr_v_mean_mean-cr_v_mean_std', xmax='cr_v_mean_mean+cr_v_mean_std'), alpha=0.5)
    + plt9.labs(x='Valence', y='Arousal', color='position', shape='test_site')
    + plt9.theme_bw(base_size=15)
    + plt9.ylim(lim_y['arousal'])
    + plt9.xlim(lim_y['valence'])
    + plt9.theme(legend_position='right', aspect_ratio=1, svg_usefonts=True)
)
# save plot
plot_pos.save(filename='av_mean_position.png', path=save_path, dpi=300, height=20, width=25, units='cm')

# %%
# plot all datapoints in arousal-valence space
plot_av = (
    plt9.ggplot(cr_av, plt9.aes(x='cr_v', y='cr_a', color='sub'))
    + plt9.geom_point(size=0.5, alpha=0.5)
    + plt9.labs(x='Valence', y='Arousal')
    + plt9.theme_bw(base_size=15)
    + plt9.ylim(lim_y['arousal'])
    + plt9.xlim(lim_y['valence'])
    + plt9.theme(legend_position='none', aspect_ratio=1, svg_usefonts=True)
)
# save plot
plot_av.save(filename='av_all.png', path=save_path, dpi=300, height=20, width=20, units='cm')

# %%
# plot mean for each participant in arousal-valence space
plot_av_sub = (
    plt9.ggplot(cr_av_mean, plt9.aes(x='cr_v_mean', y='cr_a_mean', fill='gender', color='test_site', shape='position'))
    + plt9.geom_point(size=2)
    + plt9.labs(x='Valence', y='Arousal')
    + plt9.theme_bw(base_size=15)
    + plt9.ylim(lim_y['arousal'])
    + plt9.xlim(lim_y['valence'])
    + plt9.theme(legend_position='right', aspect_ratio=1, svg_usefonts=True)
)
# save plot
plot_av_sub.save(filename='av_mean_sub.png', path=save_path, dpi=300, height=20, width=25, units='cm')

# %%
## CR-SR correlation plots
