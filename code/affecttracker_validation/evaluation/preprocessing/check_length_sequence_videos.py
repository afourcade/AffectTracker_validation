########################################################################################################################
# Script to check the length of the sequence of videos and the framerate of the CR data
# Output: ?
# Author: Antonin Fourcade
# Last version: 25.03.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import plotnine as plt9

#%matplotlib qt # plot in a separate window

# %%
# Set paths and experiment parameters
data_path = 'E:/AffectiveVR/Phase_2/Data/'
bids_folder = 'AVR' # folder with data in bids format

# experiment parameters
# blocks = ['Practice', 'Experiment']
# logging_freq = ['CR', 'SR']
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2

cr_fs = 1/0.05  # sampling frequency CR in Hz
vid_len = 90*4 + 253 + 390 + 381  # length of the sequence of videos in s: 4xScifi + Invasion + Asteroids + Underwood

# participant file
participant_file_path = 'participants.csv'
participant_file = pd.read_csv(data_path + bids_folder + '/' + participant_file_path, sep=',', header=0)

# get list of participants
sub_list = participant_file['sub'].values

# save path
save_path = 'E:/AffectiveVR/Phase_2/affectivevr/checks/'
if not os.path.exists(save_path):
    os.makedirs(save_path)

plt9.options.figure_size = (6, 6) # set default figure size in inches
plt9.options.dpi = 300 # set default dpi

# %%
# Initialize variables
rh_time_diff = []
sub = []
site = []
d = pd.DataFrame({'sub': [], 'test_site': [], 'exp_dur': [], 'rh_fs_mean': [], 
                  'neutral01_dur': [], 'neutral02_dur': [], 'neutral03_dur': [], 'neutral04_dur': [], 'invasion_dur': [], 'asteroids_dur': [], 'underwood_dur': [], 
                  'delay_neutral01_invasion': [], 'delay_invasion_neutral02': [], 'delay_neutral02_asteroids': [], 'delay_asteroids_neutral03': [], 'delay_neutral03_underwood': [], 'delay_underwood_neutral04': []})
#sub = 'sub-01'

# %%
# Read and preprocess data
for sub_idx, sub in enumerate(sub_list):

    # Read results in csv file
    sub_path = data_path + bids_folder + '/' + sub + '/'
    trial_results_filename = sub_path + 'trial_results.csv'
    trials_results = pd.read_csv(trial_results_filename)

    # Select Experiment data
    trials_exp = trials_results[trials_results['block_name'] == 'Experiment']

    # check sampling of right hand tracking (90Hz)
    rh_loc = 'righthand_movement_location_0'
    rh_path = data_path + trials_exp[rh_loc].item()
    rh = pd.read_csv(rh_path)
    rh_t = rh['time'].values
    rh_t_diff = np.append(np.nan, np.diff(rh_t))

    # get duration of the experiment and videos
    exp_dur = trials_exp['stop_neutral04'].values - trials_exp['start_neutral01'].values
    neutral01_dur = trials_exp['stop_neutral01'].values - trials_exp['start_neutral01'].values
    neutral02_dur = trials_exp['stop_neutral02'].values - trials_exp['start_neutral02'].values
    neutral03_dur = trials_exp['stop_neutral03'].values - trials_exp['start_neutral03'].values
    neutral04_dur = trials_exp['stop_neutral04'].values - trials_exp['start_neutral04'].values
    invasion_dur = trials_exp['stop_invasion'].values - trials_exp['start_invasion'].values
    asteroids_dur = trials_exp['stop_asteroids'].values - trials_exp['start_asteroids'].values
    underwood_dur = trials_exp['stop_underwood'].values - trials_exp['start_underwood'].values
    
    # get delay between videos
    delay_neutral01_invasion = trials_exp['start_invasion'].values - trials_exp['stop_neutral01'].values
    delay_invasion_neutral02 = trials_exp['start_neutral02'].values - trials_exp['stop_invasion'].values
    delay_neutral02_asteroids = trials_exp['start_asteroids'].values - trials_exp['stop_neutral02'].values
    delay_asteroids_neutral03 = trials_exp['start_neutral03'].values - trials_exp['stop_asteroids'].values
    delay_neutral03_underwood = trials_exp['start_underwood'].values - trials_exp['stop_neutral03'].values
    delay_underwood_neutral04 = trials_exp['start_neutral04'].values - trials_exp['stop_underwood'].values


    # save in dataframe d
    d.loc[sub_idx, 'sub'] = sub
    d.loc[sub_idx, 'test_site'] = participant_file['test_site'][sub_idx]
    d.loc[sub_idx, 'exp_dur'] = exp_dur
    d.loc[sub_idx, 'rh_fs_mean'] = 1 / np.nanmean(rh_t_diff)
    d.loc[sub_idx, 'neutral01_dur'] = neutral01_dur
    d.loc[sub_idx, 'neutral02_dur'] = neutral02_dur
    d.loc[sub_idx, 'neutral03_dur'] = neutral03_dur
    d.loc[sub_idx, 'neutral04_dur'] = neutral04_dur
    d.loc[sub_idx, 'invasion_dur'] = invasion_dur
    d.loc[sub_idx, 'asteroids_dur'] = asteroids_dur
    d.loc[sub_idx, 'underwood_dur'] = underwood_dur
    d.loc[sub_idx, 'delay_neutral01_invasion'] = delay_neutral01_invasion
    d.loc[sub_idx, 'delay_invasion_neutral02'] = delay_invasion_neutral02
    d.loc[sub_idx, 'delay_neutral02_asteroids'] = delay_neutral02_asteroids
    d.loc[sub_idx, 'delay_asteroids_neutral03'] = delay_asteroids_neutral03
    d.loc[sub_idx, 'delay_neutral03_underwood'] = delay_neutral03_underwood
    d.loc[sub_idx, 'delay_underwood_neutral04'] = delay_underwood_neutral04

    # show progress
    print(sub + " done")  

# %%
# save dataframe d
d.to_csv(save_path + 'exp_dur.csv', index=False)

# %%
# for plotting, create sub-id column by converting 'sub' to category and map categories to integers (starting from 1)
d.loc[:,'sub-id'] = d.loc[:,'sub'].astype('category').cat.codes + 1

# timings of interest
timings = ['exp_dur', 'neutral01_dur', 'neutral02_dur', 'neutral03_dur', 'neutral04_dur', 'invasion_dur', 'asteroids_dur', 'underwood_dur']

# plot timings of interest
for timing in timings:
    p = (
    plt9.ggplot(d, plt9.aes(x='test_site', y=timing))
    # Add the violin plot layer
    + plt9.geom_violin(plt9.aes(fill='test_site'), position=plt9.position_nudge(x=0.2), style='right', alpha=0.5)
    # Add the points layer
    + plt9.geom_point(plt9.aes(color= 'sub-id'), position=plt9.position_jitter(width=0.1, height=0), size=1.5) 
    # Add the boxplot layer
    + plt9.geom_boxplot(plt9.aes(fill='test_site'), width=0.1, alpha=0.5)
    # Set the y-axis labels and breaks
    + plt9.scale_y_continuous(name = timing + ' (s)')
    # Set the x-axis label
    + plt9.xlab(xlab='Test site')
    # Set the color scale legend
    + plt9.scale_color_continuous(name = 'sub-id')
    # Set the theme to a white background with a base font size of 14
    + plt9.theme_bw(base_size = 14)
    )
    p.save(filename=save_path + timing + '_test_site.png')

# %%
# histogram of time difference between samples for right hand tracking
plt.figure()
plt.hist(1/rh_t_diff, bins=200, range=(0, 200))
plt.title('Right hand tracking')
plt.xlabel('Framerate (Hz)')
plt.ylabel('Number of samples')  
plt.show()
#save figure
plt.savefig(save_path + 'RH_framerate.png', dpi=300)  
max_t_diff = np.nanmax(rh_t_diff)
min_t_diff = np.nanmin(rh_t_diff)
median_t_diff = np.nanmedian(rh_t_diff)
median_fs = 1 / median_t_diff
mean_fs = 1 / np.nanmean(rh_t_diff)
std_fs = np.nanstd(1/rh_t_diff)


# %%
# remove column 'sub' from dataframe d
d_stats = d.drop(columns=['sub'])
# mean and std of framerate and experiment duration grouped by test site
d_stats = d_stats.groupby('test_site').agg(['mean', 'std'])
#save dataframe
d_stats.to_csv(save_path + 'exp_dur_stats.csv', index=True)


# %%
