########################################################################################################################
# Script to get descriptive stats from demographics, CRi, SR and survey data from AVR experiment
# Need csv files containing demographics, CRi, SR and survey data, from cri_sr.py, import_PRE_questionnaire.py and import_POST_questionnaire.py
# Output: csv files (dataframe) containing descriptive stats for each data type
# Author: Antonin Fourcade
# Last version: 12.06.2025
########################################################################################################################

# %%
# import packages
import pandas as pd
import os

# %%
# Set paths and experiment parameters
phase_path = 'E:/AffectiveVR/Phase_2/'
data_path = phase_path + 'Data/'
bids_folder = 'AVR' # folder with data in bids format

# path to demographics csv file and load data
filename_dem = 'participants.csv'
dem_path = data_path + bids_folder + '/' + filename_dem
dem = pd.read_csv(dem_path)

# path to CRi-SR csv file and load data
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'} # dictionnary for rating dimensions
filename_cri = 'cri_sr_clean.csv'
cri_path = data_path + bids_folder + '/' + filename_cri
cri_sr = pd.read_csv(cri_path)
# CRis of interest
cri_select = ['cr_mean', 'cr_std', 'cr_skew', 'cr_kurtosis']

# path to survey csv files and load data
filename_pre = 'pre_survey_preprocessed.csv'
filename_post = 'post_survey_preprocessed.csv'
pre_path = data_path + bids_folder + '/' + filename_pre
post_path = data_path + bids_folder + '/' + filename_post
pre_survey = pd.read_csv(pre_path)
post_survey = pd.read_csv(post_path)

# excluded subjects
excluded_subs = ['sub-13'] # subjects to exclude from the analysis
# remove excluded subjects from demographics, CRi and SR data
dem = dem[~dem['sub'].isin(excluded_subs)]
cri_sr = cri_sr[~cri_sr['sub'].isin(excluded_subs)]
pre_survey = pre_survey[~pre_survey['sub'].isin(excluded_subs)]
post_survey = post_survey[~post_survey['sub'].isin(excluded_subs)]

# save path
save_path = 'E:/AffectiveVR/affecttracker_validation/results/evaluation/descriptive_stats/'
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# %%
# Get descriptive stats for demographics data
save_dem_path = save_path + 'demographics/'
if not os.path.exists(save_dem_path): # create folder if it doesn't exist
    os.makedirs(save_dem_path)
age_stats = dem['age'].describe()
age_stats_per_site = dem.groupby('test_site')['age'].describe().round(1)
gender_count_per_site = dem.groupby('test_site')['gender'].value_counts()
position_count_per_site = dem.groupby('test_site')['position'].value_counts()
education_count_per_site = dem.groupby('test_site')['education'].value_counts()
print(gender_count_per_site)
# save descriptive stats for demographics data
age_stats.to_csv(save_dem_path + 'age_stats.csv')
age_stats_per_site.to_csv(save_dem_path + 'age_stats_per_site.csv')
gender_count_per_site.to_csv(save_dem_path + 'gender_count_per_site.csv')
position_count_per_site.to_csv(save_dem_path + 'position_count_per_site.csv')
education_count_per_site.to_csv(save_dem_path + 'education_count_per_site.csv')

# %%
# Get descriptive stats for CRi and SR data (values rounded to 2 decimals)
save_cri_sr_path = save_path + 'cri_sr/'
if not os.path.exists(save_cri_sr_path): # create folder if it doesn't exist
    os.makedirs(save_cri_sr_path)
# descriptive stats grouped by all factors
#cri_sr.groupby(['test_site', 'position', 'gender']).describe().round(2)
# Get descriptive stats for CRi and SR data for each rating dimension
for dim in rat_dim:
    # select data only for current rating dimension
    #cri_sr_dim = cri_sr.filter(like='_' + rat_dim[dim])
    # select SR and CRis of interest
    cri_sr_select = ['sr'] + cri_select
    cri_sr_select_dim = [x + '_' + rat_dim[dim] for x in cri_sr_select]
    cri_sr_dim = cri_sr[cri_sr_select_dim]
    # get descriptive stats without grouping
    d_stats = cri_sr_dim.describe().round(2)
    d_stats.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '.csv', index=False)
    # get descriptive stats with grouping
    cri_sr_dim = pd.concat([cri_sr_dim, cri_sr[['sub', 'test_site', 'position', 'gender']]], axis=1)
    # get descriptive stats for each site
    d_stats_site = cri_sr_dim.groupby('test_site').describe().round(2)
    d_stats_site.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_site.csv', index=False)
    # get descriptive stats for each gender
    d_stats_gender = cri_sr_dim.groupby('gender').describe().round(2)
    d_stats_gender.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_gender.csv', index=False)
    # get descriptive stats for each position
    d_stats_position = cri_sr_dim.groupby('position').describe().round(2)
    d_stats_position = d_stats_position.drop('standing') # remove standing position (redundant with d_stats_site)
    d_stats_position.index = pd.MultiIndex.from_tuples([('all', 'seated')], names=['test_site', 'position']) # reformat index
    # get descriptive stats for each position and test site
    d_stats_site_position = cri_sr_dim.groupby(['test_site', 'position']).describe().round(2)
    # concatenate descriptive stats for each position and test site
    d_stats_position = pd.concat([d_stats_position, d_stats_site_position])
    d_stats_position.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_position.csv', index=False)
    
    # create another version of the table
    # describe but only mean, std, min and max
    d_stats_position = cri_sr_dim.groupby('position').describe(percentiles=None).round(2)
    d_stats_position = d_stats_position.drop('standing')
    # drop count, 25%, 50%, 75% columns and cr_std, cr_skew, cr_kurtosis columns
    d_stats_position = d_stats_position.drop(columns=['cr_std_' + rat_dim[dim], 'cr_skew_' + rat_dim[dim], 'cr_kurtosis_' + rat_dim[dim]])
    d_stats_position = d_stats_position.drop(['count', '25%', '50%', '75%'], axis=1, level=1)
    d_stats_position.index = pd.MultiIndex.from_tuples([('TOR+BER', 'seated')], names=['test_site', 'position']) # reformat index
    d_stats_site_position = cri_sr_dim.groupby(['test_site', 'position']).describe().round(2)
    d_stats_site_position = d_stats_site_position.drop(columns=['cr_std_' + rat_dim[dim], 'cr_skew_' + rat_dim[dim], 'cr_kurtosis_' + rat_dim[dim]])
    d_stats_site_position = d_stats_site_position.drop(['count', '25%', '50%', '75%'], axis=1, level=1)
    d_stats_position = pd.concat([d_stats_position, d_stats_site_position])
    # sort rows with TOR seated, BER seated, BER standing, all seated
    d_stats_position = d_stats_position.reindex(['TOR', 'BER', 'TOR+BER'], level=0)
    #save
    d_stats_position.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_position_simple.csv', index=False)

# %%
# Get descriptive stats for post survey data
save_post_survey_path = save_path + 'post_survey/'
if not os.path.exists(save_post_survey_path): # create folder if it doesn't exist
    os.makedirs(save_post_survey_path)
# get descriptive stats without grouping
post_survey_stats = post_survey.groupby('questionnaire').describe().round(2)
post_survey_stats.reset_index().to_csv(save_post_survey_path + 'post_survey.csv', index=False)
# get descriptive stats for each site
post_survey_stats_site = post_survey.groupby(['test_site', 'questionnaire']).describe().round(2)
post_survey_stats_site.reset_index().to_csv(save_post_survey_path + 'post_survey_site.csv', index=False)
# get descriptive for each gender
post_survey_stats_gender = post_survey.groupby(['gender', 'questionnaire']).describe().round(2)
post_survey_stats_gender.reset_index().to_csv(save_post_survey_path + 'post_survey_gender.csv', index=False)
# get descriptive for each position
post_survey_stats_position = post_survey.groupby(['position', 'questionnaire']).describe().round(2)
post_survey_stats_position = post_survey_stats_position.drop('standing') # remove standing position (redundant with post_survey_stats_site)
# reformat index
# get the current index values
p = post_survey_stats_position.index.get_level_values('position').unique().tolist()
q = post_survey_stats_position.index.get_level_values('questionnaire').unique().tolist()
# add the new level test_site with value 'all'
post_survey_stats_position.index = pd.MultiIndex.from_product([['all'], p, q], 
                                       names=['test_site', 'position', 'questionnaire'])
# get descriptive for each position and test site
post_survey_stats_site_position = post_survey.groupby(['test_site', 'position', 'questionnaire']).describe().round(2)
# concatenate descriptive stats for each position and test site
post_survey_stats_position = pd.concat([post_survey_stats_position, post_survey_stats_site_position])
post_survey_stats_position.reset_index().to_csv(save_post_survey_path + 'post_survey_position.csv', index=False)

# %%
# Get descriptive stats for pre survey data
save_pre_survey_path = save_path + 'pre_survey/'
if not os.path.exists(save_pre_survey_path): # create folder if it doesn't exist
    os.makedirs(save_pre_survey_path)
# get descriptive stats without grouping
pre_survey_stats = pre_survey.groupby('questionnaire').describe().round(2)
pre_survey_stats.reset_index().to_csv(save_pre_survey_path + 'pre_survey.csv', index=False)
# get descriptive stats for each site
pre_survey_stats_site = pre_survey.groupby(['test_site', 'questionnaire']).describe().round(2)
pre_survey_stats_site.reset_index().to_csv(save_pre_survey_path + 'pre_survey_site.csv', index=False)
# get descriptive stats for each gender
pre_survey_stats_gender = pre_survey.groupby(['gender', 'questionnaire']).describe().round(2)
pre_survey_stats_gender.reset_index().to_csv(save_pre_survey_path + 'pre_survey_gender.csv', index=False)
# get descriptive stats for each position
pre_survey_stats_position = pre_survey.groupby(['position', 'questionnaire']).describe().round(2)
pre_survey_stats_position = pre_survey_stats_position.drop('standing') # remove standing position (redundant with pre_survey_stats_site)
# reformat index
# get the current index values
p = pre_survey_stats_position.index.get_level_values('position').unique().tolist()
q = pre_survey_stats_position.index.get_level_values('questionnaire').unique().tolist()
# add the new level test_site with value 'all'
pre_survey_stats_position.index = pd.MultiIndex.from_product([['all'], p, q], 
                                       names=['test_site', 'position', 'questionnaire'])
# get descriptive stats for each position and test site
pre_survey_stats_site_position = pre_survey.groupby(['test_site', 'position', 'questionnaire']).describe().round(2)
# concatenate descriptive stats for each position and test site
pre_survey_stats_position = pd.concat([pre_survey_stats_position, pre_survey_stats_site_position])
pre_survey_stats_position.reset_index().to_csv(save_pre_survey_path + 'pre_survey_position.csv', index=False)


# %%
