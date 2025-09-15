########################################################################################################################
# Script to get descriptive stats from demographics, CRi, SR and survey data from AVR experiment
# Need csv files containing demographics, CRi, SR and survey data, from cri_sr.py, import_PRE_questionnaire.py and import_POST_questionnaire.py
# Output: csv files (dataframe) containing descriptive stats for each data type
# Author: Antonin Fourcade
# Last version: 22.08.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import os
from scipy.stats import ttest_ind
from statsmodels.formula.api import ols
import statsmodels.api as sm

# %%
# Set paths and experiment parameters
phase_path = 'E:/AffectiveVR/Phase_1/'
data_path = phase_path + 'Data/'
bids_folder = 'AVR' # folder with data in bids format

# experiment parameters
# blocks = ['Practice', 'Experiment']
# logging_freq = ['CR', 'SR']
# test_site = ['TOR', 'BER']  # Torino = 0, Berlin = 1
# position = ['seated', 'standing']  # Seated = 0, Standing = 1
# gender = ['male', 'female', 'non_binary'] # male = 0, female = 1, non_binary = 2
# cri_list = ['last', 'mean', 'median', 'mode', 'max', 'min', 'std', 'cv', 'range', 'iqr', 'skew', 'kurtosis', 'auc', 'cp']
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'} # dictionnary for rating dimensions
covariates = ['test_site', 'position', 'gender'] # covariates for ANOVA, note: these are categorical variables

# path to demographics csv file and load data
filename_dem = 'participants.csv'
dem_path = data_path + bids_folder + '/' + filename_dem
dem = pd.read_csv(dem_path)

# path to CRi-SR csv file and load data
filename_cri = 'cri_sr_clean.csv'
cri_path = data_path + bids_folder + '/' + filename_cri
cri_sr = pd.read_csv(cri_path)
# CRis of interest
cri_select = ['cr_mean', 'cr_std', 'cr_skew', 'cr_kurtosis']

# save path
save_path = phase_path + 'affectivevr/descriptive_stats/'
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
    cri_sr_dim = pd.concat([cri_sr_dim, cri_sr[['sj_id', 'test_site', 'rating_method', 'quadrant']]], axis=1)
    # get descriptive stats for each site
    d_stats_site = cri_sr_dim.groupby('test_site').describe().round(2)
    d_stats_site.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_site.csv', index=False)
    # get descriptive stats for each rating_method
    d_stats_rating_method = cri_sr_dim.groupby('rating_method').describe().round(2)
    d_stats_rating_method.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_rating_method.csv', index=False)
    # get descriptive stats for each quadrant
    d_stats_quadrant = cri_sr_dim.groupby('quadrant').describe().round(2)
    d_stats_quadrant.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_quadrant.csv', index=False)
    
    # get descriptive stats for each quadrant and rating_method
    d_stats_quadrant_rating_method = cri_sr_dim.groupby(['quadrant', 'rating_method'])[['cr_mean_' + rat_dim[dim],'sr_' + rat_dim[dim]]].describe().round(2)
    # drop count, 25%, 50%, 75%, columns
    d_stats_quadrant_rating_method = d_stats_quadrant_rating_method.drop(['count', '25%', '50%', '75%'], axis=1, level=1)
    # save
    d_stats_quadrant_rating_method.reset_index().to_csv(save_cri_sr_path + 'cri_sr_' + dim + '_quadrant_rating_method.csv', index=False)

    # test if differences in cr means and sr between test sites are significant
    ttest_cr_site = ttest_ind(cri_sr_dim[cri_sr_dim['test_site'] == 'Berlin']['cr_mean_' + rat_dim[dim]], cri_sr_dim[cri_sr_dim['test_site'] == 'Torino']['cr_mean_' + rat_dim[dim]], nan_policy='omit')
    ttest_sr_site = ttest_ind(cri_sr_dim[cri_sr_dim['test_site'] == 'Berlin']['sr_' + rat_dim[dim]], cri_sr_dim[cri_sr_dim['test_site'] == 'Torino']['sr_' + rat_dim[dim]], nan_policy='omit')
    # save t-test results
    with open(save_cri_sr_path + 'ttest_ind_site_' + dim + '.txt', 'w') as f:
        f.write('T-test results for cr_mean_' + rat_dim[dim] + ' between Berlin and Torino\n')
        f.write('t = ' + str(ttest_cr_site.statistic) + '\n')
        f.write('p = ' + str(ttest_cr_site.pvalue) + '\n')
        f.write('df = ' + str(ttest_cr_site.df) + '\n')
        if ttest_cr_site.pvalue < 0.05:
            f.write('Significant difference between Berlin and Torino\n')
        else:
            f.write('No significant difference between Berlin and Torino\n')
        f.write('\n')
        f.write('T-test results for sr_' + rat_dim[dim] + ' between Berlin and Torino\n')
        f.write('t = ' + str(ttest_sr_site.statistic) + '\n')
        f.write('p = ' + str(ttest_sr_site.pvalue) + '\n')
        f.write('df = ' + str(ttest_sr_site.df) + '\n')
        if ttest_sr_site.pvalue < 0.05:
            f.write('Significant difference between Berlin and Torino\n')
        else:
            f.write('No significant difference between Berlin and Torino\n')
    # test if differences in cr means and sr between test site and gender are significant with anova
    # ANOVA on CRis, with test_site, position and gender as factors
    model = ols('cr_mean_' + rat_dim[dim] + ' ~ C(test_site) + C(gender)', data=cri_sr_dim).fit()
    anova_table = sm.stats.anova_lm(model, typ=2)
    anova_table.index.name = 'cr_mean_' + dim
    anova_table.reset_index(inplace=True)




# %%
# Get descriptive stats for survey data
save_survey_path = save_path + 'assessment/'
if not os.path.exists(save_survey_path): # create folder if it doesn't exist
    os.makedirs(save_survey_path)
# load survey data
filename_survey = 'assessment_preprocessed.csv'
survey_path = data_path + bids_folder + '/' + filename_survey
survey = pd.read_csv(survey_path)
# get descriptive stats for survey data
survey_stats = survey.groupby('questionnaire')['response'].describe().round(2)
survey_stats.reset_index().to_csv(save_survey_path + 'survey_stats.csv', index=False)
# get descriptive stats for each rating_method and questionnnaire
survey_stats_rating_method = survey.groupby(['rating_method','questionnaire'])['response'].describe().round(2)
# drop count, 25%, 50%, 75%, columns
survey_stats_rating_method = survey_stats_rating_method.drop(['count', '25%', '50%', '75%'], axis=1)
# pivot table with rating_method as rows and questionnaires as columns, and mean, std, min, max as subcolumns
survey_stats_rating_method = survey_stats_rating_method.unstack()
# have mean, std, min, max as level 2 columns, for each level 1 column
survey_stats_rating_method = survey_stats_rating_method.swaplevel(axis=1)
survey_stats_rating_method = survey_stats_rating_method.stack(future_stack=True)
survey_stats_rating_method = survey_stats_rating_method.unstack()
# rename presence_score to presence, sus_score to sus and emo_rep to emotion representation
survey_stats_rating_method = survey_stats_rating_method.rename(columns={'presence_score': 'presence', 'sus_score': 'sus', 'emo_rep': 'emotion representation'})
# reorder columns with invasiveness first, then presence, sus, emo_rep and satisfaction
survey_stats_rating_method = survey_stats_rating_method[['invasiveness', 'presence', 'sus', 'emotion representation', 'satisfaction']]
# save
survey_stats_rating_method.reset_index().to_csv(save_survey_path + 'survey_stats_rating_method.csv', index=False)

# get same table but with columns and rows swapped
survey_stats_rating_method = survey.groupby(['rating_method','questionnaire'])['response'].describe().round(2)
survey_stats_rating_method = survey_stats_rating_method.drop(['count', '25%', '50%', '75%'], axis=1)
survey_stats_rating_method = survey_stats_rating_method.pivot_table(index='questionnaire', columns='rating_method').swaplevel(axis=1).stack(future_stack=True).unstack()
#rename presence_score to presence, sus_score to sus and emo_rep to emotion representation
survey_stats_rating_method = survey_stats_rating_method.rename({'presence_score': 'presence', 'sus_score': 'sus', 'emo_rep': 'emotion representation'})
# reorder rows with invasiveness first, then presence, sus, emo_rep and satisfaction
survey_stats_rating_method = survey_stats_rating_method.reindex(['invasiveness', 'presence', 'sus', 'emotion representation', 'satisfaction'])
# save
survey_stats_rating_method.to_csv(save_survey_path + 'survey_stats_rating_method_transposed.csv')


# %%
