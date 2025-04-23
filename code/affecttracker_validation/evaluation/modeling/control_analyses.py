########################################################################################################################
# Script for control analyses of demographics, CRi, SR and survey data from AVR experiment
# Need csv files containing CRi, SR and survey data, from cri_sr.py, import_PRE_questionnaire.py and import_POST_questionnaire.py
# Output: csv files (dataframe) containing ANOVA results
# Author: Antonin Fourcade
# Last version: 16.02.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import os
import statsmodels.api as sm
from statsmodels.formula.api import ols
from statsmodels.stats.mediation import Mediation
from pingouin import mediation_analysis, ttest


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
# cri_list = ['last', 'mean', 'median', 'mode', 'max', 'min', 'std', 'cv', 'range', 'iqr', 'skew', 'kurtosis', 'auc', 'cp']
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'} # dictionnary for rating dimensions
covariates = ['test_site', 'position', 'gender'] # covariates for ANOVA, note: these are categorical variables

# path to CRi-SR csv file and load data
filename_cri = 'cri_sr_clean.csv'
cri_path = data_path + bids_folder + '/' + filename_cri
cri_sr = pd.read_csv(cri_path)
# CRis of interest
cri_select = ['cr_mean', 'cr_std', 'cr_skew', 'cr_kurtosis']
# position of interest
pos_select = ['seated'] # ['seated', 'standing']

# path to survey csv files and load data
filename_pre = 'pre_survey_preprocessed.csv'
filename_post = 'post_survey_preprocessed.csv'
pre_path = data_path + bids_folder + '/' + filename_pre
post_path = data_path + bids_folder + '/' + filename_post
pre_survey = pd.read_csv(pre_path)
post_survey = pd.read_csv(post_path)

# save path
if len(pos_select) == 2:
    save_path = phase_path + 'affectivevr/control_analyses/all_positions/'
else:
    save_path = phase_path + 'affectivevr/control_analyses/' + pos_select[0] + '/'
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# %%
# Select only position of interest
cri_sr = cri_sr[cri_sr['position'].isin(pos_select)]
pre_survey = pre_survey[pre_survey['position'].isin(pos_select)]
post_survey = post_survey[post_survey['position'].isin(pos_select)]

# %%
# Survey data
# ANOVAs
# Loop through pre and post survey data
for survey in [pre_survey, post_survey]:
    # get score list
    score_list = survey['questionnaire'].unique()
    # loop through each score
    for score in score_list:
        # select score
        score_df = survey[survey['questionnaire'] == score]      
        if len(pos_select) == 2: # if both seated and standing positions
            # ANOVA on responses, with test_site, position and gender as factors
            model = ols('response ~ C(test_site) + C(gender) + C(position)', data=score_df).fit() # C() specifies categorical variables
            anova_table = sm.stats.anova_lm(model, typ=2) # Type 2 ANOVA
            anova_table.index.name = score # rename index column to score
            anova_table.reset_index(inplace=True) # turn index into column
            # ANOVA on responses, with position as factor, but only for test_site = BER
            model_BER = ols('response ~ C(position)', data=score_df[score_df['test_site'] == 'BER']).fit()
            anova_table_BER = sm.stats.anova_lm(model_BER, typ=2) # Type 2 ANOVA
            anova_table_BER.index.name = 'BER_' + score # rename index column to BER score
            anova_table_BER.reset_index(inplace=True) # turn index into column
            # combine anova tables and save
            anova_table_all = pd.concat([anova_table, anova_table_BER], axis=1)
            anova_table_all.to_csv(save_path + 'anova_' + score + '.csv', index=False, na_rep='')
            # post-hoc BER analysis with t-tests if significant effect of position
            if anova_table_BER['PR(>F)'][0] < 0.05:
                print('Significant effect in BER of position for ' + score)
                # run t-tests between seated and standing positions
                standing_score = score_df.loc[(score_df['test_site'] == 'BER') & (score_df['position'] == 'standing'), 'response']
                seated_score = score_df.loc[(score_df['test_site'] == 'BER') & (score_df['position'] == 'seated'), 'response']
                t_test = ttest(standing_score, seated_score)
                t_test.index.name = score # rename index column to score
                t_test.reset_index(inplace=True) # turn index into column
                # save results
                t_test.round(2).to_csv(save_path + 'BER_position_posthoc_ttest_' + score + '.csv', index=False)
        else:
             # ANOVA on CRis, with test_site and gender as factors
             model = ols('response ~ C(test_site) + C(gender)', data=score_df).fit()
             anova_table = sm.stats.anova_lm(model, typ=2) # Type 2 ANOVA
             anova_table.index.name = score # rename index column to score
             anova_table.reset_index(inplace=True) # turn index into column
             anova_table.to_csv(save_path + 'anova_' + score + '.csv', index=False, na_rep='')
            

# %%
# SR data
# keep only SR data
sr_all = cri_sr[['sub'] + covariates + ['sr_v', 'sr_a', 'sr_dist', 'sr_angle']]
# loop through each rating dimension
anova_table_all = pd.DataFrame()
for dim in rat_dim:
    if len(pos_select) == 2: # if both seated and standing positions
        # ANOVA on SRs, with test_site, position and gender as factors
        model = ols('sr_' + rat_dim[dim] + ' ~ C(test_site) + C(gender) + C(position)', data=sr_all).fit() # C() specifies categorical variables
        anova_table = sm.stats.anova_lm(model, typ=2) # Type 2 ANOVA
        anova_table.index.name = 'SR_' + dim # rename index column
        anova_table.reset_index(inplace=True) # turn index into column
        # ANOVA on SRs, with position as factor, but only for test_site = BER
        model_BER = ols('sr_' + rat_dim[dim] + ' ~ C(position)', data=sr_all[sr_all['test_site'] == 'BER']).fit()
        anova_table_BER = sm.stats.anova_lm(model_BER, typ=2) # Type 2 ANOVA
        anova_table_BER.index.name = 'BER_SR_' + dim # rename index column
        anova_table_BER.reset_index(inplace=True) # turn index into column
        # combine anova tables
        anova_table_dim = pd.concat([anova_table, anova_table_BER], axis=1)
        # post-hoc BER analysis with t-tests if significant effect of position
        if anova_table_BER['PR(>F)'][0] < 0.05:
                print('Significant effect in BER of position for SR_' + dim)
                # run t-tests between seated and standing positions
                standing_sr = sr_all.loc[(sr_all['test_site'] == 'BER') & (sr_all['position'] == 'standing'), 'sr_' + rat_dim[dim]]
                seated_sr = sr_all.loc[(sr_all['test_site'] == 'BER') & (sr_all['position'] == 'seated'), 'sr_' + rat_dim[dim]]
                t_test = ttest(standing_sr, seated_sr)
                t_test.index.name = 'SR_' + dim # rename index column to score
                t_test.reset_index(inplace=True) # turn index into column
                # save results
                t_test.round(2).to_csv(save_path + 'BER_position_posthoc_ttest_sr_' + dim + '.csv', index=False)
    else:
        # ANOVA on CRis, with test_site and gender as factors
        model = ols('sr_' + rat_dim[dim] + ' ~ C(test_site) + C(gender)', data=sr_all).fit()
        anova_table_dim = sm.stats.anova_lm(model, typ=2) # Type 2 ANOVA
        anova_table_dim.index.name = 'SR_' + dim # rename index column
        anova_table_dim.reset_index(inplace=True) # turn index into column
    # add to anova table for all rating dimensions
    anova_table_all = pd.concat([anova_table_all, anova_table_dim], axis=1)
        
# save anova table
anova_table_all.to_csv(save_path + 'anova_sr.csv', index=False, na_rep='')

# %%
# CRi data
# loop through each cri of interest
for cri in cri_select:
    anova_table_all = pd.DataFrame()
    # loop through each rating dimension
    for dim in rat_dim:
        # select CRi of interest for current rating dimension
        cri_select_dim = cri + '_' + rat_dim[dim]
        cri_dim = cri_sr[['sub'] + covariates + [cri_select_dim]]
        if len(pos_select) == 2:
            # ANOVA on CRis, with test_site, position and gender as factors
            model = ols(cri_select_dim + ' ~ C(test_site) + C(position) + C(gender)', data=cri_dim).fit()
            anova_table = sm.stats.anova_lm(model, typ=2)
            anova_table.index.name = cri + '_' + dim
            anova_table.reset_index(inplace=True)
            # ANOVA on CRis, with position as factor, but only for test_site = BER
            model_BER = ols(cri_select_dim + ' ~ C(position)', data=cri_dim[cri_dim['test_site'] == 'BER']).fit()
            anova_table_BER = sm.stats.anova_lm(model_BER, typ=2)
            anova_table_BER.index.name = 'BER_' + cri + '_' + dim
            anova_table_BER.reset_index(inplace=True)
            # combine anova tables
            anova_table_dim = pd.concat([anova_table, anova_table_BER], axis=1)
            # post-hoc BER analysis with t-tests if significant effect of position
            if anova_table_BER['PR(>F)'][0] < 0.05:
                    print('Significant effect in BER of position for ' + cri_select_dim)
                    # run t-tests between seated and standing positions
                    standing_cri_dim = cri_dim.loc[(cri_dim['test_site'] == 'BER') & (cri_dim['position'] == 'standing'), cri_select_dim]
                    seated_cri_dim = cri_dim.loc[(cri_dim['test_site'] == 'BER') & (cri_dim['position'] == 'seated'), cri_select_dim]
                    t_test = ttest(standing_cri_dim, seated_cri_dim)
                    t_test.index.name = cri_select_dim # rename index column to score
                    t_test.reset_index(inplace=True) # turn index into column
                    # save results
                    t_test.round(2).to_csv(save_path + 'BER_position_posthoc_ttest_' + cri_select_dim + '.csv', index=False)
        else:
             # ANOVA on CRis, with test_site and gender as factors
             model = ols(cri_select_dim + ' ~ C(test_site) + C(gender)', data=cri_dim).fit()
             anova_table_dim = sm.stats.anova_lm(model, typ=2)
             anova_table_dim.index.name = cri + '_' + dim
             anova_table_dim.reset_index(inplace=True)
        # add to anova table for all rating dimensions
        anova_table_all = pd.concat([anova_table_all, anova_table_dim], axis=1)
    
            
    # save anova table
    anova_table_all.to_csv(save_path + 'anova_' + cri + '.csv', index=False, na_rep='')
    
# %%
# Mediator analysis position on CR mean with MAIA scores as mediators
# with maia scores as multiple mediators
scores =['maia_emotional_awareness', 'maia_body_listening']
#select cri = cr_mean_a only for BER site
cr_mean_a_BER = cri_sr[cri_sr['test_site'] == 'BER'][['sub', 'position', 'gender', 'cr_mean_a']]
# select maia scores only for BER site
# loop through each score
med_data = cr_mean_a_BER
for score in scores:
    maia_score_BER = pre_survey[(pre_survey['questionnaire'] == score) & (pre_survey['test_site'] == 'BER')][['sub', 'response']]
    # rename response column to maia scores
    maia_score_BER.columns = ['sub', score]
    # merge cr_mean_a_BER and maia_score_BER
    med_data = med_data.merge(maia_score_BER, on='sub')
# convert position to numeric (seated = 0, standing = 1)
med_data['position'] = med_data['position'].map({'seated': 0, 'standing': 1})
med_data['gender'] = med_data['gender'].map({'male':0, 'female':1, 'non_binary':2})
# mediation model
med = mediation_analysis(data=med_data, x='position', m=scores, y='cr_mean_a', covar='gender', alpha=0.05, seed=42)
med.round(2).to_csv(save_path + 'mediation_BER_position_cr_mean_a.csv', index=False, na_rep='')
# another way through statsmodels
# outcome_model = sm.OLS.from_formula('cr_mean_a ~ maia_emotional_awareness + position', data=med_data)
# mediator_model = sm.OLS.from_formula('maia_emotional_awareness ~ position', data=med_data)
# med = Mediation(outcome_model, mediator_model, 'position', 'maia_emotional_awareness').fit()
# med.summary()

# for all test sites
#select cri = cr_mean_a
cr_mean_a = cri_sr[['sub', 'test_site', 'gender', 'position', 'cr_mean_a']]
# select maia scores only for BER site
# loop through each score
med_data = cr_mean_a
for score in scores:
    maia_score = pre_survey[pre_survey['questionnaire'] == score][['sub', 'response']]
    # rename response column to maia scores
    maia_score.columns = ['sub', score]
    # merge cr_mean_a_BER and maia_score_BER
    med_data = med_data.merge(maia_score, on='sub')
# convert position to numeric (seated = 0, standing = 1)
med_data['position'] = med_data['position'].map({'seated': 0, 'standing': 1})
med_data['test_site'] = med_data['test_site'].map({'TOR': 0, 'BER': 1})
med_data['gender'] = med_data['gender'].map({'male':0, 'female':1, 'non_binary':2})
# mediation model
med = mediation_analysis(data=med_data, x='position', m=scores, y='cr_mean_a', covar=['test_site', 'gender'], alpha=0.05, seed=42)
med.round(2).to_csv(save_path + 'mediation_position_cr_mean_a.csv', index=False, na_rep='')
# another way through statsmodels
# outcome_model = sm.OLS.from_formula('cr_mean_a ~ maia_emotional_awareness + position + test_site + gender', data=med_data)
# mediator_model = sm.OLS.from_formula('maia_emotional_awareness ~ position + test_site + gender', data=med_data)
# med = Mediation(outcome_model, mediator_model, 'position', 'maia_emotional_awareness').fit()
# med.summary()

# %%
# Moderator analysis position on CR mean with MAIA scores as moderators
# ANOVA on cr_mean_a, with test_site, position, gender, maia_emotional_awareness and maia_body_listening as factors
# moderation? -> interaction
model = ols('cr_mean_a ~ C(test_site) + C(position) + C(gender) + maia_emotional_awareness + maia_body_listening + C(position)*maia_emotional_awareness + C(position)*maia_body_listening', data=med_data).fit()
model = ols('cr_mean_a ~ + C(position) + maia_emotional_awareness + maia_body_listening + C(position)*maia_emotional_awareness + C(position)*maia_body_listening', data=med_data).fit()
anova_table = sm.stats.anova_lm(model, typ=2)
anova_table.index.name = 'cr_mean_a'
anova_table.reset_index(inplace=True)
anova_table.round(2).to_csv(save_path + 'moderation_position_cr_mean_a.csv', index=False, na_rep='')

# %%
