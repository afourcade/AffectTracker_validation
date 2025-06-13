########################################################################################################################
# Script for analyses of correlations between CRis and SRs from AVR experiment
# Need csv files containing CRi and SR data, from cri_sr.py
# Output: csv files (dataframe) containing results
# Author: Antonin Fourcade
# Last version: 12.06.2025
########################################################################################################################

# %%
# import packages
import pandas as pd
import os
#import statsmodels.api as sm
import numpy as np
#from statsmodels.formula.api import ols
from scipy.stats import norm, pearsonr
from pingouin import corr

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
cri_list = ['last', 'mean', 'median', 'mode', 'max', 'min', 'std', 'cv', 'range', 'iqr', 'skew', 'kurtosis', 'auc', 'cp']
rat_dim = {'valence': 'v', 'arousal': 'a', 'distance': 'dist', 'angle': 'angle'} # dictionnary for rating dimensions
covariates = ['test_site', 'position', 'gender'] # covariates for ANOVA, note: these are categorical variables

# path to CRi-SR csv file and load data
filename_cri = 'cri_sr_clean.csv'
cri_path = data_path + bids_folder + '/' + filename_cri
cri_sr = pd.read_csv(cri_path)
# CRis of interest
cri_select = ['cr_mean', 'cr_std', 'cr_skew', 'cr_kurtosis']
# position of interest
#pos_select = ['seated'] 
pos_select = ['seated', 'standing']

# excluded subjects
excluded_subs = ['sub-13']
# remove excluded subjects from cri_sr
cri_sr = cri_sr[~cri_sr['sub'].isin(excluded_subs)]

# save path
save_path = 'E:/AffectiveVR/affecttracker_validation/results/evaluation/corr_results/'
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# phase 1 data
phase1_path = 'E:/AffectiveVR/Phase_1/affectivevr/cocor_results/'
filename_cri_phase1 = 'summary_cocor_cor_results.csv'

# %% Functions
# Fisher r-to-z transformation and comparison of correlation coefficients
def independent_corr(corr1, corr2, n1, n2, twotailed=True):
    # Calculates the statistic significance between two independent correlation coefficients
    # Inputs:
    # corr1: 1st correlation coefficient
    # corr2: 2nd correlation coefficient
    # n1: number of observations with which corr1 is calculated
    # n2: number of observations with which corr2 is calculated
    # twotailed: whether to calculate a one or two tailed test
    # Outputs: z and p-val

    corr1_z = 0.5 * np.log((1 + corr1)/(1 - corr1))
    corr2_z = 0.5 * np.log((1 + corr2)/(1 - corr2))

    se_diff_r = np.sqrt(1/(n1 - 3) + 1/(n2 - 3))
    diff = corr1_z - corr2_z
    z = np.abs(diff / se_diff_r)
    p = (1 - norm.cdf(z))
    if twotailed:
        p *= 2
    return z, p

# %%
# Select only position of interest
cri_sr = cri_sr[cri_sr['position'].isin(pos_select)]
    
# %%
# Correlation between CRis and SRs, for each rating dimension
# create dataframe with correlation results, with first column the cri_list
res_df = pd.DataFrame(index=['cr_' + cri for cri in cri_list], columns=list(rat_dim.keys()))
# add level 2 columns: r and p
multi_idx = pd.MultiIndex.from_product([res_df.columns, ['r', 'p']])
res_df = pd.DataFrame(index=res_df.index, columns=multi_idx)
# loop over CRis
for cri in cri_list:
    # save path for each CRi
    save_cri_path = save_path + 'cr_' + cri + '/'
    if not os.path.exists(save_cri_path): # create folder if it doesn't exist
        os.makedirs(save_cri_path)
    # loop over rating dimensions
    for dim in rat_dim:
        cri_dim = 'cr_' + cri + '_' + rat_dim[dim]
        sr_dim = 'sr_' + rat_dim[dim]
        # correlation
        corr_res = corr(x=cri_sr[cri_dim], y=cri_sr[sr_dim], method='pearson')
        #corr_res = pearsonr(cri_sr[cri_dim], cri_sr[sr_dim])
        # degrees of freedom
        #n = len(cri_sr[cri_dim])
        #df = n - 2
        # t-statistic
        #t = corr_res[0] * np.sqrt(df / (1 - corr_res[0]**2))
        # rename index column
        corr_res.index.name = 'corr(' + cri_dim + ', ' + sr_dim + ')'
        # save to csv
        if len(pos_select) == 2:
            filename = 'corr_cr_' + cri + '_' + rat_dim[dim] + '.csv'
        else:
            filename = 'corr_cr_' + cri + '_' + rat_dim[dim] + '_' + pos_select[0] + '.csv'
        corr_res.to_csv(save_cri_path + filename)
        # save r in res_df (round to 3 decimals)
        res_df.loc['cr_' + cri, (dim,'r')] = corr_res['r'].iloc[0].round(3)
        res_df.loc['cr_' + cri, (dim,'p')] = corr_res['p-val'].iloc[0].round(3)
# rename index of res_df
res_df.index.name = 'Pearson r(CRi-SR)'
# save res_df to csv
if len(pos_select) == 2:
    save_name = 'summary_corr_cri_sr.csv'
else:
    save_name = 'summary_corr_cri_sr_' + pos_select[0] + '.csv'
res_df.to_csv(save_path + save_name, index=True, na_rep='NA')

# %% WIP - THIS IS NOT HOW IT SHOULD BE DONE, CORR ARE NOT INDEPENDENT -> NEED COCOR
# Comparison CRi-SR correlations between each CRi
# pairwise comparison of correlation coefficients
# loop over rating dimensions
for dim in rat_dim:
    # create dataframe with results
    paircomp_df = pd.DataFrame(index=['cr_' + cri for cri in cri_list], columns=['cr_' + cri for cri in cri_list])
    # add level 2 columns: z and p
    multi_idx = pd.MultiIndex.from_product([paircomp_df.columns, ['z', 'p']])
    paircomp_df = pd.DataFrame(index=paircomp_df.index, columns=multi_idx)
    # loop over CRis
    for cri1 in cri_list:
        for cri2 in cri_list:
            # get correlation coefficients
            corr1 = res_df.loc['cr_' + cri1, (dim, 'r')]
            corr2 = res_df.loc['cr_' + cri2, (dim, 'r')]
            # get number of observations
            n = len(cri_sr)
            # calculate z and p
            z, p = independent_corr(corr1, corr2, n, n)
            # fill paircomp_df
            paircomp_df.loc['cr_' + cri1, ('cr_' + cri2, 'z')] = z.round(3)
            paircomp_df.loc['cr_' + cri1, ('cr_' + cri2, 'p')] = p.round(3)
    # save to csv
    if len(pos_select) == 2:
        save_name = 'comparison_corr_cri_sr_' + rat_dim[dim] + '.csv'
    else:
        save_name = 'comparison_corr_cri_sr_' + rat_dim[dim] + '_' + pos_select[0] + '.csv'
    paircomp_df.to_csv(save_path + save_name, index=True, na_rep='NA')

# %%
# Comparison CRi-SR correlation for Flubber between phase 1 and phase 2
# load phase 1 data
res_df_phase1 = pd.read_csv(phase1_path + filename_cri_phase1, header=0)
# A bit of formatting
# keep only row 0 (for check) and 2 and remove first column
res_df_phase1 = res_df_phase1.iloc[[0, 2], 1:]
# rename index
res_df_phase1.index = ['cri', 'phase1']
# reformat columns
res_df_phase1.columns = pd.MultiIndex.from_product([['valence', 'arousal', 'distance', 'angle'], ['mean', 'std', 'skew', 'kurtosis']])
# remove first row
res_df_phase1 = res_df_phase1.iloc[1:, :]
# convert to float
res_df_phase1 = res_df_phase1.astype('float64')
# phase 2 data: in res_df, keep only rows with cr_mean, cr_std, cr_skew, cr_kurtosis
res_df_phase2 = res_df.loc[cri_select, :]
# keep only r values
res_df_phase2 = res_df_phase2.xs('r', axis=1, level=1)
# create new dataframe for comparison
res_df = res_df_phase1
# create new empty row with index phase2
res_df.loc['phase2', :] = None
# fill phase2 row of res_df with phase2 data
res_df.loc['phase2', :] = res_df_phase2.values.T.flatten()

# add empty row with index 'difference between phase1 and phase2'
res_df.loc['diff_z', :] = None
res_df.loc['diff_p_val', :] = None
# %%
# loop over columns
for dim in res_df.columns.levels[0]:
    for cri in res_df.columns.levels[1]:
        # get correlation coefficients
        corr1 = res_df.loc['phase1', (dim, cri)]
        corr2 = res_df.loc['phase2', (dim, cri)]
        # get number of observations
        n1 = 204 # phase 1 - 51 subjects * 4 blocks (rating methods)
        n2 = len(cri_sr) # phase 2
        # calculate z and p
        z, p = independent_corr(corr1, corr2, n1, n2)
        # fill res_df
        res_df.loc['diff_z', (dim, cri)] = z.round(3)
        res_df.loc['diff_p_val', (dim, cri)] = p.round(3)
# save to csv
if len(pos_select) == 2:
    save_name = 'comparison_corr_cri_sr_phase1_phase2.csv'
else:
    save_name = 'comparison_corr_cri_sr_phase1_phase2_' + pos_select[0] + '.csv'
res_df.to_csv(save_path + save_name, index=True, na_rep='NA')


# %%
# Comparison CRi-SR correlation for Flubber between phase 1 and phase 2
# Compare different CRi together
# for phase 2, seated position only:
# highest CRi-SR correlation is cr_mean for valence, and cr_std for arousal
# compare cr_mean and cr_mean between phase 1 and phase 2 for valence
# compare cr_mean and cr_std between phase 1 and phase 2 for arousal

# valence
# get correlation coefficients
corr1 = res_df.loc['phase1', ('valence', 'mean')]
corr2 = res_df.loc['phase2', ('valence', 'mean')]
# get number of observations
n1 = 204 # phase 1 - 51 subjects * 4 blocks (rating methods)
n2 = len(cri_sr) # phase 2
# calculate z and p
z, p = independent_corr(corr1, corr2, n1, n2)

# arousal
# get correlation coefficients
corr1 = res_df.loc['phase1', ('arousal', 'mean')]
corr2 = res_df.loc['phase2', ('arousal', 'std')]
# get number of observations
n1 = 204 # phase 1 - 51 subjects * 4 blocks (rating methods)
n2 = len(cri_sr) # phase 2
# calculate z and p
z, p = independent_corr(corr1, corr2, n1, n2)


# %%
