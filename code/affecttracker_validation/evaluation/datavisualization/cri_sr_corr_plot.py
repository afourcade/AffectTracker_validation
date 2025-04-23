########################################################################################################################
# Script to plot correlation CRi-SR from AVR experiment phase 1 and 2
# Need csv file containing CRi and SR from both phases
# Output: figures of CRi-SR correlations
# Author: Antonin Fourcade
# Last version: 10.10.2024
########################################################################################################################

# %%
# import packages
import pandas as pd
import numpy as np
import os
import plotnine as plt9

# %%
# Set paths and experiment parameters
phases = ['phase1', 'phase2']
phases_path = {'phase1':'E:/AffectiveVR/Phase_1/',
                'phase2':'E:/AffectiveVR/Phase_2/'}
data_path = {'phase1':'E:/AffectiveVR/Phase_1/Data/AVR/',
                'phase2':'E:/AffectiveVR/Phase_2/Data/AVR/'} # path to data

# path to CR csv file and load data
filename_cri_sr = 'cri_sr_clean.csv'

# save path
save_path = phases_path['phase2'] + 'affectivevr/corr_results/'
# create folders if not there
if not os.path.exists(save_path): # create folder if it doesn't exist
    os.makedirs(save_path)

# load data
cri_sr = {'phase1': pd.read_csv(data_path['phase1'] + filename_cri_sr),
          'phase2': pd.read_csv(data_path['phase2'] + filename_cri_sr)}

# get list of participants
sub_list = {'phase1': cri_sr['phase1']['sj_id'].unique(),
            'phase2': cri_sr['phase2']['sub'].unique()}

# CRis of interest
cri_select = {'phase1': 'cr_mean',
              'phase2': 'cr_mean'}

# dimensions of interest
dim_select = {'valence': 'v', 'arousal': 'a'}

# position of interest (phase 2)
pos_select = ['seated'] # ['seated', 'standing']

# labels for correlation values
# manually write the results here because lazy...
# be careful with order...
# cor_labels = {'valence': ['r = 0.922', 'r = 0.950', 'r = 0.920', 'r = 0.215'],
#               'arousal': ['r = 0.873', 'r = 0.938', 'r = 0.852', 'r = 0.510']} # for all positions
cor_labels = {'valence': ['r = 0.922', 'r = 0.950', 'r = 0.920', 'r = 0.186'],
              'arousal': ['r = 0.873', 'r = 0.938', 'r = 0.852', 'r = 0.559']} # for seated only

# %%
# loop over dimensions
for dim in dim_select.keys():
    # subset data for CRi and dimension of interest
    cri_sr_select = {'phase1': cri_sr['phase1'][['sj_id', 'test_site', 'quadrant', 'rating_method', cri_select['phase1'] + '_' + dim_select[dim], 'sr_' + dim_select[dim]]],
                     'phase2': cri_sr['phase2'][['sub', 'test_site', 'position', cri_select['phase1'] + '_' + dim_select[dim], 'sr_' + dim_select[dim]]]}
    # phase 1 - exclude rating method = baseline
    cri_sr_select['phase1'] = cri_sr_select['phase1'][cri_sr_select['phase1']['rating_method'] != 'Baseline']
    # phase 2 - select only position of interest 
    cri_sr_select['phase2'] = cri_sr_select['phase2'][cri_sr_select['phase2']['position'].isin(pos_select)]
    # phase 2 - add rating_method column = 'Flubber - Evaluation' and quadrant column = 'Sequence phase 2'
    cri_sr_select['phase2']['rating_method'] = 'Flubber'
    cri_sr_select['phase2']['quadrant'] = 'Evaluation Sequence'
    # add phase column
    cri_sr_select['phase1']['phase'] = 'Selection'
    cri_sr_select['phase2']['phase'] = 'Evaluation'
    
    # plot correlation between CRi and SR
    plot_cri_sr = (
        plt9.ggplot() 
        + plt9.geom_point(cri_sr_select['phase1'], plt9.aes(x=cri_select['phase1'] + '_' + dim_select[dim], y='sr_' + dim_select[dim], fill='rating_method', shape='quadrant'))
        + plt9.geom_smooth(cri_sr_select['phase1'],plt9.aes(x=cri_select['phase1'] + '_' + dim_select[dim], y='sr_' + dim_select[dim], color='rating_method', linetype='phase'), method='lm')
        + plt9.geom_point(cri_sr_select['phase2'], plt9.aes(x=cri_select['phase2'] + '_' + dim_select[dim], y='sr_' + dim_select[dim], fill='rating_method', shape='quadrant'))
        + plt9.geom_smooth(cri_sr_select['phase2'], plt9.aes(x=cri_select['phase2'] + '_' + dim_select[dim], y='sr_' + dim_select[dim], color='rating_method', linetype='phase'), method='lm')
        + plt9.annotate("text",x=-0.7,y=[0.9, 0.75, 0.6, 0.45],label=cor_labels[dim], color = ['#00BA38', '#F8766D', '#619CFF', '#00BA38'], size = 20)
        + plt9.scale_fill_manual(values = ['#00BA38', '#F8766D', '#619CFF']) 
        + plt9.scale_color_manual(values = ['#00BA38', '#F8766D', '#619CFF'])
        + plt9.scale_linetype_manual(values = ['solid', 'dashed'])
        + plt9.labs(x= 'CR mean', y= 'SR', fill = 'Feedback', shape = '360° Video', color = 'Feedback', linetype = 'Phase') 
        + plt9.ggtitle(dim.capitalize())
        + plt9.theme_bw(base_size=20)
        + plt9.theme(plot_title = plt9.element_text(hjust = 0.5), plot_background = plt9.element_blank(), legend_background = plt9.element_blank(), svg_usefonts=True)
    )
    # save plot in png and svg
    if len(pos_select) > 1:
        plot_cri_sr.save(filename=save_path + 'corr_cri_sr_' + dim + '_all_position.png', height=10, width=15, units='in', dpi=300)
        plot_cri_sr.save(filename=save_path + 'corr_cri_sr_' + dim + '_all_position.svg', height=10, width=15, units='in', dpi=600)
    else:
        plot_cri_sr.save(filename=save_path + 'corr_cri_sr_' + dim + '_' + pos_select[0] + '.png', height=10, width=15, units='in', dpi=300)
        plot_cri_sr.save(filename=save_path + 'corr_cri_sr_' + dim + '_' + pos_select[0] + '.svg', height=10, width=15, units='in', dpi=600)
    
# %%
