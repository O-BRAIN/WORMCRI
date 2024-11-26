#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 10:39:21 2024

@author: naherzog
"""
import pandas as pd
import matplotlib.pyplot as plt
import mne

#load data
file_path = '/data/p_02191/Analysis/Nadine/P300/per_trial_per_timepoint_P300.csv'
P3b = pd.read_csv(file_path)


#load random raw file to get epoch info
raw = mne.io.read_raw_eeglab('/data/p_02191/Analysis/Nadine/EEG_data/task/preprocessing/data/S088/S088post_ICA.set', preload=True)
local_triggers_dict = {'update': 'S 32', 'ignore': 'S 30', 'm1': 'S 31', 'm2': 'S 33'}
# Extract events and their ids
events, event_id = mne.events_from_annotations(raw)
# Create a dictionary with the correct mapping from string to integer
local_triggers_dict = {cond: event_id[trigger] for cond, trigger in local_triggers_dict.items()}
# create the Epochs object with the corrected dictionary
epochs = mne.Epochs(raw, events, event_id=local_triggers_dict, tmin=-0.5, tmax=1.8, baseline=(-0.5, -0.05), preload=True)


# Filter data for the 'update' condition
update_P3b = P3b[P3b['condition'] == 'update']
ignore_P3b = P3b[P3b['condition'] == 'ignore']

# Group by timepoint and calculate the mean N2_FC6_value
mean_erp_update_P3b = update_P3b.groupby('timepoint')['P300'].mean()*1000000
mean_erp_ignore_P3b = ignore_P3b.groupby('timepoint')['P300'].mean()*1000000

# Plot
plt.figure(figsize=(19, 8))
plt.plot(epochs.times[250:], mean_erp_update_P3b.values[250:], label='update', color='#8aaedcff', linewidth=4)
plt.plot(epochs.times[250:], mean_erp_ignore_P3b.values[250:], label='ignore', color='#ffc107ff', linewidth=4)
plt.fill_between(epochs.times[250:], 
                 min(mean_erp_update_P3b.values[250:])-0.05, 
                 max(mean_erp_ignore_P3b.values[250:])+0.05, 
                 where=((epochs.times[250:] >= 0.310) & (epochs.times[250:] <= 0.50)), 
                 color='#d3d3d3',  # Lighter gray color
                 alpha=0.8)  # Adjust alpha as needed

x_values = epochs.times[250::50]  
plt.xticks(x_values, fontsize=24)
y_values = [-2, -1, 0, 1, 2]  
plt.yticks(y_values, fontsize=24) 

plt.xlabel('Time (s)', fontsize=24)  # Increase font size of x-axis label
plt.ylabel('µV', fontsize=24)  # Increase font size of y-axis label
plt.title('P3b timeseries at Pz', fontsize=24, loc='left')  # Increase font size of title
plt.legend(fontsize=14)  # Increase font size of legend
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/P300_timeseries', dpi=300, format='png')
plt.show()


