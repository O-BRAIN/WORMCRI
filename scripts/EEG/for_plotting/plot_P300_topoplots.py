#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 12:30:41 2024

@author: naherzog
"""

#read data
import pandas as pd
import mne
import matplotlib.pyplot as plt

excel_file_path = '/data/p_02191/Analysis/Nadine/EEG_data/task/cleaned_data/P300/P300_per_electrode.xlsx'
df = pd.read_excel(excel_file_path)

# load random raw file to have the info
raw = mne.io.read_raw_eeglab('/data/p_02191/Analysis/Nadine/EEG_data/task/preprocessing/data/S088/S088post_ICA.set', preload=True)
local_triggers_dict = {'update': 'S 32', 'ignore': 'S 30', 'm1': 'S 31', 'm2': 'S 33'}
# Extract events and their ids
events, event_id = mne.events_from_annotations(raw)
# Create a dictionary with the correct mapping from string to integer
local_triggers_dict = {cond: event_id[trigger] for cond, trigger in local_triggers_dict.items()}
# create the Epochs object with the corrected dictionary
epochs = mne.Epochs(raw, events, event_id=local_triggers_dict, tmin=-0.5, tmax=1.8, baseline=(-0.5, -0.05), preload=True)
# Assuming 'info' is your EEG data layout from the MNE object
info = epochs.info  # Replace 'epochs' with your actual epochs or raw data object



#which data
ignore = df['ignore_p3b_erp_avg']*1000000
update = df['update_p3b_erp_avg']*1000000
#data = data1-data2
cmap = plt.cm.get_cmap('RdBu_r', 15)

#########################################################     topoplot IGNORE      ######################################
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(ignore, epochs.info, cmap=cmap,
                             sensors='k.', contours=0,
                             mask_params=dict(marker='o', markerfacecolor='w',
                                              markeredgecolor='k', linewidth=0, markersize=10),
                             axes=ax, show=False)

im.set_clim(vmin=min(ignore), vmax=max(ignore))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('µV', fontsize=38)   # Adjust the label font size here
cbar.ax.tick_params(labelsize=38)   # Adjust the tick font size here 
ax.set_title('ignore', fontsize=38) # add title
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/ignore_P300.png', dpi=300, format='png')   #save


#########################################################     topoplot UPDATE      ######################################
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(update, epochs.info, cmap=cmap,
                             sensors='k.', contours=0,
                             mask_params=dict(marker='o', markerfacecolor='w',
                                              markeredgecolor='k', linewidth=0, markersize=10),
                             axes=ax, show=False)

im.set_clim(vmin=min(ignore), vmax=max(ignore))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('µV', fontsize=38)  # Adjust the label font size here
cbar.ax.tick_params(labelsize=38)  # Adjust the tick font size here 

# Add a title indicating the time window
ax.set_title('update', fontsize=38)
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/update_P300.png', dpi=300, format='png')
plt.show()

#########################################################     topoplot difference     ######################################
data = ignore-update
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(data, epochs.info, cmap=cmap,
                             sensors='k.', contours=0,
                             mask_params=dict(marker='o', markerfacecolor='w',
                                              markeredgecolor='k', linewidth=0, markersize=10),
                             axes=ax, show=False)

im.set_clim(vmin=min(data), vmax=max(data))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('µV', fontsize=38)  # Adjust the label font size here
cbar.ax.tick_params(labelsize=38)  # Adjust the tick font size here 

# Add a title indicating the time window
ax.set_title('ignore-update', fontsize=38)

