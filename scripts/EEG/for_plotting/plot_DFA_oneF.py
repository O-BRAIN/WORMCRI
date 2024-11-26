#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 29 10:43:50 2023

@author: naherzog
"""
import matplotlib.pyplot as plt
import mne
import scipy.io
import pandas as pd



# load random raw file to have the info
raw = mne.io.read_raw_eeglab('/data/p_02191/Analysis/Nadine/EEG/task/preprocessing/data/S088/S088post_ICA.set', preload=True)
local_triggers_dict = {'update': 'S 32', 'ignore': 'S 30', 'm1': 'S 31', 'm2': 'S 33'}
# Extract events and their ids
events, event_id = mne.events_from_annotations(raw)
# Create a dictionary with the correct mapping from string to integer
local_triggers_dict = {cond: event_id[trigger] for cond, trigger in local_triggers_dict.items()}
# create the Epochs object with the corrected dictionary
epochs = mne.Epochs(raw, events, event_id=local_triggers_dict, tmin=-0.5, tmax=1.8, baseline=(-0.5, -0.05), preload=True)
# Assuming 'info' is your EEG data layout from the MNE object
info = epochs.info  # Replace 'epochs' with your actual epochs or raw data object

##################################################  plot avg DFA   ########################################################
# load avg. DFA data
mat = scipy.io.loadmat('/data/p_02191/Analysis/Nadine/EEG/rest/DFA/avgDFA.mat')
DFA = mat['DFA']
DFA.shape
# Reshape the DFA data to have shape (61,)
DFA_reshaped = DFA.reshape(61,)


# Create the topomap plot
import matplotlib.colors as mcolors
colors = colors = [(0.9, 0.9, 0.9), (0.4, 0.4, 0.4)]  # Light gray to dark gray
cmap_name = 'custom_gray'
cmap = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=8)  # Specify the number of levels

# Create the topomap plot for DFA
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(DFA_reshaped, info,cmap=cmap,
                             mask_params=dict(marker='o', markerfacecolor='none', markeredgecolor='w', markeredgewidth=4, markersize=15),
                             sensors='k.',
                             contours=0,
                             axes=ax, show=False)

im.set_clim(vmin=min(DFA_reshaped), vmax=max(DFA_reshaped))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('exponent', fontsize=24)  # Adjust the label font size here
cbar.ax.tick_params(labelsize=24)  # Adjust the tick font size here

# Add a title indicating the time window
ax.set_title('average DFA alpha', fontsize=24)
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/DFA_avg.png', dpi=300, bbox_inches='tight')
plt.show()


##################################################  plot sig DFA t-val  ########################################################
# load T-val DFA data
data = pd.read_excel("/data/p_02191/Analysis/Nadine/result_matrix_DFA.xlsx")
DFA = data['T_val_upt_ign']

# significance mask
mask = data["sig"] == "***"


# Create the topomap plot
import matplotlib.colors as mcolors
colors = [
    (0.9, 0.9, 0.9),  # gray 
    (0.95, 0.85, 0.2),        # yellow
    (1, 0.6, 0),        # orange
    (1, 0, 0),          # red
    (0.8, 0, 0.15)      # dark red
]

# Create the custom colormap
cmap_name = 'custom_jet_gray'
cmap = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N = 10)

# Create the topomap plot for DFA
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(DFA, info,mask = mask, cmap = cmap,
                             mask_params=dict(marker='o', markerfacecolor='white', markeredgecolor='black', markeredgewidth=3, markersize=15),
                             sensors='k.',
                             contours=0,
                             axes=ax, show=False)

im.set_clim(vmin=1.5, vmax=max(DFA))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('T-value', fontsize=24)  # Adjust the label font size here
cbar.ax.tick_params(labelsize=24)  # Adjust the tick font size here
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/DFA_tvals.png', dpi=300, bbox_inches='tight')
plt.show()




#########################################   average 1/f slope over subjects   ########################
# load 1/f slope data
mat = scipy.io.loadmat('/data/p_02191/Analysis/Nadine/EEG/rest/1overf/avg_1f.mat')
oneF = mat['avg1f1']
oneF = oneF[:, 0]
# Reshape the DFA data to have shape (61,)
oneF_reshaped = oneF.reshape(61,)

colors = colors = [(0.9, 0.9, 0.9), (0.2, 0.2, 0.2)]  # Light gray to dark gray
cmap_name = 'custom_gray'
cmap = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N=8)  # Specify the number of levels


# Create the topomap plot
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(oneF_reshaped, info, cmap=cmap,
                             mask_params=dict(marker='o', markerfacecolor='none', markeredgecolor='w', markeredgewidth=4, markersize=15),
                             sensors='k.',
                             contours=0,
                             axes=ax, show=False)


im.set_clim(vmin=min(oneF_reshaped), vmax=max(oneF_reshaped))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('exponent', fontsize=24)  # Adjust the label font size here
cbar.ax.tick_params(labelsize=24)  # Adjust the tick font size here
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/oneF_avg.png', dpi=300, bbox_inches='tight')
plt.show()



############################################          tvals 1/f slope       ##########################
data = pd.read_excel("/data/p_02191/Analysis/Nadine/results/result_matrix_1f_AA.xlsx")
oneF = data['t_val']
# significance mask
mask = data["sig"] == "***"

colors = [
    (0.9, 0.9, 0.9),  # gray 
    (0.95, 0.85, 0.2),        # yellow
    (1, 0.6, 0),        # orange
    (1, 0, 0),          # red
    (0.8, 0, 0.15)      # dark red
]

# Create the custom colormap
cmap_name = 'custom_jet_gray'
cmap = mcolors.LinearSegmentedColormap.from_list(cmap_name, colors, N = 10)


# Create the topomap plot
fig, ax = plt.subplots(figsize=(7, 5))
im, _ = mne.viz.plot_topomap(oneF, info, mask=mask, cmap=cmap,
                             mask_params=dict(marker='o', markerfacecolor='white', markeredgecolor='black', markeredgewidth=3, markersize=15),
                             sensors='k.',
                             contours=0,
                             axes=ax, show=False)


im.set_clim(vmin=1.75, vmax=max(oneF))
# Add colorbar with the correct mappable
cbar = plt.colorbar(im, ax=ax)
cbar.set_label('T-value', fontsize=24)  # Adjust the label font size here
cbar.ax.tick_params(labelsize=24)  # Adjust the tick font size here
# Add a title indicating the time window
#ax.set_title('average 1/f slope', fontsize=24)
plt.savefig('/data/p_02191/Analysis/Nadine/results/figures/png/oneF_AA_tvals.png', dpi=300, bbox_inches='tight')
plt.show()



