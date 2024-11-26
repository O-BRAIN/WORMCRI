#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 19 09:39:33 2023

@author: naherzog
"""
import mne
import os
import numpy as np
import glob
import pandas as pd

#create subject list
subject_dir = "/data/p_02191/Analysis/Nadine/EEG_data/task/preprocessing/data/SUBJECT/SUBJECTpost_ICA.set"
data_dir = "/data/p_02191/Analysis/Nadine/EEG_data/task/preprocessing/data/"
subjects = [folder for folder in os.listdir(data_dir) if folder.startswith("S")] #what is with S091 task EEG?

# exclude specific subjects from the subjects list (they had VERY noisy data)
exclude_subjects = ["S056", "S060", "S096", "S099", "S081"]
subjects = [subject for subject in subjects if subject not in exclude_subjects]

#EEG tiggers for update and ignore events
triggers_dict = {'update': 'S 32', 'ignore': 'S 30'}#, 'm1': 'S 31', 'm2': 'S 33'}

all_data = [] #initiallize the data frame to store the per_trial data in

#loop over all subjects to get at P300 amplitude per trial
for subject in subjects:
    local_triggers_dict = triggers_dict.copy()
    path_pattern = os.path.join(subject_dir.replace("SUBJECT", subject))  #create subject-specific data path
    matched_files = glob.glob(path_pattern)

    if matched_files:
        file_path = matched_files[0]
        raw = mne.io.read_raw_eeglab(file_path, preload=True) #load raw eeg data (it was preprocessed already, i.e. filtered and ICA components rejected)
    else:
        print(f"No .set file found for {subject} at {path_pattern}")
    #continue

    events, event_id = mne.events_from_annotations(raw)
    for key, value in local_triggers_dict.items():
        local_triggers_dict[key] = event_id[f"{value}"]
    epochs = mne.Epochs(raw, events, event_id=local_triggers_dict, tmin=-0.5, tmax=1.8, baseline=(-0.5, -0.05), preload=True) #epoch data based on condition events (ignore and update)

    correctness_dict = {'correct': 'S 61', 'incorrect': 'S 60'}   #initialize EEG triggers for when the trial was correct vs. when it was incorrect
    for key, value in correctness_dict.items():
        correctness_dict[key] = event_id[value]
    conditions = [] 
    correctness_list = []

    #loop over correct vs. incorrect trigger events to find if trial was correct or incorrect
    sampling_rate = raw.info['sfreq']
    time_limit_samples = int(6 * sampling_rate)
    for i, event in enumerate(events[:-1]):
        if event[2] in local_triggers_dict.values():
            condition = [k for k, v in local_triggers_dict.items() if v == event[2]][0]
            conditions.append(condition)
            for j in range(i+1, len(events)):
                if (events[j][0] - event[0]) <= time_limit_samples:
                    if events[j][2] == correctness_dict['correct']:
                        correctness_list.append('correct')
                        break
                    elif events[j][2] == correctness_dict['incorrect']:
                        correctness_list.append('incorrect')
                        break
                else:
                    correctness_list.append('unknown')
                    break

    assert len(conditions) == len(epochs), "Conditions list length doesn't match number of epochs!"
    assert len(correctness_list) == len(epochs), "Correctness list length doesn't match number of epochs!"
    epochs.metadata = pd.DataFrame({'condition': conditions, 'correctness': correctness_list}) #append if trial was correct or incorrect to epoch info metadata

    #custom epoch rejection: 
        #does more than 30% of the samples in a single channel exceed an amplitude threshold of 15microvolt?
        #does this happen in more than 2 channels?
    threshold = 15e-6
    percentage_limit = 0.3
    channel_limit = 2
    bad_epochs = []
    n_samples = len(epochs.times)

    for i, epoch in enumerate(epochs):
        n_over_threshold_per_channel = np.sum(np.abs(epoch) > threshold, axis=1)
        bad_channels = np.where(n_over_threshold_per_channel > percentage_limit * n_samples)[0]
        if len(bad_channels) > channel_limit:
            bad_epochs.append(i)
    epochs.drop(indices=bad_epochs, reason='Custom amplitude rejection')
    print(f"Dropped {len(bad_epochs)} epochs based on custom criteria.")
    
    epochs.apply_baseline(baseline= (-0.5, -0.05)) #apply baseline correction to epoch
    
    #loop over epochs to append subID, condition and correctness to each epoch (which is a single trial)
    #also extract average amplitude at 310 to 500 msec per epoch (trial)
    for i, metadata in enumerate(epochs.metadata.iterrows()):
        this_epoch = epochs[i:i+1].copy().pick_channels(['Pz']).get_data()[0, :]
        
        row_data = {
            'subject_id': subject,
            'condition': metadata[1]['condition'],
            'correctness': metadata[1]['correctness']
        }

        #extract average P300 amplitude at 310 to 500 msec
        start_index = np.searchsorted(epochs.times, 0.310) # Find the indices in epochs.times that correspond to the start and end times
        end_index = np.searchsorted(epochs.times, 500)
        selected_data = this_epoch[0, start_index:end_index+1] # Extract the data from this_epoch within the time range
        row_data['P300'] = float(np.mean(selected_data))  # Store the peak values in your row data

        all_data.append(row_data)

# Convert the accumulated results into a DataFrame and save
data_rows = []
for trial in all_data:
    condition = trial['condition']
    correctness = trial['correctness']
    subject_id = trial['subject_id']
    P300 = trial['P300']
    data_rows.append([subject_id, condition, correctness, P300])

# Create a DataFrame from the data
df = pd.DataFrame(data_rows, columns=['subject_id', 'condition', 'correctness', 'P300'])

# Convert 'correctness' to a binary variable (assuming 'correct' is the positive class)
df['correctness'] = df['correctness'].map({'correct': 1, 'incorrect': 0})
# Optionally, to save as a CSV
df.to_csv('/data/p_02191/Analysis/Nadine/per_trial_per_timepoint_P300.csv')




