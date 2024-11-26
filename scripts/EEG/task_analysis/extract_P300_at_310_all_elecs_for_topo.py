import mne
import os
import numpy as np
import glob
import pandas as pd

subject_dir = "/data/p_02191/Analysis/Nadine/EEG_data/task/preprocessing/data/SUBJECT/SUBJECTpost_ICA.set"
data_dir = "/data/p_02191/Analysis/Nadine/EEG_data/task/preprocessing/data/"
out_dir = "/data/p_02191/Analysis/Nadine/EEG_data/task/cleaned_data/python/"
subjects = [folder for folder in os.listdir(data_dir) if folder.startswith("S")]
#subjects = ['S088', 'S011', 'S066', 'S054', 'S023', 'S029', 'S008',]
# Exclude specific subjects from the subjects list (they had VERY noisy data)
exclude_subjects = ["S056", "S060", "S096", "S099", "S081"]
subjects = [subject for subject in subjects if subject not in exclude_subjects]


triggers_dict = {'update': 'S 32', 'ignore': 'S 30'}

all_data = []
counter = 0

for subject in subjects:
    counter += 1
    local_triggers_dict = triggers_dict.copy()
    path_pattern = os.path.join(subject_dir.replace("SUBJECT", subject))
    matched_files = glob.glob(path_pattern)

    if matched_files:
        file_path = matched_files[0]
        raw = mne.io.read_raw_eeglab(file_path, preload=True)
    else:
        print(f"No .set file found for {subject} at {path_pattern}")
    #continue

    events, event_id = mne.events_from_annotations(raw)
    for key, value in local_triggers_dict.items():
        local_triggers_dict[key] = event_id[f"{value}"]
        
    probe_triggers = [event_id['S 50'], event_id['S 51'], event_id['S 52']]
    response_triggers = [event_id['S 58'], event_id['S 59']]
        
    epochs = mne.Epochs(raw, events, event_id=local_triggers_dict, tmin=-0.5, tmax=1.8, baseline=(-0.5, -0.05), preload=True)
   
    correctness_dict = {'correct': 'S 61', 'incorrect': 'S 60'}
   
    for key, value in correctness_dict.items():
        correctness_dict[key] = event_id[value]
    conditions = []
    correctness_list = []
    reaction_times = []
    trial_count_list = []
    trial_count = 0
    sampling_rate = raw.info['sfreq']
    time_limit_samples = int(6 * sampling_rate)
    
    for i, event in enumerate(events[:-1]):
        if event[2] in local_triggers_dict.values():
            condition = [k for k, v in local_triggers_dict.items() if v == event[2]][0]
            conditions.append(condition)
            probe_time = None
            response_time = None
            trial_count = trial_count +1
            trial_count_list.append(trial_count)
            

            for j in range(i+1, len(events)):
                if (events[j][0] - event[0]) <= time_limit_samples:
                    if events[j][2] == correctness_dict['correct']:
                        correctness_list.append('correct')
                        break
                    elif events[j][2] == correctness_dict['incorrect']:
                        correctness_list.append('incorrect')
                        break
                   
                    # Check for probe and response triggers to calculate reaction time
                    if events[j][2] in probe_triggers:
                        probe_time = events[j][0]
                    elif events[j][2] in response_triggers:
                        response_time = events[j][0]
                   
                else:
                    correctness_list.append('unknown')
                    break
           
            # Calculate and store reaction time if both probe and response are detected
            if probe_time and response_time:
                reaction_time = (response_time - probe_time) / sampling_rate  # Convert to seconds
                reaction_times.append(reaction_time)
            else:
                reaction_times.append(None)  # Placeholder for missing data
            

    epochs.metadata = pd.DataFrame({'condition': conditions, 'correct': correctness_list, 'reaction_time': reaction_times,"trial_numb":trial_count_list})

    threshold = 15e-6
    percentage_limit = 0.3
    channel_limit = 2
    bad_epochs = []
    n_samples = len(epochs.times)
    epoch_quality = []  # Initialize a list to store quality labels
    
    for i, epoch in enumerate(epochs):
        n_over_threshold_per_channel = np.sum(np.abs(epoch) > threshold, axis=1)
        bad_channels = np.where(n_over_threshold_per_channel > percentage_limit * n_samples)[0]
        if len(bad_channels) > channel_limit:
            bad_epochs.append(i)
    epochs.drop(indices=bad_epochs, reason='Custom amplitude rejection')
    print(f"Dropped {len(bad_epochs)} epochs based on custom criteria.")

    epochs = epochs.apply_baseline((-0.5, -0.05))  # Baseline correction from -0.5 to -0.05 seconds
     
    # Define the time ranges 
    p3b_start_index = np.argmin(np.abs(epochs.times - 0.310))
    p3b_end_index = np.argmin(np.abs(epochs.times - 0.500))

    # Compute average ERP for condition 'update' 
    update_epochs = epochs['update']
    p3b_update_erp = np.mean(update_epochs.get_data()[:, :, p3b_start_index:p3b_end_index], axis=(0, 2))

    # Compute average ERP for condition 'ignore' for P3b
    ignore_epochs = epochs['ignore']
    p3b_ignore_erp = np.mean(ignore_epochs.get_data()[:, :, p3b_start_index:p3b_end_index], axis=(0, 2))

    all_data.append({
        'subject_id': subject,
        'update_p3b_erp': p3b_update_erp.tolist(),
        'ignore_p3b_erp': p3b_ignore_erp.tolist()})


# Average ERPs across subjects for each condition
update_p3b_erp_avg = np.mean([subject['update_p3b_erp'] for subject in all_data], axis=0)
ignore_p3b_erp_avg = np.mean([subject['ignore_p3b_erp'] for subject in all_data], axis=0)

### SAVE

import pandas as pd

# Create a DataFrame with the ERP averages
data = {
    'update_p3b_erp_avg': update_p3b_erp_avg.flatten(),
    'ignore_p3b_erp_avg': ignore_p3b_erp_avg.flatten()
}

df = pd.DataFrame(data)

# Save DataFrame to Excel
excel_file_path = '/data/p_02191/Analysis/Nadine/EEG_data/task/cleaned_data/P300/P300_per_electrode.xlsx'
df.to_excel(excel_file_path, index=False)








