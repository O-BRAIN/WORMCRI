# Import the FOOOF object
import fooof
from fooof import FOOOF
import mne
import os
import numpy as np
import pandas as pd


root = "path/to/where/data/is"
dirlist = [f for f in os.listdir(root) if f.endswith('.set')]
dirlist[1]


#initialize obejct to store 1/f in
#load signle subject to get electrode names
#raw = mne.io.read_raw_eeglab('/data/p_02191/Analysis/Nadine/EEG/rest/preprocessing/data/S069/S069post_ICA_rejected_data.set', preload=True)
#elec_list = raw.ch_names
all_subject_data = []

for file_name in dirlist:
    sub_data = []
    # Construct the full path to the .set file
    subfile = os.path.join(root, file_name)
    
    # Extract the subject identifier from the filename if needed
    sub = file_name.split(".set")[0]
    
    # Load the data
    print(f"Loading file: {subfile}")
    raw = mne.io.read_epochs_eeglab(subfile)
    elec_list = raw.ch_names
    #calculate power spectrum
    spectrum = raw.compute_psd(method="welch")
    data, freqs = spectrum.get_data(return_freqs=True)

    # Initialize FOOOF object
    fm = FOOOF()
    # Define frequency range across which to model the spectrum
    freq_range = [3, 40]

    # Model the 1/f power spectrum with FOOOF for each electrode
    for elec in range(len(elec_list)):#61
        # Average power spectrum across all segments for the current electrode
        avg_power_spectrum = data[:, elec, :].mean(axis=0)  # Shape: (1025,)
        fm.fit(freqs, avg_power_spectrum, freq_range)     # Fit the FOOOF model on the averaged power spectrum
        slope = fm.aperiodic_params_[1]     # Extract the 1/f slope for this electrode
        sub_data.append(slope)
        
        # Prepend subject ID to the sub_data list
    sub_data = [sub] + sub_data
    
    # Append data for each subject to the list `all_subject_data`
    all_subject_data.append(sub_data) 

# Create the DataFrame without specifying 'ID' initially
df2 = pd.DataFrame(all_subject_data)  # Only electrode names as columns for now
elec_list = ["oneF" + elec for elec in elec_list]
df2.columns = ['ID'] + elec_list      
df2.to_excel('/path/to/where/excel/should/be/saved.xlsx',index = False)


