function [paths, p, v] = source_space_config()
    %% Set up paths
    % Path to the behavioral data
    paths.beh = '/data/p_02191/Analysis/Nadine/final_scripts_and_data/data/for_analysis/all_data_behavioral.xlsx';

    % Path to the folder with preprocessed EEG data
    paths.preproc = '/data/p_02191/Analysis/Nadine/EEG_data/rest/preprocessing/data/';
    
    % Where to save the eLORETA inverse operator
    paths.eloreta = '/data/p_02490/Scripts/working_memory_Nadine/data_minimal/sa_eLoreta.mat';

    % Path to python (fooof should be installed)
    paths.python = '/data/u_kapralov_software/miniconda3/envs/fooof-for-matlab/bin/python3';

    % Where to save the intermediate results
    paths.data = '/data/p_02490/Scripts/working_memory_Nadine/data_minimal/';

    % Where to save the plots
    paths.plot = '/data/p_02490/Scripts/working_memory_Nadine/results_minimal/';               

    % Location of the scripts for DFA
    paths.toolbox.dfa = '/data/p_02191/Analysis/Nadine/final_scripts_and_data/scripts/EEG/EEG_rest_analysis/DFA/';

    % EEGLAB toolbox
    paths.toolbox.eeglab = '/data/p_02490/Toolboxes/eeglab2021.0/';

    % MATLAB wrapper for FOOOF
    paths.toolbox.fooof = '/data/p_02490/Scripts/working_memory_Nadine/toolboxes/fooof_mat/';
    
    % Scripts for source space analysis 
    % (authored by Guido Nolte and Stefan Haufe)
    paths.toolbox.haufe = '/data/p_02490/Scripts/working_memory_Nadine/toolboxes_minimal/haufe/';

    %% Set up parameters
    p = struct();

    % DFA
    p.dfa = struct();
    p.dfa.Fs = 500;                     % sampling frequency.
    p.dfa.DFA_SmallTime = 0.5; 		    % Smallest time window (in seconds) to be computed in DFA.
    p.dfa.DFA_LargeTime = 180;          % Largest time window (in seconds) to be computed in DFA.
    p.dfa.DFA_SmallTimeFit = 2; 		% Smallest time window (in seconds) to be include in the DFA fit.
    p.dfa.DFA_LargeTimeFit = 25;		% Largest time window (in seconds) to be include in the DFA fit.
    p.dfa.DFA_Overlap = 0.5;		    % Overlap between windows in DFA.
    p.dfa.DFA_Plot = 1;    
    
    % Filtering
    p.fir = struct();
    p.fir.hp = 8;                       % highpass frequency (Hz)
    p.fir.lp = 14;                      % lowpass frequency (Hz)
    p.fir.fir_order = 2 / p.fir.hp;     % filter order
    p.fir.Fs = p.dfa.Fs;                % sampling frequency
    
    % PSD (might be different from MNE)
    p.psd = struct();
    p.psd.nfft = 4 * p.dfa.Fs;          % 4 s windows -> 0.25 Hz frequency resolution
    p.psd.noverlap = 2 * p.dfa.Fs;      % 50 % overlap
    
    % FOOOF (should be the same as in Python)
    p.fooof = struct();
    p.fooof.fit_range = [3 40];         % in Hz
    p.fooof.fooof_settings = struct();  % use default
    p.fooof.return_model = 0;
    
    % Views to display for source space plots
    v = struct();
    v.empty = 0;
    v.leftlateral = 1;
    v.leftmedial = 2;
    v.rightlateral = 3;
    v.rightmedial = 4;
    v.dorsal = 5;
    v.dorsalhorizontal = 6;
    v.ventral = 7; 
    v.ventralhorizontal = 8;
    v.colorbar = 9;
end

