%% Initialize the environment
% Paths and analysis-specific parameters
[paths, p, v] = source_space_config();

% MATLAB scripts and toolboxes
addpath(paths.toolbox.dfa);             % Scripts for calculating DFA
addpath(paths.toolbox.eeglab);          % EEGLAB
addpath(genpath(paths.toolbox.fooof));  % MATLAB wrapper for FOOOF
addpath(paths.toolbox.haufe);           % Scripts for source space analysis
eeglab

% Activate Python within MATLAB
pe = pyenv;
if strcmp(pe.Executable, paths.python) && pe.Status == 0
    terminate(pe);
end

% Executing Python code in a separate process protects from crashing
% the whole MATLAB if Python crashes
pyenv('Version', paths.python, 'ExecutionMode', 'OutOfProcess');


%% Prepare the eLORETA inverse operator
if ~exist(paths.eloreta, "file")
    prepare_eLORETA_inverse_operator;
end


%% Run source space analysis for each subject
num_workers = 8;   % parallelize
subject_ids = dir([paths.preproc 'S*']);
n_subjects = numel(subject_ids);
parfor (i = 1:n_subjects, num_workers)
    subject_id = subject_ids(i).name;
    source_space_analysis(subject_id, paths, p);
end


%% Combine the results for all subjects
source_space_ga;