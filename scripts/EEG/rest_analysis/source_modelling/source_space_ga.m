%% Load the results from individual subjects
% Load subject IDs from behavioral results
% 5 subjects were excluded + 1 pilot participant
T = readtable(paths.beh);
subject_ids = unique(T.ID(:))';
n_subjects = numel(subject_ids);

load(paths.eloreta);

%% Allocate space for combining the results
dfa_alpha_sensor = zeros(n_sensors, n_subjects);
power_alpha_sensor = zeros(n_sensors, n_subjects);
slope_sensor = zeros(n_sensors, n_subjects);

dfa_alpha_source = zeros(n_sources, n_subjects);
power_alpha_source = zeros(n_sources, n_subjects);
slope_source = zeros(n_sources, n_subjects);

%% Loop over subjects
n_processed = 0;
for i_subject = 1:n_subjects
    subject_id = subject_ids{i_subject};
    disp(subject_id);

    % Load single subject results
    filepath = [paths.data subject_id '_DFA_power_slope.mat'];
    if ~exist(filepath, 'file')
        warning([subject_id ': Could not open ' filepath]);
        continue
    end
    subject_data = load(filepath);
    
    % Store in the common matrix
    dfa_alpha_sensor(:, i_subject) = subject_data.dfa_alpha_sensor;
    power_alpha_sensor(:, i_subject) = subject_data.power_alpha_sensor;
    slope_sensor(:, i_subject) = subject_data.slope_sensor;

    dfa_alpha_source(:, i_subject) = subject_data.dfa_alpha_source;
    power_alpha_source(:, i_subject) = subject_data.power_alpha_source;
    slope_source(:, i_subject) = subject_data.slope_source;

    n_processed = n_processed + 1;
end
if n_processed < n_subjects
    warning('Some subjects were not processed');
end


%% Save the results
savefile = [paths.data 'results_DFA_power_slope.mat'];
fprintf('Saving the results to %s\n', savefile);
save(savefile, 'subject_ids', 'n_subjects', ...
    "dfa_alpha_sensor", "power_alpha_sensor", "slope_sensor", ...
    "dfa_alpha_source", "power_alpha_source", "slope_source");


%% Calculate grand-average
dfa_alpha_sensor_ga = mean(dfa_alpha_sensor, 2);
power_alpha_sensor_ga = mean(power_alpha_sensor, 2);
slope_sensor_ga = mean(slope_sensor, 2);

dfa_alpha_source_ga = mean(dfa_alpha_source, 2);
power_alpha_source_ga = mean(power_alpha_source, 2);
slope_source_ga = mean(slope_source, 2);


%% Save the grand average results
savefile = [paths.data 'results_ga_DFA_power_slope.mat'];
fprintf('Saving the results to %s\n', savefile);
save(savefile, 'subject_ids', 'n_subjects', ...
    "dfa_alpha_sensor_ga", "power_alpha_sensor_ga", "slope_sensor_ga", ...
    "dfa_alpha_source_ga", "power_alpha_source_ga", "slope_source_ga");


%% Export the single-subject source-space results to R
% Load behavioral and sensor space results
T = readtable(paths.beh);
columns_to_pick = [1:15 138:139];
T_beh = T(:, columns_to_pick);
dfa_names = cellfun(@(x) sprintf('DFA_%d', x), num2cell(1:n_sources), ...
    'UniformOutput', false);
slope_names = cellfun(@(x) sprintf('oneF_%d', x), num2cell(1:n_sources), ...
    'UniformOutput', false);
T_dfa = array2table(dfa_alpha_source', 'VariableNames', dfa_names);
T_dfa.ID = subject_ids';
T_slope = array2table(slope_source', 'VariableNames', slope_names);
T_slope.ID = subject_ids';
T_beh_dfa = innerjoin(T_beh, T_dfa, 'Keys', 'ID');
T_beh_dfa_slope = innerjoin(T_beh_dfa, T_slope, 'Keys', 'ID');

savefile = [paths.data 'source_dfa_slope.csv'];
fprintf('Saving the results for R to %s\n', savefile);
writetable(T_beh_dfa_slope, savefile);
