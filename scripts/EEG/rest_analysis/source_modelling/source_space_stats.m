%% Initialize the environment
% Paths and analysis-specific parameters
[paths, p, v] = source_space_config();

% Scripts and data for source space visualization
addpath(paths.toolbox.haufe);
load(paths.eloreta);
load cm17;

%% Significant results - DFA * condition

% DFA * condition (before FDR)
T_dfa_source = readtable([paths.data 'results_matrix_DFA.csv']);

h = plot_stats(T_dfa_source.p_DFA < 0.01, sa, v, [0 1], cm17a);
sgtitle('eLORETA | fixed | DFA * condition | no FDR');
exportgraphics(h, [paths.plot 'significant_sources_DFA_condition.png'], ...
    "Resolution", 500);
close(h);

num_voxels_roi_DFA = count_significant(T_dfa_source.p_DFA < 0.01, sa);
fileID = fopen([paths.plot 'significant_sources_DFA_before_FDR.txt'], 'w+');
fprintf(fileID, 'DFA * condition - significant sources before FDR correction:\n%s\n', ...
    list_significant(num_voxels_roi_DFA, sa));
fclose(fileID);

% DFA * condition (after FDR)
h = plot_stats(T_dfa_source.p_DFA_FDR < 0.05, sa, v, [0 1], cm17a);
sgtitle('eLORETA | fixed | DFA * condition | after FDR');
exportgraphics(h, [paths.plot 'significant_sources_DFA_condition_FDR.png'], ...
    "Resolution", 500);
close(h);

num_voxels_roi_DFA_FDR = count_significant(T_dfa_source.p_DFA_FDR < 0.05, sa);
fileID = fopen([paths.plot 'significant_sources_DFA_after_FDR.txt'], 'w+');
fprintf(fileID, 'DFA * condition - significant sources after FDR correction:\n%s\n', ...
    list_significant(num_voxels_roi_DFA_FDR, sa));
fclose(fileID);

% DFA * condition (t-values)
lim = max(abs(T_dfa_source.t_update_DFA));
h = plot_stats(T_dfa_source.t_update_DFA, sa, v, [-lim lim], cm17);
sgtitle('eLORETA | fixed | DFA * condition | t-value, update');
exportgraphics(h, [paths.plot 'tvalue_update_DFA_condition.png'], ...
    "Resolution", 500);
close(h);


%% Significant results - slope * AA ratio * condition

% slope * AA * condition (before FDR)
T_slope_source = readtable([paths.data 'results_matrix_oneF.csv']);

h = plot_stats(T_slope_source.p_oneF_AA < 0.01, sa, v, [0 1], cm17a);
sgtitle('eLORETA | fixed | slope * AAratio * condition | no FDR');
exportgraphics(h, [paths.plot 'significant_sources_slope_condition_AAratio.png'], ...
    "Resolution", 500);
close(h);

num_voxels_roi_oneF_AA = count_significant(T_slope_source.p_oneF_AA < 0.01, sa);
fileID = fopen([paths.plot 'significant_sources_oneF_AA_before_FDR.txt'], 'w+');
fprintf(fileID, 'slope * AAratio * condition - significant sources before FDR correction:\n%s\n', ...
    list_significant(num_voxels_roi_oneF_AA, sa));
fclose(fileID);

% slope * AA * condition (after FDR)
h = plot_stats(T_slope_source.p_oneF_AA_FDR < 0.05, sa, v, [0 1], cm17a);
sgtitle('eLORETA | fixed | slope * AAratio * condition | after FDR');
exportgraphics(h, [paths.plot 'significant_sources_slope_condition_AAratio_FDR.png'], ...
    "Resolution", 500);
close(h);

num_voxels_roi_oneF_AA_FDR = count_significant(T_slope_source.p_oneF_AA_FDR < 0.05, sa);
fileID = fopen([paths.plot 'significant_sources_oneF_AA_after_FDR.txt'], 'w+');
fprintf(fileID, 'slope * AAratio * condition - significant sources after FDR correction:\n%s\n', ...
    list_significant(num_voxels_roi_oneF_AA_FDR, sa));
fclose(fileID);

% slope * AA * condition (t-values)
lim = max(abs(T_slope_source.t_update_oneF_AA));
h = plot_stats(T_slope_source.t_update_oneF_AA, sa, v, [-lim lim], cm17);
sgtitle('eLORETA | fixed | slope * AAratio * condition | t-value, update');
exportgraphics(h, [paths.plot 'tvalue_update_slope_condition_AAratio.png'], ...
    'Resolution', 500);
close(h);
