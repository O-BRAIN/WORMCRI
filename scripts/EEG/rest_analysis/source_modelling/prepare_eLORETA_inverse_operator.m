%%% 
% Pre-compute eLORETA inverse operator
%
% Outputs:
% 1. sa_eLoreta.mat with the following variables:
%
%   clab - labels of all channels
%   sa - source analysis structure
%   L_free, L_fixed - lead field (free and fixed orientation)
%   n_sensors - number of channels
%   n_voxels - number of voxels
%   A_eloreta_free, A_eloreta_fixed - eLORETAinverse operator (free / fixed
%     dipole orientations)
%%%

%% Paths
% Take channel locations from any participant
filepath = [paths.preproc 'S001/'];
filename = 'S001post_ICA_rejected_data.set';

% Where to save the results
output_path = paths.eloreta;

%% Prepare the leadfield matrix
% Load channel locations
EEG = pop_loadset('filename', filename, 'filepath', filepath);
all_chanlocs = EEG.chanlocs;
clab = {EEG.chanlocs(:).labels};
n_sensors = numel(all_chanlocs);

% Prepare the source analysis structure
sa = prepare_sourceanalysis(clab, 'nyhead', struct('newlocs', 1));

% Lead field for 5K vertices restricted to the cortical surface with either
% free (L_free) or fixed (L_fixed) dipole orientations
sa.voxels_5K_cort = sa.cortex5K.in_from_cortex75K(sa.cortex5K.in_cort);
L_free = sa.cortex75K.V_fem(:, sa.voxels_5K_cort, :);
L_fixed = sa.cortex75K.V_fem_normal(:, sa.voxels_5K_cort);
[~, n_sources, ~] = size(L_free);

% Prepare the common average reference operator
M = length(sa.clab_electrodes);
H = eye(M) - ones(M) / M;

% Apply common average reference to match preprocessing
L_free = reshape(H * reshape(L_free, M, []), M, n_sources, 3);
L_fixed = H * L_fixed;

%% Prepare the eLoreta inverse operator
lambda = 0.05;    % regularization parameter

fprintf('eLORETA for free dipole orientations...  ');
A_eloreta_free = mkfilt_eloreta2(L_free, lambda);
fprintf('Done\n');

fprintf('eLORETA for fixed dipole orientations... ');
A_eloreta_fixed = mkfilt_eloreta2(L_fixed, lambda);
fprintf('Done\n');

fprintf('Saving the results...                    ');
save(output_path, ...
     'clab', 'sa', 'n_sensors', 'n_sources', ...
     'all_chanlocs', 'L_free', 'L_fixed', ...
     'A_eloreta_free', 'A_eloreta_fixed');
fprintf('Done\n');