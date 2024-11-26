function [num_voxels_roi] = count_significant(data, sa)
%COUNT_SIGNIFICANT Summary of this function goes here
%   Detailed explanation goes here
    n_rois_HO = numel(sa.HO_labels) - 1; % ignore Subcortical
    num_voxels_roi = zeros(n_rois_HO, 1);
    significance_mask = data > 0;
    for roi_ind = 1:n_rois_HO
        num_voxels_roi(roi_ind) = numel(get_voxels_roi(sa, roi_ind, significance_mask));
    end
end

