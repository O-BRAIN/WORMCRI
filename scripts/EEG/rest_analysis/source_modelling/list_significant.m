function [result] = list_significant(num_voxels_roi, sa)
%LIST_SIGNIFICANT Summary of this function goes here
%   Detailed explanation goes here
    n_rois_HO = numel(sa.HO_labels) - 1; % ignore Subcortical
    result = "";
    
    [~, sort_idx] = sort(num_voxels_roi, 'descend');
    for i = 1:n_rois_HO
        if num_voxels_roi(sort_idx(i)) == 0
            break
        end
    
        newStr = sprintf('%d\t%s\n', ...
            num_voxels_roi(sort_idx(i)), sa.HO_labels{sort_idx(i)});
        result = result + newStr;
    end

    if strlength(result) == 0
        result = "No significant sources were found.";
    end
end

