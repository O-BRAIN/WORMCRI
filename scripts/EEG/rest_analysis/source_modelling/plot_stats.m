function [h] = plot_stats(data, sa, v, lims, cm17a)
%PLOT_STATS Summary of this function goes here
%   Detailed explanation goes here

    % Open the figure
    h = figure('Position', [10 10 1000 600]);
    
    % Setup the layout
    views = [
        v.leftlateral, v.dorsal, v.rightlateral; ...
        v.rightmedial, v.colorbar, v.leftmedial;
    ];
    
    allplots_cortex_subplots(sa, cort5K2full(data, sa), ...
        lims, cm17a, '', 1, 'views', views, 'cb_location', 'south');

end

