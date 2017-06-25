clc
clear
close all

% put the folder in your MATLAB path
loc = pwd;
addpath(genpath(loc));

% load the file
file = cytof_data_min('filename.fcs');

% three markers of interest:
marker1 = 'marker1';
marker2 = 'marker2';
marker3 = 'marker3';

% get the N x 3 matrix of the three markers
data = file.get_data({marker1, marker2, marker3});

% parameters:
noise_threshold = 0.9;
number_of_tides_steps = 257; % to get TIDES at 256 places, please provide the value 257

minx = min(data(:, 1));
miny = min(data(:, 2));  
minz = min(data(:, 3));

maxx = max(data(:, 1));
maxy = max(data(:, 2));  
maxz = max(data(:, 3));  

% If slices of DREVI needed, it is a vector of locations where DREVI is
% needed. The following will produce a DREVI slice at x = 0.1 only.
get_drevi_slices_at = [0.1];

% get threeD-DREVI:
% Can turn off some of the inputs:
% display_plots will display the plots
% smooth will smooth the surfaces
% slice plots the DREVI slice   

file.threeD_edge_visualize(marker1, marker2, marker3, get_drevi_slices_at, 'display_plots', 'smooth', 'slice', ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);
    
% get TIDES. The TIDES curve may need to be smoothed.
[tides_scores, trajectory_points] = file.compute_windowed_DREMI_interpolate(marker1, marker2, marker3, number_of_tides_steps, noise_threshold, ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);

% get threeD-dremi
dremi_score = file.compute_threeD_dremi(marker1, marker2, marker3, noise_threshold, ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);               
