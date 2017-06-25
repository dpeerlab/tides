clc
clear
close all

loc = pwd;
addpath(genpath(loc));

file = cytof_data_min('tgfb_48h_1_wanderlust.fcs');
data = file.get_data({'wanderlust', 'ps6', 'pgsk3b'});

marker1 = 'wanderlust';
marker2 = 'ps6';
marker3 = 'pgsk3b';

% parameters:
noise_threshold = 0.9;
number_of_tides_steps = 257; % to get TIDES at 256 places, please provide the value 257

minx = min(data(:, 1));
miny = min(data(:, 2));
minz = min(data(:, 3));

maxx = max(data(:, 1));
maxy = max(data(:, 2));
maxz = max(data(:, 3));

get_drevi_slices_at = [0.1, 0.9];

figure;
% get the threeD-DREVI:
file.threeD_edge_visualize(marker1, marker2, marker3, get_drevi_slices_at, 'display_plots', 'smooth', 'slice', ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);
