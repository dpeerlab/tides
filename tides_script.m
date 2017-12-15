clc
clear
close all

% put the folder in your MATLAB path
loc = pwd;
addpath(genpath(loc));

% load the fcs file
file = cytof_data_min('EM141482.fcs');

% arch-sinh transform the data with a co-factor of 5:
file = file.transform_data(5);

% three markers of interest:
marker1 = 'il-7ra';
marker2 = 'perforin';
marker3 = 't-bet';

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
get_drevi_slices_at = linspace(minx, maxx, 40);

% get threeD-DREVI:
% Can turn off some of the inputs:
% display_plots will display the plots
% smooth will smooth the surfaces
% slice plots the DREVI slice   
figure;
file.threeD_edge_visualize(marker1, marker2, marker3, get_drevi_slices_at, 'display_plots', 'smooth', 'slice', ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);
xlabel(marker1)
ylabel(marker2)
zlabel(marker3)
    
% If you want to make a video of the slices, please use the following code:
clearvars M;
name_of_video = 'example.avi';
v = VideoWriter(name_of_video);
v.FrameRate = 3;
N=length(get_drevi_slices_at); % Number of frames
j = linspecer(256);
for i = 1:N
    %figure();
    file.threeD_edge_visualize(marker1, marker2, marker3, get_drevi_slices_at(i), j, 'smooth', 'slice');       
    xlabel(marker1)
    ylabel(marker2)
    zlabel(marker3)
    M(i)=getframe(gcf); % leaving gcf out crops the frame in the movie.
end
%movie2avi(M, name_of_video, 'fps', 4.5);
open(v)
writeVideo(v, M);
close(v);
    
    
% get TIDES. The TIDES curve may need to be smoothed.
[tides_scores, trajectory_points] = file.compute_windowed_DREMI_interpolate(marker1, marker2, marker3, number_of_tides_steps, noise_threshold, ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);
figure; 
plot(trajectory_points, tides_scores, 'LineWidth', 2)
axis tight
box off
xlabel(marker1)
ylabel('DREMI')
title([marker2, ' -> ', marker3])
set(gca, 'FontSize', 14)

% if you want to smooth the TIDES curve:
[X, Y] = make_trend_plots(trajectory_points(:), tides_scores(:), 'smoothing_factor', 0.5);
figure; 
plot(X, Y, 'LineWidth', 2)
axis tight
box off
xlabel(marker1)
title([marker2, ' -> ', marker3])
ylabel('DREMI')
set(gca, 'FontSize', 14)
    
% get threeD-dremi
dremi_score = file.compute_threeD_dremi(marker1, marker2, marker3, noise_threshold, ...
        'minx', minx, 'miny', miny, 'minz', minz, 'maxx', maxx, 'maxy', maxy, 'maxz', maxz);
    
disp(['The 3D-DREMI score is: ' num2str(dremi_score, 3)])
