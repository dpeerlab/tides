function [x_plot, y_plot] = make_trend_plots(x, y, varargin)

set_normalize = 0;
smoothing_factor = 0.5;
if ~isempty(varargin)
    for j = 1:length(varargin)
        if strcmp(varargin{j}, 'normalize')
            set_normalize = 1;
        end
        
        if strcmp(varargin{j}, 'smoothing_factor')
            smoothing_factor = varargin{j+1};
        end
    end
end
%smoothing_factor


num_bins = 256;
weights = zeros(num_bins, length(x));

for i=1:num_bins
	weights(i, :) = compute_weights(x, (i/num_bins)*range(x)+min(x), smoothing_factor);
end


x_plot = linspace(min(x), max(x), num_bins);
y_plot = weights*y./repmat(sum(weights, 2), 1, size(y, 2));


if set_normalize == 1
    disp('doing')
    y_plot = normalize_unit(y_plot);
end
end

function weights = compute_weights(points, loc, smoothing_factor)
    
    %range = quantile(points, .98) - quantile(points, .02);
    %min_std_dev = factor*.18*range; % minimum std_dev for dense regions
    %max_std_dev = .19*range; % max std_dev for sparse regions  
    %linear_slope = 10/range;
    
    
    min_std_dev = std(points) * 1.34 * length(points)^(-1/5) * smoothing_factor;
    
    weights = ((2*pi*(min_std_dev)^2)^(-0.5))*exp(-.5*((points - loc)/min_std_dev).^2);
    
% 	weights = mynormalize(weights(:), 100);
end

function W_new = normalize_unit(W)
    W_new = (W - min(W))/(max(W) - min(W));
end