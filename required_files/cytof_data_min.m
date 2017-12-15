%How do I extend this to time-series data?
%Need help extending this to time series. I need to make a struct of data
%qualities of time series data is that there arent equally as many cells
%in each time point so the data unless its chopped can't be in one array
%how do I initialize the marker channels?

classdef cytof_data_min
    
    properties
        
        %basic data
        data
        
        data_backup
        
        discretized_data
        channel_levels
        
        
        
        
        header
        %doublets or debris
        %event information
        debris
        wells
        event_nums
        invalid
        well_2_names
        names_2_well
        beads
        event_gate
        name_channel_map
        channel_name_map
        neg_control_well
        apoptosis_well
        pos_control_well
        empty_well
        
        DREMI_adjacency
        DREMI_sig_adjacency
        DREMI_sig_edges
        DREMI_channels
        activity_matrix
        
        
        dist_graph
        dist_clusters
        dist_hierarchy
        
        ord_clusters
        time_pred_clusters
        
        barcodes
        barcode_discretized
        bc_clusters
        bc_nd_gmm_means
        bc_nd_gmm_covars
        bc_nlogl
        bc_post
        bc_hierarchy
        bc_mahalad
        bead_nlogl
        bead_post
        
        gmfit
        bead_gmfit
        dna_gmfit
        %should I remove the doublets? yes probably.
        
        
        plates
        well_key %key that maps from barcodes to wells
        num_time_points
        cell_length_col
        
        %stims
        %drugs
        %doses
        time_points
        
        %channel information
        markers
        channels
        proteins
        surface_channels
        internal_channels
        barcoding_channels
        marker_channels
        dna_channels
        bead_channels
        eventnum_channel
        time_col
        
        %computed stats per channel
        means
        corrs
        vars
        skews
        kurtosis
        peaks
        maxes
        mins
        %peak_densities
        
        well_means
        well_corrs
        well_vars
        well_mins
        well_maxes
        well_distances
        well_knockdowns
        well_knockups
        well_L2s
        well_counts
        well_spearman
        well_skews
        well_kurtosis
        well_peaks
        
        %well_peak_densities
        
        well_KS_pvals
        well_KS_stats
        well_KS_qvals
        
        well_densities
        channel_densities
        channel_density_xvals
        
        
        %individual types of distances
        well_mean_distances
        well_corr_distances
        well_var_distances
        %well_kurt_distances
        %well_skew_distances
        %well_peak_difference
        %well_peak_density_distances
        
        %noise and measurement effects
        channel_decays
        mean_decay
        total_decay
        bc_gmm_means
        bc_gmm_sigmas
        debris_gmm_means
        debris_gmm_sigmas
        
        
        interpolated_bead_X;
        interpolated_bead_Y;
        
        
        %stuff that used to be in bnorm
        bead_data
        acquisition_times
        acquisition_rate
        acquisition_rate_index
        
        channel_averages
        channel_intercepts
        %if I use the channel averages then I have to compute them.
        %I should do it on the fly
        
        data_recon
        fluctuation_vector
        dna_low
        
        start_times
        end_times
        
        initial_y
        initial_x
        
        shift_x
        shift_y
        
        
        smooth_decays
        
        
        smooth_corrs
        recon_corrs
        
        smooth_index
        smooth_data
        smooth_diffs
        smooth_data_recon
        smooth_fluctuation_vector
        
        
        
        %tsne
        tsne_mapped_data
        sampled_events
        tsne_gated_indices
        
        
        %louvain
        COMTY
        edge_distances
        
        %wanderlust
        wanderlust_channel
        trajectories
        trajectory_times
        
        %ctbn
        ctbn_parent_map
        
        
        
    end
    
    methods
        %input output
        function obj = cytof_data_min( filename, name2, well_key)
            
            
            
            if(strcmp(filename,'empty'))
                
                sprintf('no file specified: %s', filename)
                %have to do everything else manually
                return;
            end
            [obj.data, obj.header, ~] = fca_readfcs(filename); %using arg1
            use_name = 0;
            if (nargin > 1)
                
                if(strcmp(name2,'name'))
                    use_name = 1;
                end
                
            end
            
            
            if(use_name ==0)
                obj = obj.compute_name_channel_maps();
            else
                obj = obj.compute_name_channel_maps('name');
            end
            
            
            [~,c] = size(obj.data);
            obj.marker_channels = linspace(1,c,c);
            
            
            obj.time_col = obj.find_channel_prefix('time');
            if(length(obj.time_col)>0)
                
                obj.marker_channels(find(obj.marker_channels == obj.time_col)) = [];
            end
            
            
            obj.cell_length_col = obj.find_channel_prefix('cell_length');
            if(length(obj.cell_length_col)>0)
                
                obj.marker_channels(find(obj.marker_channels == obj.cell_length_col)) = [];
            end
            
            
            obj.dna_channels = obj.find_channel_prefix('dna');
            [~, num_dna]=size(obj.dna_channels);
            for j=1:num_dna
                
                obj.marker_channels(find(obj.marker_channels == obj.dna_channels(j))) = [];
                
            end
            
            obj.barcoding_channels = obj.find_channel_prefix('bc');
            [~, num_barcodes]=size(obj.barcoding_channels);
            
            for j=1:num_barcodes
                
                obj.marker_channels(find(obj.marker_channels == obj.barcoding_channels(j))) = [];
                
            end
            
            obj.bead_channels = obj.find_channel_prefix('bead');
            [~, num_bead]=size(obj.bead_channels);
            for j=1:num_bead
                
                obj.marker_channels(find(obj.marker_channels == obj.bead_channels(j))) = [];
                
            end
            
            trash_channel = obj.find_channel_prefix('barcode');
            [~, num_trash]=size(trash_channel);
            for j=1:num_trash
                
                obj.marker_channels(find(obj.marker_channels == trash_channel(j))) = [];
                
            end
            
            %using arg4
            
            %specifyed barcoding channels (in ORDER)
            if (nargin > 2)
                
                obj.well_key = well_key;
                
            end
            
            %             %event number column is already computed and specified
            %             if (nargin > 2)
            %
            %                obj.event_no_col = event_no_col;
            %                obj.marker_channels(find(obj.marker_channels == obj.event_no_col)) = [];
            %
            %
            %             end
            %
            obj.time_points = ones(size(obj.data,1),1);
        end
        
        function channel_indices =  find_channel_prefix(obj, prefix_string)
            
            num_channels = size(obj.data,2);
            channel_indices = [];
            for i=1:num_channels
                channel_name = obj.channel_name_map(i);
                index = strfind(channel_name,prefix_string);
                if(isempty(index{1}))
                    continue;
                end
                
                if(index{1}(1)~=1)
                    continue;
                end
                channel_indices = [channel_indices i] ;
                
            end
            
        end
        
        function write_data(obj, filename, varargin)
            
            %have to figure out what marker_names and channel names are
            optargin = size(varargin, 2);
            [num_events, num_markers] = size(obj.data);
            marker_names = cell(num_markers,1);
            channel_names = cell(num_markers,1);
            
            for i=1:num_markers
                
                marker_names(i) = obj.channel_name_map(i);
                channel_names(i) = obj.channel_name_map(i);
                
            end
            
            if(optargin >0)
                
                data = obj.get_well_data(varargin{1});
                
            else
                
                data = obj.data;
            end
            
            fca_writefcs(filename, data, marker_names, channel_names);
            
            
            
        end
        
        function obj = recompute_name_channel_map(obj)
            
            obj.name_channel_map = containers.Map();
            for i=1:length(obj.channel_name_map)
                obj.name_channel_map(obj.channel_name_map{i}) = i;
            end
            
        end
        
        
        function obj = compute_name_channel_maps(obj, varargin)
            
            %s = sprintf('in name channel map! \n');
            
            use_name = 0;
            optargin = length(varargin);
            for i=1:optargin
                
                if(strcmp('name',varargin{i}))
                    use_name = 1
                end
            end
            
            header_size = length(obj.header.par);
            obj.name_channel_map = containers.Map();
            obj.channel_name_map = cell(header_size,1);
            for i=1:header_size
                
                hname = obj.header.par(i).name2;
                if(use_name)
                    hname = obj.header.par(i).name;
                end
                
                if(length(strfind(hname, '('))==0)
                    
                    
                    %stuff like the DNA channels do not have the isotope
                    %names
                    
                    if(length(hname)==0)
                        hname = sprintf('channel%d',i);
                        %if its still null for some reason just use the
                        %channel number
                    end
                    hname = lower(hname);
                    obj.name_channel_map(hname) = i;
                    obj.channel_name_map{i} = hname;
                    
                else
                    n=strfind(hname, '(');
                    
                    if(n==1)
                        %using isotope name
                        
                        n_end = strfind(hname, ')');
                        real_name = hname(n+1:n_end-1);
                        
                    else
                        
                        real_name = hname(1:n-1);
                        
                    end
                    
                    
                    if(length(real_name)==0)
                        real_name = sprintf('channel%d',i);
                        %if its still null for some reason just use the
                        %channel number
                    end
                    real_name = lower(real_name);
                    obj.name_channel_map(real_name) = i;
                    obj.channel_name_map{i} = real_name;
                end
                
            end
            
            
        end
        
        function obj = compute_name_channel_map_from_array(obj, varargin)
            
            
            
            header_size = length(obj.header);
            obj.name_channel_map = containers.Map();
            obj.channel_name_map = cell(header_size,1);
            for i=1:header_size
                
                hname = obj.header{i};
                
                
                if(length(strfind(hname, '('))==0)
                    
                    
                    %stuff like the DNA channels do not have the isotope
                    %names
                    
                    if(length(hname)==0)
                        hname = sprintf('channel%d',i);
                        %if its still null for some reason just use the
                        %channel number
                    end
                    hname = lower(hname);
                    obj.name_channel_map(hname) = i;
                    obj.channel_name_map{i} = hname;
                    
                else
                    n=strfind(hname, '(');
                    
                    if(n==1)
                        %using isotope name
                        
                        n_end = strfind(hname, ')');
                        real_name = hname(n+1:n_end-1);
                        
                    else
                        
                        real_name = hname(1:n-1);
                        
                    end
                    
                    
                    if(length(real_name)==0)
                        real_name = sprintf('channel%d',i);
                        %if its still null for some reason just use the
                        %channel number
                    end
                    real_name = lower(real_name);
                    obj.name_channel_map(real_name) = i;
                    obj.channel_name_map{i} = real_name;
                end
                
            end
            
            
        end
        
        
        function [name] = get_name_from_channel(obj, channel)
            
            name = obj.channel_name_map(channel);
            
        end
        
        function [channel] = get_channel_from_name(obj, name)
            
            name = lower(name);
            channel = obj.name_channel_map(name);
        end
        
        function obj = add_data_matrix_header(obj, data, header, wells)
            
            obj.data = data;
            obj.header = header;
            %obj = obj.compute_name_channel_maps();
            obj = obj.compute_name_channel_map_from_array();
            
            [~, c] = size(obj.data);
            for i=1:c
                
                if i==obj.cell_length_col
                    continue;
                end
                
                if i==obj.time_col
                    continue;
                end
                
                obj.marker_channels = [obj.marker_channels i];
            end
            
            [~, num_dna]=size(obj.dna_channels);
            for j=1:num_dna
                
                obj.marker_channels(find(obj.marker_channels == obj.dna_channels(j))) = [];
                
            end
            
            [~, num_barcodes]=size(obj.barcoding_channels);
            
            for j=1:num_barcodes
                
                obj.marker_channels(find(obj.marker_channels == obj.barcoding_channels(j))) = [];
                
            end
            
            
            %obj.marker_channels(find(obj.marker_channels == obj.event_no_col)) = [];
            obj.wells = wells;
            %debarcoding info so that I do not have to repeat.
            
        end
        
        
        function obj = add_dna_channels(obj,dna_channels)
            
            %specifying the DNA channels
            obj.dna_channels = dna_channels;
            [~, num_dna]=size(dna_channels);
            for j=1:num_dna
                
                obj.marker_channels(find(obj.marker_channels == dna_channels(j))) = [];
                
            end
            
        end
        
        
        function obj = add_bead_channels(obj,bead_channels)
            
            obj.bead_channels = bead_channels;
            [~, num_beads]=size(bead_channels);
            
            %                 for j=1:num_beads
            %
            %                     obj.marker_channels(find(obj.marker_channels == bead_channels(j))) = [];
            %
            %                 end
            
        end
        
        %specify eventnum channel
        
        function obj = add_eventnum_channel(obj, eventnum_channel)
            
            obj.eventnum_channel = eventnum_channel;
            obj.marker_channels(find(obj.marker_channels == eventnum_channel)) = [];
            
        end
        
        
        
        function obj = add_data(obj, filename, varargin)
            
            %assuming same header for this data otherwise use a higher
            %level construct
            %ayyo I guess the panels are not the same for these
            %oh well i guess I'll have to make sure I know what the panels
            %are at some point
            filename
            
            [new_data, ~ , ~] = fca_readfcs(filename);
            
            
            
            optargin = size(varargin, 2);
            
            if(optargin == 1)
                
                if(strcmp(varargin{1}, 'fix_times'))
                    
                    [old_data_size, ~] = size(obj.data)
                    final_acquisition_time = obj.data(old_data_size,1)
                    
                    [new_data_size, ~] = size(new_data)
                    
                    for i=1:new_data_size
                        
                        new_data(i,1) = new_data(i,1)+final_acquisition_time;
                        
                    end
                    
                    obj.data = [obj.data; new_data];
                end
                
                if(strcmp(varargin{1}, 'add timepoints'))
                    
                    
                    
                    obj.data = [obj.data; new_data];
                    
                    time_point_end = length(obj.time_points);
                    new_time_point = obj.time_points(time_point_end)+1;
                    time_points = ones(size(new_data,1),1);
                    time_points = time_points .* new_time_point;
                    obj.time_points = [obj.time_points; time_points];
                    
                    
                end
            else
                
                obj.data = [obj.data; new_data];
                
            end
            
            
            
        end
        
        %figure out barcode for each event
        function obj = compute_wells(obj)
            
            % obj = obj.compute_barcode_gmm();
            [r,~] = size(obj.data);
            obj.wells = zeros(r,1);
            % barcodes = zeros(100,1);
            obj.barcode_discretized = zeros(r, 7);
            
            for i=1:r
                
                [barcode, ~] = obj.compute_event_barcode(i);
                obj.barcodes(i) = barcode;
                %obj.barcode_discretized(i,:) = barcode_array;
                
                %barcodes(i) = barcode;
                if(barcode == -1 || barcode==0)
                    obj.wells(i) = 0;
                    obj.barcodes(i) = 0;
                else
                    
                    obj.wells(i) = obj.well_key(barcode);
                    %already has the +1
                end
                
                
            end
            
            %obj.wells = wells(i);
            %lets see if we return the proper wells
            %delete the events that cant be barcoded properly
            %optional if this would be easier
            %obj.wells
            
        end
        function plot_scatter(obj, channel1_name, channel2_name, varargin)
            c_by = 0;
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{i}, 'color')
                        c_by = 1;
                        color_by = varargin{i+1};
                    end
                end
            end
            channel1 = obj.data(:, obj.name_channel_map(channel1_name));
            channel2 = obj.data(:, obj.name_channel_map(channel2_name));
            if c_by == 0;
                %colormap(linspecer(256))
                plot(channel1, channel2, '.', 'MarkerSize', 10)
            else
                plot(channel1, channel2, '.', 'MarkerSize', 10, 'col', color_by)
            end
            xlabel(channel1_name, 'FontSize', 20)
            ylabel(channel2_name, 'FontSize', 20)
%             xlim([0, Inf])
%             ylim([0, Inf])
        end
        
        function [barcode, barcode_array] = compute_event_barcode(obj,eventnum)
            
            [~,num_barcodes] = size(obj.barcoding_channels);
            barcode = 0;
            
            barcode_array = zeros(1,7);
            
            for i=1:num_barcodes
                
                value = obj.barcode_value(i, eventnum);
                barcode_array(i) = value;
                
                if value == -1
                    %need to throw away event
                    barcode = -1;
                    break
                end
                if value == 1
                    barcode = barcode + 2^(num_barcodes-i);
                    %barcodes go the other way according to the
                    %script I got from Rachel.
                end
            end
            if(barcode >= 0)
                barcode = barcode + 1;
                %the range of the keys is 1-2^7 and not 0-(2^7-1)
            end
            sprintf('event %d barcode %d ', eventnum, barcode);
        end
        
        function bc_val = barcode_value(obj,bc_index,eventnum)
            
            bc_channel = obj.barcoding_channels(bc_index);
            
            prob0 = normpdf(obj.data(eventnum,bc_channel), obj.bc_gmm_means(bc_index,1), obj.bc_gmm_sigmas(bc_index,1));
            prob1 = normpdf(obj.data(eventnum,bc_channel), obj.bc_gmm_means(bc_index,2), obj.bc_gmm_sigmas(bc_index,2));
            
            if(prob0 > prob1)
                bc_val = 0;
            elseif (prob1 > prob0)
                bc_val = 1;
            else
                bc_val = -1;
            end
            
        end
        
        function dna_val = debris_cluster_value(obj,dna_index,eventnum)
            
            dna_channel = obj.dna_channels(dna_index);
            
            prob0 = normpdf(obj.data(eventnum,dna_channel), obj.debris_gmm_means(dna_index,1), obj.debris_gmm_sigmas(dna_index,1));
            prob1 = normpdf(obj.data(eventnum,dna_channel), obj.debris_gmm_means(dna_index,2), obj.debris_gmm_sigmas(dna_index,2));
            
            if(prob0 > prob1)
                dna_val = 0;
            elseif (prob1 > prob0)
                dna_val = 1;
            else
                dna_val = -1;
            end
            
        end
        
        
        function obj = compute_barcode_gmm(obj)
            
            clear X;
            
            [~,num_barcodes] = size(obj.barcoding_channels);
            obj.bc_gmm_means = zeros(num_barcodes,2);
            
            for i=1:num_barcodes
                
                barcode_col_data = obj.data(:,obj.barcoding_channels(i));
                [~, M, V, ~] = EM_GM_fast(barcode_col_data,2);
                if (M(1) < M(2))
                    obj.bc_gmm_means(i,1) = M(1);
                    obj.bc_gmm_means(i,2) = M(2);
                    obj.bc_gmm_sigmas(i,1) = V(1,1,1);
                    obj.bc_gmm_sigmas(i,2) = V(1,1,2);
                else
                    obj.bc_gmm_means(i,1) = M(2);
                    obj.bc_gmm_means(i,2) = M(1);
                    obj.bc_gmm_sigmas(i,1) = V(1,1,2);
                    obj.bc_gmm_sigmas(i,2) = V(1,1,1);
                    
                end
                
                
            end
            
            
            
            
        end
        
        function [mi_matrix ] = channel_pair_mi_compute(obj, xchannel, channel_names, noise_threshold, varargin)
            
            ctrl_specified = -1;
            ctrl_data = [];
            
            mi_matrix = zeros(1, length(channel_names));
            
            
            
            flip_channel = 0;
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'ctrl_data'))
                    
                    ctrl_specified = 1;
                    ctrl_data = varargin{i+1};
                    
                    
                end
                
                if(strcmp(varargin{i}, 'y'))
                    
                    
                    flip_channel = 1;
                    
                    
                end
                
                
            end
            
            
            
            
            for j=1:length(channel_names)
                
                
                if(flip_channel==1)
                    channel2_name = xchannel
                    channel1_name = channel_names{j}
                    'flipped channel'
                    %flipped it.
                else
                    channel1_name = xchannel
                    channel2_name = channel_names{j}
                    
                end
                
                
                if(ctrl_specified<0)
                    [mi_matrix(j),~] = obj.compute_dremi(channel1_name, channel2_name, noise_threshold);
                else
                    [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(obj, channel1_name, channel2_name, 25, 255);
                    [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(ctrl_data, channel1_name, channel2_name, 25, 255);
                    maxy = max(maxy1,maxy2);
                    [mi_matrix(j),~] = obj.compute_twoD_dremi(channel1_name, channel2_name, noise_threshold, 'maxy', maxy);
                end
                
                
                
                
                
            end
            
            
            
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'plot'))
                    
                    
                    CLIM = [0 1];
                    imagesc(mi_matrix);
                    set(gca,'ytick',1);
                    set(gca,'yticklabel',{xchannel});
                    xticklabel_rotate([1:length(channel_names)],45,channel_names);
                    colorbar
                    
                end
                
            end
            
            
        end
        
        
        function obj = compute_debris_gmm(obj)
            clear X;
            
            [~,num_DNA] = size(obj.dna_channels);
            obj.debris_gmm_means = zeros(num_DNA,2);
            
            for i=1:num_DNA
                
                DNA_col_data = obj.data(:,obj.dna_channels(i));
                [~, M, V, ~] = EM_GM_fast(DNA_col_data,2);
                
                if (M(1) < M(2))
                    obj.debris_gmm_means(i,1) = M(1);
                    obj.debris_gmm_means(i,2) = M(2);
                    obj.debris_gmm_sigmas(i,1) = V(1,1,1);
                    obj.debris_gmm_sigmas(i,2) = V(1,1,2);
                else
                    obj.debris_gmm_means(i,1) = M(2);
                    obj.debris_gmm_means(i,2) = M(1);
                    obj.debris_gmm_sigmas(i,1) = V(1,1,2);
                    obj.debris_gmm_sigmas(i,2) = V(1,1,1);
                    
                end
                
                
            end
            
            
        end
        
        
        
        function obj = compute_stats(obj)
            
            obj = obj.compute_means();
            obj = obj.compute_corrcoefs();
            obj = obj.compute_vars();
            obj = obj.compute_distro_stats();
            obj = obj.compute_maxes_mins();
        end
        
        
        
        function obj = compute_well_stats(obj, num_wells)
            
            %how many wells do we have, 96
            means = zeros(num_wells,length(obj.marker_channels));
            vars = zeros(num_wells, length(obj.marker_channels));
            corrs = zeros(num_wells,length(obj.marker_channels), length(obj.marker_channels));
            %also including dna channels and cell length
            mins = zeros(num_wells, length(obj.marker_channels));
            maxes = zeros(num_wells, length(obj.marker_channels));
            well_skews = zeros(num_wells, length(obj.marker_channels));
            well_kurtosis = zeros(num_wells, length(obj.marker_channels));
            well_peaks = zeros(num_wells, length(obj.marker_channels));
            
            for i=1:num_wells
                well_no = i
                well_indices = find(obj.wells==i);
                if(length(well_indices)<10)
                    
                    continue;
                end
                means(i,:) = obj.compute_well_means(i);
                corrs(i,:,:) = obj.compute_well_corrcoefs(i);
                vars(i,:) = obj.compute_well_vars(i);
                vars(i,:) = obj.compute_well_vars(i);
                [mins(i,:), maxes(i,:)] = obj.compute_well_min_max(i);
                [well_skews(i,:), well_kurtosis(i,:), well_peaks(i,:)] = obj.compute_well_distro_stats(i);
                
            end
            
            obj.well_means = means;
            obj.well_vars = vars;
            obj.well_corrs = corrs;
            obj.well_mins = mins;
            obj.well_maxes = maxes;
            obj.well_skews = well_skews;
            obj.well_peaks = well_peaks;
            obj.well_kurtosis = well_kurtosis;
            
        end
        
        function obj = compute_well_distance_stats(obj)
            
            %how many wells do we have, 96
            mean_distances = zeros(96,length(obj.marker_channels));
            var_distances = zeros(96, length(obj.marker_channels));
            corr_distances = zeros(96,length(obj.marker_channels), length(obj.marker_channels));
            
            
            for i=1:96
                
                mean_distances(i,:) = obj.compute_well_mean_distances(i);
                corr_distances(i,:,:) = obj.compute_well_corrcoef_distances(i);
                var_distances(i,:) = obj.compute_well_var_distances(i);
                
            end
            
            obj.well_mean_distances = mean_distances;
            obj.well_var_distances = var_distances;
            obj.well_corr_distances = corr_distances;
            %need to add kurtosis and skew to this.
            
            
        end
        
        
        function obj = compute_means(obj)
            
            
            [~,num_markers] = size(obj.marker_channels);
            obj.means = zeros(1,num_markers);
            for i=1:num_markers
                c=obj.marker_channels(i);
                obj.means(i) = median(obj.data(:, c));
            end
        end
        
        function obj = compute_maxes_mins(obj)
            
            
            [~,num_markers] = size(obj.marker_channels);
            obj.mins = zeros(1,num_markers);
            obj.maxes = zeros(1,num_markers);
            for i=1:num_markers
                c=obj.marker_channels(i);
                obj.mins(i) = min(obj.data(:, c));
                obj.maxes(i) = max(obj.data(:,c));
                
            end
        end
        
        function obj = compute_distro_stats(obj)
            
            
            [~,num_markers] = size(obj.marker_channels);
            
            obj.skews = zeros(1, num_markers);
            obj.peaks = zeros(1, num_markers);
            obj.kurtosis = zeros(1, num_markers);
            for i=1:num_markers
                c=obj.marker_channels(i);
                obj.skews(i) = skewness(obj.data(:, c))
                obj.kurtosis(i) = kurtosis(obj.data(:,c));
                [f, xi] = ksdensity(obj.data(:,c));
                obj.peaks(i)= length(findpeaks(f));
            end
        end
        
        
        
        
        function channel_mean = get_channel_mean(obj, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            
            channel_mean = median(obj.data(:,channel));
            
        end
        
        function channel_var = get_channel_var(obj, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            
            channel_mean = var(obj.data(:,channel));
            
        end
        
        
        function means = compute_well_means(obj, well_no)
            
            
            %find the markers that belong to that well
            well_indices = find(obj.wells == well_no);
            
            [~,num_markers] = size(obj.marker_channels);
            means = zeros(1, num_markers);
            
            for i=1:num_markers
                
                c=obj.marker_channels(i);
                data_subset = obj.data(:,c);
                %means(i) = mean(data_subset(well_indices));
                means(i) = median(data_subset(well_indices));
            end
            
        end
        
        
        
        function [well_skews, well_kurtosis, well_peaks] = compute_well_distro_stats(obj, well_no)
            
            well_data = obj.get_well_data(well_no);
            [~,num_markers] = size(obj.marker_channels);
            well_skews = zeros(1, num_markers);
            well_peaks = zeros(1, num_markers);
            well_kurtosis = zeros(1, num_markers);
            
            for i=1:num_markers
                
                c = obj.marker_channels(i);
                well_skews(i) = skewness(well_data(:, c));
                well_kurtosis(i) = kurtosis(well_data(:,c));
                [f, xi] = ksdensity(well_data(:,c));
                well_peaks(i) = length(findpeaks(f));
            end
        end
        
        
        
        function mean_distances = compute_well_mean_distances(obj, well_no)
            
            [~,num_markers] = size(obj.marker_channels);
            mean_distances = zeros(1, num_markers);
            
            for i=1:num_markers
                
                
                mean_distances(i) = obj.well_means(well_no, i)-obj.means(i);
            end
            
            
        end
        
        
        
        
        function [mins, maxes] = compute_well_min_max(obj, well_no)
            
            
            %find the markers that belong to that well
            well_indices = find(obj.wells == well_no);
            
            [~,num_markers] = size(obj.marker_channels);
            mins = zeros(1, num_markers);
            maxes = zeros(1, num_markers);
            
            for i=1:num_markers
                
                c=obj.marker_channels(i);
                data_subset = obj.data(:,c);
                mins(i) = min(data_subset(well_indices));
                maxes(i) = max(data_subset(well_indices));
            end
            
        end
        
        
        
        function scatter_and_color(obj, chn1, chn2, chn3)
            x = obj.data(:, obj.name_channel_map(chn1));
            y = obj.data(:, obj.name_channel_map(chn2));
            col_by = obj.data(:, obj.name_channel_map(chn3));
            colormap('jet');
            scatter(x, y, 40, col_by, 'filled');
            colorbar()
            xlabel(chn1);
            ylabel(chn2);
            set(gca, 'FontSize', 15);
        end
        
        function [obj, corrs, range] = compute_corrcoefs(obj)
            
            [~,num_markers] = size(obj.marker_channels);
            range = zeros(num_markers,num_markers);
            
            for i=1:num_markers
                for j = 1: num_markers
                    
                    x = obj.marker_channels(i);
                    y = obj.marker_channels(j);
                    R = corrcoef(obj.data(:, x), obj.data(:, y));
                    obj.corrs(i,j) = R(1,2);
                    range1 = max(obj.data(:,x))-min(obj.data(:,x));
                    range2 = max(obj.data(:,y))-min(obj.data(:,y));
                    range(i,j) = (range1+range2)/2;
                end
            end
            
            corrs = obj.corrs;
            
            
        end
        
        function corrs = compute_well_corrcoefs(obj, well_no)
            
            
            %find the markers that belong to that well
            well_indices = find(obj.wells == well_no);
            
            
            [~,num_markers] = size(obj.marker_channels);
            corrs = zeros(num_markers,num_markers);
            %I want correlations computed with markers, DNA and
            %cell length
            
            for i=1:num_markers
                for j=1:num_markers
                    if(i==j)
                        continue;
                    end
                    x = obj.marker_channels(i);
                    y = obj.marker_channels(j);
                    
                    data_subset_x = obj.data(:, x);
                    data_subset_y = obj.data(:, y);
                    
                    R = corrcoef(data_subset_x(well_indices), data_subset_y(well_indices));
                    corrs(i,j) = R(1,2);
                end
            end
            
        end
        
        function corr_distances = compute_well_corrcoef_distances(obj, well_no)
            
            
            [~,num_markers] = size(obj.marker_channels);
            corrs_distances = zeros(num_markers,num_markers);
            %I want correlations computed with markers, DNA and
            %cell length
            
            for i=1:num_markers
                for j=1:num_markers
                    if(i==j)
                        continue;
                    end
                    
                    corr_distances(i,j) = obj.well_corrs(well_no, i, j)-obj.corrs(i,j);
                end
            end
            
        end
        
        function obj = compute_vars(obj)
            
            [~,num_markers] = size(obj.marker_channels);
            for i=1:num_markers
                c=obj.marker_channels(i);
                x = obj.data(:,c);
                obj.vars(i) = var(x);
                
            end
            
        end
        
        function vars = compute_well_vars(obj, well_no)
            
            
            %find the markers that belong to that well
            well_indices = find(obj.wells == well_no);
            
            
            [~,num_markers] = size(obj.marker_channels);
            vars = zeros(1,num_markers);
            
            for i=1:num_markers
                
                c=obj.marker_channels(i);
                data_subset = obj.data(:,c);
                vars(i) = var(data_subset(well_indices));
            end
            
        end
        
        function var_distances = compute_well_var_distances(obj, well_no)
            
            
            %find the markers that belong to that well
            
            [~,num_markers] = size(obj.marker_channels);
            var_distances = zeros(1,num_markers);
            
            for i=1:num_markers
                var_distances(i) = obj.well_vars(well_no, i)-obj.vars(i);
            end
            
        end
        
        function obj = compute_all_KS_stats(obj)
            
            [~,num_markers] = size(obj.marker_channels);
            obj.well_KS_pvals = zeros(96, num_markers);
            obj.well_KS_stats = zeros(96, num_markers);
            obj.well_KS_qvals = zeros(96, num_markers);
            
            for i =1:96
                for j = 1:num_markers
                    
                    mc = obj.marker_channels(j);
                    [obj.well_KS_pvals(i,j), obj.well_KS_stats(i,j)] = obj.compute_KS_test(i,mc);
                    
                end
            end
            
            %computing the qvalues based on the pvalues
            for i=1:num_markers
                
                obj.well_KS_qvals(:,i) = NewFDR(obj.well_KS_pvals(:,i),0)
                
            end
            
            % qValues=NewFDR(obj._well_KS_pvals(i,j), 0);
            % compute the qvalues here after all the pvalues are computed
        end
        
        function [pval, stat] = compute_KS_test(obj, well_no, channel_name, varargin)
            
            
            
            well_data = obj.get_well_data(well_no, channelname);
            
            optargin = size(varargin,2);
            if (optargin>0)
                channel_data = varargin{1};
            else
                channel_data  = obj.get_nonedge_data(channelname);
            end
            
            %[h, KS_pval, KS_stat] = kstest2(well_data, channel_data);
            [pval,~, stat] = RANKSUM(well_data,channel_data)
            
        end
        
        %need a function to compare two signals.
        function [dist, dist_vector] = compute_distance(obj, obj2)
            
            [~,num_markers] = size(obj.marker_channels);
            dist_vector = zeros(num_markers);
            for i=1:num_markers
                
                dist_vector(i)=abs(obj.means(i)- obj2.means(i));
                
            end
            
            offset = num_markers;
            
            for i=1:num_markers
                
                for j=1:num_markers
                    
                    dist_vector(offset+j)=abs(obj.corrs(i,j)- obj2.corrs(i,j));
                    
                end
                
                offset = offset + num_markers;
                
            end
            
            dist=sum(dist_vector);
        end
        
        
        %what should the distance be?
        %distance between the means can be one dimension of the distance
        %L2 norm of the difference can be another dimension
        %what is the one that Dana wanted me to try?
        function [dist, dist_vector] = compute_distance_from_control(obj, well_no, varargin)
            
            [~,num_markers] = size(obj.marker_channels);
            dist_vector = zeros(1,num_markers+(num_markers^2)+num_markers);
            %trying z-score for now
            
            
            optargin = size(varargin,2);
            if(optargin==1)
                control_well = varargin{1};
            else
                control_well = 0;
            end
            
            
            well_means_zscore = zscore(obj.well_means);
            well_vars_zscore = zscore(obj.well_vars);
            
            for i=1:num_markers
                
                
                if(optargin==1)
                    
                    dist_vector(i) = well_means_zscore(well_no, i)- well_means_zscore(control_well, i);
                else
                    dist_vector(i) = well_means_zscore(well_no, i);
                    
                end
                
                
                
            end
            
            offset = num_markers;
            
            %I need the
            well_corrs_zscore = zscore(obj.well_corrs(well_no,:,:));
            
            if optargin==1
                control_corrs_zscore = zscore(obj.well_corrs(control_well,:, :));
            end
            
            
            for i=1:num_markers
                
                for j=1:num_markers
                    
                    if optargin==1
                        
                        dist_vector(offset+j)= well_corrs_zscore(1, i,j)- control_corrs_zscore(1, i, j);
                    else
                        dist_vector(offset+j)= well_corrs_zscore(1, i,j);
                        
                    end
                    
                end
                
                offset = offset + num_markers;
                
            end
            
            for i=1:num_markers
                
                if optargin==1
                    
                    dist_vector(offset+i) = well_vars_zscore(well_no, i) - well_vars_zscore(control_well, i);
                else
                    dist_vector(offset+i) = well_vars_zscore(well_no, i);
                    
                end
                
            end
            
            dist=norm(dist_vector);
            
        end
        
        
        
        
        function obj = transform_data(obj,c)
            
            obj.data_backup = obj.data;
            
            am_being_transformed = sprintf('I am being transformed! \n')
            
            
            %should I just transform all of this data?
            %well I don't understand this transform well
            %maybe I should just transform everything and then check
            
            
            [~,num_markers] = size(obj.marker_channels);
            for i=1:num_markers
                d=obj.marker_channels(i);
                
                x = obj.data(:,d);
                %obj.data(:,d) = log(1/c*x + sqrt((1/c*x).^2+1));
                obj.data(:,d) = asinh(x/c);
                %isn't this right?
                %I should be transforming the marker channels and putting
                %them back where they belong.
            end
            [~,bc_channels] = size(obj.barcoding_channels);
            for i=1:bc_channels
                d=obj.barcoding_channels(i);
                x = obj.data(:,d);
                %obj.data(:,d) = log(1/c*x + sqrt((1/c*x).^2+1));
                obj.data(:,d) = asinh(x/c);
                %isn't this right?
                %I should be transforming the marker channels and putting
                %them back where they belong.
            end
            
            [~,dn_channels] = size(obj.dna_channels);
            for i=1:dn_channels
                d=obj.dna_channels(i);
                x = obj.data(:,d);
                %obj.data(:,d) = log(1/c*x + sqrt((1/c*x).^2+1));
                obj.data(:,d) = asinh(x/c);
                %isn't this right?
                %I should be transforming the marker channels and putting
                %them back where they belong.
            end
            
        end
        
        function obj = add_computed_channel(obj, channel_name, channel_data)
            
            num_channels = size(obj.data,2);
            obj.data = [obj.data channel_data];
            obj.marker_channels = [obj.marker_channels num_channels+1];
            channel_num = num_channels+1;
            obj.name_channel_map(channel_name) = channel_num;
            obj.channel_name_map{channel_num} = channel_name;
            
        end
        
        function obj = untransform_data(obj)
            
            obj.data = obj.data_backup;
            obj.data_backup = [];
            
            
        end
        
        function obj = remove_low_dna(obj)
            
            obj = obj.compute_dna_gmm_clusters(2, obj.dna_channels);
            obj.debris = obj.dna_low;
            obj = obj.remove_debris();
            
        end
        
        
        
        function obj = clean_beads(obj)
            
            %basically take the beads you get from identify_beads and discard anything with anything with other channels active
            num_events = length(obj.beads);
            for j = 1:num_events
                
                if(obj.beads(j) == 0)
                    continue;
                end
                is_bead = 1;
                
                for i=1:length(obj.marker_channels)
                    
                    marker_channel = obj.marker_channels(i);
                    
                    is_bead_channel = find(obj.bead_channels==marker_channel);
                    if(length(is_bead_channel)>0)
                        continue;
                    end
                    
                    marker_channel_value = obj.data(j,marker_channel);
                    if(marker_channel_value>0)
                        is_bead = 0;
                        break;
                    end
                    
                    
                end
                
                %do it for dna channels too
                for i=1:length(obj.dna_channels)
                    
                    dna_channel = obj.dna_channels(i);
                    dna_channel_value = obj.data(j,dna_channel);
                    if(dna_channel_value>0)
                        is_bead = 0;
                        break;
                    end
                    
                    
                end
                
                obj.beads(j) = is_bead;
            end
            
            
            
        end
        
        
        
        
        function obj = identify_beads(obj)
            
            
            found_bead_cluster = 0;
            
            
            filter_channels = [obj.bead_channels obj.dna_channels];
            %filter_channels = union(filter_channels, obj.marker_channels)
            
            %dna channels and bead channels, the dna channels should be very
            %low or negative in fact, adding other channels also
            
            bead_cluster_index = 0;
            
            
            %doublets should be eliminated if I take only cells with
            %bead-bead
            sprintf(' in bead identify ')
            
            obj = obj.compute_dna_gmm_clusters(2, obj.dna_channels);
            obj.invalid = ~obj.dna_low;
            
            
            %obj = obj.compute_dna_gmm_clusters(2, 4);
            %obj.invalid = obj.invalid | obj.dna_low
            sprintf('done doing the dna clusters \n')
            
            
            %figure out an approximate percentage of beads automatically
            
            potential_bead_indices = find(obj.dna_low == 1);
            num_beads = 0;
            %obj.beads = obj.dna_low;
            
            size(potential_bead_indices)
            
            for i=1:length(potential_bead_indices)
                
                index = potential_bead_indices(i);
                is_bead = 1;
                for j=1:length(obj.bead_channels)
                    
                    if(obj.data(index,obj.bead_channels(j))<0)
                        
                        is_bead = 0;
                        break;
                    end
                end
                
                if(is_bead)
                    
                    num_beads = num_beads + 1;
                    %obj.beads(index) = 1;
                    
                end
            end
            [num_events, ~] = size(obj.data);
            
            approx_percentage = num_beads/ num_events
            approx_num_beads = num_beads
            
            
            for i=2:10
                
                %precluster the DNA channels and only do it on those.
                
                
                obj = obj.compute_bead_gmm_clusters(i, filter_channels);
                
                sprintf('done with bead gmm clusters for k= %d \n', i)
                
                if(length(obj.bead_gmfit.mu)==0)
                    %I think that this means it did not converge
                    continue;
                end
                
                %these were sorted this time
                %so now I have to find the bead channels and the crappy
                %nonbead channels
                bead_indices =zeros(1, length(obj.bead_channels));
                
                for j=1:length(obj.bead_channels)
                    
                    bead_indices(j) = find(filter_channels == obj.bead_channels(j))
                    
                end
                
                nonbead_channels = setdiff(filter_channels, obj.bead_channels)
                nonbead_indices=zeros(1, length(nonbead_channels));
                
                for j=1:length(nonbead_channels)
                    
                    nonbead_indices(j) = find(filter_channels == nonbead_channels(j))
                end
                
                mean_beadsum = sum(obj.bead_gmfit.mu(:,bead_indices),2);
                mean_nonbeadsum = sum(obj.bead_gmfit.mu(:,nonbead_indices),2);
                
                [~,max_bead_index] = max(mean_beadsum);
                
                
                [~,min_nonbead_index] = min(mean_nonbeadsum);
                
                
                if(max_bead_index ~= min_nonbead_index)
                    
                    %sprintf('WARNING max bead and min non beads are not on the same cluster. \n')
                    
                end
                
                bead_cluster_index = max_bead_index;
                %this one is better because the other stuff should mostly be
                %debris
                
                
                cluster_sizes = obj.compute_bc_cluster_sizes(obj.beads)
                [num_events, ~] = size(obj.data);
                percentage = cluster_sizes(bead_cluster_index)/num_events;
                
                if(percentage < approx_percentage)
                    found_bead_cluster = 1;
                    s = sprintf('found a bead cluster! %d components %f percent of events \n', i, percentage)
                    break;
                end
                
            end
            
            if(found_bead_cluster == 0)
                
                s = sprintf('Did not find a bead cluster ! \n')
            else
                
                s = sprintf('Found a bead cluster ! \n')
                num_events = length(obj.beads);
                
                for i=1:num_events
                    
                    if(obj.beads(i) == bead_cluster_index)
                        %just add a criteria about how strong you want that
                        %cluster potential to be
                        if(obj.bead_post(i,bead_cluster_index)>.7)
                            if(obj.dna_low(i) == 1)
                                obj.beads(i) = 1;
                                sprintf('found a bead');
                            else
                                obj.beads(i) = 0;
                            end
                        else
                            obj.beads(i) = 0;
                        end
                    else
                        obj.beads(i) = 0;
                    end
                end
                
                
                %                 bead_data = obj.data(bead_indices,obj.bead_channels);
                %                 all_data = obj.data(:, obj.bead_channels);
                %                 within_cluster_distances = mahal(bead_data, bead_data);
                %                 distance_threshold = 2*std(within_cluster_distances);
                %                 distances_to_beads = mahal(all_data, bead_data);
                %
                %
                %
                %
                %                 new_bead_indices = find(distances_to_beads<distance_threshold);
                %
                %                 for m=1:length(new_bead_indices)
                %
                %                     index = new_bead_indices(m);
                %
                %                     if(obj.beads(index)==1)
                %                         continue;
                %                     end
                %
                %                     if(obj.dna_low(index)==1)
                %                         obj.beads(index)=1;
                %                     end
                %                 end
                
                
            end
            
            obj.clean_beads();
            obj.debris = obj.dna_low;
            
        end
        
        
        
        
        function bead_indices = get_bead_indices(obj)
            
            bead_indices = find(obj.beads==1);
            
        end
        
        function obj = correct_channel_by_marker(obj, channel, windowsize)
            
            
            smooth_data = zeros(1,num_events-window_size);
            smooth_index = zeros(1,num_events-window_size);
            
            for i=1:num_events-window_size
                
                smooth_data(i) = median(event_data(i:i+window_size,channel));
                smooth_index(i) = median(event_data(i:i+window_size,1));
                
            end
            
            
        end
        
        function obj = plot_bead_channels_vs_dna(obj)
            
            dna_channel = obj.dna_channels(1);
            
            for i=1:length(obj.bead_channels)
                
                subplot(1,5,i)
                bead_channel = obj.bead_channels(i);
                data = [obj.data(:,bead_channel) obj.data(:,dna_channel)];
                
                [bandwidth,density,X,Y]=kde2d(data);
                
                
                contour(X,Y,density,30), hold on,
                plot(data(:,1),data(:,2),'b.','MarkerSize',5)
                
                bead_indices = find(obj.beads==1);
                
                data = [obj.data(bead_indices,bead_channel) obj.data(bead_indices,dna_channel)];
                
                [bandwidth,density,X,Y]=kde2d(data);
                
                
                contour(X,Y,density,30), hold on,
                plot(data(:,1),data(:,2),'r.','MarkerSize',5)
                
                
                
                
            end
            
            
        end
        
        function obj = plot_bead_channel_density(obj, channel1, channel2)
            
            
            
            data = [obj.data(:,channel1) obj.data(:,channel2)];
            
            [bandwidth,density,X,Y]=kde2d(data);
            
            
            contour(X,Y,density,30), hold on,
            plot(data(:,1),data(:,2),'b.','MarkerSize',5)
            
            bead_indices = find(obj.beads==1);
            
            data = [obj.data(bead_indices,channel1) obj.data(bead_indices,channel2)];
            
            [bandwidth,density,X,Y]=kde2d(data);
            
            
            contour(X,Y,density,30), hold on,
            plot(data(:,1),data(:,2),'r.','MarkerSize',5)
            
            
        end
        
        function obj = manual_bead_identify(obj, filter_channel1, filter_channel2)
            
            fhandle = figure;
            
            obj.plot_2d_channel_density(filter_channel1, filter_channel2);
            
            rect = getrect(fhandle)
            left = rect(1);
            bottom = rect(2);
            right = rect(3);
            top = rect(4);
            
            %format for the rectangle is left bottom right top
            
            channel1 = obj.name_channel_map(filter_channel1);
            channel2 = obj.name_channel_map(filter_channel2);
            
            [num_events, ~] = size(obj.data);
            obj.beads = zeros(num_events,1);
            
            for i=1:num_events
                
                channel1_data = obj.data(i,channel1);
                channel2_data = obj.data(i, channel2);
                if ((left<channel1_data) && (channel1_data<left+right))
                    
                    if ((bottom<channel2_data) && (channel1_data<bottom+top))
                        obj.beads(i) = 1;
                    end
                end
                
            end
            
            num_beads = size(find(obj.beads == 1))
            bead_percentage = num_beads/num_events
            
        end
        
        function obj = compute_dna_gmm_clusters(obj, num_clusters, filter_channels)
            
            
            clustering_data = obj.data(:,filter_channels);
            options = statset('MaxIter',100);
            
            gmfit = gmdistribution.fit(clustering_data,num_clusters,'options',options);
            
            
            obj.dna_gmfit = gmfit;
            
            [dna_low,nlogl,post] = gmfit.cluster(clustering_data);
            
            
            %wait i don't know which cluster is which
            %otherwise I may have to switch 0's and ones
            
            dna_cluster = 0;
            min_sum = 100000;
            
            for i=1:length(gmfit.mu)
                
                if(sum(gmfit.mu(i))<min_sum)
                    
                    min_sum = sum(gmfit.mu(i));
                    dna_cluster = i;
                end
                
            end
            
            dna_cluster
            
            for i=1:length(dna_low)
                
                if(dna_low(i) == dna_cluster)
                    
                    dna_low(i) = 1;
                else
                    
                    dna_low(i) = 0;
                end
                
            end
            
            
            %             for i=1:length(dna_low)
            %
            %                  if(dna_low(i) == 0)
            %                      continue;
            %
            %                  end
            %
            %                 if(post(i) < .97)
            %
            %
            %                     dna_low(i) = 0;
            %
            %                end
            %
            %             end
            
            
            obj.dna_low = dna_low;
            
        end
        
        
        function obj=compute_bead_gmm_clusters(obj, num_clusters, filter_channels)
            
            
            %n-dimensional gaussian mixture model
            % Inputs:
            %   X(n,d) - input data, n=number of observations, d=dimension of variable
            %   k - maximum number of Gaussian components allowed
            %   ltol - percentage of the log likelihood difference between 2 iterations ([] for none)
            %   maxiter - maximum number of iteration allowed ([] for none)
            %   pflag - 1 for plotting GM for 1D or 2D cases only, 0 otherwise ([] for none)
            %   Init - structure of initial W, M, V: Init.W, Init.M, Init.V ([] for none)
            obj.beads = [];
            
            bead_data=obj.data(:,filter_channels);
            %gchannels = [3 4];
            %bead_data=obj.data(:,gchannels);
            invalids = find(obj.invalid==1);
            valids = find(obj.invalid==0);
            bead_data(invalids,:)=[];
            
            
            %learning without the obviously invalid ones
            
            %d = number of barcoding channels here.
            
            % Ouputs:
            %   W(1,k) - estimated weights of GM
            %   M(d,k) - estimated mean vectors of GM
            %   V(d,d,k) - estimated covariance matrices of GM
            %   L - log likelihood of estimates
            
            options = statset('MaxIter',400);
            
            gmfit = gmdistribution.fit(bead_data,num_clusters,'options',options);
            %one should be like a crap-cluster as far as I can tell
            
            % mu           A K-by-D matrix of component means.
            % Sigma        An array or a matrix containing the component covariance
            %        matrices.  Sigma is one of the following
            %           * A D-by-D-by-K
            %[W,obj.bc_nd_gmm_means,obj.bc_nd_gmm_covars,L] = EM_GM(bc_data,20)
            
            
            obj.bead_gmfit = gmfit;
            
            %bead_data = obj.data(:,filter_channels);
            
            [beads,obj.bead_nlogl,bead_post] = gmfit.cluster(bead_data);
            
            num_valids = length(valids);
            [num_total_events,~] = size(obj.data);
            obj.beads = zeros(num_total_events,1);
            [~,cols] = size(bead_post);
            obj.bead_post = zeros(num_total_events, cols);
            
            for i=1:num_valids
                obj.beads(valids(i)) = beads(i);
                obj.bead_post(valids(i),:) = bead_post(i,:);
            end
            obj.beads(invalids) = 0; %reset the invalids to 0
            num_invalids = length(invalids);
            for i=1:num_invalids
                obj.bead_post(i,:) = bead_post(1,:);
            end
            %filling it with some random crap
            
        end
        
        function obj = remove_beads(obj)
            
            bead_indices = find(obj.beads==1);
            
            bead_data = obj.data(bead_indices, obj.bead_channels);
            
            obj.data(bead_indices,:)=[];
            obj.beads(bead_indices)=[];
            %               obj = obj.remove_debris();
            %               %currently removing all debris including beads
            %
            %
            %               distance_from_beads = mahal(obj.data(:,obj.bead_channels), bead_data);
            %
            %               distance_from_itself = mahal(bead_data, bead_data);
            %
            %               threshold = max(distance_from_itself);
            %
            %               bead_doublet_indices = find(distance_from_beads<threshold);
            %
            %               obj.data(bead_doublet_indices,:)=[];
            
        end
        
        function obj = debarcode_general(obj, well_key, num_wells, post_threshold, mahalad_threshold, varargin)
            
            %bc_estimates = obj.compute_bc_estimates_general(well_key, num_wells);
            bc_estimates = obj.compute_bc_estimates_extrema();
            
            num_events = size(obj.data,1);
            num_bc_channels = length(obj.barcoding_channels)
            obj.invalid = zeros(num_events,1);
            valid = 0;
            invalid = 0;
            obj.wells = zeros(num_events,1);
            optargin = size(varargin,2);
            bc_estimates_user = [];
            
            if(optargin>0)
                bc_estimates_user = varargin{1};
            end
            
            
            for i=1:num_events
                
                
                ambiguous = 0;
                barcode_value = 0;
                
                for j=1:num_bc_channels
                    
                    channel = obj.barcoding_channels(j);
                    if(length(bc_estimates_user)>0)
                        cutoff_low = bc_estimates_user(j,1);
                        cutoff_high = bc_estimates_user(j,2);
                    else
                        
                        cutoff_low = bc_estimates(j,1);
                        cutoff_high = bc_estimates(j,2);
                    end
                    
                    if (cutoff_low < obj.data(i,channel)) && (obj.data(i,channel) < cutoff_high)
                        sprintf('ambiguous event! \n');
                        
                        ambiguous = 1;
                        obj.invalid(i) =1;
                        
                    end
                    
                    
                    
                    if(obj.data(i,channel)>cutoff_high)
                        
                        barcode_value = barcode_value + 2^(num_bc_channels-j);
                        
                    end
                    
                    
                    
                end
                
                if(ambiguous==0)
                    valid = valid+1;
                    sprintf('valid event! \n');
                    if(barcode_value == 0)
                        obj.wells(i)  = 0;
                    else
                        obj.wells(i) = well_key(barcode_value+1);
                    end
                end
                
                
                
                
            end
            
            
            
            invalids = find(obj.wells==0);
            
            num_invalids = size(invalids);
            
            obj = obj.compute_barcoding_nd_gmm_clusters(num_wells, post_threshold, mahalad_threshold);
            
            invalids_after_clustering = find(obj.wells==0);
            
            num_invalids_after_clustering = size(invalids_after_clustering);
            
        end
        
        
        
        
        
        function obj = debarcode_aml_single_cell_estimate(obj, well_key_aml,high_low_threshold, mahalad_threshold)
            
            %bc_estimates = obj.compute_bc_estimates();
            
            obj = obj.compute_3one_bc_normalize(well_key_aml, high_low_threshold);
            
            invalids = find(obj.wells==0);
            
            num_invalids = size(invalids)
            
            %obj = obj.compute_barcoding_nd_gmm_clusters(20, post_threshold, mahalad_threshold);
            obj = obj.mahalanobis_correction_well_clusters(20, mahalad_threshold);
            
            invalids_after_clustering = find(obj.wells==0);
            
            num_invalids_after_clustering = size(invalids_after_clustering)
            
            
        end
        
        
        function obj = debarcode_aml(obj, well_key_aml, post_threshold, mahalad_threshold)
            
            bc_estimates = obj.compute_bc_estimates();
            
            obj = obj.compute_3one_bc(bc_estimates, well_key_aml);
            
            invalids = find(obj.wells==0);
            
            num_invalids = size(invalids)
            
            obj = obj.compute_barcoding_nd_gmm_clusters(20, post_threshold, mahalad_threshold);
            
            invalids_after_clustering = find(obj.wells==0);
            
            num_invalids_after_clustering = size(invalids_after_clustering)
            
            
        end
        
        function obj = mahalanobis_correction_well_clusters(obj, num_wells, mahal_threshold)
            
            
            for i=1:num_wells
                well_indices = find(obj.wells == i);
                num_cells = length(well_indices);
                well_matrix = obj.data(well_indices, obj.barcoding_channels);
                D2 = mahal(well_matrix, well_matrix);
                
                for j=1:length(well_indices)
                    well_index = well_indices(j);
                    
                    if(D2(j)<mahal_threshold)
                        
                        obj.wells(well_index) = i;
                    else
                        obj.wells(well_index) = 0;
                    end
                end
                
            end
            
            
            
            
            
            
        end
        
        
        function obj=compute_barcoding_nd_gmm_clusters(obj, num_wells, post_threshold, mahalad_threshold)
            
            
            %n-dimensional gaussian mixture model
            % Inputs:
            %   X(n,d) - input data, n=number of observations, d=dimension of variable
            %   k - maximum number of Gaussian components allowed
            %   ltol - percentage of the log likelihood difference between 2 iterations ([] for none)
            %   maxiter - maximum number of iteration allowed ([] for none)
            %   pflag - 1 for plotting GM for 1D or 2D cases only, 0 otherwise ([] for none)
            %   Init - structure of initial W, M, V: Init.W, Init.M, Init.V ([] for none)
            bc_data=obj.data(:,obj.barcoding_channels);
            invalids = find(obj.wells==0);
            bc_data(invalids,:)=[];
            presorted_wells = obj.wells;
            presorted_wells(invalids,:)=[];
            %learning without the obviously invalid ones
            
            %d = number of barcoding channels here.
            
            % Ouputs:
            %   W(1,k) - estimated weights of GM
            %   M(d,k) - estimated mean vectors of GM
            %   V(d,d,k) - estimated covariance matrices of GM
            %   L - log likelihood of estimates
            
            options = statset('MaxIter',400);
            
            %I need to put in a start parameter with the clusters in
            
            gmfit = gmdistribution.fit(bc_data,num_wells,'options',options, 'start', presorted_wells);
            
            
            % mu           A K-by-D matrix of component means.
            % Sigma        An array or a matrix containing the component covariance
            %        matrices.  Sigma is one of the following
            %           * A D-by-D-by-K
            %[W,obj.bc_nd_gmm_means,obj.bc_nd_gmm_covars,L] = EM_GM(bc_data,20)
            
            obj.bc_nd_gmm_means = gmfit.mu;
            obj.bc_nd_gmm_covars = gmfit.Sigma;
            obj.gmfit = gmfit;
            %should I just put this in the wells again because this time it
            %should have which well its in correct
            cluster_data=obj.data(:,obj.barcoding_channels);
            [obj.wells,obj.bc_nlogl,obj.bc_post, logpdf, obj.bc_mahalad] = gmfit.cluster(cluster_data);
            
            num_events = length(obj.wells)
            
            
            
            
            for i=1:num_events
                
                well_no = obj.wells(i);
                if(obj.bc_post(i,well_no)<post_threshold)
                    
                    obj.wells(i) = 0;
                end
                
                %for p-value = .05 the mahalanobis distance is approx 13
                %for a variable with 6-degrees of freedom
                if(obj.bc_mahalad(i,well_no)>mahalad_threshold)
                    
                    obj.wells(i) = 0;
                end
                
            end
            
            for i=1:length(invalids)
                
                index = invalids(i);
                obj.wells(index) = 0;
                
            end
            
            %not allowing any extra assignments
            
        end
        
        %now I need something that makes clusters out of these
        function obj=compute_linkage_clusters(obj, dist_channels, num_clusters)
            
            dist_data=obj.data(:,dist_channels);
            Y = pdist(dist_data,'euclid');
            Z = linkage(Y, 'average');
            obj.dist_clusters = cluster(Z, 'cutoff',1.1547005);
            %cutoff = 1.1547005 gives 32 clustesr for beads
            
            %obj.dist_clusters = clusterdata(dist_data, 1.17); %1 is the 'cutoff'
            
        end
        
        
        
        function write_dist_graph(obj, filename)
            
            file = fopen(filename,'w');
            
            for i=1:length(obj.dist_graph)
                fprintf(file,'%d %d \n', obj.dist_graph(i,1),obj.dist_graph(i,2));
            end
            fclose(file);
            
        end
        
        
        function write_data_matrix(obj, filename, channels)
            
            
            M = obj.data(:,channels);
            size_vector = size(M);
            dlmwrite(filename, size_vector, 'delimiter', '\t');
            dlmwrite(filename, M,'delimiter','\t', '-append');
            
            
            
        end
        
        function write_data_matrix_csv(obj, filename)
            
            %Do I want to write out all the channels, I guess so right
            [num_elements, num_channels] = size(obj.data);
            
            fid = fopen(filename,'w');
            for i=1:num_channels-1
                
                
                fprintf(fid,'%s,',obj.channel_name_map{i});
                
                
                
            end
            
            fprintf(fid,'%s \n',obj.channel_name_map{num_channels});
            fclose(fid);
            
            for i=1:num_channels
                
                
                dlmwrite(filename, obj.data,'delimiter',',', '-append');
                
            end
            
            
        end
        
        function obj = compute_distance_graph(obj, dist_channels, n, k)
            %k is for kNN clustering so how many neighbors you want
            %n should probably be more than k its how many candidates you
            %examine
            obj.dist_graph = [];
            num_dist_channels = length(dist_channels);
            num_events = size(obj.data,1);
            bc_sorts = zeros(num_events, num_bc_channels);
            
            
            
            for i=1:num_dist_channels
                
                [~, bc_sorts(:,i)] = sort(obj.data(:,dist_channels(i)));
                
            end
            
            obj.dist_graph = zeros(num_events*k,2);
            graph_index = 1;
            
            for i=1:num_events
                
                candidates = [];
                for j=1:num_dist_channels
                    
                    event_index = find(bc_sorts(:,j)==i);
                    start_pos = max(event_index-n,1);
                    end_pos = min(event_index+n, num_events);
                    candidates = [candidates; bc_sorts(start_pos:end_pos,j)];
                end
                
                candidates = unique(candidates, 'rows');
                self_node = find(candidates==i);
                candidates(self_node)=[];
                candidate_distances = zeros(1,length(candidates));
                for j=1:length(candidates)
                    
                    candidate_distances(j) = obj.compute_distance(i,candidates(j),dist_channels);
                    
                end
                
                
                [~, neighbors] = sort(candidate_distances);
                
                num_neighbors = min(k,length(neighbors));
                
                for j=1:num_neighbors
                    
                    neighboring_event = candidates(neighbors(j));
                    
                    obj.dist_graph(graph_index,1) = i;
                    obj.dist_graph(graph_index,2) = neighboring_event;
                    graph_index = graph_index + 1;
                end
                
                
                
            end
            
        end
        
        function dist = compute_bc_distance(obj, event1, event2, dist_channels)
            
            dist = 0;
            for i=1:length(dist_channels)
                
                channel = dist_channels(i);
                val1 = obj.data(event1,channel);
                val2 = obj.data(event2,channel);
                dist = dist+(val1-val2)^2;
            end
            
            dist = sqrt(dist);
            
        end
        
        function obj = import_bc_clusters(obj, cluster_file)
            
            obj.bc_clusters = importdata(cluster_file);
            
        end
        
        
        %ok if i do this first and then do a 20 element cluster what do I
        %get?
        
        
        function bc_estimates = compute_bc_estimates(obj)
            
            sprintf('in compute bc estimates \n')
            %take each of the barcoding channels
            %compute their cdf and then see when the cdf is about 50%
            %is this a good estimate? give some leeway for the estimates
            %later;
            num_bc_channels = length(obj.barcoding_channels);
            bc_estimates=zeros(length(num_bc_channels), 2);
            %I will change the bc_estimates so that it does both low and
            %high things
            
            yi = linspace(.01,.99,99);
            
            for i=1:num_bc_channels
                
                bc_channel = obj.barcoding_channels(i);
                g = ksdensity(obj.data(:,bc_channel),yi,'function','icdf');
                
                %first column will have the low cutoffs
                %second column will give me the high cutoffs
                
                bc_estimates(i,1) = g(49);
                bc_estimates(i,2) = g(51);
                %should give the value thats at the 50th percentile
                
            end
            
            
            
        end
        
        
        function bc_estimates = compute_bc_estimates_extrema(obj)
            
            sprintf('in compute bc estimates extrema \n')
            %take each of the barcoding channels
            %compute their cdf and then see when the cdf is about 50%
            %is this a good estimate? give some leeway for the estimates
            %later;
            
            num_bc_channels = length(obj.barcoding_channels);
            
            bc_estimates=zeros(length(num_bc_channels), 2);
            
            %I will change the bc_estimates so that it does both low and
            %high things
            
            yi = linspace(.01,.99,99);
            
            for i=1:num_bc_channels
                
                bc_channel = obj.barcoding_channels(i);
                [f,xi] = ksdensity(obj.data(:,bc_channel));
                [xmax,imax,xmin,imin] = extrema(f)
                imax_vals = xi(imax);
                imin_vals = xi(imin);
                
                g = ksdensity(obj.data(:,bc_channel),yi,'function','icdf');
                c = ksdensity(obj.data(:,bc_channel),xi,'function','cdf');
                
                if(length(imin)==3)
                    
                    sprintf('only 3 mins\n')
                    [sorted_imin, sorted_imin_index] = sort(imin, 'ascend');
                    %should be second one
                    second_imin_index = sorted_imin_index(2)
                    second_imin_value = xi(imin(second_imin_index))
                    cdf_index = find(xi==second_imin_value);
                    cdf_value = c(cdf_index)
                    cdf_value_low = floor((cdf_value*100))
                    cdf_value_high = ceil((cdf_value*100))
                    
                    
                    
                    bc_estimates(i,1) = g(cdf_value_low);
                    bc_estimates(i,2) = g(cdf_value_high);
                    
                    
                    continue;
                end
                
                %if min and max are too close get rid of them from the list but which one?
                %both,
                
                combined_imin_imax = [imin imax]
                combined_xmin_xmax = [xmin xmax]
                imin_or_imax_binary = [zeros(1,length(imin)) ones(1,length(imax))]
                
                [sorted_combined_imin_imax, sorted_combined_indices] = sort(combined_imin_imax,'ascend')
                sorted_combined_xmin_xmax = combined_xmin_xmax(sorted_combined_indices)
                sorted_imin_or_imax_binary = imin_or_imax_binary(sorted_combined_indices)
                remove_list = [];
                blip_threshold = .02;
                
                for j=1:length(combined_imin_imax)-1
                    
                    extrema1 = sorted_combined_xmin_xmax(j);
                    extrema2 = sorted_combined_xmin_xmax(j+1);
                    if(abs(extrema1-extrema2)<blip_threshold)
                        remove_list = [j j+1];
                    end
                end
                
                remove_list
                sorted_combined_imin_imax(remove_list)=[];
                sorted_combined_xmin_xmax(remove_list)=[];
                sorted_imin_or_imax_binary(remove_list) =[];
                
                imin_indices = find(sorted_imin_or_imax_binary==0);
                imax_indices = find(sorted_imin_or_imax_binary==1);
                imax = sorted_combined_imin_imax(imax_indices);
                imin = sorted_combined_imin_imax(imin_indices);
                xmax = sorted_combined_xmin_xmax(imax_indices);
                xmin = sorted_combined_xmin_xmax(imin_indices);
                
                
                %criteria1: closest one to this --does it really have to be
                %close to 50%?? i guess the next criteria may subsume this
                %one.
                [sorted_imin, sorted_imin_index] = sort(imin, 'ascend');
                second_imin_index = sorted_imin_index(2);
                second_imin_value = xi(imin(second_imin_index));
                
                
                
                
                cdf_index = find(xi==second_imin_value);
                cdf_value = c(cdf_index)
                cdf_value_low = floor((cdf_value*100))
                cdf_value_high = ceil((cdf_value*100))
                
                
                
                bc_estimates(i,1) = g(cdf_value_low);
                bc_estimates(i,2) = g(cdf_value_high);
                
                
            end
            
            
        end
        
        
        function bc_estimates = compute_bc_estimates_general(obj, well_key, num_wells)
            
            num_bc_channels = length(obj.barcoding_channels);
            bc_counts = zeros(1,num_bc_channels);
            for i=1:length(well_key)
                
                if(well_key(i) == 0)
                    continue;
                end
                
                for j=1:num_bc_channels
                    
                    bc_channel = num_bc_channels - j +1;
                    
                    if(bitget(i, bc_channel)==1)
                        
                        bc_counts(j) = bc_counts(j)+1;
                        
                    end
                    
                end
                
                
                
            end
            
            
            
            bc_counts = bc_counts ./ num_wells;
            
            bc_counts = bc_counts .* 100
            
            
            %the well key that I usually give it is the computed value
            %but I can take the computed value and change it into binary
            
            
            sprintf('in compute bc estimates \n')
            %take each of the barcoding channels
            %compute their cdf and then see when the cdf is about 50%
            %is this a good estimate? give some leeway for the estimates
            %later;
            
            bc_estimates=zeros(length(num_bc_channels), 2);
            %I will change the bc_estimates so that it does both low and
            %high things
            
            yi = linspace(.01,.99,99);
            
            for i=1:num_bc_channels
                
                bc_channel = obj.barcoding_channels(i);
                g = ksdensity(obj.data(:,bc_channel),yi,'function','icdf');
                
                %first column will have the low cutoffs
                %second column will give me the high cutoffs
                
                low_cutoff = floor(bc_counts(i))-1;
                high_cutoff = ceil(bc_counts(i))+1;
                
                bc_estimates(i,1) = g(low_cutoff);
                bc_estimates(i,2) = g(high_cutoff);
                %should give the value thats at the 50th percentile
                
            end
            
            
            
        end
        
        
        function obj = compute_3one_bc_normalize(obj, well_key_aml, threshold)
            
            %added some leeway to the barcoding experiments
            num_events = size(obj.data,1);
            num_bc_channels = length(obj.barcoding_channels);
            obj.invalid = zeros(num_events,1);
            
            obj.wells = zeros(num_events,1);
            bc_data = obj.data(:, obj.barcoding_channels);
            normalized_bc_data = zscore(bc_data);
            
            
            for i=1:num_events
                
                [sorted_bc_data, sorted_bc_indices] = sort(normalized_bc_data(i,:),'descend');
                top_3_inds = sorted_bc_indices(1:3);
                
                barcode_value = 2^(num_bc_channels-sorted_bc_indices(1));
                barcode_value = barcode_value + 2^(num_bc_channels-sorted_bc_indices(2));
                barcode_value = barcode_value + 2^(num_bc_channels-sorted_bc_indices(3));
                
                if ( (sorted_bc_data(3)-sorted_bc_data(4))<threshold)
                    obj.invalid(i) = 1;
                    obj.wells(i) = 0;
                else
                    
                    obj.wells(i) = well_key_aml(barcode_value+1);
                    if(obj.wells(i)==0)
                        sprintf('alert! invalid bc value %d \n',barcode_value)
                    end
                end
                
                
            end
            
        end
        
        function obj = compute_3one_bc(obj, bc_estimates, well_key_aml)
            
            %added some leeway to the barcoding experiments
            num_events = size(obj.data,1);
            num_bc_channels = length(obj.barcoding_channels);
            obj.invalid = zeros(num_events,1);
            valid = 0;
            invalid = 0;
            obj.wells = zeros(num_events,1);
            
            for i=1:num_events
                
                count_zeros = 0;
                count_ones = 0;
                barcode_value = 0;
                
                for j=1:num_bc_channels
                    
                    channel = obj.barcoding_channels(j);
                    cutoff = bc_estimates(j,1);
                    
                    if(obj.data(i,channel)<cutoff)
                        
                        count_zeros = count_zeros+1;
                        
                    end
                    
                    cutoff = bc_estimates(j,2);
                    
                    
                    
                    if(obj.data(i,channel)>cutoff)
                        count_ones = count_ones+1;
                        barcode_value = barcode_value + 2^(num_bc_channels-j);
                        
                    end
                    
                end
                
                
                
                if(count_zeros > 3)
                    obj.invalid(i) = 1;
                    invalid = invalid+1;
                end
                if(count_ones > 3)
                    
                    obj.invalid(i) = 1;
                    invalid = invalid+1;
                    
                end
                
                if(count_ones == 3)
                    
                    valid = valid+1;
                    sprintf('valid event ! \n');
                    
                    obj.wells(i) = well_key_aml(barcode_value+1);
                end
                
                
            end
            
            valid
            
            invalids = find(obj.invalid==1);
            
            num_invalids = size(invalids)
            
            %obj.data(invalids,:)=[];
            
        end
        
        %need to compute cluster sizes, cluster averages
        function cluster_sizes = compute_bc_cluster_sizes(obj,bc_clusters)
            
            
            %num_clusters = max(obj.bc_clusters(:,2))+1;
            %I think it numbers the clusters from 0 to cluster-1
            num_clusters = max(bc_clusters);
            cluster_sizes = zeros(1,num_clusters);
            
            for i=1:num_clusters
                
                cluster_sizes(i) = length(find(bc_clusters==(i)));
            end
            
        end
        
        function cluster = get_cluster(obj, cluster_no)
            
            cluster = find(obj.bc_clusters==(i-1));
            
        end
        
        function obj = add_eventnum(obj)
            %assuming in order of time
            %fcs files seem to have ordered events anyways
            [r,c] = size(obj.data);
            eventnum = linspace(1,r,r);
            obj.data(:,c+1) = transpose(eventnum);
            obj.eventnum_channel = c+1;
            obj.marker_channels(find(obj.marker_channels == obj.eventnum_channel)) = [];
            
        end
        
        %need a function to compare two signals.
        function obj = compute_decay(obj)
            
            [~,n] = size(obj.data);
            for i=1:n
                channel = i;
                p = polyfit(obj.data(:, obj.eventnum_channel), obj.data(:,channel),1);
                obj.channel_decays(i) = p(1);
            end
            
            obj.mean_decay = mean(obj.channel_decays);
            obj.total_decay = sum(obj.channel_decays);
        end
        
        function P = compute_channel_decay(obj, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            
            P = polyfit(obj.data(:, obj.eventnum_channel), obj.data(:,channel),1);
            
            
        end
        
        
        function [smooth_index, smooth_data] =  plot_channel_smooth(obj, channel_name, window_size, varargin)
            
            channel = obj.name_channel_map(channel_name);
            
            [num_events, n] = size(obj.data);
            smooth_data = zeros(1,num_events-window_size);
            for i=1:num_events-window_size
                smooth_data(i) = median(obj.data(i:i+window_size,channel));
            end
            smooth_index=linspace(1,num_events-window_size,num_events-window_size);
            
            optargin = size(varargin,2);
            if(optargin==1)
                plot(smooth_index,smooth_data,varargin{1});
            else
                plot(smooth_index,smooth_data);
            end
            
            
        end
        
        function [smooth_index, smooth_data] =  plot_channel_vs_time_smooth(obj, channel_name, window_size, varargin)
            
            
            channel = obj.name_channel_map(channel_name);
            [num_events, n] = size(obj.data);
            smooth_data = zeros(1,num_events-window_size);
            smooth_index = zeros(1,num_events-window_size);
            for i=1:num_events-window_size
                smooth_data(i) = median(obj.data(i:i+window_size,channel));
                smooth_index(i) = median(obj.data(i:i+window_size,1));
            end
            
            optargin = size(varargin,2);
            if(optargin==1)
                plot(smooth_index,smooth_data,varargin{1});
            else
                plot(smooth_index,smooth_data);
            end
            
            
        end
        
        function obj = compute_invalid(obj)
            
            obj.invalid = find(obj.wells==0);
            
        end
        
        function obj = remove_invalid(obj)
            
            obj.invalid = find(obj.wells==0);
            obj.wells(obj.invalid) = [];
            obj.data(obj.invalid,:) = [];
            
        end
        
        
        function obj = remove_nonbeads(obj)
            
            non_bead_events = find(obj.beads==0);
            bead_events = find(obj.beads==1);
            
            sprintf('num bead events %d num nonbead events %d ', length(bead_events), length(non_bead_events))
            
            obj.data(non_bead_events,:) = [];
            
        end
        
        function obj = remove_nongated(obj)
            
            non_gated_events = find(obj.event_gate==0);
            gated_events = find(obj.event_gate==1);
            
            sprintf('num gated events %d num nongated events %d ', length(gated_events), length(non_gated_events))
            
            obj.data(non_gated_events,:) = [];
            
        end
        
        function obj = remove_debris(obj)
            
            
            %the low dna clusters that are not beads
            debris_indices = find(obj.debris==1);
            obj.data(debris_indices,:) = [];
        end
        
        
        function plot_channel_means(obj, channel_name, window_size)
            
            channel = obj.name_channel_map(channel_name);
            [num_events, ~] = size(obj.data);
            data_means=zeros(num_events-window_size,1);
            for i=1:num_events-window_size
                
                data_means(i) = mean(obj.data(i:i+window_size,channel));
            end
            
            plot(data_means);
            
        end
        
        function obj = compute_debris(obj)
            
            %For now lets try removing events with -DNA and -length
            [num_events, ~] = size(obj.data);
            
            obj.debris = zeros(1, num_events);
            for i=1:num_events
                
                debris_cluster1 = obj.debris_cluster_value(1,i);
                debris_cluster2 = obj.debris_cluster_value(2,i);
                
                if (debris_cluster1 == 0) || (debris_cluster2 == 0)
                    obj.debris(i) = 1;
                else
                    obj.debris(i) = 0;
                end
            end
            
            %If it seemed like it was generated from the first rather
            %than the second cluster then discarded as the debris
            
        end
        
        function plot_2d_channel_density(obj, channel1_name, channel2_name, varargin)
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            
            data = [obj.data(:,channel1) obj.data(:,channel2)];
            limits = [];
            for i=1:length(varargin)
                
                
                if(strcmp(varargin{i}, 'limits'))
                    
                    limits = varargin{i+1};
                    
                end
            end
            
            
            %[bandwidth,density,X,Y]=kde2d(data,256, [min(obj.data(:,channel1)) min(obj.data(:,channel2))], [max(obj.data(:,channel1)) max(obj.data(:,channel2))]);
            maxx = max(obj.data(:,channel1));
            
            if(length(limits)>0)
                
                [bandwidth,density,X,Y]=kde2d(data,256, [0 0], [limits(3) limits(4)]);
            else
                [bandwidth,density,X,Y]=kde2d(data,256, [0 0], [max(obj.data(:,channel1)) max(obj.data(:,channel2))]);
                %[bandwidth,density,X,Y]=kde2d(data,256);
            end
            
            %http://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
            %Dani recommends this method of density estimation
            %http://code.google.com/p/danapeerlab/source/browse/trunk/freecell/src/axes.py
            %python for doing this.
            %Dani says that the contour3 stuff shoudl work.
            
            %leaving out optional arguments
            %[bandwidth,density,X,Y]=kde2d(data,n,MIN_XY,MAX_XY)
            
            slices_size = maxx/8;
            
            
            
            %surf(X,Y,density)
            %view([0,70])
            %colormap hot
            %hold on
            %alpha(.8)
            %set(gca, 'color', 'blue');
            %plot(data(:,1),data(:,2),'w.','MarkerSize',5)
            
            optargin = length(varargin);
            if(optargin == 0)
                
                plot(data(:,1),data(:,2),'b.','MarkerSize',5)
                hold on,
                contour(X,Y,density,30),
            end
            %
            %plot(data(:,1),data(:,2),'b.','MarkerSize',5)
            
            for i=1:length(varargin)
                
                
                if(strcmp(varargin{i}, 'imagesc'))
                    
                    %                   density_filtered = density>.15;
                    %                   density = density.*density_filtered;
                    %                   j = linspecer(256);
                    %                   j(1,:) = [ 1 1 1 ];
                    %                   colormap(j);
                    colormap(linspecer(256));
                    imagesc(X(1,:),Y(:,1),density);
                    set(gca,'YDir','normal');
                    set(gca,'XTick',[]);
                    set(gca,'YTick',[]);
                    %plot_as_vertical_lines([1 2 3 4 5 6 7].*slices_size, 'w');
                    xlim([min(X(1, :)), max(X(1, :))])
                    ylim([min(Y(:, 1)), max(Y(:, 1))])
                end
                if(strcmp(varargin{i}, 'contour'))
                    
                    plot(data(:,1),data(:,2),'b.','MarkerSize',5)
                    hold on,
                    contour(X,Y,density,30),
                    set(gca,'XLim',[0 max(data(:,1))]);
                    set(gca,'YLim',[0 max(data(:,2))]);
                end
                
            end
            
            
            
        end
        
        
        
        
        function obj =  tsne_map_data(obj,  varargin)
            
            optargin = length(varargin);
            
            num_events = size(obj.data,1);
            
            
            channels = [];
            if(optargin == 0)
                channels = obj.marker_channels;
            else
                
                
                channel_names = varargin{1};
                for i=1:length(channel_names)
                    
                    channels = [channels obj.name_channel_map(channel_names{i})];
                    
                end
                %                channels = varargin{1};
                
            end
            channels
            obj.tsne_mapped_data = fast_tsne(obj.data(:,channels));
            
            
            
        end
        
        function draw_tsne(obj,varargin)
            
            
            if(length(varargin)>0)
                
                coloring_channel = varargin{1};
                channel = obj.name_channel_map(coloring_channel);
                %color_channel_data=varargin{1};
                
                color_channel_data = obj.data(:,channel);
                scatter(obj.tsne_mapped_data(:,1),obj.tsne_mapped_data(:,2),14 ,color_channel_data,'fill');
            else
                
                scatter(obj.tsne_mapped_data(:,1),obj.tsne_mapped_data(:,2),14);
                
                
            end
            
            
            
            
            %
            %
            %             [~, density, x, y] = kde2d([obj.tsne_mapped_data(:,1) obj.tsne_mapped_data(:,2)], 256);
            %
            %             contour(x, y, density, 5);
            %
            %             hold
            %
            %             if(length(varargin)>0)
            %
            % %                coloring_channel = varargin{1};
            % %                channel = obj.name_channel_map(coloring_channel);
            %
            % %                color_channel_data = obj.data(:,channel);
            %                 color_channel_data = varargin{1};
            %
            %                rx = linspace( 0, 1 ); % bins used in color coding
            %                colors = linspecer(256)( length( rx ) - 1 ); % coloring scheme
            %                color_channel_data = color_channel_data/max(color_channel_data);
            %
            %
            %
            %             	for i = 1:length(rx)-1
            %                     rx_indices = find( rx(i) <= color_channel_data & color_channel_data <= rx( i + 1 ) );
            %               		plot( obj.tsne_mapped_data( rx_indices, 1 ), obj.tsne_mapped_data( rx_indices, 2 ), '.', 'Color', colors( i, : ) );
            %
            %                 end
            %
            %             end
            
        end
        
        
        
        
        function compute_and_graph_render_edges_dremi(obj, edges, nodes, varargin)
            
            ctrl_specificed = 0;
            ctrl_data = [];
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'ctrl_data'))
                    
                    ctrl_specified = 1;
                    ctrl_data = varargin{i+1};
                    
                    
                end
            end
            
            adjMatrix = zeros(length(nodes), length(nodes));
            edges_indexed = zeros(length(edges),2);
            for i = 1:length(edges)
                
                node_index1 = find(strcmp(edges{i,1},nodes));
                node_index2 = find(strcmp(edges{i,2},nodes));
                adjMatrix(node_index1, node_index2) = 1;
                edges_indexed(i,1) = node_index1;
                edges_indexed(i,2) = node_index2;
                
                
            end
            gObj = biograph(adjMatrix,nodes);
            
            
            for i=1:length(edges)
                
                
                if(ctrl_specified==0)
                    dremi_values(i) = obj.compute_dremi(edges{i,1},edges{i,2},.80);
                else
                    
                    channel1_name = edges{i,1};
                    channel2_name = edges{i,2};
                    [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(obj, channel1_name, channel2_name, 25, 255);
                    [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(ctrl_data, channel1_name, channel2_name, 25, 255);
                    maxy = max(maxy1,maxy2);
                    [dremi_values(i),~] = obj.compute_dremi(channel1_name, channel2_name, noise_threshold, 'maxy', maxy);
                end
            end
            
            range = .65;
            load 'continuous_BuPu9.mat';
            colormap(continuous_BuPu9);
            
            
            
            for i = 1:length(dremi_values)
                
                [color_value] = get_color_value(dremi_values(i), .65, continuous_BuPu9);
                set(gObj.edges(i),'LineColor',color_value);
                set(gObj.edges(i),'LineWidth',2.0);
            end
            
            all_nodes = 1:length(nodes);
            set(gObj.nodes(all_nodes),'LineColor',[.4 .4 .4]);
            set(gObj.nodes(all_nodes),'Color',[1 1 1]);
            set(gObj.nodes(all_nodes),'Shape','ellipse');
            set(gObj.nodes(all_nodes),'LineWidth',1.1);
            set(gObj.nodes(all_nodes),'fontSize',14);
            view(gObj);
            dremi_values
        end
        
        
        
        function [channel_names] = get_marker_channel_names(obj)
            
            channel_names = cell(1,length(obj.marker_channels));
            
            for i=1:length(channel_names)
                
                channel_names{i} = obj.channel_name_map{obj.marker_channels(i)};
            end
            
        end
        
        
        function [weight_matrix] = compute_lingam(obj, channel_names)
            
            for i=1:length(channel_names)
                channel_numbers(i) = obj.name_channel_map(channel_names{i});
            end
            data_for_lingam = obj.data(:,channel_numbers);
            data_for_lingam = transpose(data_for_lingam);
            [B stde ci k W] = estimate(data_for_lingam);
            Bpruned = prune(data_for_lingam, k, 'method', 'olsboot', 'B', B);
            weight_matrix = transpose(Bpruned);
            
            %plot_causal_graph(channel_names, weight_matrix);
        end
        
        function obj = tsne_gate(obj, num_rects)
            
            fhandle = figure;
            mapped_data = obj.tsne_mapped_data;
            
            [~, density, x, y] = kde2d([mapped_data(:,1) mapped_data(:,2)], 256);
            
            contour(x, y, density, 12);
            
            
            gated_indices = [];
            num_events = size(mapped_data,1);
            
            
            for j=1:num_rects
                
                rect = getrect(fhandle)
                left = rect(1);
                bottom = rect(2);
                width = rect(3);
                height = rect(4);
                
                for i=1:num_events
                    
                    if((mapped_data(i,1)>left)&&(mapped_data(i,1)<left+width))
                        if((mapped_data(i,2)>bottom)&&(mapped_data(i,2)<bottom+height))
                            
                            gated_indices = union(gated_indices ,[i]);
                            
                        end
                    end
                    
                end
                
            end
            obj.tsne_gated_indices = gated_indices
            
            percentage = length(obj.tsne_gated_indices)/length(obj.sampled_events)
            
        end
        
        function visualize_markers(obj,varargin)
            
            marker_channels = obj.marker_channels;
            
            optargin = length(varargin);
            
            if(optargin>0)
                
                draw_channel_names = varargin{1};
                
                for i=1:length(draw_channel_names)
                    
                    new_channel = obj.name_channel_map(draw_channel_names(i))
                    marker_channels = [marker_channels new_channel];
                    
                end
                
            end
            
            ncols = 5;
            nrows = ceil(length(marker_channels)/ncols);
            
            for i=1:length(marker_channels)
                
                channel_name = obj.channel_name_map{marker_channels(i)};
                subplot(nrows,ncols,i);
                obj.plot_channel_density(channel_name,'b');
                
                xlabel(channel_name, 'fontsize', 11);
                
                
                
            end
            
        end
        
        
        function [signaling_direction_vector]  = visualize_cluster_markers(obj, cluster_indices, varargin)
            
            num_events = size(obj.data, 1);
            
            non_cluster_indices = linspace(1, num_events, num_events);
            non_cluster_indices(cluster_indices) = [];
            
            signaling_direction_vector = [];
            
            optargin = length(varargin);
            channels = [];
            if(optargin == 0)
                
                channels = obj.marker_channels;
            else
                
                channel_names = varargin{1};
                for i=1:length(channel_names)
                    
                    channel_no = obj.name_channel_map(channel_names{i});
                    channels = [channels channel_no];
                    
                end
                
                
            end
            
            ncols = 5;
            nrows = ceil(length(channels)/5);
            
            for i=1:length(channels)
                
                subplot(nrows,ncols,i);
                channel = channels(i);
                
                channel_name = obj.channel_name_map{channel};
                marker_data = obj.data(:,channel);
                
                non_cluster_data = marker_data(non_cluster_indices);
                cluster_data = marker_data(cluster_indices);
                
                %signaling_direction_vector{1,i} = channel_name;
                population_median = median(obj.data(:,channel));
                %signaling_direction_vector{2,i} = median(cluster_data)-median(population_median)/median(population_median);
                signaling_direction_vector(i) = median(cluster_data)-median(population_median)/median(population_median);
                
                [f, xi] = ksdensity(non_cluster_data);
                plot(xi,f,'LineWidth',1.3);
                hold
                f2 = ksdensity(cluster_data,xi);
                plot(xi,f2, 'r','LineWidth',1.3);
                xlabel(channel_name, 'fontsize', 11);
                
                
            end
            
            
            
        end
        
        function plot_2d_channel_scatter(obj, channel1_name, channel2_name,varargin)
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            channel1_data = obj.data(:,channel1);
            channel2_data = obj.data(:,channel2);
            low_x_indices = find(channel1_data<0);
            channel1_data(low_x_indices)=[];
            channel2_data(low_x_indices)=[];
            
            low_x_indices = find(channel2_data<0);
            channel1_data(low_x_indices)=[];
            channel2_data(low_x_indices)=[];
            
            optargin = length(varargin);
            if(optargin>0)
                plot(channel1_data,channel2_data,varargin{1},'*');
            else
                
                plot(channel1_data, channel2_data,'*');
            end
            minx = min(obj.data(:,channel1));
            miny = min(obj.data(:,channel2));
            maxx = max(obj.data(:,channel1));
            maxy = max(obj.data(:,channel2));
            xlim([minx maxx]);
            ylim([miny maxy]);
            
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            %box on
            
        end
        
        
        function pairwise_visual_with_density(obj, channel1_name, channel2_name, varargin)
            
            subplot(2,1,1);
            [points_x, points_y, normalized_density] = obj.pairwise_visualize(channel1_name, channel2_name, varargin);
            
            channel1 = obj.name_channel_map(channel1_name);
            
            f = ksdensity(obj.data(:,channel1), points_x);
            
            maxf = max(f);
            f = f/maxf;
            
            subplot(2,1,2);
            
            plot(points_x, f);
            
            
        end
        
        function [ p ] = compute_anovan(obj, categorical_variable_names, continuous_variable_name, num_conditions)
            %p = anovan(y,{g1 g2 g3});
            num_categorical = length(categorical_variable_names);
            cgroup = cell(1, num_categorical);
            for i = 1:num_categorical
                
                cgroup{i} = obj.categorize_variable(categorical_variable_names{i}, num_conditions);
            end
            channel = obj.name_channel_map(continuous_variable_name);
            y = obj.data(:,channel);
            
            p = anovan(y, cgroup, 'model','interaction');
            
            
        end
        
        function [categorical_value] = categorize_variable(obj, channel_name, num_conditions)
            
            categorical_value = zeros(size(obj.data,1),1);
            
            condition_increment = 1/num_conditions;
            
            for j = 1:num_conditions
                
                
                
                
                percent_start = 0+ (condition_increment*(j-1));
                percent_end = percent_start+condition_increment;
                [~, gated_indices] = obj.high_low_gate_events(channel_name, percent_start, percent_end);
                
                categorical_value(gated_indices) = j;
                
            end
            
            
            
            
            
        end
        
        function [conditional_dremi] = compute_conditional_dremi(obj, channel1_name, channel2_name, channel3_name, noise_threshold, varargin)
            
            %channel 1 is the additional variable that I am conditioning on
            %here, using compute_dremi as a blackbox as much as possible
            %density "reweighting" channel1 as well, but it has to be gross
            %slices
            
            
            num_slices = 4;
            
            
            for i=1:length(varargin)
                
                
                if(strcmp(varargin{i}, 'num_slices'))
                    num_slices = varargin{i+1};
                end
                
                
            end
            
            [minx, miny, maxx, maxy] = find_data_cutoffs(obj, channel1_name, channel2_name, 50, 255);
            
            
            xincrement = (maxx-minx)/num_slices;
            
            edges_x = [minx];
            
            for i=1:num_slices-1;
                
                edges_x = [edges_x edges_x(i)+xincrement];
                
                
            end
            
            edges_x = [edges_x maxx];
            
            
            
            conditional_dremi = 0;
            total_events = size(obj.data,1);
            
            for i=1:(length(edges_x)-1)
                
                
                cdata_thresholded = obj.high_low_gate_events_interval(channel1_name, edges_x(i), edges_x(i+1));
                num_events = size(cdata_thresholded.data,1);
                %prob = num_events/total_events;
                prob = 1/num_slices;
                dremi = cdata_thresholded.compute_dremi(channel2_name, channel3_name, noise_threshold);
                conditional_dremi = conditional_dremi + (dremi*prob);
            end
            
            
        end
        
        
        
        
        
        function [PT] = compute_PT(obj, y_channel)
            
            y_num = obj.name_channel_map(y_channel);
            y_data = obj.discretized_data(:,y_num);
            y_states = max(y_data);
            data_size = length(y_data);
            
            PT = zeros(y_states,1);
            
            for i=1:y_states
                
                y_indices = find(y_data ==i);
                PT(i) = length(y_indices)/data_size;
            end
            
            
        end
        
        function [CPT] = compute_CPT(obj, y_channel, x_channels)
            
            %dimension of this matrix should be
            %channel code biinary to decimal
            channel_nums = [];
            for i=1:length(x_channels)
                
                channel_nums = [channel_nums obj.name_channel_map(x_channels{i})];
            end
            
            
            y_num = obj.name_channel_map(y_channel);
            
            num_x_channels = length(x_channels);
            xchannel_data = obj.discretized_data(:,channel_nums);
            y_data = obj.discretized_data(:,y_num);
            
            
            xchannel_states = max(xchannel_data);
            num_CPT_rows = prod(xchannel_states);
            y_states = max(y_data);
            
            CPT = zeros(num_CPT_rows,y_states);
            
            counter = [];
            
            for i=1:num_x_channels
                svector = transpose(1: xchannel_states(i));
                if(length(counter)==0)
                    counter = svector;
                else
                    prev_counter_length = size(counter,1);
                    counter = repmat(counter,xchannel_states(i),1);
                    rep_counter = ones(prev_counter_length,1);
                    new_counter = kron(svector,rep_counter);
                    counter = [new_counter counter];
                end
            end
            
            
            
            for i=1:num_CPT_rows
                
                
                loa = ismember(xchannel_data, counter(i,:), 'rows');
                indices = find(loa==1);
                valid_indices = find(y_data>0);
                indices = intersect(indices,valid_indices);
                
                for j=1:y_states
                    
                    y_indices = find(y_data(indices) == j);
                    CPT(i,j) = length(y_indices)/length(indices);
                    
                end
            end
            
            
            
        end
        
        
        function [mi_CPT] = compute_CPT_mi(obj, y_channel, x_channels)
            
            x_channels
            CPT = obj.compute_CPT(y_channel, x_channels)
            inputPT = [];
            for i=1:length(x_channels);
                
                if(i==1)
                    inputPT = obj.compute_PT(x_channels{1});
                else
                    
                    PT = obj.compute_PT(x_channels{i});
                    inputPT = kron(PT,inputPT);
                    
                end
            end
            
            
            outputPT = obj.compute_PT(y_channel)
            
            [num_rows, num_cols] = size(CPT);
            
            entropy = -1 * sum(log2(outputPT).* outputPT);
            centropy = 0;
            
            for i=1:num_rows
                
                row_entropy = -1 * (sum(log2(CPT(i,:)).* CPT(i,:)));
                centropy = centropy+(row_entropy*inputPT(i));
                
            end
            
            
            
            mi_CPT = entropy-centropy;
            
        end
        
        function [obj] = discretize_channel_from_edge(obj, channel1_name, channel2_name, num_levels, noise_threshold, varargin)
            
            
            transpose_edge = 0;
            
            
            for i=1:length(varargin)
                
                
                
                if(strcmp(varargin{i}, 'transpose_edge'))
                    transpose_edge = 1;
                end
                
                
            end
            
            
            [opt_partition, noise] = obj.partition_y(channel1_name, channel2_name, num_levels, noise_threshold, varargin);
            
            if(transpose_edge==0)
                obj = obj.discretize_channel(channel2_name, opt_partition, noise);
            else
                transpose_string = sprintf('transpose_edge')
                obj = obj.discretize_channel(channel1_name, opt_partition, noise);
            end
            
            
        end
        
        
        
        function [obj] = discretize_channel(obj, channel_name, opt_partition, noise)
            
            [num_events, num_channels] = size(obj.data);
            
            if(length(obj.discretized_data)==0)
                
                obj.discretized_data = zeros(num_events, num_channels);
            end
            
            channel_num = obj.name_channel_map(channel_name);
            
            discrete_vector = zeros(num_events,1);
            
            minval = min(obj.data(:,channel_num));
            maxval = max(obj.data(:,channel_num));
            opt_partition = [opt_partition maxval];
            channel_data = obj.data(:, channel_num);
            for i=1:length(opt_partition)
                
                indices = find((channel_data>=minval)&(channel_data<opt_partition(i)));
                discrete_vector(indices) = i;
                minval = opt_partition(i);
            end
            discrete_vector = discrete_vector .* noise;
            obj.discretized_data(:,channel_num) = discrete_vector;
            
        end
        
        
        
        
        
        function [opt_partition] = partition_y_from_edges(obj, edge_list, channel_name, num_levels, merge_threshold)
            % function [new_opt_partition] = discretize_y_from_edges(obj, opt_partition)
            % merge_threshold = .2;
            num_levels
            
            opt_partition = [];
            for i = 1:length(edge_list)
                
                channel1_name = edge_list{i,1};
                channel2_name = edge_list{i,2};
                num_levels
                
                if(strcmp(channel1_name, channel_name))
                    
                    new_partition = obj.partition_y(channel1_name, channel2_name, num_levels, 'transpose_edge', 'no_plot');
                else
                    
                    new_partition = obj.partition_y(channel1_name, channel2_name, num_levels, 'no_plot');
                end
                
                opt_partition = [opt_partition new_partition];
                
            end
            
            
            opt_partition = sort(opt_partition, 'ascend');
            
            %only one partition under 1 allowed
            
            low_partition_indices = find(opt_partition<1);
            
            opt_partition(low_partition_indices(1:end-1))  = [];
            
            
            
            current_index =1;
            
            
            while(current_index <= length(opt_partition))
                current_index
                remaining_vector = opt_partition(current_index+1:end);
                current_val = opt_partition(current_index);
                close_indices = find(abs(remaining_vector-current_val)<merge_threshold);
                if(length(close_indices)==0)
                    current_index = current_index+1;
                    continue;
                end
                
                close_indices = close_indices+current_index;
                average_val = mean([opt_partition(current_index) opt_partition(close_indices)]);
                opt_partition(current_index) = average_val;
                opt_partition(close_indices) = [];
                close_indices = [];
                
            end
            
        end
        
        
        
        function [opt_partition, noise] = partition_y(obj, channel1_name, channel2_name, num_levels, noise_threshold, varargin)
            
            
            max_resolution = 50;
            num_slices = 8;
            num_partitions = num_levels - 1;
            transpose_edge = 0;
            no_plot = 0;
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'max_resolution'))
                    max_resolution = varargin{i+1};
                end
                
                
                if(strcmp(varargin{i}, 'num_slices'))
                    num_slices = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'transpose_edge'))
                    transpose_edge = 1;
                end
                
                if(strcmp(varargin{i}, 'no_plot'))
                    no_plot = 1;
                end
            end
            
            
            [minx, miny, maxx, maxy] = find_data_cutoffs(obj, channel1_name, channel2_name, 50, 255);
            
            if(no_plot == 0)
                [points_x, points_y, point_weights, ~, normalized_density, xaxis, yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', .9);
            else
                [points_x, points_y, point_weights, ~, normalized_density, xaxis, yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', .9, 'no_plot');
                
            end
            if(transpose_edge == 1)
                temp = points_x;
                points_x = points_y;
                points_y = temp;
                temp = channel2_name;
                channel2_name = channel1_name;
                channel1_name = temp;
            end
            [opt_partition] = optimal_partition(transpose([points_x; points_y]), channel1_name, channel2_name, num_partitions, max_resolution, num_slices, minx, miny, maxx, maxy);
            
            if(no_plot == 0)
                
                if(transpose_edge == 1)
                    
                    plot_as_vertical_lines(opt_partition, 'w');
                    
                else
                    
                    plot_as_horizontal_lines(opt_partition, 'w');
                    
                end
                
            end
            
            channel1_num = obj.name_channel_map(channel1_name);
            
            channel1_data = obj.data(:,channel1_num);
            
            channel2_num = obj.name_channel_map(channel2_name);
            
            channel2_data = obj.data(:,channel2_num);
            noise = ones(length(channel2_data),1);
            if(transpose_edge==0)
                for i=1:length(channel1_data)
                    
                    xaxis = abs(xaxis - channel1_data(i));
                    yaxis = abs(yaxis - channel2_data(i));
                    [~,x_index] = min(xaxis);
                    [~,y_index] = min(yaxis);
                    
                    if(normalized_density(y_index, x_index)<noise_threshold)
                        noise(i) = 0;
                    end
                end
            end
            
        end
        
        function [obj] = prune_transitive_edges(obj)
            %this is a really limited case of dismissing edges
            
            
            
            
            for i=1:length(obj.DREMI_channels)
                for j=1:length(obj.DREMI_channels)
                    
                    
                    if(obj.DREMI_sig_adjacency(i,j)==0)
                        continue;
                    end
                    
                    %found an edge
                    
                    for k=1:length(obj.DREMI_channels)
                        
                        if(obj.DREMI_sig_adjacency(j,k)==0)
                            continue;
                        end
                        
                        if(obj.DREMI_sig_adjacency(i,k)==0)
                            continue;
                        end
                        
                        %found a potential case of dismissable transitivity
                        
                        if(obj.DREMI_sig_adjacency(i,k)<obj.DREMI_sig_adjacency(j,k))
                            
                            obj.DREMI_sig_adjacency(i,k) = 0;
                            
                        end
                        
                        
                        
                    end
                end
                
                
            end
            
            
            
            
        end
        
        
        
        function [reduced_channels, reduced_sig_adjacency_dremi_values, color_values] = draw_denovo_DREMI_graph(obj)
            
            
            binary_adjacency = (obj.DREMI_sig_adjacency>0);
            colsum = sum(binary_adjacency,1);
            rowsum = sum(binary_adjacency,2);
            index = 0;
            node_map = zeros(1,length(colsum));
            
            nodes_to_keep=[];
            nodes_to_delete=[];
            for i=1:length(colsum)
                if((colsum(i)>0)|(rowsum(i)>0))
                    index = index+1;
                    node_map(i)=index;
                    nodes_to_keep = [nodes_to_keep i];
                else
                    nodes_to_delete = [nodes_to_delete i];
                end
            end
            
            
            
            reduced_adjacency = binary_adjacency(nodes_to_keep, nodes_to_keep);
            num_nodes = length(reduced_adjacency);
            reduced_channels = obj.DREMI_channels(nodes_to_keep);
            reduced_sig_adjacency = obj.DREMI_sig_adjacency(nodes_to_keep, nodes_to_keep);
            
            gObj = biograph(reduced_adjacency,reduced_channels);
            
            
            range = .7;
            
            load 'continuous_BuPu9.mat';
            colormap(continuous_BuPu9);
            all_nodes = 1:num_nodes;
            set(gObj.nodes(all_nodes),'LineColor',[.4 .4 .4]);
            set(gObj.nodes(all_nodes),'Color',[1 1 1]);
            set(gObj.nodes(all_nodes),'Shape','ellipse');
            set(gObj.nodes(all_nodes),'LineWidth',1.1);
            set(gObj.nodes(all_nodes),'fontSize',10);
            set(gObj,'ShowArrows','off');
            %            all_nodes = 1:length(obj.DREMI_channels);
            %            set(gObj.nodes(all_nodes),'LineColor',[.4 .4 .4]);
            %            set(gObj.nodes(all_nodes),'Color',[1 1 1]);
            %            set(gObj.nodes(all_nodes),'Shape','ellipse');
            %            set(gObj.nodes(all_nodes),'LineWidth',1.1);
            %            set(gObj.nodes(all_nodes),'fontSize',10);
            
            
            tot_length_color_value = sum(size(reduced_adjacency, 1), size(reduced_adjacency, 2));
            color_values = zeros(tot_length_color_value, 3);
            count_count = 0;
            
            reduced_sig_adjacency_dremi_values = reduced_sig_adjacency;
            
            
            reduced_sig_adjacency_vec = reduced_sig_adjacency(:);
            reduced_sig_adjacency_vec_non_zero = reduced_sig_adjacency_vec(reduced_sig_adjacency_vec ~= 0);
            linewidth_mat = reduced_sig_adjacency_vec_non_zero;
            min_use = min(reduced_sig_adjacency_vec_non_zero);
            max_use = max(reduced_sig_adjacency_vec_non_zero);
            count = 0;
            for ctct_row = 1:size(reduced_sig_adjacency, 1)
                for ctct_col = 1:size(reduced_sig_adjacency, 2)
                    if reduced_sig_adjacency(ctct_row, ctct_col) ~= 0
                        count = count + 1;
                        element = reduced_sig_adjacency(ctct_row, ctct_col);
                        reduced_sig_adjacency(ctct_row, ctct_col) = ((element - min_use)/(max_use - min_use))*0.6 + 0.3;
                        linewidth_mat(count) = ((element - min_use)/(max_use - min_use))*5 + 1;
                    end
                end
            end
            
            
            %             reduced_sig_adjacency_vec = (reduced_sig_adjacency_vec - min(reduced_sig_adjacency_vec))/(max(reduced_sig_adjacency_vec) - min(reduced_sig_adjacency_vec));
            %             reduced_sig_adjacency_new = reshape(reduced_sig_adjacency_vec, size(reduced_sig_adjacency));
            %             reduced_sig_adjacency = reduced_sig_adjacency_new;
            
            for i = 1:size(reduced_adjacency,1)
                for j = 1:size(reduced_adjacency,2)
                    
                    
                    if(reduced_adjacency(i,j)==0)
                        continue;
                    end
                    count_count  = count_count + 1;
                    
                    edgeID = getedgesbynodeid(gObj,get(gObj.Nodes([i j]),'ID'));
                    
                    [color_value] = get_color_value(reduced_sig_adjacency(i,j), 1.1, continuous_BuPu9);
                    color_values(count_count, :) = color_value;
                    
                    set(edgeID,'LineColor',color_value);
                    set(edgeID,'LineWidth',linewidth_mat(count_count));
                    
                end
            end
            
            dolayout(gObj);
            view(gObj);
            
            
        end
        
        
        
        
        
        function [obj, sig_edges, sig_mis, sig_edge_matrix, channel_names, val] = write_mi_graph(obj, channel_names, mi_matrix, threshold, filename, varargin)
            
            num_edges = length(find(mi_matrix>threshold));
            sig_edges = {num_edges, 2};
            sig_mis = zeros(num_edges,1);
            file = fopen(filename,'w');
            undirected = 0;
            choose_direction = 0;
            sig_edge_matrix = zeros(length(channel_names), length(channel_names));
            
            if(~isempty(varargin))
                if(strcmp(varargin{1},'undirected'))
                    undirected = 1;
                end
                if(strcmp(varargin{1},'choose'))
                    choose_direction = 1;
                end
            end
            
            
            current = 1;
            for j=1:length(channel_names)
                
                for k=j:length(channel_names)
                    %                       if k == 15
                    %                           disp(j)
                    %                       end
                    if j ~= k
                        %                         if(undirected ==1)
                        %                             if(j==k)
                        %                                 break;
                        %                             end
                        %                         end
                        
                        val = mi_matrix(j,k);
                        %                         if j == 5 && k == 15
                        %                             disp(mi_matrix(j, k))
                        %                         end
                        if(undirected == 1)
                            if mi_matrix(j, k) > mi_matrix(k, j)
                                val = mi_matrix(j, k);
                                ind_first = j;
                                ind_second = k;
                            else
                                val = mi_matrix(k, j);
                                ind_first = k;
                                ind_second = j;
                            end
                            
                        end
                        
                        
                        %                          if k == 20 && j == 18
                        %                             disp(mi_matrix(j,k))
                        %
                        %                             disp(mi_matrix(k,j))
                        %
                        %                             disp(val)
                        %                          end
                        %
                        if(val>threshold)
                            
                            if(choose_direction ==1)
                                if(val<mi_matrix(ind_first,ind_second)) %will only enter one direction
                                    continue;
                                end
                            end
                            %channel_names{j}
                            
                            fprintf(file,'%s %s %f\n', channel_names{ind_first}, channel_names{ind_second}, val);
                            
                            sig_edges{current,1} = channel_names{ind_first};
                            sig_edges{current,2} = channel_names{ind_second};
                            sig_mis(current) = val;
                            current = current+1;
                            sig_edge_matrix(ind_first,ind_second) = val;
                        end
                    end
                end
            end
            obj.DREMI_sig_adjacency = sig_edge_matrix;
            obj.DREMI_channels = channel_names;
            obj.DREMI_sig_edges = sig_edges;
            %obj.DREMI_sig_edge_mi = sig_mis;
            
        end
        
        function [molecule_edges] = get_edges_for_y_molecule(obj, channel_name)
            
            edge_indices = [];
            num_edges = length(obj.DREMI_sig_edges);
            for i=1:num_edges
                
                edge = obj.DREMI_sig_edges(i,:);
                if(strcmp(edge{2}, channel_name))
                    edge_indices = [edge_indices i];
                end
                
            end
            
            molecule_edges = obj.DREMI_sig_edges(edge_indices,:);
        end
        
        
        function [auc, eval_points_x, eval_points_y] = compute_edge_auc(obj, channel1_name, channel2_name, varargin)
            
            
            
            eval_points_x = [];
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'eval_points'))
                    
                    
                    eval_points_x = varargin{i+1};
                    
                    
                end
                
                
            end
            
            
            
            
            [points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
            
            if(length(eval_points_x)==0)
                
                eval_points_x = points_x;
                
            end
            
            eval_points_y = interp1(points_x,points_y,eval_points_x, 'linear', 'extrap');
            
            %area under curve
            auc = sum(eval_points_y .* abs(eval_points_x(2) - eval_points_x(1)));
            
            %normalize the area under the curve.
            %auc
            xrange = max(eval_points_x)-min(eval_points_x);
            yrange = max(eval_points_y)-min(eval_points_y);
            auc = auc/(xrange*yrange);
            
            
        end
        
        function [obj, activity_matrix] = pairwise_auc_compute(obj, channel_names, varargin)
            
            activity_matrix = zeros(length(channel_names), length(channel_names));
            ctrl_specified = 0;
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'ctrl_data'))
                    
                    ctrl_specified = 1;
                    ctrl_data = varargin{i+1};
                    
                    
                end
                
            end
            
            
            for i=1:length(channel_names)
                for j=1:length(channel_names)
                    
                    if i==j
                        continue;
                    end
                    
                    
                    if(ctrl_specified==1)
                        
                        %sprintf('computing activity matrices, control specified')
                        
                        [~, points_x1, points_y1] = obj.compute_edge_auc(channel_names{i},channel_names{j});
                        [~, points_x2, points_y2] = ctrl_data.compute_edge_auc(channel_names{i},channel_names{j});
                        
                        if(max(points_x1) < max(points_x2))
                            points_x = points_x1;
                        else
                            points_x = points_x2;
                        end
                        %i dont think union or intersection works.
                        
                        [auc1, points_x, points_y1] = obj.compute_edge_auc(channel_names{i},channel_names{j}, 'eval_points', points_x);
                        [auc2, points_x, points_y2] = ctrl_data.compute_edge_auc(channel_names{i},channel_names{j}, 'eval_points', points_x);
                        
                        points_y = abs(points_y1-points_y2);
                        %                         auc = sum(points_y);
                        
                        auc = sum(points_y .* abs(points_x(2) - points_x(1)));
                        xrange = max(points_x)-min(points_x);
                        yrange = max(points_y)-min(points_y);
                        auc = auc/(xrange*yrange);
                        dremi1 = obj.compute_dremi(channel_names{i}, channel_names{j}, 0.8);
                        dremi2 = ctrl_data.compute_dremi(channel_names{i}, channel_names{j}, 0.8);
                        dremi = max(dremi1, dremi2);
                        
                        if(dremi>.15)
                            activity_matrix(i,j) = auc;
                        end
                        
                    else
                        
                        activity_matrix(i,j) = obj.compute_edge_auc(channel_names{i},channel_names{j});
                        
                        
                        
                    end
                    
                end
            end
            
            obj.activity_matrix = activity_matrix;
            sprintf('set object activity matrix to activity matrix')
            
            
            
        end
        
        function [causal_graph, causal_score, runtimes ] = pairwise_causality(obj, channel_names, method_to_use, method_options, varargin)
            
            causal_graph = zeros(length(channel_names), length(channel_names));
            causal_score = zeros(length(channel_names), length(channel_names));
            runtimes = [];
            
            for i=1:length(channel_names)
                for j=1:length(channel_names)
                    
                    if(j>=i)
                        break;
                    end
                    
                    %dremi_score = max(obj.DREMI_sig_adjacency(i,j), obj.DREMI_sig_adjacency(j,i));
                    dremi_score = max(obj.DREMI_adjacency(i,j), obj.DREMI_adjacency(j,i));
                    
                    if(dremi_score <= 0.1)
                        
                        
                        continue;
                    end
                    
                    channel1_name = channel_names{i}
                    channel2_name = channel_names{j}
                    %cdata_subsample = obj.subsample_data(1000);
                    tStart = tic;
                    [C_xy] = compute_causal_scores(obj, channel1_name, channel2_name, method_to_use, method_options);
                    tElapsed = toc(tStart)
                    runtimes = [runtimes tElapsed];
                    
                    if(C_xy.decision==1)
                        causal_graph(i,j) = 1;
                        causal_score(i,j) = C_xy.conf;
                    else
                        causal_graph(j,i) = 1;
                        causal_score(j,i) = C_xy.conf;
                    end
                    
                    
                end
            end
            
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'plot'))
                    
                    
                    colormap(linspecer(256))
                    CLIM = [-1 1];
                    imagesc(causal_graph);
                    set(gca,'ytick',1:length(channel_names));
                    set(gca,'yticklabel',channel_names);
                    xticklabel_rotate([1:length(channel_names)],45,channel_names);
                    colorbar
                    
                end
                
                
                
            end
            
            
        end
        
        function [obj, mi_matrix ] = pairwise_mi_compute(obj, channel_names, noise_threshold, varargin)
            
            
            ctrl_specified = -1;
            ctrl_data = [];
            
            mi_matrix = zeros(length(channel_names), length(channel_names));
            
            skiplist = [];
            
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'ctrl_data'))
                    
                    ctrl_specified = 1;
                    ctrl_data = varargin{i+1};
                    
                    
                end
                if(strcmp(varargin{i}, 'skip_list'))
                    
                    
                    skiplist = varargin{i+1};
                    
                    
                end
                
                
                
            end
            
            
            
            for i=1:length(channel_names)
                for j=1:length(channel_names)
                    if(i==j)
                        continue;
                    end
                    channel1_name = channel_names{i}
                    channel2_name = channel_names{j}
                    
                    skip=0;
                    for k=1:size(skiplist,1)
                        
                        if(strcmp(skiplist(k,1), channel1_name))
                            if(strcmp(skiplist(k,2),channel2_name))
                                
                                skip = 1;
                                break;
                            end
                        end
                    end
                    
                    if(skip == 1)
                        continue;
                        %skip this computation
                    end
                    
                    
                    if(ctrl_specified<0)
                        [mi_matrix(i,j),~] = obj.compute_dremi(channel1_name, channel2_name, noise_threshold);
                    else
                        [minx1, miny1, maxx1, maxy1] = find_data_cutoffs(obj, channel1_name, channel2_name, 25, 255);
                        [minx2, miny2, maxx2, maxy2] = find_data_cutoffs(ctrl_data, channel1_name, channel2_name, 25, 255);
                        maxy = max(maxy1,maxy2);
                        [mi_matrix(i,j),~] = obj.compute_dremi(channel1_name, channel2_name, noise_threshold, 'maxy', maxy);
                    end
                    
                    
                    
                    
                    
                end
            end
            
            obj.DREMI_adjacency = mi_matrix;
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'plot'))
                    
                    
                    colormap(linspecer(256))
                    CLIM = [0 1];
                    imagesc(mi_matrix);
                    set(gca,'ytick',1:length(channel_names));
                    set(gca,'yticklabel',channel_names);
                    xticklabel_rotate([1:length(channel_names)],45,channel_names);
                    colorbar
                    
                end
                
                
                
                
            end
            
            
        end
        
        
        
        
        function [obj] = subsample_data(obj, number_sampled)
            
            data_length = size(obj.data,1);
            new_sample = randsample(data_length,number_sampled);
            obj.data = obj.data(new_sample,:);
            
            
        end
        
        
        
        function [dremi, pvalue, samples_x, samples_y] = compute_dremi(obj, channel1_name, channel2_name, noise_threshold , varargin)
            
            
            compute_pvalue = 0;
            set_maxy  = 0;
            num_permutations = 0;
            max_yval = 0;
            num_slices = 8;
            non_kde_style = 0;
            min_yval = 0;
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i}, 'compute_pvalue'))
                    compute_pvalue = 1;
                    num_permutations = varargin{i+1};
                end
                if(strcmp(varargin{i}, 'maxy'))
                    set_maxy = 1;
                    
                    maxy_val = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'num_slices'))
                    num_slices = varargin{i+1};
                end
                if(strcmp(varargin{i}, 'no_kde'))
                    non_kde_style = 1;
                end
                
                
            end
            
            if(non_kde_style==0)
                
                if(set_maxy==0);
                    [points_x, points_y, point_weights, ~, normalized_density, ~,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', noise_threshold,'no_plot');
                else
                    
                    [points_x, points_y, point_weights, ~, normalized_density, ~,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'non_averaged_pts', noise_threshold,'no_plot', 'MinMaxY', 0, maxy_val);
                end
                
                total_slice_samples = sum(point_weights) * 1000;
                samples_x=[];
                samples_y=[];
                
                for i = 1:length(points_x)
                    
                    num_repeats = floor(point_weights(i) * 1000);
                    new_samples_x = ones(num_repeats,1).*points_x(i);
                    new_samples_y = ones(num_repeats,1).*points_y(i);
                    
                    samples_x = [samples_x; new_samples_x];
                    samples_y = [samples_y; new_samples_y];
                end
                
                
                
                data = [samples_x samples_y];
                
                %             [points_x, points_y, point_weights, ~, normalized_density, xaxis,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
                %             data = [transpose(points_x) transpose(points_y)];
                
                
                %[dremi, entropy_y, cond_entropy_y] = dremi_resampled(data, 8, 8);
                minx = min(data(:,1));
                miny = min(data(:,2));
                maxx = max(data(:,1));
                maxy = max(data(:,2));
                
            else
                
                [minx, miny, maxx, maxy] = find_data_cutoffs(obj, channel1_name, channel2_name, 50, 255);
                
                
            end
            
            
            
            if(set_maxy==1)
                maxy = maxy_val;
            end
            
            if(non_kde_style == 0)
                
                dremi = delta_entropyreweight_rawdata(data, minx, miny, maxx, maxy, num_slices,num_slices);
                
            else
                
                dremi = delta_entropyreweight(obj, channel1_name, channel2_name, num_slices, num_slices, minx, miny, maxx, maxy);
                
            end
            
            pvalue = 0;
            if(num_permutations == 0)
                compute_pvalue = 0;
            end
            
            if(compute_pvalue ==1)
                s= sprintf('computing pvalue!\n');
                dremi_random = zeros(1,num_permutations);
                total_perms_performed = num_permutations;
                num_greater = 0;
                for i=1:num_permutations
                    
                    if(non_kde_style==0)
                        
                        dremi_random(i) = delta_entropyreweight_rawdata(data, minx, miny, maxx, maxy, num_slices, num_slices, 'permute_y');
                        
                        
                    else
                        
                        dremi_random(i) = delta_entropyreweight(obj, channel1_name, channel2_name, num_slices, num_slices, minx, miny, maxx, maxy);
                        
                    end
                    
                    
                    random_dremi = dremi_random(i);
                    if(dremi_random(i)>dremi)
                        s=sprintf('dremi random greater than dremi');
                        
                        num_greater = num_greater+1;
                        if(num_greater==10)
                            total_perms_performed = i;
                            break;
                        end
                        
                    end
                    
                end
                
                
                
                above_dremi = find(dremi_random(1:total_perms_performed)>dremi);
                if(total_perms_performed<num_permutations)
                    pvalue = 10/total_perms_performed;
                else
                    pvalue = (length(above_dremi)+1)/(num_permutations+1);
                    
                end
            end
            
            
            
            
            
        end
        
        
        
        
        
        function [responding_fractions] = pairwise_visualize_responding_fraction(obj, channel1_name, channel2_name, num_slices, threshold, varargin)
            
            maxy_val = -1;
            limits = [];
            
            for i=1:length(varargin)
                
                
                if(strcmp(varargin{i}, 'maxy'))
                    set_maxy = 1;
                    limits = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'Limits'))
                    
                    limits = varargin{i+1};
                end
            end
            
            if(length(limits)==0)
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.8, 'no_plot');
            else
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.8,'Limits', limits, 'no_plot');
            end
            responding_fractions = [];
            
            minx = min(x);
            
            maxx = max(x);
            slice_width = (maxx-minx)/num_slices;
            slice_lines = [minx];
            
            for i=1:num_slices
                
                next_val = slice_lines(end) + slice_width;
                slice_lines = [slice_lines next_val];
                
            end
            
            slices_lines = [slice_lines maxx];
            
            for i=1:length(slice_lines)-1
                
                slice_end = slice_lines(i+1)
                slice_start = slice_lines(i);
                
                slice_indices = find((x>=slice_start)&(x<=slice_end));
                
                slice_data = y(slice_indices);
                
                responding_indices = find(slice_data>threshold);
                
                slice_fraction = length(responding_indices)/length(slice_data);
                responding_fractions = [responding_fractions slice_fraction];
                
            end
            plot_as_horizontal_lines([threshold],'w');
            plot_as_vertical_lines(slice_lines(2:num_slices),'w');
            maxresponse = max(responding_fractions);
            %figure;
            bar(responding_fractions,.3)
            
            set(gca,'YLim',[0 1]);
            set(gca,'XTick',[]);
            set(gca,'YTick',[]);
            
            
        end
        
        
        function [F_fitted, MSE, x, y_fit] = double_sigmoid_fit_edge(obj, channel1_name, channel2_name, varargin)
            
            [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.95, 'no_plot');
            %            [x,y] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
            %             channel1 = obj.name_channel_map(channel1_name);
            %                 channel2 = obj.name_channel_map(channel2_name);
            %                 x = transpose(obj.data(:,channel1));
            %                 y =  transpose(obj.data(:,channel2));
            hold on;
            
            for i=1:length(varargin)
                if(strcmp(varargin{i},'cutoff_data'))
                    
                    x_cutoff = varargin{i+1};
                    index = find(x>x_cutoff,1,'first');
                    x(index:end) = [];
                    y(index:end) = [];
                    
                end
            end
            
            lower_intercept = 0;
            mid_intercept = 1;
            upper_intercept = 2;
            
            fix_intercepts = 0;
            for i=1:length(varargin)
                if(strcmp(varargin{i},'fix_intercepts'))
                    
                    fix_intercepts = 1;
                    
                    lower_intercept = varargin{i+1};
                    mid_intercept = varargin{i+2};
                    upper_intercept = varargin{i+3};
                    %slope = varargin{4};
                end
            end
            
            if(fix_intercepts==1)
                
                f = @(p,x) (((lower_intercept - mid_intercept) ./ (1 + (x/p(1)).^p(2)))+mid_intercept) + (((mid_intercept - upper_intercept) ./ (1 + (x/p(3)).^p(4)))+upper_intercept);
                [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[1 1 1 1]);
                
            else
                
                
                f = @(p,x) (((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)))+p(4))+(((p(4) - p(7)) ./ (1 + (x/p(6)).^p(5)))+p(7));
                %p1 = lower, p4 = mid, p7 = upper p3 = inflection1 p2 =
                %slope 1 p6 = infection p5 = slope 2
                
                
                [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[lower_intercept 1 1 mid_intercept 1 1 upper_intercept]);
                
                
                
                
            end
            
            disp(['F = ',num2str(F_fitted)])
            % Plot the data and fit
            %figure(1)
            y_fit = f(F_fitted,x);
            plot(x,y,'*',x,f(F_fitted,x),'r');
            %plot(x,y_fit,'w','LineWidth',2.0);
            
            
        end
        
        
        function [x, y] = conditional_mean_edge(obj, channel1_name, channel2_name, varargin)
            
            Limits = [];
            for i=1:length(varargin)
                if(strcmp(varargin{i},'Limits'))
                    Limits = varargin{i+1};
                end
            end
            
            if(length(Limits)==0)
                
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name);
                
            else
                
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'Limits', Limits);
                
            end
            
            plot(x,y,'w','LineWidth',3.0);
            
        end
        
        function [noise_DREMI] = compute_noise_dremi(obj, channel1_name, channel2_name, val)
            
            %subtract the sigmoid fit from the data. Then compute DREMI on
            %the remainder.
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            dataX = obj.data(:,channel1);
            dataY = obj.data(:,channel2);
            
            
            %             figure
            %              [~,~, points_x_unsorted,points_y_unsorted] = sigmoid_fit_edge(obj, channel1_name, channel2_name,'nearest','extrap');
            %              [points_x,sorted_indices]=sort(points_x_unsorted);
            %              points_y = points_y_unsorted(sorted_indices);
            %              points_x = transpose(points_x);
            %              points_y = transpose(points_y);
            %              [points_x, unique_indices] = unique(points_x);
            %              points_y = points_y(unique_indices);
            %             [points_x, points_y]=obj.pairwise_visualize(channel1_name, channel2_name);
            %              noise_DREMI =0;
            %              Y_causal = interp1(points_x,points_y,transpose(dataX));
            %
            %              plot(points_x, points_y,'w','LineWidth',3.0);
            %
            %               dataY_shifted = dataY-transpose(Y_causal);
            % %              sizedataY = size(dataY)
            %               minyval = min(dataY_shifted);
            %              if(minyval<0)
            %                  dataY_shifted = dataY_shifted+abs(minyval);
            %              end
            %              cdata = obj;
            %
            %              cdata.data(:,channel1) = dataX;
            %              cdata.data(:,channel2) = dataY_shifted;
            %              noise_DREMI = cdata.compute_dremi(channel1_name, channel2_name,0.5)
            %              figure
            %              [points_x2, points_y2] = cdata.pairwise_visualize(channel1_name,channel2_name);
            
            %              plot(points_x2, points_y2, 'w','LineWidth', 3.0);
            [shifted_normalized_density, xaxis, yaxis] = mean_shifted_density(obj, channel1_name, channel2_name);
            
            [noise_DREMI, samples_x, samples_y] = compute_dremi_densitymatrix(shifted_normalized_density,xaxis, yaxis, val);
            
            
        end
        
        
        function [F_fitted, MSE, x,y_fit, Kd ] = sigmoid_fit_edge(obj, channel1_name, channel2_name, varargin)
            
            
            Limits = [];
            for i=1:length(varargin)
                if(strcmp(varargin{i},'Limits'))
                    Limits = varargin{i+1};
                end
            end
            
            if(length(Limits)==0)
                
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.95);
                %plot(x,y)
                %[x,y] = obj.pairwise_visualize(channel1_name, channel2_name);
                %                  plot(x,y)
                %                  channel1 = obj.name_channel_map(channel1_name);
                %                  channel2 = obj.name_channel_map(channel2_name);
                %                  x = transpose(obj.data(:,channel1));
                %                  y =  transpose(obj.data(:,channel2));
            else
                
                [x,y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts', .9,  'Limits', Limits);
                
            end
            
            
            
            for i=1:length(varargin)
                if(strcmp(varargin{i},'cutoff_data'))
                    
                    x_cutoff = varargin{i+1};
                    index = find(x>x_cutoff,1,'first');
                    x(index:end) = [];
                    y(index:end) = [];
                    
                end
            end
            
            
            
            %f = @(F,x) y = F(1)./(1+exp(-F(2)*(x-F(3))));
            %f = @(p,x) ((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)).^p(5))+p(4);
            
            %p1 = lower asymptote
            %p2 = hill slope
            %p3 = disassociatin constant
            %p4 = upper asymptote
            ysmooth = zeros(1,length(x));
            ysmooth(1) = (y(1)+y(2))/2;
            for i = 2:length(x)-1
                
                ysmooth(i) = (y(i-1)+y(i)+y(i+1))/3;
                
            end
            ysmooth(length(x)) = (y(length(x))+y(length(x)-1))/2;
            
            smooth_diff = diff(ysmooth);
            min_diff = min(smooth_diff);
            std_diff = std(smooth_diff)*2;
            
            index_first = find(smooth_diff>min_diff+std_diff,1, 'first');
            
            index_last = find(fliplr(smooth_diff)>min_diff+std_diff,1,'first');
            index_last = length(x) - index_last;
            lower_intercept = mean(y(1:index_first))
            upper_intercept = mean(y(index_last:length(x)))
            x_index_first = x(index_first);
            x_index_last = x(index_last);
            
            R = corrcoef(x,y);
            
            if(R<0)
                temp = upper_intercept;
                upper_intercept = lower_intercept;
                lower_intercept = temp;
            end
            
            
            fix_intercepts = 0;
            for i=1:length(varargin)
                if(strcmp(varargin{i},'fix_intercepts'))
                    fix_intercepts = 1;
                    lower_intercept = varargin{i+1};
                    upper_intercept = varargin{i+2};
                    %slope = varargin{4};
                end
            end
            
            if(fix_intercepts==1)
                
                f = @(p,x) ((lower_intercept - upper_intercept) ./ (1 + (x/p(1)).^p(2)))+upper_intercept;
                [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[1 1]);
                %f = @(p,x) ((lower_intercept - upper_intercept) ./ (1 + (x/p(1)).^slope))+upper_intercept;
                %[F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[1]);
                Kd = F_fitted(1)^F_fitted(2);
                Hc = F_fitted(2);
                inflection = F_fitted(1);
                hill_slope = F_fitted(2);
            else
                
                %f = @(p,x) ((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)))+p(4);
                %f = @(p,x) ((p(1) - p(4)) ./ (1 + (x.*p(3)).^p(2)))+p(4);
                f = @(p,x) ((p(1) - p(4)) ./ (1 + (x/p(3)).^p(2)))+p(4);
                %f = @(p,x) ((x.^p(2) ./ (p(3)^p(2) + x.^p(2))).* p(1))+p(4);
                [F_fitted,R,J,COVB,MSE] = nlinfit(x,y,f,[lower_intercept 1 1 upper_intercept]);
                Kd = F_fitted(3) ^ F_fitted(2);
                Hc = F_fitted(2);
                inflection = F_fitted(3);
                hill_slope = F_fitted(2);
            end
            %slope = F_fitted(2)
            %Ka = F_fitted(3)
            %Kd = F_fitted(3)^F_fitted(2)
            % Display fitted coefficients
            %disp(['F = ',num2str(F_fitted)])
            % Plot the data and fit
            %figure(1)
            y_fit = f(F_fitted,x);
            %plot(x,y,'*',x,y_fit,'r');
            
            
            disp('done');
            
            plot(x,y_fit,'w','LineWidth',3.0);
            %             xlabel_string = sprintf('inflection: %.2f hill slope: %.2f ', inflection, hill_slope);
            %             xlabel(xlabel_string);
            %area(x,y_fit,'FaceColor',[0.65,0.65,0.65])
            %set(gca,'XLim',[Limits(1) Limits(3)]);
            %set(gca,'YLim',[Limits(2) Limits(4)]);
            
            %legend('data','fit')
            
            
        end
        
        function [point_density] = compute_point_density(obj, normalized_density, xaxis, yaxis, point_xval, point_yval)
            
            minx = min(xaxis);
            maxx = max(xaxis);
            
            if(point_xval<minx)
                point_density=0;
                return;
            end
            if(point_xval>maxx)
                point_density=0;
                return;
            end
            
            miny = min(yaxis);
            maxy = max(yaxis);
            
            if(point_yval<miny)
                point_density=0;
                return;
            end
            if(point_yval>maxy)
                point_density=0;
                return;
            end
            
            xdistance_vector = abs(xaxis-point_xval);
            [~,closex_index] = min(xdistance_vector);
            
            ydistance_vector = abs(yaxis-point_yval);
            [~,closey_index] = min(ydistance_vector);
            
            point_density = normalized_density(closey_index,closex_index);
        end
        
        function [C_xy] = compute_causal_scores(obj, channel1_name, channel2_name, method_to_use, method_options)
            
            
            %          [data_x, data_y, point_weights] = obj.pairwise_visualize(channel1_name,channel2_name,'no_plot','non_averaged_pts',0.8);
            %           total_slice_samples = sum(point_weights) * 1000;
            %                   samples_x=[];
            %                   samples_y=[];
            %
            %                   for i = 1:length(data_x)
            %
            %                       num_repeats = floor(point_weights(i) * 1000);
            %                       new_samples_x = ones(num_repeats,1).*data_x(i);
            %                       new_samples_y = ones(num_repeats,1).*data_y(i);
            %
            %                       samples_x = [samples_x; new_samples_x];
            %                       samples_y = [samples_y; new_samples_y];
            %                   end
            %
            %
            %
            %                  data = [samples_x samples_y];
            %
            %           pointsx = data(:,1);
            %           pointsy = data(:,2);
            %[~, ~, ~, ~, normalized_density, xaxis,yaxis] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            
            %          pointsx = obj.data(:,channel1);
            %          pointsy = obj.data(:,channel2);
            %          throw = [];
            %          for i=1:length(pointsx)
            %             point_density = obj.compute_point_density(normalized_density, xaxis, yaxis, pointsx(i), pointsy(i));
            %             if(point_density<.2)
            %                 throw = [throw i];
            %             end
            %          end
            %
            %          points_thrown = length(throw)
            %          pointsx(throw)=[];
            %          pointsy(throw)=[];
            %
            %           [~, ~, ~, ~, normalized_density, xaxis,yaxis] = obj.pairwise_visualize(channel2_name, channel1_name,'no_plot');
            %
            %
            %          throw = [];
            %          for i=1:length(pointsx)
            %             point_density = obj.compute_point_density(normalized_density, xaxis, yaxis, pointsy(i), pointsx(i));
            %             if(point_density<.2)
            %                 throw = [throw i];
            %             end
            %          end
            %
            %          points_thrown = length(throw)
            %          pointsx(throw)=[];
            %          pointsy(throw)=[];
            
            
            
            datax = obj.data(:,channel1);
            datay = obj.data(:,channel2);
            
            %Calling Joris code now
            
            C_xy = 0;
            
            if(strcmp(method_to_use,'IGCI'))
                
                %refMeasure 1 = uniform, 2 = Gaussian
                methodpars = struct('refMeasure', 2,'estimator',method_options,'entest','KSG');
                %C_xy = cep_igci(pointsx, pointsy, methodpars);
                
                C_xy = cep_igci(datax, datay, methodpars)
            end
            
            if(strcmp(method_to_use,'ANM'))
                %          methodpars       Struct containing the parameters for the method
                %            .nrperm          number of permutations to use with HSIC test,
                %                             0 means gamma approximation (default: 0)
                %            .gaussianize     whether to gaussianize the input before fitting the GP (default: 0)
                %            .FITC            if nonzero, uses this amount of points for FITC approximation (using linear grid)
                %            .splitdata       whether to split the data into a training set (for regression) and a test set (for independence testing)
                %            .evaluation      model selection criterion: 'pHSIC' (default), 'HSIC', 'entropy', 'FN', 'Gauss', or 'MML'
                %            .bandwidths      bandwidths for HSIC kernel (default: [0,0])
                %            .meanf           GP pars mean function (default: 'meanAffine')
                %            .minimize        GP pars minimization function: either 'minimize' (default) or 'minimize_lbfgsb'
                %            .entest          Entropy estimation method to use (see entropy.m for options)
                methodpars = struct('nrperm',0,'gaussianize',0,'FITC',0,'splitdata',0,'evaluation', method_options,'bandwidths',[0,0],'meanf','meanAffine', 'minimize','minimize','entest','KSG');
                C_xy = cep_anm(datax,datay, methodpars);
            end
            
            if(strcmp(method_to_use,'DREMI_Residual'))
                
                dremi1 = obj.compute_dremi(channel1_name, channel2_name,.8);
                dremi2 = obj.compute_dremi(channel2_name, channel1_name,.8);
                %actually use DREMI score under certain conditions.
                if(dremi1>dremi2)
                    
                    if(dremi2<=.15)
                        
                        if((dremi1-dremi2)>0.5)
                            C_xy.decision =1;
                            C_xy.conf = (dremi1-dremi2)/dremi1;
                            return;
                        end
                    end
                end
                
                if(dremi2>dremi1)
                    
                    if(dremi1<=.15)
                        
                        if((dremi2-dremi1)>0.5)
                            C_xy.decision =-1;
                            C_xy.conf = (dremi2-dremi1)/dremi2;
                            return;
                        end
                    end
                end
                
                noise_dremi1 =  obj.compute_noise_dremi(channel1_name, channel2_name, method_options);
                noise_dremi2 = obj.compute_noise_dremi(channel2_name, channel1_name, method_options);
                %remember decision is backwards
                if(noise_dremi1>noise_dremi2)
                    C_xy.decision = -1;
                    C_xy.conf = (noise_dremi1-noise_dremi2)/noise_dremi1;
                else
                    C_xy.decision = 1;
                    C_xy.conf = (noise_dremi2-noise_dremi1)/noise_dremi2;
                end
            end
            
            if(strcmp(method_to_use,'DREMI'))
                
                dremi1 =  obj.compute_dremi(channel1_name, channel2_name,0.8);
                dremi2 =  obj.compute_dremi(channel2_name, channel1_name,0.8);
                if(dremi1>dremi2)
                    C_xy.decision = 1;
                    C_xy.conf = dremi1;
                else
                    C_xy.decision = -1;
                    C_xy.conf = dremi2
                end
            end
            
            
            
        end
        
        function [dremi_diff_matrix] = compute_max_DREMI_diff(obj, window_channel_name, marker_channels, window_size)
            
            
            dremi_diff_matrix = zeros(length(marker_channels),length(marker_channels));
            
            for i=1:length(marker_channels)
                i
                for j=1:length(marker_channels)
                    j
                    if (i==j)
                        continue;
                    end
                    
                    DREMI_values = obj.compute_windowed_DREMI(window_channel_name,marker_channels{i}, marker_channels{j}, window_size);
                    
                    dremi_diff_matrix(i,j) = max(DREMI_values)-min(DREMI_values);
                end
            end
            
            
        end
        
        function [DREMI_values, Traj_time, Vqq, meshX2d_ex, meshY2d, meshZ2d, meshX] = compute_windowed_DREMI_interpolate(obj, window_channel, channel1_name, channel2_name, steps, prob_threshold, varargin)
            num_partitions = 8;
            
            call_at = 0.1; set_min = 1; num_bins_val = 128;
            if ~isempty(varargin)
                for i = 1:length(varargin)
                    if strcmp(varargin{i}, 'call_at')
                        call_at = varargin{i+1};
                    end
                    
                    if strcmp(varargin{i}, 'set_min')
                        set_min = 1;
                    end
                    
                    if strcmp(varargin{i}, 'num_bins')
                        num_bins_val = varargin{i+1};
                    end
                end
            end
            
            
            channel1 = obj.name_channel_map(window_channel);
            channel2 = obj.name_channel_map(channel1_name);
            channel3 = obj.name_channel_map(channel2_name);
            
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            
%             ind_for_now = find(edge_data(:, 1) > 0.1 & edge_data(:, 1) < 0.9);
%             edge_data = edge_data(ind_for_now, :);
            
            
            [x_chn1, y_chn1] = obj.pairwise_visualize_new(channel1_name, channel2_name, 'no_plot');
            %[x_chn1, y_chn1] = obj.pairwise_visualize(channel1_name, channel2_name, 'no_plot');
            %max(x_chn1)
            
            if set_min == 1
                %[ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128, 'maxy', max(x_chn1), 'set_min');%, 'maxz', max(y_chn1));
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional_min(edge_data,num_bins_val,  'set_min', varargin{:});%, 'maxz', max(y_chn1));
            else
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional_min(edge_data,num_bins_val, 'maxy', max(x_chn1));%, 'maxz', max(y_chn1));
            end
            
            
            % [~,  meshX, meshY, meshZ, estimated_points] =  obj.threeD_edge_visualize(window_channel, channel1_name, channel2_name, 0.1);
            
            
            
            %[ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128);%, 'maxy', max(x_chn1));
            
            %
            maxX = max(max(max(meshX)));
            minX = min(min(min(meshX)));
            xvals = meshX(1,:,1);
            currentX = (xvals(1)+xvals(2))/2;
            DREMI_values_raw = [];
            Traj_times_raw = [];
            step_size = (maxX-minX)/(steps);
            
            
            [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
            [size1,size2] = size(meshY2d);
            
            while(currentX<maxX)
                
                meshX2d = ones(size1,size2).*currentX;
                
                
                Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d);
                
                if currentX < call_at
                    Vqq = Vq;
                    meshX2d_ex = meshX2d;
                end
                
                %                  figure
                %                   slice(meshX, meshY, meshZ, estimated_points, currentX, [], []);
                %                  figure
                %                   imagesc(Vq)
                %                   set(gca,'YDir','normal');
                
                
                Traj_times_raw = [Traj_times_raw currentX];
                DREMI_value = obj.compute_twoD_dremi(channel1_name, channel2_name, prob_threshold, num_partitions, 'estimated_pts', Vq, meshY2d, meshZ2d);
                %[DREMI_value] = compute_dremi_matrix(Vq, meshY2d, meshZ2d, .80);
                currentX = currentX+step_size;
                if DREMI_value < 0
                    DREMI_value = 0;
                end
                DREMI_values_raw = [DREMI_values_raw DREMI_value];
                
                %                  titlestring = sprintf('DREMI: %.2f',DREMI_value);
                %                  title(titlestring);
            end
            
            %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
            
            
             DREMI_values = DREMI_values_raw;
             Traj_time = Traj_times_raw;
            
%             for i=1:(length(DREMI_values_raw)-4)
%                 
%                 DREMI_values(i) = mean(DREMI_values_raw(i:i+4));
%                 Traj_time(i) = mean(Traj_times_raw(i:i+4));
%             end
%             




            % plot(Traj_time, DREMI_values);
            % DREMI_values = DREMI_values_raw;
            % Traj_time = Traj_times_raw;
        end
        
        
        %         function [DREMI_values, Traj_time] = compute_windowed_DREMI_interpolate_threeD(obj, window_channel, channel1_name, channel2_name, steps, prob_threshold, num_partitions)
        %
        %            channel1 = obj.name_channel_map(window_channel);
        %            channel2 = obj.name_channel_map(channel1_name);
        %            channel3 = obj.name_channel_map(channel2_name);
        %
        %            edge_data = obj.data(:,[channel1 channel2 channel3]);
        %            [x_chn1, ~] = obj.pairwise_visualize_new(channel1_name, channel2_name, 'no_plot');
        %            [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128, 'maxy', max(x_chn1));
        %
        %             %[estimated_points, meshX, meshY, meshZ] =  obj.threeD_edge_visualize(window_channel, channel1_name, channel2_name);
        %             maxX = max(max(max(meshX)));
        %             minX = min(min(min(meshX)));
        %             xvals = meshX(1,:,1);
        %             currentX = (xvals(1)+xvals(2))/2;
        %             DREMI_values_raw = [];
        %             Traj_times_raw = [];
        %             step_size = (maxX-minX)/steps;
        %
        %
        %             [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
        %              [size1,size2] = size(meshY2d);
        %
        %             while(currentX<maxX)
        %
        %
        %
        %                  meshX2d = ones(size1,size2).*currentX;
        %
        %
        %                  Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d);
        %
        %
        % %                  figure
        % %                   slice(meshX, meshY, meshZ, estimated_points, currentX, [], []);
        % %                  figure
        % %                   imagesc(Vq)
        % %                   set(gca,'YDir','normal');
        %
        %
        %                  Traj_times_raw = [Traj_times_raw currentX];
        %                  DREMI_value = obj.compute_threeD_dremi(channel1_name, channel2_name, prob_threshold, num_partitions, 'estimated_pts', Vq, meshY2d, meshZ2d);
        %                  %[DREMI_value] = compute_dremi_matrix(Vq, meshY2d, meshZ2d, .80);
        %                  currentX = currentX+step_size;
        %                  DREMI_values_raw = [DREMI_values_raw DREMI_value];
        %
        % %                  titlestring = sprintf('DREMI: %.2f',DREMI_value);
        % %                  title(titlestring);
        %             end
        %
        %             %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
        %
        %             for i=1:(length(DREMI_values_raw)-4)
        %
        %                DREMI_values(i) = mean(DREMI_values_raw(i:i+4));
        %                Traj_time(i) = median(Traj_times_raw(i:i+4));
        %            end
        %
        %           % plot(Traj_time, DREMI_values);
        %          % DREMI_values = DREMI_values_raw;
        %          % Traj_time = Traj_times_raw;
        %         end
        
        
        function [dremi_raw, traj_raw] = compute_windowed_threeD_DREMI(obj, window_channel, channel1_name, channel2_name, prob_threshold, num_partitions)
            
            channel1 = obj.name_channel_map(window_channel);
            channel2 = obj.name_channel_map(channel1_name);
            channel3 = obj.name_channel_map(channel2_name);
            
            
            
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            
            [~, ind] = sort(edge_data(:, 1));
            
            size(ind)
            %            shift_by = 500;
            %            count = 3000;
            start = 1;
            dremi_raw = [];
            traj_raw = [];
            for i = 1:124
                
                edge_data_new = edge_data(ind(start:start+3000), :);
                fca_writefcs('xyz.fcs', edge_data_new, {window_channel, channel1_name, channel2_name}, {window_channel, channel1_name, channel2_name});
                re = cytof_data_new('xyz.fcs');
                dremi = re.compute_threeD_dremi(window_channel, channel1_name, channel2_name, prob_threshold, num_partitions);
                dremi_raw = [dremi_raw, dremi];
                traj_raw = [traj_raw, mean(edge_data(ind(start:start+3000), 1))];
                start = start + 74;
                
            end
            
            
            %for i = 1:128
            
            %            [x_chn1, ~] = obj.pairwise_visualize_new(channel1_name, channel2_name, 'no_plot');
            %            [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128, 'maxy', max(x_chn1));
            %
            %             %[estimated_points, meshX, meshY, meshZ] =  obj.threeD_edge_visualize(window_channel, channel1_name, channel2_name);
            %             maxX = max(max(max(meshX)));
            %             minX = min(min(min(meshX)));
            %             xvals = meshX(1,:,1);
            %             currentX = (xvals(1)+xvals(2))/2;
            %             DREMI_values_raw = [];
            %             Traj_times_raw = [];
            %             step_size = (maxX-minX)/steps;
            %
            %
            %             [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
            %              [size1,size2] = size(meshY2d);
            %
            %             while(currentX<maxX)
            %
            %
            %
            %                  meshX2d = ones(size1,size2).*currentX;
            %
            %
            %                  Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d);
            %
            %
            % %                  figure
            % %                   slice(meshX, meshY, meshZ, estimated_points, currentX, [], []);
            % %                  figure
            % %                   imagesc(Vq)
            % %                   set(gca,'YDir','normal');
            %
            %
            %                  Traj_times_raw = [Traj_times_raw currentX];
            %                  DREMI_value = obj.compute_twoD_dremi(channel1_name, channel2_name, prob_threshold, num_partitions, 'estimated_pts', Vq, meshY2d, meshZ2d);
            %                  %[DREMI_value] = compute_dremi_matrix(Vq, meshY2d, meshZ2d, .80);
            %                  currentX = currentX+step_size;
            %                  DREMI_values_raw = [DREMI_values_raw DREMI_value];
            %
            % %                  titlestring = sprintf('DREMI: %.2f',DREMI_value);
            % %                  title(titlestring);
            %             end
            %
            %             %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
            %
            %             for i=1:(length(DREMI_values_raw)-4)
            %
            %                DREMI_values(i) = mean(DREMI_values_raw(i:i+4));
            %                Traj_time(i) = median(Traj_times_raw(i:i+4));
            %            end
            
            % plot(Traj_time, DREMI_values);
            % DREMI_values = DREMI_values_raw;
            % Traj_time = Traj_times_raw;
        end
        
        function [parent_dremis, child_trajectory] = compute_windowed_twoD_DREMI(obj, window_channel, names, num_dremis)
            
            channel1 = obj.get_data(window_channel);
            
            [sorted_window, index_sort] = sort(channel1);
            
            window_grid = linspace(min(channel1), max(channel1), num_dremis+2);
            
            %edge_data = obj.data(:,[channel1 channel2 channel3]);
            dremi_raw = zeros(128, 1 );
            
            
            parent_dremis = zeros((length(names) * length(names)), num_dremis+1);
            child_trajectory = zeros((length(names) * length(names)), num_dremis+1);
            %[~, ind] = sort(edge_data(:, 1));
            for i = 2:length(window_grid)-1
                grid_left = sum(sorted_window < window_grid(i));%, 1, 'last');
                grid_right = sum(sorted_window > window_grid(i));
                num = 900;
                tot_num = 1800;
                if grid_left < num && grid_right > num
                    ind_left = find(sorted_window < window_grid(i), 1, 'last');
                    left_half = obj.data(index_sort(1:ind_left), :);
                    right_half = obj.data(index_sort(ind_left+1:tot_num), :);
                    full_data = [left_half; right_half];
                elseif grid_left > num && grid_right < num
                    ind_right = find(sorted_window > window_grid(i), 1, 'first');
                    right_half = obj.data(index_sort(ind_right:end), :);
                    number_of_cells_needed = tot_num - grid_right;
                    left_half = obj.data(index_sort(ind_right-number_of_cells_needed:ind_right-1), :);
                    full_data = [left_half; right_half];
                elseif grid_left > num && grid_right > num
                    %disp('yes');
                    ind_center = find(sorted_window > window_grid(i), 1, 'first');
                    ind_center = ind_center - 1;
                    left_half = obj.data(index_sort(ind_center-(num-1):ind_center), :);
                    %size(left_half, 1)
                    right_half = obj.data(index_sort(ind_center+1:ind_center+num), :);
                    %size(right_half, 1)
                    full_data = [left_half; right_half];
                end
                fca_writefcs('xyz.fcs', full_data, obj.get_all_channel_names(), obj.get_all_channel_names());
                re = cytof_data_new('xyz.fcs');
                %figure();
                %re.plot_scatter('ecadherin', 'vimentin');
                %count = count +1
                %dremi_raw(i-1) = size(full_data, 1);
                %disp(dremi_raw(i-1))
                count = 0;
                for ct = 1:length(names)
                    for ct2 = 1:length(names)
                        count = count + 1
                        if i == 2
                            parent_dremis(count, 1) = ct;
                            child_trajectory(count, 1) = ct2;
                        end
                        if ct ~= ct2
                            parent_dremis(count, i) = re.compute_dremi(names{ct}, names{ct2}, 0.8);
                            child_trajectory(count, i) = window_grid(i);
                        else
                            parent_dremis(count, i) = 0;
                            child_trajectory(count, i) = 0;
                        end
                    end
                end
            end
            
            
            %
            %            start = 1;
            %            dremi_raw = [];
            %            traj_raw = [];
            %            for i = 1:124
            %
            %                 edge_data_new = edge_data(ind(start:start+3000), :);
            %
            %
            %                 dremi = re.compute_threeD_dremi(window_channel, channel1_name, channel2_name, prob_threshold, num_partitions);
            %                 dremi_raw = [dremi_raw, dremi];
            %                 traj_raw = [traj_raw, mean(edge_data(ind(start:start+3000), 1))];
            %                 start = start + 74;
            %
            %            end
            
            
            %for i = 1:128
            
            %            [x_chn1, ~] = obj.pairwise_visualize_new(channel1_name, channel2_name, 'no_plot');
            %            [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128, 'maxy', max(x_chn1));
            %
            %             %[estimated_points, meshX, meshY, meshZ] =  obj.threeD_edge_visualize(window_channel, channel1_name, channel2_name);
            %             maxX = max(max(max(meshX)));
            %             minX = min(min(min(meshX)));
            %             xvals = meshX(1,:,1);
            %             currentX = (xvals(1)+xvals(2))/2;
            %             DREMI_values_raw = [];
            %             Traj_times_raw = [];
            %             step_size = (maxX-minX)/steps;
            %
            %
            %             [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
            %              [size1,size2] = size(meshY2d);
            %
            %             while(currentX<maxX)
            %
            %
            %
            %                  meshX2d = ones(size1,size2).*currentX;
            %
            %
            %                  Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d);
            %
            %
            % %                  figure
            % %                   slice(meshX, meshY, meshZ, estimated_points, currentX, [], []);
            % %                  figure
            % %                   imagesc(Vq)
            % %                   set(gca,'YDir','normal');
            %
            %
            %                  Traj_times_raw = [Traj_times_raw currentX];
            %                  DREMI_value = obj.compute_twoD_dremi(channel1_name, channel2_name, prob_threshold, num_partitions, 'estimated_pts', Vq, meshY2d, meshZ2d);
            %                  %[DREMI_value] = compute_dremi_matrix(Vq, meshY2d, meshZ2d, .80);
            %                  currentX = currentX+step_size;
            %                  DREMI_values_raw = [DREMI_values_raw DREMI_value];
            %
            % %                  titlestring = sprintf('DREMI: %.2f',DREMI_value);
            % %                  title(titlestring);
            %             end
            %
            %             %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
            %
            %             for i=1:(length(DREMI_values_raw)-4)
            %
            %                DREMI_values(i) = mean(DREMI_values_raw(i:i+4));
            %                Traj_time(i) = median(Traj_times_raw(i:i+4));
            %            end
            
            % plot(Traj_time, DREMI_values);
            % DREMI_values = DREMI_values_raw;
            % Traj_time = Traj_times_raw;
        end
        
        
        
        
        
        function [parent_dremis, child_trajectory] = compute_sweeping_window_twoD_DREMI(obj, window_channel, names, num_dremis)
            
            channel1 = obj.get_data(window_channel);
            
            parent_dremis = zeros((length(names) * length(names)), num_dremis+1);
            child_trajectory = zeros((length(names) * length(names)), num_dremis+1);
            
            [~, index_sort] = sort(channel1);
            window_size = 700;
            
            increment = floor((size(obj.data, 1) - window_size)/(num_dremis + 1));
            
            for i = 1:num_dremis
                i
                if i ~= num_dremis
                    ind_use = index_sort(1+(increment*(i-1)):window_size+(increment*(i-1)));
                    data_use = obj.data(ind_use, :);
                else
                    ind_use = index_sort(window_size+(increment*(i-1))+1:end);
                    data_use = obj.data(ind_use, :);
                    size(data_use, 1)
                end
                fca_writefcs('xyz.fcs', data_use, obj.get_all_channel_names(), obj.get_all_channel_names());
                re = cytof_data_new('xyz.fcs');
                size(re.data)
                
                count = 0;
                for ct = 1:length(names)
                    for ct2 = 1:length(names)
                        count = count + 1;
                        if i == 1
                            parent_dremis(count, 1) = ct;
                            child_trajectory(count, 1) = ct2;
                        end
                        if ct ~= ct2
                            parent_dremis(count, i+1) = re.compute_dremi(names{ct}, names{ct2}, 0.8);
                            child_trajectory(count, i+1) = mean(re.get_data(window_channel));
                        else
                            parent_dremis(count, i+1) = 0;
                            child_trajectory(count, i+1) = 0;
                        end
                    end
                end
            end            
        end
                    
        
        
        
        function [DREMI_value] = compute_windowed_DREMI_slices(obj, window_channel, channel1_name, channel2_name, currentX, prob_threshold, num_partitions)
            
            channel1 = obj.name_channel_map(window_channel);
            channel2 = obj.name_channel_map(channel1_name);
            channel3 = obj.name_channel_map(channel2_name);
            
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            [x_chn1, ~] = obj.pairwise_visualize_new(channel1_name, channel2_name, 'no_plot');
            [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128, 'maxy', max(x_chn1));
            
            %[estimated_points, meshX, meshY, meshZ] =  obj.threeD_edge_visualize(window_channel, channel1_name, channel2_name);
            % maxX = max(max(max(meshX)));
            % minX = min(min(min(meshX)));
            %xvals = meshX(1,:,1);
            %currentX = (xvals(1)+xvals(2))/2;
            %DREMI_values_raw = [];
            %Traj_times_raw = [];
            %step_size = (maxX-minX)/steps;
            
            
            [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
            [size1,size2] = size(meshY2d);
            
            %while(currentX<maxX)
            
            
            
            meshX2d = ones(size1,size2).*currentX;
            
            
            Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d);
            
            
            %                  figure
            %                   slice(meshX, meshY, meshZ, estimated_points, currentX, [], []);
            %                  figure
            %                   imagesc(Vq)
            %                   set(gca,'YDir','normal');
            
            
            %  Traj_times_raw = [Traj_times_raw currentX];
            DREMI_value = obj.compute_twoD_dremi(channel1_name, channel2_name, prob_threshold, num_partitions, 'estimated_pts', Vq, meshY2d, meshZ2d);
            %[DREMI_value] = compute_dremi_matrix(Vq, meshY2d, meshZ2d, .80);
            % currentX = currentX+step_size;
            % DREMI_values_raw = [DREMI_values_raw DREMI_value];
            
            %                  titlestring = sprintf('DREMI: %.2f',DREMI_value);
            %                  title(titlestring);
            %end
            
            %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
            
            % for i=1:(length(DREMI_values_raw)-4)
            
            %    DREMI_values(i) = mean(DREMI_values_raw(i:i+4));
            %    Traj_time(i) = median(Traj_times_raw(i:i+4));
            %end
            
            % plot(Traj_time, DREMI_values);
            % DREMI_values = DREMI_values_raw;
            % Traj_time = Traj_times_raw;
        end
        
        function [points_x, points_y, density, xaxis_to_plot,yaxis_to_plot, smoothed_normalized_density,  point_weights, xaxis,yaxis, aa, normalized_density, num_slices] = pairwise_visualize_new(obj, channel1_name, channel2_name, varargin)
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            X = obj.data(:,channel1);
            Y = obj.data(:,channel2);
            total_cells = size(obj.data,1);
            
            
            num_slices = 256;
            minxval = max(0, min(X));
            minyval = max(0, min(Y));
            cutoff = 50;
            draw_contour = 0;
            show_density = 0;
            draw_plot = 1;
            avg_pts = 1;
            avg_pts_threshold = .9;
            fix_limits = 0;
            maxyval = max(Y);
            maxxval = max(X);
            fixy = 0;
            visual_threshold = 0;
            
            %disp(length(varargin))
            
            for i=1:length(varargin)-1
                
                if(strcmp(varargin{i},'Slices'))
                    num_slices = varargin{i+1};
                    disp(num_slices)
                end
                if(strcmp(varargin{i},'MinMaxY'))
                    
                    minyval = varargin{i+1};
                    maxyval = varargin{i+2};
                    
                    fixy = 1;
                    
                end
                if(strcmp(varargin{i},'MinMaxX'))
                    
                    minxval = varargin{i+1};
                    maxxval = varargin{i+2};
                    fixy = 1;
                    
                end
                if(strcmp(varargin{i},'MinMaxFitX'))
                    
                    minxval = min(X);
                end
                if(strcmp(varargin{i},'MinMaxFitY'))
                    
                    minyval = min(Y);
                end
                
                if(strcmp(varargin{i},'Cutoff'))
                    
                    cutoff = varargin{i+1};
                    
                end
                
                if(strcmp(varargin{i},'Minval'))
                    
                    minxval = varargin{i+1};
                    minyval = minxval;
                end
                if(strcmp(varargin{i},'Limits'))
                    % uncomment the line below
                    %fix_limits = 1;
                    limitvector = varargin{i+1};
                    %minx = limitvector(1);
                    %miny = limitvector(2);
                    %maxx = limitvector(3);
                    %maxy = limitvector(4);
                    
                    % Comment all the lines below
                    minxval = limitvector(1);
                    maxxval = limitvector(3);
                    minyval = limitvector(2);
                    maxyval = limitvector(4);
                    
                    
                end
                if(strcmp(varargin{i},'non_averaged_pts'))
                    
                    avg_pts = 0;
                    avg_pts_threshold = varargin{i+1};
                    
                end
                if(strcmp(varargin{i}, 'visual_threshold'))
                    visual_threshold = varargin{i+1};
                end
                
                
            end
            
            
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i},'draw_contour'))
                    draw_contour = 1;
                end
                if(strcmp(varargin{i},'show_density'))
                    show_density = 1;
                end
                if(strcmp(varargin{i},'no_plot'))
                    draw_plot = 0;
                end
                
            end
            
            
            if (fix_limits == 0)
                Xmin = minxval;
                Xmax = maxxval;
                Ymin = minyval;
                Ymax = maxyval;
            else
                Xmin = minx;
                Xmax = maxx;
                Ymin = miny;
                Ymax = maxy;
            end
            
            %                  R = Xmax - Xmin;
            %                  dx = R/(num_slices - 1);
            %                  xmesh1 = Xmin + [0:dx:R];
            %                  [init_X, ind_hist] = histc(X, xmesh1);
            %                  z_init_X = (init_X - mean(init_X))/std(init_X);
            %                  ind_high = find(z_init_X > 6);
            %                  ind_to_remove = [];
            %                  if ~isempty(ind_high)
            %                     hist_sorted = sort(init_X, 'descend');
            %                     for i = 1:length(ind_high)
            %                         to_add = hist_sorted(length(ind_high) + 1);
            %                         ind_orig = find(ind_hist == ind_high(i));
            %                         ind_to_remove = [ind_to_remove; randsample(ind_orig, length(ind_orig) - to_add)];
            %                     end
            %                     X(ind_to_remove) = [];
            %                     Y(ind_to_remove) = [];
            %                  end
            %                 %sum(X == 0)
            %
            
            %Xmin = 0;
            %Xmax = max(X);
            R = Xmax - Xmin;
            
            dx = R/(num_slices - 1);
            xmesh1 = Xmin + [0:dx:R];
            [init_X, ind_hist] = histc(X, xmesh1);
            z_init_X = (init_X - mean(init_X))/std(init_X);
            ind_high = find(z_init_X > 3);
            ind_to_remove = [];
            if ~isempty(ind_high)
                hist_sorted = sort(init_X, 'descend');
                for i = 1:length(ind_high)
                    to_add = ceil(mean(init_X) + mean(init_X)/5); %hist_sorted(length(ind_high) + 1);
                    ind_orig = find(ind_hist == ind_high(i));
                    ind_to_remove = [ind_to_remove; randsample(ind_orig, length(ind_orig) - to_add)];
                end
                X(ind_to_remove) = [];
                Y(ind_to_remove) = [];
            end
            
            
            
            %Xmin = 0;
            %Xmax = max(X);
            %Ymin = 0;
            %Ymax = max(Y);
            Ry = Ymax - Ymin;
            dy = Ry/(num_slices - 1);
            ymesh1 = Ymin + [0:dy:Ry];
            [init_Y, ind_hist_y] = histc(Y, ymesh1);
            z_init_Y = (init_Y - mean(init_Y))/std(init_Y);
            ind_high_y = find(z_init_Y > 3);
            ind_to_remove_y = [];
            if ~isempty(ind_high_y)
                hist_sorted_y = sort(init_Y, 'descend');
                for i = 1:length(ind_high_y)
                    %to_add_y = hist_sorted_y(length(ind_high_y) + 1);
                    to_add_y = ceil(mean(init_Y) + mean(init_Y)/5);
                    ind_orig_y = find(ind_hist_y == ind_high_y(i));
                    %ind_to_remove_y = [ind_to_remove_y; randsample(ind_orig_y, ceil(mean(init_X) + 1.5*std(init_X)))];
                    ind_to_remove_y = [ind_to_remove_y; randsample(ind_orig_y, length(ind_orig_y) - to_add_y)];
                end
                X(ind_to_remove_y) = [];
                Y(ind_to_remove_y) = [];
            end
            
            
            
            if(fix_limits == 0)
                
                [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[Xmin Ymin],[Xmax Ymax]);
                
            else
                
                [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[minx miny],[maxx maxy]);
                
            end
            
            %rounding down small values
            
            for i=1:num_slices
                for j=1:num_slices
                    
                    if(density(i,j)<0.00001)
                        density(i,j) = 0;
                    end
                end
            end
            
            
            %weeding out sparse ends
            
            
            
            row=0;
            
            total_cells_so_far = 0;
            while(total_cells_so_far<cutoff)
                
                
                row = row+1;
                total_cells_so_far = length(find(Y>Grid_Y(num_slices-row,1)));
                
                
            end
            
            maxy = Grid_Y(num_slices-row,1);
            
            
            total_cells_so_far = 0;
            start_row = 0;
            while(total_cells_so_far<cutoff)
                
                
                
                start_row = start_row + 1;
                total_cells_so_far = length(find(Y<Grid_Y(start_row,1)));
                
            end
            
            
            miny = Grid_Y(start_row,1);
            
            %row = 0;
            %start_row = 1;
            
            total_cells_so_far = 0;
            col=0;
            while(total_cells_so_far<cutoff)
                
                
                
                col = col+1;
                total_cells_so_far = length(find(X>Grid_X(1,num_slices-col)));
                
            end
            
            maxx = Grid_X(1,num_slices-col);
            
            
            total_cells_so_far = 0;
            start_col=0;
            while(total_cells_so_far<cutoff)
                
                
                
                start_col = start_col+1;
                total_cells_so_far = length(find(X<Grid_X(1,start_col)));
                
            end
            
            minx = Grid_X(1,start_col);
            
            
            
            if(fix_limits == 1)
                
                start_row = 1;
                start_col = 1;
                row = 0;
                col = 0;
            end
            
            if(fixy==1)
                
                start_row = 1;
                row = 0;
                
            end
            density = density(start_row:num_slices-row,start_col:num_slices-col);
            num_cols = size(density,2);
            num_rows = size(density,1);
            xaxis = Grid_X(1,start_col:num_slices-col);
            yaxis = Grid_Y(start_row:num_slices-row,1);
            
            normalized_density = zeros(num_rows,num_cols);
            prob_normalized_density = zeros(num_rows,num_cols);
            %normalized by column for plotting the data
            for i=1:num_cols
                
                %normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                normalized_density(:,i) = density(:,i)/max(density(:,i));
                prob_normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                
            end
            
            
            
            %now create the side bars
            
            
            colsum = sum(density,1);
            normalized_colsum = colsum./max(colsum);
            
            rowsum = sum(density,2);
            normalized_rowsum = rowsum./max(rowsum);
            
            
            
            
            %the corner is a fudge
            
            
            %blueval = min(normalized_colsum);
            blueval = 0;
            corner = ones(11,11).*blueval;
            
            %make the top bar
            
            %yaxis_increment = abs(yaxis(2)-yaxis(1,1));
            yaxis_increment = .01;
            yaxis_top_bar = [];
            top_bar = [];
            zero_vector = zeros(1,length(normalized_colsum));
            for i=1:1
                top_bar = [top_bar; zero_vector];
                yaxis_top_bar = [yaxis_top_bar; max(yaxis)+(yaxis_increment*i)];
                
            end
            for i=1:10
                top_bar = [top_bar; normalized_colsum];
                yaxis_top_bar = [yaxis_top_bar; max(yaxis)+(yaxis_increment*i)];
            end
            
            
            %make the side bar
            %xaxis_increment = abs(xaxis(2)-xaxis(1));
            xaxis_increment = .01;
            xaxis_side_bar = [];
            side_bar = [];
            zero_vector = zeros(length(normalized_rowsum),1);
            
            for i=1:1
                side_bar = [side_bar zero_vector];
                xaxis_side_bar = [xaxis_side_bar max(xaxis)+(xaxis_increment*i)];
                
            end
            
            for i=1:10
                side_bar = [side_bar normalized_rowsum];
                xaxis_side_bar = [xaxis_side_bar max(xaxis)+(xaxis_increment*i)];
            end
            
            
            
            %find the trace through the peak regions for the return value
            points_x = [];
            points_y = [];
            point_weights = [];
            if(avg_pts==1)
                for i=1:num_cols
                    
                    
                    
                    max_indices = find(normalized_density(:,i)>= avg_pts_threshold);
                    % points_y = [points_y mean(Grid_Y(max_indices,i))];
                    
                    points_x = [points_x xaxis(i)];
                    
                    %new_point_y = dot(Grid_Y(start_row+max_indices,start_col+i),normalized_density(max_indices,i));
                    points_y = [points_y mean(yaxis(max_indices))];
                    %points_y = [points_y new_point_y];
                    
                end
                point_weights = ones(1,length(points_y));
            else
                
                for i=1:num_cols
                    
                    %instead of referring to the grid maybe just take all the points
                    %in the high density squares ??
                    
                    
                    
                    max_indices = find(normalized_density(:,i)>= avg_pts_threshold);
                    % points_y = [points_y mean(Grid_Y(max_indices,i))];
                    new_points = ones(1,length(max_indices)).*xaxis(i);
                    new_point_weights = transpose(normalized_density(max_indices,i));
                    new_point_weights = new_point_weights ./ (sum(new_point_weights));
                    points_x = [points_x new_points];
                    
                    %points_y(i) = dot(Grid_Y(start_row:255-row,i),orig_normalized_density(:,i));
                    y_indices = max_indices;
                    new_points_y = transpose(yaxis(y_indices));
                    points_y = [points_y new_points_y];
                    point_weights = [point_weights new_point_weights];
                    
                end
                
            end
            
            
            
            smoothed_normalized_density = zeros(num_rows, num_cols);
            
            for i=1:num_rows
                for j = 2:num_cols-1
                    smoothed_normalized_density(i,j) = (normalized_density(i,j-1)+normalized_density(i,j)+normalized_density(i,j+1))/3;
                end
                
            end
            %imagesc(flipud(Grid_X(:,1)), flipud(transpose(Grid_Y(1,:))), normalized_density);
            
            if(visual_threshold>0)
                smoothed_normalized_density = (smoothed_normalized_density>visual_threshold).*smoothed_normalized_density;
                
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % To get the side bars remove % from the next two %-ed lines
            % And, comment out the next line.
            %colormap(linspecer(256));
            %density_filtered = smoothed_normalized_density>.6;
            %smoothed_normalized_density = smoothed_normalized_density.*density_filtered;
            
            
            
            %matrix_to_plot = smoothed_normalized_density;
            %matrix_to_plot = [smoothed_normalized_density side_bar];
            %top_bar = [top_bar corner];
            %matrix_to_plot = [matrix_to_plot; top_bar];
            
            
            
            xaxis_to_plot = [xaxis xaxis_side_bar];
            yaxis_to_plot = [yaxis; yaxis_top_bar];
            
            
            
            
            
            
            %imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%% FOR BARS ON THE SIDE AND ABOVE
            %smoothed_normalized_density = [smoothed_normalized_density side_bar];
            %top_bar = [top_bar corner];
            %smoothed_normalized_density = [smoothed_normalized_density; top_bar];
            density_filtered = smoothed_normalized_density>.6;
            smoothed_normalized_density = smoothed_normalized_density.*density_filtered;
            
            density_filtered = smoothed_normalized_density>.6;
            smoothed_normalized_density = smoothed_normalized_density.*density_filtered;
            
            %%%%%%%%%%%%
            % The following is the original
            %density_filtered = matrix_to_plot>0.6;
            
            
            
            %density_filtered = matrix_to_plot>0.6;
            %density_filtered = matrix_to_plot>0.1;
            %%%%%%%%%%%%
            
            
            %density_filtered = matrix_to_plot;
            %  matrix_to_plot = matrix_to_plot.*density_filtered;
            %aa = matrix_to_plot;
            %
            %
            if(draw_plot)
                colormap(linspecer(256));
                %%%%%%%%%%%%
                j = linspecer(256);
                j(1,:) = [ 1 1 1 ];
                colormap(j);
                %colormap(gray);
                %%%%%%%%%%%%%%%
                %    colormap(linspecer(256));
                
                
                
                
                %imagesc(xaxis_to_plot, yaxis_to_plot, smoothed_normalized_density);
                %figure();
                
                
                
                %imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                imagesc(xaxis_to_plot,yaxis_to_plot, smoothed_normalized_density);
                set(gca,'YDir','normal');
                xlabel(channel1_name, 'FontSize', 12, 'FontWeight', 'bold')
                ylabel(channel2_name, 'FontSize', 12, 'FontWeight', 'bold')
                %set(gca,'XTick',[]);
                %set(gca,'YTick',[]);
                %set(gca, 'XTickLabel','');
                %set(gca, 'YTickLabel','');
                
                
                %         set(gca, 'FontSize',16);
                %         xlabel(channel1_name);
                %         ylabel(channel2_name);
                
                hold
            end
            
            if(draw_contour)
                [bandwidth,rdensity,rGrid_X,rGrid_Y]=kde2d([X Y],num_slices+1,[minx miny],[maxx maxy]);
                
                contour(rGrid_X, rGrid_Y, rdensity, 12);
                
            end
            
            if(show_density)
                
                f = ksdensity(X, points_x);
                plot(points_x,f, 'w', 'LineWidth',1.3);
                
            end
            
            
            
            
            for i=1:length(varargin)-1
                
                if(strcmp(varargin{i},'Title'))
                    
                    
                    title(varargin{i+1});
                    
                    
                end
            end
            
            
            
            
        end
        
        
        
        
        function [DREMI_values, Traj_time] = compute_windowed_DREMI(obj, window_channel_name, channel1_name, channel2_name, window_size)
            
            window_channel = obj.name_channel_map(window_channel_name);
            [num_events, num_channels] = size(obj.data);
            [sorted_vals,sorted_id] = sort(obj.data(:, window_channel));
            i=1;
            DREMI_values_raw=[];
            Traj_time_raw = [];
            channel2 = obj.name_channel_map(channel2_name);
            
            maxy = max(obj.data(:,channel2));
            
            while i< (num_events-window_size)
                
                window_indices = sorted_id(i:i+window_size);
                cdata_window = obj.cluster_gate_events(window_indices);
                size(cdata_window.data)
                DREMI_values_raw(end+1) = cdata_window.compute_dremi(channel1_name,channel2_name,0.80, 'maxy',maxy);
                %DREMI_values_raw(end+1) = cdata_window.compute_twoD_dremi(channel1_name, channel2_name,.6,8);
                
                Traj_time_raw(end+1) = median(sorted_vals(i:i+window_size));
                i=i+50;
            end
            size_dremi_values_raw = size(DREMI_values_raw);
            
            for i=1:(length(DREMI_values_raw)-10)
                
                DREMI_values(i) = mean(DREMI_values_raw(i:i+10));
                Traj_time(i) = median(Traj_time_raw(i:i+10));
            end
            
            plot(Traj_time, DREMI_values);
            
        end
        
        
        function [obj] = normalize_channel(obj, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            channel_data = obj.data(:,channel);
            channel_data = channel_data - min(channel_data);
            channel_data = channel_data ./ max(channel_data);
            obj.data(:,channel) = channel_data;
        end
        
        
        function [Difference_values, Traj_time] = compute_windowed_difference(obj, obj2, window_channel_name, channel_name, window_size)
            
            window_channel = obj.name_channel_map(window_channel_name);
            [num_events, num_channels] = size(obj.data);
            [sorted_vals,sorted_id] = sort(obj.data(:, window_channel));
            i=1;
            Difference_values_raw=[];
            Traj_time_raw = [];
            channel = obj.name_channel_map(channel_name);
            
            
            
            while i< (num_events-window_size)
                
                window_indices = sorted_id(i:i+window_size);
                cdata_window = obj.cluster_gate_events(window_indices);
                Difference_values_raw(end+1) = cdata_window.get_channel_mean(channel_name);
                Traj_time_raw(end+1) = median(sorted_vals(i:i+window_size));
                i=i+10;
            end
            
            window_channel = obj2.name_channel_map(window_channel_name);
            [num_events, num_channels] = size(obj2.data);
            [sorted_vals2,sorted_id2] = sort(obj2.data(:, window_channel));
            i=1;
            Difference_values_raw2=[];
            Traj_time_raw2 = [];
            channel = obj2.name_channel_map(channel_name);
            
            
            while i< (num_events-window_size)
                
                window_indices = sorted_id2(i:i+window_size);
                cdata_window = obj2.cluster_gate_events(window_indices);
                Difference_values_raw2(end+1) = cdata_window.get_channel_mean(channel_name);
                Traj_time_raw2(end+1) = median(sorted_vals2(i:i+window_size));
                i=i+10;
            end
            
            traj_time_max = 0.8;
            traj_time_min = 0.2;
            size_traj_time_raw = size(Traj_time_raw)
            size_Difference_values_raw = size(Difference_values_raw)
            points = linspace(traj_time_min,traj_time_max,200);
            
            
            
            diff_values_interp1 = interp1(Traj_time_raw, Difference_values_raw, points,'spline');
            diff_values_interp2 = interp1(Traj_time_raw2, Difference_values_raw2, points,'spline');
            
            
            
            diff_values = abs(diff_values_interp1-diff_values_interp2);
            
            
            
            for i=1:(length(diff_values)-10)
                
                Difference_values(i) = mean(diff_values(i:i+10));
                Traj_time(i) = median(points(i:i+10));
            end
            
            plot(Traj_time, Difference_values);
            
        end
        
        function scatter_plot_indices(obj, channel1, channel2, ind_use)
            figure();
            xy = obj.get_data({channel1, channel2});
            plot(xy(ind_use, 1), xy(ind_use, 2), '.', 'MarkerSize', 20, 'col', 'blue')
            xlabel(channel1)
            ylabel(channel2)
            set(gca, 'FontSize', 16);
        end
        
        function [DREMI, DREMIXZ, DREMIYZ ] = compute_threeD_dremi_bin_first(obj, channel1_name, channel2_name,channel3_name, prob_threshold, num_partitions,varargin)
            
            num_partitions
            threshold = 50;
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            channel3 = obj.name_channel_map(channel3_name);
            
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            widthx = max(edge_data(:,1))-min(edge_data(:,1));
            widthx = widthx/num_partitions;
            minedgedata = min(edge_data(:,1))
            maxedgedata = max(edge_data(:,1))
            xbin_edges = min(edge_data(:,1)):widthx:max(edge_data(:,1));
            
            widthy = max(edge_data(:,2))-min(edge_data(:,2));
            widthy = widthy/num_partitions;
            
            ybin_edges = min(edge_data(:,2)):widthy:max(edge_data(:,2));
            
            edges = cell(1,2);
            edges{1} = xbin_edges;
            edges{2} = ybin_edges;
            cellhist = hist3(edge_data(:,[1 2]), 'Edges', edges)
            
            %just to get the XI
            [~, XI ] = ksdensity(edge_data(:,3),'npoints', 256);
            
            partition_entropies = zeros(num_partitions,num_partitions);
            partition_density_estimates = cell(num_partitions,num_partitions);
            valid_partitions = 0;
            total_density_estimate = zeros(1,256);
            
            for i=1:num_partitions
                for j=1:num_partitions
                    
                    if(cellhist(i,j) > threshold)
                        
                        
                        valid_partitions = valid_partitions + 1;
                        %binned data only
                        bin_indices = find((edge_data(:,1)>=xbin_edges(i))&(edge_data(:,1)<xbin_edges(i+1))&(edge_data(:,2)>=ybin_edges(j))&(edge_data(:,2)<ybin_edges(j+1)));
                        bin_data = edge_data(bin_indices,3);
                        
                        
                        density_estimate = ksdensity(bin_data,XI);
                        density_estimate = density_estimate / max(density_estimate);
                        density_estimate_filter  = density_estimate>prob_threshold;
                        density_estimate_refined = density_estimate .* density_estimate_filter;
                        
                        partition_density_estimates{i,j} = density_estimate_refined;
                        partition_entropies(i,j) = compute_density_entropy(density_estimate_refined, XI, num_partitions);
                        total_density_estimate = total_density_estimate + density_estimate_refined;
                        
                    else
                        
                        partition_density_estimates{i,j} = zeros(1,length(XI));
                        
                    end
                    
                end
            end
            
            
            total_entropy = compute_density_entropy(total_density_estimate, XI, num_partitions);
            
            avg_partition_entropy = sum(sum(partition_entropies))/valid_partitions;
            DREMI = total_entropy - avg_partition_entropy;
            
            x_partition_entropies = zeros(1,num_partitions);
            valid_x_partitions = 0;
            
            for i=1:num_partitions
                
                x_partition_density_estimate = zeros(1,length(XI));
                valid_bit = 0;
                
                for j=1:num_partitions
                    if(cellhist(i,j)>100)
                        
                        x_partition_density_estimate = x_partition_density_estimate + partition_density_estimates{i,j};
                        valid_bit = 1;
                    end
                end
                
                if(valid_bit>0)
                    x_partition_entropies(i) = compute_density_entropy(x_partition_density_estimate, XI, num_partitions);
                end
                valid_x_partitions = valid_x_partitions + valid_bit;
                
            end
            
            avg_partition_entropyX = sum(x_partition_entropies)/valid_x_partitions;
            DREMIXZ = total_entropy - avg_partition_entropyX;
            
            y_partition_entropies = zeros(1,num_partitions);
            valid_y_partitions = 0;
            
            for j=1:num_partitions
                
                y_partition_density_estimate = zeros(1,length(XI));
                valid_bit = 0;
                
                for i=1:num_partitions
                    if(cellhist(i,j)>100)
                        
                        y_partition_density_estimate = y_partition_density_estimate + partition_density_estimates{i,j};
                        valid_bit = 1;
                    end
                end
                
                if(valid_bit>0)
                    y_partition_entropies(j) = compute_density_entropy(y_partition_density_estimate, XI, num_partitions);
                end
                valid_y_partitions = valid_y_partitions + valid_bit;
                
            end
            
            avg_partition_entropyY = sum(y_partition_entropies)/valid_y_partitions;
            DREMIYZ = total_entropy - avg_partition_entropyY;
            
            
            
        end
        
        
        function [hyperedge_scores] = find_initial_network_edges(obj, marker_channels)
            
            num_channels = length(marker_channels);
            combination_vector = 1:num_channels;
            all_combinations = nchoosek(combination_vector,3);
            hyperedge_scores = cell(size(all_combinations,1)*3, 2)
            next_score = 1;
            for i=1:size(all_combinations,1)
                
                channel1_name = marker_channels{all_combinations(i,1)};
                channel2_name = marker_channels{all_combinations(i,2)};
                channel3_name = marker_channels{all_combinations(i,3)};
                
                
                %[DREMI, DREMIXZ, DREMIYZ ] = obj.compute_threeD_dremi_bin_first(channel1_name,channel2_name,channel3_name, .5 , 16);
                [DREMI, DREMIXZ, DREMIYZ ] = obj.compute_threeD_dremi(channel1_name,channel2_name,channel3_name, .5 , 8);
                channel_array = {channel1_name, channel2_name, channel3_name};
                scores_array = [DREMI, DREMIXZ, DREMIYZ];
                hyperedge_scores{next_score,1} = channel_array;
                hyperedge_scores{next_score,2} = scores_array;
                next_score = next_score+1;
                
                %[DREMI, DREMIXZ, DREMIYZ ] = obj.compute_threeD_dremi_bin_first(channel2_name,channel3_name,channel1_name, .5 , 16);
                [DREMI, DREMIXZ, DREMIYZ ] = obj.compute_threeD_dremi(channel2_name,channel3_name,channel1_name, .5 , 8);
                channel_array = {channel2_name, channel3_name, channel1_name};
                scores_array = [DREMI, DREMIXZ, DREMIYZ];
                hyperedge_scores{next_score,1} = channel_array;
                hyperedge_scores{next_score,2} = scores_array;
                next_score = next_score+1;
                
                %[DREMI, DREMIXZ, DREMIYZ ] = obj.compute_threeD_dremi_bin_first(channel1_name,channel3_name,channel2_name, .5 , 16);
                [DREMI, DREMIXZ, DREMIYZ ] = obj.compute_threeD_dremi(channel1_name,channel3_name,channel2_name, .5 , 8);
                channel_array = {channel1_name, channel3_name, channel2_name};
                scores_array = [DREMI, DREMIXZ, DREMIYZ];
                hyperedge_scores{next_score,1} = channel_array;
                hyperedge_scores{next_score,2} = scores_array;
                next_score = next_score+1;
                
            end
            
            %for now i probably just want to see it right?
            
            
        end
        
        
        
        function [DREMI, DREMIXZ, DREMIYZ ] = compute_threeD_dremi(obj, channel1_name, channel2_name,channel3_name, prob_threshold, varargin)
            num_partitions = 8;
%             maxx = max(obj.get_data(channel1_name));
%             maxy = max(obj.get_data(channel2_name));
%             maxz = max(obj.get_data(channel3_name));
            
            edge_data = obj.get_data({channel1_name, channel2_name, channel3_name});
            [~, y0, ~, ~, ~, x1, y1] = pairwise_visualize_matrix_min(edge_data(:, 1:2), 'no_plot');
            max_vals = [1, max(y1), max(edge_data(:, 3))];
            %max_vals
            maxx = max_vals(1);
            maxy = max_vals(2);
            maxz = max_vals(3);


            if(length(varargin)>0)
                for i=1:length(varargin)
                    if(strcmp(varargin{i},'maxy'))
                        maxy = varargin{i+1};
                    end
                    if(strcmp(varargin{i},'maxx'))
                        maxx = varargin{i+1};
                    end
                    if(strcmp(varargin{i},'maxz'))
                        maxz = varargin{i+1};
                    end
                end
            end
            
            
            [~, ~, ~, ~, estimated_points] =  obj.threeD_edge_visualize( channel1_name, channel2_name, channel3_name,[0.1], 'set_min', varargin{:});%, 'maxx', maxx, 'maxy', maxy, 'maxz',maxz, varargin{:});
            
            estimated_points_joint_flat = sum(estimated_points, 3);
            
            %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
            %so I'd just call this repeatedly and then use a subset of
            %the grid
            
            %instead of doing this with a histogram I should do it
            %with a fraction of the density
            
            %  [numy, numx, numz] = size(estimated_points);
            %zero_vector = zeros(numz,1);
            
            [y_length, x_length, ~] = size(estimated_points);
            
            x_partition_ends = 0:((x_length)/num_partitions):x_length;
            
            y_partition_ends = 0:((y_length)/num_partitions):y_length;
            
            
            
            total_entropy = compute_sample_entropy(estimated_points, num_partitions, prob_threshold);
            
            partition_entropies = zeros(num_partitions, num_partitions);
            valid_partitions = 0;
            
            for i = 1: num_partitions
                for j = 1:num_partitions
                    
                    xpartbegin = x_partition_ends(i)+1;
                    xpartend = x_partition_ends(i+1);
                    
                    ypartbegin = y_partition_ends(j)+1;
                    ypartend = y_partition_ends(j+1);
                    
                    estimated_points_submatrix = estimated_points(ypartbegin:ypartend, xpartbegin:xpartend,:);
                    
                    %cellhist_submatrix = cellhist(ypartbegin:ypartend, xpartbegin:xpartend);
                    
                    [part_entropy, valid_bit] = compute_sample_entropy(estimated_points_submatrix, num_partitions, prob_threshold);
                    valid_partitions = valid_partitions + valid_bit;
                    partition_entropies(i,j) = 0;
                    if(valid_bit>0)
                        partition_weight = sigmoid_4pl(estimated_points_joint_flat(i,j), 1, 0.5, .05);
                        %sigmoid_4pl(x, max_asymptote, slope, shift)
                        %y = max_asymptote./(1+exp(-1*slope.*(x-shift)));
                        
                        partition_entropies(i,j) =  partition_weight*part_entropy;
                    end
                end
            end
            
            
            avg_partition_entropies = sum(sum(partition_entropies))*(1/valid_partitions);
            DREMI = total_entropy - avg_partition_entropies;
            
            partition_entropiesX = zeros(1,num_partitions);
            valid_partitionsX=0;
            for i = 1: num_partitions
                
                
                xpartbegin = x_partition_ends(i)+1;
                xpartend = x_partition_ends(i+1);
                
                
                
                estimated_points_submatrix = estimated_points(:, xpartbegin:xpartend, :);
                
                
                
                [part_entropy, valid_bit] = compute_sample_entropy(estimated_points_submatrix, num_partitions, prob_threshold);
                valid_partitionsX = valid_partitionsX + valid_bit;
                partition_entropiesX(i) = part_entropy;
                
            end
            
            avg_partition_entropiesX = sum(partition_entropiesX)*(1/valid_partitionsX);
            DREMIXZ = total_entropy - avg_partition_entropiesX;
            
            partition_entropiesY = zeros(1,num_partitions);
            valid_partitionsY = 0;
            
            for i = 1: num_partitions
                
                ypartbegin = y_partition_ends(i)+1;
                ypartend = y_partition_ends(i+1);
                
                estimated_points_submatrix = estimated_points(ypartbegin:ypartend,:, :);
                
                
                
                [part_entropy, valid_bit] = compute_sample_entropy(estimated_points_submatrix, num_partitions, prob_threshold);
                valid_partitionsY = valid_partitionsY + valid_bit;
                partition_entropiesY(i) = part_entropy;
                
            end
            
            avg_partition_entropiesY = sum(partition_entropiesY)*(1/valid_partitionsY);
            DREMIYZ = total_entropy - avg_partition_entropiesY;
            
        end
        
        %         function [DREMI, DREMIXZ, DREMIYZ ] = compute_threeD_dremi(obj, channel1_name, channel2_name,channel3_name, prob_threshold, num_partitions,varargin)
        %                 maxx = max(obj.get_data(channel1_name));
        %                 maxy = max(obj.get_data(channel2_name));
        %                 maxz = max(obj.get_data(channel3_name));
        %
        %
        %             if(length(varargin)>0)
        %                 for i=1:length(varargin)
        %                     if(strcmp(varargin{i},'maxy'))
        %                         maxy = varargin{i+1};
        %                     end
        %                     if(strcmp(varargin{i},'maxx'))
        %                         maxx = varargin{i+1};
        %                     end
        %                     if(strcmp(varargin{i},'maxz'))
        %                         maxz = varargin{i+1};
        %                     end
        %                 end
        %             end
        %
        %
        %             [estimated_points, meshX, meshY, meshZ, ~, ~, ~, cellhist, estimated_points_joint] =  obj.threeD_edge_visualize( channel1_name, channel2_name, channel3_name,'maxx', maxx, 'maxy', maxy, 'maxz',maxz);
        %
        %             estimated_points_joint_flat = sum(estimated_points, 3);zeros(size(estimated_points,1),size(estimated_points,2));
        %
        %                  for i = 1:size(estimated_points,1);
        %                     for j = 1:size(estimated_points,2);
        %
        %                         estimated_points_joint_flat(i,j) = sum(estimated_points_grid(i,j,:));
        %
        %                     end
        %                  end
        %
        %                  %Vq = interp3(X,Y,Z,V,Xq,Yq,Zq)
        %                  %so I'd just call this repeatedly and then use a subset of
        %                  %the grid
        %
        %                  %instead of doing this with a histogram I should do it
        %                  %with a fraction of the density
        %
        %
        %
        %                  [numx, numy, numz] = size(estimated_points);
        %                  zero_vector = zeros(numz,1);
        %
        %
        %
        %                  [x_length, y_length, ~] = size(estimated_points);
        %
        %                  x_partition_ends = 0:((x_length)/num_partitions):x_length;
        %
        %                  y_partition_ends = 0:((y_length)/num_partitions):y_length;
        %
        %
        %
        %                  total_entropy = compute_sample_entropy(estimated_points, num_partitions, prob_threshold);
        %
        %                  partition_entropies = zeros(num_partitions, num_partitions);
        %                  valid_partitions = 0;
        %
        %                  for i = 1: num_partitions
        %                     for j = 1:num_partitions
        %
        %                         xpartbegin = x_partition_ends(i)+1;
        %                         xpartend = x_partition_ends(i+1);
        %
        %                         ypartbegin = y_partition_ends(j)+1;
        %                         ypartend = y_partition_ends(j+1);
        %
        %                         estimated_points_submatrix = estimated_points(xpartbegin:xpartend, ypartbegin:ypartend,:);
        %
        %                         cellhist_submatrix = cellhist(xpartbegin:xpartend, ypartbegin:ypartend);
        %
        %                         [part_entropy, valid_bit] = compute_sample_entropy(estimated_points_submatrix, num_partitions, prob_threshold);
        %                         valid_partitions = valid_partitions + valid_bit;
        %                         partition_entropies(i,j) = 0;
        %                         if(valid_bit>0)
        %                             partition_weight = sigmoid_4pl(estimated_points_joint_flat(i,j), 1, 0.5, .05);
        %                             %sigmoid_4pl(x, max_asymptote, slope, shift)
        %                             % y = max_asymptote./(1+exp(-1*slope.*(x-shift)));
        %
        %                           partition_entropies(i,j) =  partition_weight*part_entropy;
        %                         end
        %                       end
        %                  end
        %
        %
        %                  avg_partition_entropies = sum(sum(partition_entropies))*(1/valid_partitions)
        %                  DREMI = total_entropy - avg_partition_entropies;
        %
        %                 partition_entropiesX = zeros(1,num_partitions);
        %                 valid_partitionsX=0;
        %                 for i = 1: num_partitions
        %
        %
        %                         xpartbegin = x_partition_ends(i)+1;
        %                         xpartend = x_partition_ends(i+1);
        %
        %
        %
        %                         estimated_points_submatrix = estimated_points(xpartbegin:xpartend, :,:);
        %
        %
        %
        %                        [part_entropy, valid_bit] = compute_sample_entropy(estimated_points_submatrix, num_partitions, prob_threshold);
        %                        valid_partitionsX = valid_partitionsX + valid_bit;
        %                        partition_entropiesX(i) = part_entropy;
        %
        %                 end
        %
        %                 avg_partition_entropiesX = sum(partition_entropiesX)*(1/valid_partitionsX);
        %                 DREMIXZ = total_entropy - avg_partition_entropiesX;
        %
        %                 partition_entropiesY = zeros(1,num_partitions);
        %                 valid_partitionsY = 0;
        %
        %                 for i = 1: num_partitions
        %
        %                         ypartbegin = y_partition_ends(i)+1;
        %                         ypartend = y_partition_ends(i+1);
        %
        %                         estimated_points_submatrix = estimated_points(:, ypartbegin:ypartend,:);
        %
        %
        %
        %                        [part_entropy, valid_bit] = compute_sample_entropy(estimated_points_submatrix, num_partitions, prob_threshold);
        %                        valid_partitionsY = valid_partitionsY + valid_bit;
        %                        partition_entropiesY(i) = part_entropy;
        %
        %                 end
        %
        %                 avg_partition_entropiesY = sum(partition_entropiesY)*(1/valid_partitionsY);
        %                 DREMIYZ = total_entropy - avg_partition_entropiesY;
        %
        %         end
        %
        %         function [DREMI ] = compute_twoD_dremi(obj, channel1_name, channel2_name, prob_threshold, num_partitions, varargin)
        %
        %                  threshold = 0; %set_minx = 0; set_maxx = 0; set_miny = 0; set_maxy = 0;
        %
        %                  [estimated_points, meshX, meshY, ~, ~, cellhist, ~] =  obj.twoD_edge_visualize( channel1_name, channel2_name,varargin);
        %
        %                  for i=1:length(varargin)
        %
        %                      if(strcmp(varargin{i},'estimated_pts'))
        %
        %                          %set_estimated_points = 1;
        %                          estimated_points = varargin{i+1};
        %                          meshX = varargin{i+2};
        %                          meshY = varargin{i+3};
        %                      end
        %
        % %                      if(strcmp(varargin{i}, 'minx'))
        % %                          set_minx = 1;
        % %                          minX = varargin{i+1};
        % %                      end
        % %
        % %                      if(strcmp(varargin{i}, 'maxx'))
        % %                          set_maxx = 1;
        % %                          maxX = varargin{i+1};
        % %                      end
        % %
        % %                      if(strcmp(varargin{i}, 'miny'))
        % %                          set_miny = 1;
        % %                          minY = varargin{i+1};
        % %                      end
        % %
        % %                      if(strcmp(varargin{i}, 'maxy'))
        % %                          set_maxy = 1;
        % %                          maxY = varargin{i+1};
        % %                      end
        %                  end
        %
        %
        %
        % %                  [xindex] = find(cellhist>=threshold);
        % %                  %zero out the entries where this isn't true.
        %                   [numx, numy] = size(estimated_points);
        % %                  zero_vector = zeros(numy,1);
        % %
        % %                  for i=1:length(xindex)
        % %
        % %                     estimated_points(xindex(i),:) = zero_vector;
        % %
        % %                  end
        %
        %                  [x_length, ~] = size(estimated_points);
        %
        %                  x_partition_ends = 0:((x_length)/num_partitions):x_length;
        %
        %                  total_entropy = compute_sample_entropy_2d(estimated_points, num_partitions, prob_threshold);
        %
        %                  partition_entropies = zeros(num_partitions,1);
        %
        %                  valid_partitions = 0;
        %
        %                  avg_partition_entropies =0;
        %                  for i = 1: num_partitions
        %
        %
        %                         xpartbegin = x_partition_ends(i)+1;
        %                         xpartend = x_partition_ends(i+1);
        %
        %
        %
        %                         estimated_points_submatrix = estimated_points(xpartbegin:xpartend, :);
        %                         %estimated_points_submatrix = estimated_points(:, xpartbegin:xpartend);
        %
        %                         %cellhist_submatrix = cellhist(xpartbegin:xpartend);
        %
        %                        [part_entropy, valid_bit] = compute_sample_entropy_2d(estimated_points_submatrix, num_partitions, prob_threshold);
        %                        valid_bit;
        %                        partition_entropies(i) = part_entropy;
        %                        valid_partitions = valid_partitions+valid_bit;
        %                      %  if(valid_bit>0)
        %                            avg_partition_entropies = avg_partition_entropies+part_entropy;
        %                      %  end
        %                  end
        %
        %                 partition_entropies;
        %                 valid_partitions;
        %                 avg_partition_entropies = avg_partition_entropies/valid_partitions;
        %                 DREMI = total_entropy - avg_partition_entropies;
        %                 DREMI = DREMI/log2(num_partitions);
        %
        %         end
        
        function [DREMI, estimated_points, meshX, meshY] = compute_twoD_dremi(obj, channel1_name, channel2_name, prob_threshold, varargin)
            threshold = 0; set_minx = 0; set_maxx = 0; set_miny = 0; set_maxy = 0; set_estimated_points = 0;
            num_partitions = 8;
            
            
            %[estimated_points, meshX, meshY, ~, ~, cellhist, ~] =  obj.twoD_edge_visualize( channel1_name, channel2_name,varargin);
            
            for i=1:length(varargin)
                
                if(strcmp(varargin{i},'estimated_pts'))
                    
                    set_estimated_points = 1;
                    estimated_points_int = varargin{i+1};
                    meshX = varargin{i+2};
                    meshY = varargin{i+3};
                end
                
                if(strcmp(varargin{i}, 'minx'))
                    set_minx = 1;
                    minX = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'maxx'))
                    set_maxx = 1;
                    maxX = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'miny'))
                    set_miny = 1;
                    minY = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'maxy'))
                    set_maxy = 1;
                    maxY = varargin{i+1};
                end
            end
            
            %set_maxx
            %set_maxy
            
            if set_estimated_points == 1
                estimated_points = estimated_points_int;
            else
%                 if set_maxx == 1 && set_maxy == 0 
%                     [estimated_points, meshX, meshY] =  obj.twoD_edge_visualize( channel1_name, channel2_name, 'maxx', maxX);
%                 elseif set_maxx == 0 && set_maxy == 1
%                     [estimated_points, meshX, meshY] =  obj.twoD_edge_visualize( channel1_name, channel2_name, 'maxy', maxY);
%                 elseif set_maxx == 1 && set_maxy == 1
%                     [estimated_points, meshX, meshY] =  obj.twoD_edge_visualize( channel1_name, channel2_name, 'maxx', maxX, 'maxy', maxY);
%                 else
%                     [estimated_points, meshX, meshY] = obj.twoD_edge_visualize(channel1_name, channel2_name);                    
%                 end

                [estimated_points, meshX, meshY] = obj.twoD_edge_visualize(channel1_name, channel2_name, varargin{:});                    
            end
            
            
            
            %                  [xindex] = find(cellhist>=threshold);
            %                  %zero out the entries where this isn't true.
            [numx, numy] = size(estimated_points);
            %                  zero_vector = zeros(numy,1);
            %
            %                  for i=1:length(xindex)
            %
            %                     estimated_points(xindex(i),:) = zero_vector;
            %
            %                  end
            
            [~, x_length] = size(estimated_points);
            
            x_partition_ends = 0:((x_length)/num_partitions):x_length;
            
            total_entropy = compute_sample_entropy_2d_min(estimated_points, num_partitions, prob_threshold);
            
            partition_entropies = zeros(num_partitions,1);
            
            valid_partitions = 0;
            
            avg_partition_entropies = 0;
            for i = 1: (num_partitions-1)
                
                xpartbegin = x_partition_ends(i)+1;
                xpartend = x_partition_ends(i+1);
                
                estimated_points_submatrix = estimated_points(:, xpartbegin:xpartend);
                
                %estimated_points_submatrix = estimated_points(:, xpartbegin:xpartend);
                
                %cellhist_submatrix = cellhist(xpartbegin:xpartend);
                
                [part_entropy, valid_bit] = compute_sample_entropy_2d_min(estimated_points_submatrix, num_partitions, prob_threshold);
                %part_entropy
                ypart = meshY(:, xpartbegin:xpartend);
                %min(meshY(:)) - min(ypart(:))
                partition_entropies(i) = part_entropy;
                valid_partitions = valid_partitions+valid_bit;
                if(valid_bit>0)
                    avg_partition_entropies = avg_partition_entropies+part_entropy;
                end
            end
            
            partition_entropies;
            valid_partitions;
            avg_partition_entropies = avg_partition_entropies/valid_partitions;
            DREMI = total_entropy - avg_partition_entropies;
            DREMI = DREMI;%/log2(num_partitions);
            range(meshX(:));
            %DREMI = DREMI/(max(meshX(:)) - min(meshX(:)));
            %DREMI = DREMI * (range(meshY(:)));
            
        end
        
        
        function [estimated_points, meshX, meshY, conditional_mean, meshX1d, cellhist, edgesx] = twoD_edge_visualize(obj, channel1_name, channel2_name,varargin)
            set_maxx = 0; set_maxy = 0;
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            
            edge_data = obj.data(:,[channel1 channel2]);
            
            for i=1:length(varargin)
                if(strcmp(varargin{i}, 'maxx'))
                    set_maxx = 1;
                    maxX = varargin{i+1};
                end
                
                if(strcmp(varargin{i}, 'maxy'))
                    set_maxy = 1;
                    maxY = varargin{i+1};
                end
                
                
            end
            
            
%             if set_maxx == 1 && set_maxy == 0
%                 [  estimated_points, meshX, meshY, conditional_mean, meshX1d] = conditional_kde_min(edge_data,256, 'maxx', maxX);          
%             elseif set_maxx == 0 && set_maxy == 1
%                 [ estimated_points, meshX, meshY, conditional_mean, meshX1d] = conditional_kde_min(edge_data,256, 'maxy', maxY);
%             elseif set_maxx == 1 && set_maxy == 1
%                 [ estimated_points, meshX, meshY, conditional_mean, meshX1d] = conditional_kde_min(edge_data,256, 'maxx', maxX, 'maxy', maxY);
%             else
%                 [ estimated_points, meshX, meshY, conditional_mean, meshX1d] = conditional_kde_min(edge_data,256);
%             end
          
            [  estimated_points, meshX, meshY, conditional_mean, meshX1d] = conditional_kde_min_0(edge_data,256, varargin{:});          
            
            
            
            %calculate it with the meshgrid centers
            
            edgesx = meshX1d(1,:);
            widthx = edgesx(2)-edgesx(1);
            edgesx = edgesx - (widthx/2);
            edgesx(end+1) = edgesx(end)+widthx;
            edgesx = transpose(edgesx);
            
            
            
            cellhist = histc(edge_data(:,1),edgesx);
            cellhist = cellhist(1:size(estimated_points,1));
            %take off the last row from this histogram or no need?
            %      figure();
            %            imagesc(edgesx, edgesy, cellhist);%
            %     set(gca,'YDir','normal');
            %           bar(cellhist);
            %visualizing the conditional mean
            
            
            %figure
            %plot(meshX1d, conditional_mean);
            
            
        end
        
        
        function chn_names = get_all_channel_names(obj)
            len = size(obj.data, 2);
            chn_names = cell(1, len);
            for i= 1:len
                chn_names(i) = obj.channel_name_map(i);
            end
        end
        
        function [dremi_values] = compute_3d_interpolated_slice_DREMI(obj, window_channel, channel1_name, channel2_name, prob_threshold, slice_values)
            
            channel1 = obj.name_channel_map(window_channel);
            channel2 = obj.name_channel_map(channel1_name);
            channel3 = obj.name_channel_map(channel2_name);
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional_kde_3d(edge_data,128);
            for i=1:length(slice_values)
                
                [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
                [size1,size2] = size(meshY2d);
                meshX2d = ones(size1,size2).*slice_values(i);
                Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d,'spline');
                dremi_values(i) = obj.compute_twoD_dremi(channel1_name, channel2_name, prob_threshold, 8, 'estimated_pts', Vq, meshY2d, meshZ2d)
            end
            
        end
        
        function [estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = threeD_edge_visualize_new(obj, channel1_name, channel2_name, channel3_name, varargin)
            smooth_yes = 0; set_maxx = 0; set_maxy = 0; set_maxz = 0; display = 0;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'smooth')
                    smooth_yes = 1;
                end
                
                if strcmp(varargin{i}, 'set_maxx')
                    set_maxx = 1;
                    maxx = varargin{i+1};
                end
                
                if strcmp(varargin{i}, 'set_maxy')
                    set_maxy = 1;
                    maxy = varargin{i+1};
                end
                
                if strcmp(varargin{i}, 'set_maxz')
                    set_maxz = 1;
                    maxz = varargin{i+1};
                end
                
                if strcmp(varargin{i}, 'display_plots')
                    display = 1;
                end
            end
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            channel3 = obj.name_channel_map(channel3_name);
            
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            
            if (set_maxx == 1 && set_maxy == 0 && set_maxz == 0)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,80,'maxx', maxx);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128,'maxx', maxx);
                disp('Heat. Eq. Method done!')
            elseif (set_maxx == 0 && set_maxy == 1 && set_maxz == 0)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,80,'maxy', maxy);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128,'maxy', maxy);
                disp('Heat. Eq. Method done!')
            elseif (set_maxx == 0 && set_maxy == 0 && set_maxz == 1)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,80,'maxz', maxz);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128,'maxz', maxz);
                disp('Heat. Eq. Method done!')
            elseif (set_maxx == 1 && set_maxy == 1 && set_maxz == 0)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,80,'maxy', maxy, 'maxx', maxx);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128,'maxx', maxx, 'maxy', maxy);
                disp('Heat. Eq. Method done!')
            elseif (set_maxx == 0 && set_maxy == 1 && set_maxz == 1)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,80,'maxy', maxy, 'maxz', maxz);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128,'maxy', maxy, 'maxz', maxz);
                disp('Heat. Eq. Method done!')
            elseif (set_maxx == 1 && set_maxy == 0 && set_maxz == 1)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,80,'maxx', maxx, 'maxz', maxz);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128,'maxx', maxx, 'maxz', maxz);
                disp('Heat. Eq. Method done!')
            elseif (set_maxx == 0 && set_maxy == 0 && set_maxz == 0)
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,50);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128);
                disp('Heat. Eq. Method done!')
            else
                %[ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,50, 'maxx', maxx, 'maxy', maxy,'maxz', maxz);
                %disp('Comp. Geo. Method done!')
                [ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional(edge_data,128, 'maxx', maxx, 'maxy', maxy,'maxz', maxz);
                disp('Heat. Eq. Method done!')
            end
            
            %Edit it to take off sparse areas. What is a rigorious way of
            %doing this?
            
            %[ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional_kde3d_changed(edge_data,128,varargin);
            
            
            %plot_density_isosurfaces(meshX, meshY, meshZ, estimated_points, {channel1_name, channel2_name, channel3_name});
            
            if smooth_yes == 1
                for i = 1:size(conditional_mean, 1)
                    conditional_mean(i, :) = smooth(conditional_mean(i, :), 20);
                end
                
                
                for j = 1:size(conditional_mean, 2)
                    conditional_mean( :, j) = smooth(conditional_mean(:, j), 20);
                end
            end
            %calculate it with the meshgrid centers
            
            %            edgesx = meshX2d(1,:);
            %            widthx = edgesx(2)-edgesx(1);
            %            edgesx = edgesx - (widthx/2);
            %            edgesx(end+1) = edgesx(end)+widthx;
            %            edgesx = transpose(edgesx);
            %
            %            edgesy = meshY2d(:,1);
            %            widthy = edgesy(2)-edgesy(1);
            %            edgesy = edgesy - (widthy/2);
            %            edgesy(end+1) = edgesy(end)+widthy;
            %
            %            edges = cell(1,2);
            %            edges{1} = edgesx;
            %            edges{2} = edgesy;
            %            cellhist = hist3(edge_data(:,[1 2]),'Edges',edges);
            %            cellhist = cellhist(1:size(estimated_points,1), 1:size(estimated_points,2));
            
            %           imagesc(edgesx, edgesy, cellhist);
            
            
            
            if display == 1
                figure
                %set(gca,'YDir','normal');
                %visualizing the conditional mean
                
                
                %figure
                use_meshx = meshX(:,:,1);
                use_meshy = meshY(:,:,1);
                
                surf(use_meshx, use_meshy, conditional_mean,'FaceColor','red','EdgeColor','none')
                %surf(meshX2d, meshY2d, conditional_mean,'FaceColor','red','EdgeColor','none')
                %surf(conditional_mean, 'FaceColor','red','EdgeColor','none');
                camlight left; lighting phong
                xlabel(channel1_name, 'FontSize', 16);
                ylabel(channel2_name, 'FontSize', 16);
                zlabel(channel3_name, 'FontSize', 16);
                
                zmax = ceil(max(max(conditional_mean)));
                title('Heat Eq. based method', 'FontSize', 16);
                
                figure
                colormap(linspecer(256));
                view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},1,  .1);
                zlim([min(meshZ(:)) zmax]);
                %set(gca, 'edgecolor', 'none');
                shading flat
                
                
                figure
                %subplot(1, 3, 2)
                colormap(linspecer(256));
                view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},1, .4);
                zlim([0 zmax]);
                %set(gca, 'edgecolor', 'none');
                shading flat;
                
                figure
                colormap(linspecer(256));
                
                view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},1, .75);
                zlim([0 zmax]);
                %set(gca, 'edgecolor', 'none');
                shading interp;
            end
            
            
            %             figure
            %            %set(gca,'YDir','normal');
            %            %visualizing the conditional mean
            %
            %
            %            %figure
            %            use_meshx_comp = meshX_comp(:,:,1);
            %            use_meshy_comp = meshY_comp(:,:,1);
            %
            %            surf(use_meshx_comp, use_meshy_comp, conditional_mean_comp,'FaceColor','red','EdgeColor','none')
            %            %surf(meshX2d, meshY2d, conditional_mean,'FaceColor','red','EdgeColor','none')
            %            %surf(conditional_mean, 'FaceColor','red','EdgeColor','none');
            %            camlight left; lighting phong
            %            xlabel(channel1_name, 'FontSize', 16);
            %            ylabel(channel2_name, 'FontSize', 16);
            %            zlabel(channel3_name, 'FontSize', 16);
            %            title('Comp. Geo. based method', 'FontSize', 16);
            %zmax = ceil(max(max(conditional_mean)));
        end
        
        
        
        
        function [conditional_mean, meshX, meshY, meshZ, estimated_points, Movie_ko_lagi, F_fitted, MSE, x,y_fit, Kd, conditional_estimated_points_grid, meshX2d, meshY2d, cellhist, edgesx, edgesy] = threeD_edge_visualize(obj, channel1_name, channel2_name, channel3_name, points_to_slice_at, varargin)
            F_fitted = 0; MSE = 0; x = 0; y_fit = 0;  Kd = 0;
            draw_response = 0; set_min = 1; num_bins_val = 128;
            smooth_yes = 0; show_comp_geo = 0; slice_or_not = 'no'; display_plots = 0; x_or_y = 1; yes_bandwidth = 0;
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'smooth')
                    smooth_yes = 1;
                end
                if strcmp(varargin{i}, 'show_comp_geo')
                    show_comp_geo  = 1;
                end
                if strcmp(varargin{i}, 'slice')
                    slice_or_not = 'slice';                    
                end
                
                if strcmp(varargin{i}, 'display_plots')
                    display_plots = 1;
                end
                
                if strcmp(varargin{i}, 'xy')
                    x_or_y = varargin{i+1};
                end
                if strcmp(varargin{i}, 'bandwidth')
                    yes_bandwidth = 1;
                    bandwidth = varargin{i+1};
                end
                   
                if strcmp(varargin{i}, 'display_response')
                    draw_response = 1;
                    F_fitted = zeros(length(points_to_slice_at), 4);
                end
                
                if strcmp(varargin{i}, 'set_min')
                    set_min = 1;
                end
                
                if strcmp(varargin{i}, 'num_bins')
                    num_bins_val = varargin{i+1};
                end
                
            end
            
            
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            channel3 = obj.name_channel_map(channel3_name);
            
            edge_data = obj.data(:,[channel1 channel2 channel3]);
            if show_comp_geo == 1
                [ estimated_points_comp, meshX_comp, meshY_comp, meshZ_comp, conditional_mean_comp, meshX2d_comp, meshY2d_comp] = conditional_kde_3d(edge_data,50,varargin);
                disp('Comp. Geo. Method done!')
            end
            
            %Edit it to take off sparse areas. What is a rigorious way of
            %doing this?
            
            %[ estimated_points, meshX, meshY, meshZ, conditional_mean, meshX2d, meshY2d] = conditional_kde3d_changed(edge_data,128,varargin);
            if yes_bandwidth == 1
                if set_min == 1
                    [ estimated_points, meshX, meshY, meshZ, conditional_mean, conditional_estimated_points_grid, meshX2d, meshY2d] = conditional_min(edge_data,num_bins_val, 'set_min', varargin{:});
                else
                    [ estimated_points, meshX, meshY, meshZ, conditional_mean, conditional_estimated_points_grid, meshX2d, meshY2d] = conditional(edge_data,128,'bandwidth', bandwidth);
                end
            else
                if set_min ==1                    
                    [ estimated_points, meshX, meshY, meshZ, conditional_mean, conditional_estimated_points_grid, meshX2d, meshY2d] = conditional_min(edge_data,num_bins_val, 'set_min', varargin{:});
                else
                    [ estimated_points, meshX, meshY, meshZ, conditional_mean, conditional_estimated_points_grid, meshX2d, meshY2d] = conditional_min(edge_data,num_bins_val, varargin{:});
                end
            end
            %          [estimated_points_grid, meshX, meshY, meshZ, conditional_mean, conditional_estimated_points_grid
            disp('Heat. Eq. Method done!')
%             min(meshX(:))
%             min(meshY(:))
%             min(meshZ(:))

            %plot_density_isosurfaces(meshX, meshY, meshZ, estimated_points, {channel1_name, channel2_name, channel3_name});
            
%             for i = 1:size(estimated
%             estimated_points(:,  = 
            
            
            
            
            
            
            if smooth_yes == 1
                for i = 1:size(conditional_mean, 1)
                    conditional_mean(i, :) = smooth(conditional_mean(i, :), 20);
                end
                
                
                for j = 1:size(conditional_mean, 2)
                    conditional_mean( :, j) = smooth(conditional_mean(:, j), 20);
                end
                
                %calculate it with the meshgrid centers
                
                %            edgesx = meshX2d(1,:);
                %            widthx = edgesx(2)-edgesx(1);
                %            edgesx = edgesx - (widthx/2);
                %            edgesx(end+1) = edgesx(end)+widthx;
                %            edgesx = transpose(edgesx);
                %
                %            edgesy = meshY2d(:,1);
                %            widthy = edgesy(2)-edgesy(1);
                %            edgesy = edgesy - (widthy/2);
                %            edgesy(end+1) = edgesy(end)+widthy;
                %
                %            edges = cell(1,2);
                %            edges{1} = edgesx;
                %            edges{2} = edgesy;
                %            cellhist = hist3(edge_data(:,[1 2]),'Edges',edges);
                %            cellhist = cellhist(1:size(estimated_points,1), 1:size(estimated_points,2));
                
            end
            
            
            
            %           imagesc(edgesx, edgesy, cellhist);
            
            xmax1 = max(meshX(:));
            ymax1 = max(meshY(:));
            zmax = max(meshZ(:));
            
            if display_plots == 1
                
                %figure
                %set(gca,'YDir','normal');
                %visualizing the conditional mean
                
                use_col = [100	149	237]/255;
                %figure
                use_meshx = meshX(:,:,1);
                use_meshy = meshY(:,:,1);
                
                %surf(use_meshx, use_meshy, conditional_mean,'FaceColor','red','EdgeColor','none')
                temp_edge_data = edge_data(:, 1:2);
                
                %conditional_mean  = conditional_mean + repmat(10, size(conditional_mean));
                [~, density_to_color_by] = kde2d(temp_edge_data, 128, min(temp_edge_data), max(temp_edge_data));
                %surf(use_meshx, use_meshy, conditional_mean, density_to_color_by, 'EdgeColor','none');%, 'FaceColor', 'red')
                surf(use_meshx, use_meshy, conditional_mean, 'EdgeColor','none', 'FaceColor', 'red', 'FaceAlpha', 1)
                %freezeColors();
                colormap(linspecer(256));
                
                camlight right;
                lighting phong
                %lighting flat
                
                
                
                %surf(meshX2d, meshY2d, conditional_mean,'FaceColor','red','EdgeColor','none')
                %surf(conditional_mean, 'FaceColor','red','EdgeColor','none');
                
                if strcmp(channel1_name, 'cycler')
                    use_x_name = 'EMT Time';
                else
                    use_x_name = channel1_name;
                end
                %title('Heat Eq. based method', 'FontSize', 16);
                xlabel('X', 'FontSize', 16);
                ylabel('Y', 'FontSize', 16);
                zlabel('Z', 'FontSize', 16);
                set(gca, 'FontSize', 16);
            end
            
            if strcmp(slice_or_not, 'slice')
                
                %title('Heat Eq. based method', 'FontSize', 16);
                
                clearvars Movie_ko_lagi;
                %figure;
                for slice_at = 1:length(points_to_slice_at)
                    figure;
                    %whitebg([230, 230, 250]/256);
                    j = linspecer(256);
                    j(1:100, :) = 1;
                    %min(meshY(:))
                    %min(meshZ(:))
                    
                    colormap(j)
                    view_2dslice_min(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name}, x_or_y,  points_to_slice_at(slice_at));
                    %hold on
                    %zlim([0 zmax+1]);             
                    %xlabel('EMT-time')
                    %ylabel('pERK')
                    %zlabel('Snail1')
                    %set(gca, 'XTick', [], 'YTick', [], 'ZTick', [], 'LineWidth', 2)
                    box on
                    axes_orig = axis;
                    hold on
                    slice_line_width = 2;
                    sl1 = points_to_slice_at(slice_at);
                    plot3([sl1, sl1], [axes_orig(3), axes_orig(4)], [axes_orig(5), axes_orig(5)], 'col', 'k', 'LineWidth', slice_line_width);
                    plot3([sl1, sl1], [axes_orig(4), axes_orig(4)], [axes_orig(5), axes_orig(6)], 'col', 'k', 'LineWidth', slice_line_width);
                    plot3([sl1, sl1], [axes_orig(4), axes_orig(3)], [axes_orig(6), axes_orig(6)], 'col', 'k', 'LineWidth', slice_line_width);
                    plot3([sl1, sl1], [axes_orig(3), axes_orig(3)], [axes_orig(6), axes_orig(5)], 'col', 'k', 'LineWidth', slice_line_width);

                    
                    if draw_response == 1
                        
                        maxX = max(max(max(meshX)));
                        minX = min(min(min(meshX)));
                        xvals = meshX(1,:,1);
                        currentX = (xvals(1)+xvals(2))/2;
                        steps = length(meshX);
                        step_size = (maxX-minX)/steps;
            
            
                        [meshY2d,meshZ2d]= meshgrid(meshY(:,1,1),meshZ(1,1,:));
                        [size1, size2] = size(meshX(:, :, 1));
            
                        %while(currentX<maxX)                
                         meshX2d = ones(size1,size2).*points_to_slice_at(slice_at);
                         Vq = interp3(meshX, meshY, meshZ, estimated_points, meshX2d, meshY2d, meshZ2d);
                         [F_fitted(slice_at, :), MSE, x,y_fit, Kd ] = sigmoid_fit_edge_standalone(Vq, meshX2d, meshY2d, meshZ2d);
                        
                         

                    end
                    
                    
%                     if x_or_y == 2  
%                         
%                         plot3([0, xmax1, xmax1, 0, 0], [pts, pts, pts, pts, pts], [0, 0, zmax, zmax, 0], 'LineWidth', 2, 'col', 'black')
%                     else
%                         
%                         plot3([pts, pts, pts, pts, pts], [0, ymax1, ymax1, 0, 0], [0, 0, zmax, zmax, 0], 'LineWidth', 2, 'col', 'black')                        
%                     end
                    %set(gca, 'grid', 'off')
                    %set(gca, 'edgecolor', 'none');
                    %set(gca, 'GridLineStyle', 'none')
                    shading flat
                    Movie_ko_lagi(slice_at) = getframe(gcf);
                    %xlabel(channel1_name, 'FontSize', 16);
                    %ylabel(channel2_name, 'FontSize', 16);
                    %zlabel(channel3_name, 'FontSize', 16);
                    %set(gca, 'FontSize', 16);
                    
                    
                end
            end
            
            %                 figure
            %                 %subplot(1, 3, 2)
            %                 colormap(pmkmp(256));
            %                 view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},1, points_to_slice_at(2));
            %                 zlim([0 zmax]);
            %                 %set(gca, 'edgecolor', 'none');
            %                 shading flat;
            %                 xlabel(channel1_name, 'FontSize', 16);
            %                 ylabel(channel2_name, 'FontSize', 16);
            %                 zlabel(channel3_name, 'FontSize', 16);
            %                 set(gca, 'FontSize', 16);
            %
            %                 figure
            %                 colormap(pmkmp(256));
            %                 view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},1, points_to_slice_at(3));
            %                 zlim([0 zmax]);
            %                 %set(gca, 'edgecolor', 'none');
            %                 shading flat;
            %                 xlabel(channel1_name, 'FontSize', 16);
            %                 ylabel(channel2_name, 'FontSize', 16);
            %                 zlabel(channel3_name, 'FontSize', 16);
            %                 set(gca, 'FontSize', 16);
            
            
            
            
            %                 figure
            %                 colormap(pmkmp(256));
            %                 view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},2,  points_to_slice_at(1));
            %                 zlim([0 zmax]);
            %                 %set(gca, 'edgecolor', 'none');
            %                 shading flat;
            %                 xlabel(channel1_name, 'FontSize', 16);
            %                 ylabel(channel2_name, 'FontSize', 16);
            %                 zlabel(channel3_name, 'FontSize', 16);
            %                 set(gca, 'FontSize', 16);
            %
            %
            %
            %
            %                 figure
            %                 %subplot(1, 3, 2)
            %                 colormap(pmkmp(256));
            %                 view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},2, points_to_slice_at(2));
            %                 zlim([0 zmax]);
            %                 %set(gca, 'edgecolor', 'none');
            %                 shading flat;
            %
            %                 figure
            %                 colormap(pmkmp(256));
            %                 view_2dslice(estimated_points, meshX, meshY, meshZ, {channel1_name, channel2_name, channel3_name},2, points_to_slice_at(3));
            %                 zlim([0 zmax]);
            %                 %set(gca, 'edgecolor', 'none');
            %                 shading flat;
            %
            %                 xlabel(channel1_name, 'FontSize', 16);
            %                 ylabel(channel2_name, 'FontSize', 16);
            %                 zlabel(channel3_name, 'FontSize', 16);
            %                 set(gca, 'FontSize', 16);
            %end
            
            
            if show_comp_geo == 1
                figure
                %set(gca,'YDir','normal');
                %visualizing the conditional mean
                
                
                %figure
                use_meshx_comp = meshX_comp(:,:,1);
                use_meshy_comp = meshY_comp(:,:,1);
                
                surf(use_meshx_comp, use_meshy_comp, conditional_mean_comp,'FaceColor','red','EdgeColor','none')
                %surf(meshX2d, meshY2d, conditional_mean,'FaceColor','red','EdgeColor','none')
                %surf(conditional_mean, 'FaceColor','red','EdgeColor','none');
                camlight left; lighting phong
                %xlabel(channel1_name, 'FontSize', 16);
                %ylabel(channel2_name, 'FontSize', 16);
                %zlabel(channel3_name, 'FontSize', 16);
                title('Comp. Geo. based method', 'FontSize', 16);
                %zmax = ceil(max(max(conditional_mean)));
            end
       end
        
        
        function [Pfit,Delta, points_x, y_fit] = linear_fit_edge(obj, channel1_name, channel2_name, varargin)
            
            %
            Limits = [];
            for i=1:length(varargin)
                if(strcmp(varargin{i},'Limits'))
                    Limits = varargin{i+1};
                end
            end
            
            if(length(Limits)==0)
                
                %                [points_x,points_y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.90, 'no_plot');
                %                  channel1 = obj.name_channel_map(channel1_name);
                %                  channel2 = obj.name_channel_map(channel2_name);
                %                  points_x = transpose(obj.data(:,channel1));
                %                  points_y =  transpose(obj.data(:,channel2));
                [points_x,points_y] = obj.pairwise_visualize(channel1_name, channel2_name);
                
            else
                
                [points_x,points_y] = obj.pairwise_visualize(channel1_name, channel2_name, 'non_averaged_pts',.9, 'Limits', Limits);
                
            end
            
            hold on;
            
            %[points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name);
            
            %need to search for a fit
            
            [Pfit, S1] = polyfit(points_x,points_y,1);
            [y_fit, Delta] = polyval(Pfit, points_x, S1);
            %hold on;
            
            %y_fit = f(F_fitted,x);
            %plot(points_x,points_y,'*',points_x,y_fit,'r');
            plot(points_x,y_fit,'w','LineWidth',3.0);
            slope = Pfit(1)
            intercept = Pfit(2)
            %             area(points_x,y_fit,'FaceColor',[0.65,0.65,0.65]);
            %set(gca,'XLim',[Limits(1) Limits(3)]);
            %set(gca,'YLim',[Limits(2) Limits(4)]);
            %[r2 rmse] = rsquare(points_y,y_fit)
            
            
        end
        
        function [R ] = corrcoef_edge(obj, channel1_name, channel2_name)
            
            %[points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name,'no_plot');
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            points_x = obj.data(:,channel1);
            points_y = obj.data(:,channel2);
            
            R = corrcoef(points_x, points_y);
        end
        
        function [split, delta, Y1, Y2, Pfit1, Pfit2] = two_line_fit_edge(obj, channel1_name, channel2_name)
            
            [points_x, points_y] = obj.pairwise_visualize(channel1_name, channel2_name);
            
            %need to search for a fit
            delta = 1000;
            split = 0;
            Pfit1 = [];
            Pfit2 = [];
            Y1 = [];
            Y2 = [];
            
            
            
            num_points = length(points_x);
            
            for i=2:num_points-2
                
                
                
                [P1, S1] = polyfit(points_x(1:i),points_y(1:i),1);
                [P2, S2] = polyfit(points_x(i+1:num_points),points_y(i+1:num_points),1);
                
                [x_eval,Delta1] = polyval(P1, points_x(1:i), S1);
                [y_eval,Delta2] = polyval(P2, points_x(i+1:num_points),S2);
                d1 = mean(Delta1);
                d2 = mean(Delta2);
                
                if(max(d1,d2)<delta)
                    
                    delta = max(d1,d2);
                    split = i;
                    Pfit1 = P1;
                    Pfit2 = P2;
                    Y1 = x_eval;
                    Y2 = y_eval;
                end
                
                
            end
            
            hold on;
            %plot(points_x, points_y , '*');
            plot(points_x(1:split),Y1,'w','LineWidth',2.0);
            plot(points_x(split+1:num_points),Y2,'w','LineWidth',2.0);
            
        end
        
        function [shifted_normalized_density, xaxis, yaxis] = mean_shifted_density(obj, channel1_name, channel2_name)
            
            
            [points_x, points_y, ~, ~, normalized_density, xaxis,yaxis] = pairwise_visualize(obj, channel1_name, channel2_name,'no_plot');
            
            
            
            noise_DREMI = 0;
            [rows, cols] = size(normalized_density);
            
            
            grid_step = abs(yaxis(2)-yaxis(1));
            max_shift = floor(max(points_y)/grid_step);
            %now just make a bigger matrix
            
            shifted_normalized_density = zeros(rows+max_shift,cols);
            
            for i=1:cols
                
                
                column_shift = floor(points_y(i)/grid_step);
                start_position = max_shift - column_shift+1;
                shifted_normalized_density(start_position:start_position+rows-1,i)=normalized_density(:,i);
            end
            
            %                 plot(points_x, points_y,'w','LineWidth',3.0);
            %
            additional_points = (1:max_shift).*grid_step;
            additional_points = min(yaxis) - additional_points;
            yaxis_old = yaxis;
            yaxis = [transpose(additional_points); yaxis];
            %                I have to get the gridded values
            %                shifted_normalized_density_size = size(shifted_normalized_density);
            %                 for i=1:cols
            %
            %                     grid_position_shift = floor(points_y(i)/grid_step);
            %                     insertion_position = max_shift-grid_position_shift;
            %
            %                     shifted_normalized_density((insertion_position+1):(insertion_position+rows),i)=normalized_density(1:rows,i);
            %
            %                 end
            
            %                shifted_sum = sum(shifted_normalized_density,2);
            %                nonzero_rows = find(shifted_sum>0);
            %                yaxis = yaxis(nonzero_rows);
            %                shifted_normalized_density = shifted_normalized_density(nonzero_rows,:);
            
            colormap(linspecer(256));
            imagesc(xaxis,yaxis, shifted_normalized_density);
            set(gca,'YDir','normal');
            
            
            
        end
        
        
        
        
        function d = get_data(obj, channel_names)
            if iscell(channel_names)
                d = zeros(size(obj.data, 1), length(channel_names));
                for i = 1:length(channel_names)
                    d(:, i) = obj.data(:, obj.name_channel_map(channel_names{i}));
                end
            else
                d = obj.data(:, obj.name_channel_map(channel_names));
            end
        end
        
        function [points_x, points_y, point_weights, density, normalized_density, xaxis,yaxis, matrix_to_plot_no_side_bar] = pairwise_visualize(obj, channel1_name, channel2_name, varargin)
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            X = obj.data(:,channel1);
            Y = obj.data(:,channel2);
            total_cells = size(obj.data,1);
            
            
            num_slices = 256;
            minxval = min(X);
            minyval = min(Y);
            cutoff = 50;
            draw_contour = 0;
            show_density = 0;
            draw_plot = 1;
            avg_pts = 1;
            avg_pts_threshold = .9;
            fix_limits = 0;
            maxyval = max(Y);
            maxxval = max(X);
            fixy = 0;
            visual_threshold = 0;
            double_normalize = 0;
            show_side_bar = 1;
            
            for i=1:length(varargin)-1
                
                
                if(strcmp(varargin{i},'Slices'))
                    num_slices = varargin{i+1};
                end
                if(strcmp(varargin{i},'MinMaxY'))
                    
                    minyval = varargin{i+1};
                    maxyval = varargin{i+2};
                    
                    fixy = 1;
                    
                end
                if(strcmp(varargin{i},'MinMaxX'))
                    
                    minxval = varargin{i+1};
                    maxxval = varargin{i+2};
                    fixy = 1;
                    
                end
                if(strcmp(varargin{i},'MinMaxFitX'))
                    
                    minxval = min(X);
                    
                end
                if(strcmp(varargin{i},'MinMaxFitY'))
                    
                    minyval = min(Y);
                    
                    
                end
                
                if(strcmp(varargin{i},'Cutoff'))
                    
                    cutoff = varargin{i+1};
                    
                end
                
                if(strcmp(varargin{i},'Minval'))
                    
                    minxval = varargin{i+1};
                    minyval = minxval;
                end
                if(strcmp(varargin{i},'Limits'))
                    
                    fix_limits = 1;
                    limitvector = varargin{i+1};
                    minx = limitvector(1);
                    miny = limitvector(2);
                    maxx = limitvector(3);
                    maxy = limitvector(4);
                    
                    
                end
                if(strcmp(varargin{i},'non_averaged_pts'))
                    
                    avg_pts = 0;
                    avg_pts_threshold = varargin{i+1};
                    
                end
                if(strcmp(varargin{i}, 'visual_threshold'))
                    visual_threshold = varargin{i+1};
                end
                
                if (strcmp(varargin{i}, 'no_plot'))
                    draw_plot = 0;
                end
                
                
            end
            
            
            
            for i=1:length(varargin)
                if(strcmp(varargin{1},'double_normalize'))
                    double_normalize = 1;
                end
                if(strcmp(varargin{i},'draw_contour'))
                    draw_contour = 1;
                end
                if(strcmp(varargin{i},'show_density'))
                    show_density = 1;
                end
                if(strcmp(varargin{i},'no_plot'))
                    draw_plot = 0;
                end
                
                
                if strcmp(varargin{i}, 'no_side_bar')
                    show_side_bar = 0;
                end
                
            end
            
            
            
            
            
            
            
            
            if(fix_limits == 0)
                
                [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[minxval minyval],[maxxval maxyval]);
                
            else
                
                [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[minx miny],[maxx maxy]);
                %[minx, maxx, miny, maxy]
                
                
            end
            bandwidth;
            %rounding down small values
            
            for i=1:num_slices
                for j=1:num_slices
                    
                    if(density(i,j)<0.00001)
                        density(i,j) = 0;
                    end
                end
            end
            
            
            %weeding out sparse ends           
            row=0;
            
            total_cells_so_far = 0;
            while(total_cells_so_far<cutoff)
                
                
                row = row+1;
                total_cells_so_far = length(find(Y>Grid_Y(num_slices-row,1)));
                
                
            end
            
            maxy = Grid_Y(num_slices-row,1);
            
            
            total_cells_so_far = 0;
            start_row = 0;
            while(total_cells_so_far<cutoff)
                
                
                
                start_row = start_row + 1;
                total_cells_so_far = length(find(Y<Grid_Y(start_row,1)));
                
            end
            
            
            miny = Grid_Y(start_row,1);
            
            %row = 0;
            %start_row = 1;
            
            total_cells_so_far = 0;
            col=0;
            while(total_cells_so_far<cutoff)
                
                
                
                col = col+1;
                total_cells_so_far = length(find(X>Grid_X(1,num_slices-col)));
                
            end
            
            maxx = Grid_X(1,num_slices-col);
            
            
            total_cells_so_far = 0;
            start_col=0;
            while(total_cells_so_far<cutoff)
                
                
                
                start_col = start_col+1;
                total_cells_so_far = length(find(X<Grid_X(1,start_col)));
                
            end
            
            minx = Grid_X(1,start_col);
            
            
            
            if(fix_limits == 1)
                
                start_row = 1;
                start_col = 1;
                row = 0;
                col = 0;
            end
            
            if(fixy==1)
                
                start_row = 1;
                row = 0;
                
            end
            density = density(start_row:num_slices-row,start_col:num_slices-col);
            num_cols = size(density,2);
            num_rows = size(density,1);
            xaxis = Grid_X(1,start_col:num_slices-col);
            yaxis = Grid_Y(start_row:num_slices-row,1);
            
            normalized_density = zeros(num_rows,num_cols);
            prob_normalized_density = zeros(num_rows,num_cols);
            %normalized by column for plotting the data
            if(double_normalize == 0)
                for i=1:num_cols
                    
                    %normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                    normalized_density(:,i) = density(:,i)/max(density(:,i));
                    prob_normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                    
                end
                
            end
            if(double_normalize == 1)
                sprintf('double_normalize')
                for i=1:num_cols
                    for j = 1:num_rows
                        %normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                        %normalized_density(j,i) = density(j,i)/min(max(density(j,:)),max(density(:,i)));
                        normalized_density(j,i) = log(density(j,i)/(sum(density(j,:))*sum(density(:,i))));
                        
                        
                    end
                    prob_normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                end
                
                %           for j = 1:num_rows
                %                 %normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                %                 normalized_density(j,:) = normalized_density(j,:)/max(normalized_density(j,:));
                %
                %           end
                %
                %           for i=1:num_cols
                %
                %             %normalized_density(:,i) = normalized_density(:,i)/norm(normalized_density(:,i),1);
                %             normalized_density(:,i) = normalized_density(:,i)/max(normalized_density(:,i));
                %             %prob_normalized_density(:,i) = density(:,i)/norm(density(:,i),1);
                %
                %           end
                
                
            end
            
            
            
            
            
            %now create the side bars
            
            
            colsum = sum(density,1);
            normalized_colsum = colsum./max(colsum);
            
            rowsum = sum(density,2);
            normalized_rowsum = rowsum./max(rowsum);
            
            
            
            
            %the corner is a fudge
            
            
            %blueval = min(normalized_colsum);
            blueval = 0;
            corner = ones(11,11).*blueval;
            
            %make the top bar
            
            %yaxis_increment = abs(yaxis(2)-yaxis(1,1));
            yaxis_increment = .01;
            yaxis_top_bar = [];
            top_bar = [];
            zero_vector = zeros(1,length(normalized_colsum));
            for i=1:1
                top_bar = [top_bar; zero_vector];
                yaxis_top_bar = [yaxis_top_bar; max(yaxis)+(yaxis_increment*i)];
                
            end
            for i=1:10
                top_bar = [top_bar; normalized_colsum];
                yaxis_top_bar = [yaxis_top_bar; max(yaxis)+(yaxis_increment*i)];
            end
            
            
            %make the side bar
            %xaxis_increment = abs(xaxis(2)-xaxis(1));
            xaxis_increment = .01;
            xaxis_side_bar = [];
            side_bar = [];
            zero_vector = zeros(length(normalized_rowsum),1);
            
            for i=1:1
                side_bar = [side_bar zero_vector];
                xaxis_side_bar = [xaxis_side_bar max(xaxis)+(xaxis_increment*i)];
                
            end
            
            for i=1:10
                side_bar = [side_bar normalized_rowsum];
                xaxis_side_bar = [xaxis_side_bar max(xaxis)+(xaxis_increment*i)];
            end
            
            
            
            %find the trace through the peak regions for the return value
            points_x = [];
            points_y = [];
            point_weights = [];
            if(avg_pts==1)
                for i=1:num_cols
                    
                    
                    
                    max_indices = find(normalized_density(:,i)>= avg_pts_threshold);
                    % points_y = [points_y mean(Grid_Y(max_indices,i))];
                    
                    points_x = [points_x xaxis(i)];
                    
                    new_point_y = dot(yaxis(max_indices),normalized_density(max_indices,i))/sum(normalized_density(max_indices,i));
                    if(isnan(new_point_y))
                        new_point_y = 0;
                    end
                    %points_y = [points_y mean(yaxis(max_indices))];
                    points_y = [points_y new_point_y];
                    
                end
                point_weights = ones(1,length(points_y));
            else
                
                for i=1:num_cols
                    
                    %instead of referring to the grid maybe just take all the points
                    %in the high density squares ??
                    
                    max_indices = find(normalized_density(:,i)>= avg_pts_threshold);
                    % points_y = [points_y mean(Grid_Y(max_indices,i))];
                    new_points = ones(1,length(max_indices)).*xaxis(i);
                    new_point_weights = transpose(normalized_density(max_indices,i));
                    new_point_weights = new_point_weights ./ (sum(new_point_weights));
                    points_x = [points_x new_points];
                    
                    %points_y(i) = dot(Grid_Y(start_row:255-row,i),orig_normalized_density(:,i));
                    y_indices = max_indices;
                    new_points_y = transpose(yaxis(y_indices));
                    points_y = [points_y new_points_y];
                    point_weights = [point_weights new_point_weights];
                    
                end
                
            end
            
            
            
            %    smoothed_normalized_density = zeros(num_rows, num_cols);
            %
            %    for i=1:num_rows
            %        for j = 2:num_cols-1
            %             smoothed_normalized_density(i,j) = (normalized_density(i,j-1)+normalized_density(i,j)+normalized_density(i,j+1))/3;
            %        end
            %
            %    end
            %imagesc(flipud(Grid_X(:,1)), flipud(transpose(Grid_Y(1,:))), normalized_density);
            smoothed_normalized_density=normalized_density;
            
            if(visual_threshold>0)
                smoothed_normalized_density = (smoothed_normalized_density>visual_threshold).*smoothed_normalized_density;
                
            end
            matrix_to_plot_no_side_bar = smoothed_normalized_density;
            matrix_to_plot = [smoothed_normalized_density side_bar];
            top_bar = [top_bar corner];
            matrix_to_plot = [matrix_to_plot; top_bar];
            
            xaxis_to_plot = [xaxis xaxis_side_bar];
            xaxis_to_plot_no_side = xaxis;
            yaxis_to_plot = [yaxis; yaxis_top_bar];
            yaxis_to_plot_no_side = yaxis;
            
            
            
            
            
            if(draw_plot)
                %       colormap(linspecer(256));
                %        imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                %         density_filtered = smoothed_normalized_density>.6;
                %         smoothed_normalized_density = smoothed_normalized_density.*density_filtered;
                
                %      density_filtered = matrix_to_plot>.6;
                %      matrix_to_plot = matrix_to_plot.*density_filtered;
                %
                %
                %    j = linspecer(256);
                %  j(1,:) = [ 1 1 1 ];
                %  colormap(j);
                colormap(linspecer(256));
                if(double_normalize==1)
                    imagesc(xaxis, yaxis, smoothed_normalized_density,[-16 -6]);
                else
                    
                    if show_side_bar == 1
                    %imagesc(xaxis, yaxis, smoothed_normalized_density);
                        imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                    else
                        imagesc(xaxis_to_plot_no_side,yaxis_to_plot_no_side, matrix_to_plot_no_side_bar);
                    end
                    %imagesc(xaxis, yaxis, smoothed_normalized_density);
                    %imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                end
                
                
                set(gca,'YDir','normal');
                %      set(gca,'XTick',[]);
                %      set(gca,'YTick',[]);
                %     set(gca, 'XTickLabel','');
                %     set(gca, 'YTickLabel','');
                
                
                %         set(gca, 'FontSize',16);
                %         xlabel(channel1_name);
                %         ylabel(channel2_name);
                xlabel(channel1_name, 'FontSize', 16)
                ylabel(channel2_name, 'FontSize', 16)
                %hold
            end
            
            if(draw_contour)
                [bandwidth,rdensity,rGrid_X,rGrid_Y]=kde2d([X Y],num_slices+1,[minx miny],[maxx maxy]);
                
                contour(rGrid_X, rGrid_Y, rdensity, 12);
                
            end
            
            if(show_density)
                
                f = ksdensity(X, points_x);
                plot(points_x,f, 'w', 'LineWidth',1.3);
                
            end
            
            
            
            
            for i=1:length(varargin)-1
                
                if(strcmp(varargin{i},'Title'))
                    
                    
                    title(varargin{i+1});
                    
                    
                end
            end
            
            
        end
        
        function cluster_data_objects = split_to_clusters(obj, cluster_indices)
            
            num_indices = max(cluster_indices);
            cluster_data_objects = cell(1,num_indices);
            
            for i=1:num_indices
                cluster_data_objects{i} = obj;
                new_data = obj.data(cluster_indices,:);
                cluster_data_objects{i}.data = new_data;
                
                
            end
            
            
        end
        
        
        function obj = cluster_gate_events(obj, cluster_indices)
            
            new_data = obj.data(cluster_indices,:);
            obj.data = new_data;
        end
        
        function [obj, gated_indices] = threshold_gate_events(obj, channel_name, thresh, greater)
            
            channel = obj.name_channel_map(channel_name);
            
            channel_data = obj.data(:,channel);
            
            if(strcmp(greater,'gt'))
                gated_indices = find(channel_data>thresh);
            end
            if(strcmp(greater,'lt'))
                gated_indices = find(channel_data<thresh);
            end
            
            
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            
        end
        
        function [obj, gated_indices] = high_low_gate_events(obj, channel_name, percentage_start, percentage_end)
            
            channel1 = obj.name_channel_map(channel_name);
            
            
            channel_data = obj.data(:,channel1);
            
            
            
            num_percentage_start = floor(percentage_start*length(channel_data))+1;
            num_percentage_end = floor(percentage_end*length(channel_data));
            
            
            
            [~, sorted_indices] = sort(channel_data,'ascend');
            
            
            
            
            gated_indices = sorted_indices(num_percentage_start:num_percentage_end);
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            
            
            
        end
        
        function [obj, gated_indices] = high_low_gate_events_interval(obj, channel_name, interval_start, interval_end)
            
            channel1 = obj.name_channel_map(channel_name);
            
            
            channel_data = obj.data(:,channel1);
            
            high_indices = find(channel_data>interval_start);
            low_indices = find(channel_data(high_indices)<interval_end);
            gated_indices = high_indices(low_indices);
            
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            
            
            
        end
        
        function obj = well_gate_events(obj, well_no)
            
            well_indices = find(obj.wells == well_no);
            new_data = obj.data(well_indices,:);
            obj.data = new_data;
            
            
        end
        
        function [obj] = manual_gate_events_box(obj, channel_name1, channel_name2, left, bottom, width, height)
            
            channel1 = obj.name_channel_map(channel_name1);
            channel2 = obj.name_channel_map(channel_name2);
            channel_data = [obj.data(:,channel1) obj.data(:,channel2)];
            candidate_indices = find((channel_data(:,1)>left)&(channel_data(:,1)<left+width));
            
            shorter_data_set = [channel_data(candidate_indices,1) channel_data(candidate_indices,2)];
            
            nested_indices = find((shorter_data_set(:,2)>bottom)&(shorter_data_set(:,2)<bottom+height));
            
            gated_indices = candidate_indices(nested_indices);
            
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            if(length(obj.time_points)>0)
                new_time_points = obj.time_points(gated_indices,:);
                obj.time_points = new_time_points;
            end
        end
        
        function obj = manual_gate_events(obj, channel_name1, channel_name2, num_rects)
            
            fhandle = figure;
            ax = gca;
            channel1 = obj.name_channel_map(channel_name1);
            channel2 = obj.name_channel_map(channel_name2);
            channel_data = [obj.data(:,channel1) obj.data(:,channel2)];
            [~, density, x, y] = kde2d([channel_data(:,1) channel_data(:,2)], 256);
            %hold on
            %plot(channel_data(:,1),channel_data(:,2),'b.','MarkerSize',5)
            
            contour(x, y, density, 30);
            hold
            
            gated_indices = [];
            num_events = size(obj.data,1);
            
            
            for j=1:num_rects
                
                rect = getrect(ax)
                left = rect(1)
                bottom = rect(2)
                width = rect(3)
                height = rect(4)
                
                candidate_indices = find((channel_data(:,1)>left)&(channel_data(:,1)<left+width) & (channel_data(:,2)>bottom)&(channel_data(:,2)<bottom+height));
                
                
                if(length(gated_indices)==0)
                    gated_indices = candidate_indices;
                else
                    gated_indices = union(gated_indices,candidate_indices);
                end
                
                
            end
            
            new_data = obj.data(gated_indices,:);
            obj.data = new_data;
            
            
        end
        
        function [F, XI] = plot_channel_density(obj, channel_name, varargin)
            col_use = [0, 0, 1];
            XI_given = 0;
            XI = [];
            for i=1:length(varargin)
                if(strcmp(varargin{i},'XI'))
                    XI_given = 1;
                    XI = varargin{i+1};
                end
                
                if strcmp(varargin{i}, 'col')
                    col_use = varargin{i+1};
                end
            end
            
            
            if iscell(channel_name)
                for i = 1:length(channel_name)
                    channel = obj.name_channel_map(channel_name{i});
                
                    if(XI_given==0)
                        [F, XI] = ksdensity(obj.data(:,channel));  
                        plot(XI,F,'LineWidth',3, 'col', col_use(i, :));
                        legend(channel_name)
                        set(gca, 'FontSize', 16)
                    end
                end
            else
                    channel = obj.name_channel_map(channel_name);
                    [F, XI] = ksdensity(obj.data(:,channel));
                    plot(XI,F,'LineWidth',3, 'col', col_use);
                    title(channel_name);
                    set(gca,  'FontSize', 16)
            end       
                
                
                %             optargin = size(varargin,2);
                %             if(optargin==1)
                %                 plot(XI,F,varargin{1},'LineWidth',5);
                %             else
                %figure;
%                 plot(XI,F,'LineWidth',2);
%                 title(channel_name, 'FontSize', 16)
                
                %             end            
            
            %             data = obj.data(:,channel);
            %             sorted_data = sort(data);
            %             num_points = length(data);
            %             low_index = floor(num_points*.01);
            %             high_index = floor(num_points-(num_points*.01));
            %             low_value = sorted_data(low_index);
            %             high_value = sorted_data(high_index);
            %             plot_as_vertical_lines([low_value high_value],'r');
            
            %             box on;
            %             set(gca,'XTick', []);
            %             set(gca, 'YTick', []);
            
        end
        
        function plot_channel_hist(obj, channel_name, num_bins, varargin)
            
            channel = obj.name_channel_map(channel_name);
            data = obj.data(:,channel);
            optargin = size(varargin,2);
            h = findobj(gca,'Type','patch');
            
            if(optargin==1)
                set(h,'FaceColor',varargin{1},'EdgeColor','w');
                hist(data, num_bins);
            else
                set(h,'FaceColor', 'b','EdgeColor','w');
                hist(data, num_bins);
            end
            
        end
        
        function [gmfit, cluster_ids, cluster_percents] = compute_channel_gm_mixture(obj, channel_name, K)
            
            channel = obj.name_channel_map(channel_name);
            gmfit = gmdistribution.fit(obj.data(:,channel),K);
            cluster_ids = cluster(gmfit,obj.data(:,channel));
            cluster_percents = zeros(K,1);
            total = size(obj.data,1);
            for i=1:K
                
                cluster_percents(i) = length(find(cluster_ids==i))/total;
            end
            
            %bar(cluster_sizes);
            
        end
        
        function plot_well_stat_density(obj, channel_name, mode)
            
            channel = obj.name_channel_map(channel_name);
            
            
            if strcmp(mode,'mean')
                
                [F, XI] = ksdensity(obj.well_means(:,channel));
                
            elseif strcmp(mode, 'var')
                
                [F, XI] = ksdensity(obj.well_vars(:,channel));
                
            elseif strcmp(mode, 'L2')
                
                [F, XI] = ksdensity(obj.well_L2s(:,channel));
                
            end
            
            plot(XI, F);
            
        end
        
        function  plot_well_channel_density(obj, well_no, channel_name, varargin)
            
            
            channel = obj.name_channel_map(channel_name)
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            well_no
            
            [F, XI] = ksdensity(well_data(:,channel));
            
            optargin = size(varargin,2);
            if(optargin==1)
                plot(XI,F,varargin{1},'LineWidth',2);
            else
                plot(XI,F,'LineWidth',2);
            end
            
            
            
        end
        
        function [high_pop_indices, low_pop_indices, gmfit] = find_high_low_well_populations(obj,well_no, channel_name, threshold)
            
            well_no
            channel = obj.name_channel_map(channel_name);
            well_data = obj.get_well_data(well_no, channel_name);
            initialization = zeros(length(well_data),1);
            
            for i=1:length(well_data)
                
                if(well_data(i)<threshold)
                    initialization(i) = 1;
                else
                    initialization(i) = 2;
                end
                
            end
            
            %do it a different way right here
            %gmdistribution.fit(bc_data,num_wells,'options',options, 'start', presorted_wells);
            options = statset('MaxIter',100);
            gmfit = gmdistribution.fit(well_data,2, 'options', options,'start', initialization);
            gmfit.mu
            idx = cluster(gmfit,well_data);
            
            if(gmfit.mu(1)>gmfit.mu(2))
                
                high_cluster_index = 1;
                low_cluster_index = 2;
                
            else
                
                high_cluster_index = 2;
                low_cluster_index = 1;
                
            end
            
            high_pop_indices = find(idx==high_cluster_index);
            low_pop_indices = find(idx==low_cluster_index);
            
            num_high = length(high_pop_indices);
            num_low = length(low_pop_indices);
            
            %             obj.plot_well_channel_density(well_no, channel_name);
            %             hold
            %             X1=[mvnrnd(gmfit.mu(1),gmfit.Sigma(1),500)];
            %             X2=[mvnrnd(gmfit.mu(2),gmfit.Sigma(2),500)];
            %             [f1,xi1] = ksdensity(X1);
            %             [f2,xi2] = ksdensity(X2);
            %             plot(xi1,f1,'r');
            %             plot(xi2,f2,'g');
            
            
        end
        
        function [mean_dist, var_dist, kurtosis_dist, skew_dist, num_peaks_before, num_peaks_after, peak_densities_before, peak_densities_after] = plot_well_channel_density_verbose(obj, well_no, channel_name, varargin)
            
            channel = obj.name_channel_map(channel_name)
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            well_no
            
            [F, XI] = ksdensity(well_data(:,channel));
            
            optargin = size(varargin,2);
            if(optargin==1)
                plot(XI,F,varargin{1});
            else
                plot(XI,F);
            end
            
            %now plot the control
            hold on
            %             %I am going to make the other one the other apoptosis control
            %            well_no2 = obj.get_well_no(8,12);
            %             well_indices2 = find(obj.wells == well_no2);
            %             well_data2 = obj.data(well_indices2,:);
            
            non_edge_plate_data = obj.get_nonedge_data(channel_name);
            
            F2 = ksdensity(non_edge_plate_data, XI);
            %             [F2, XI2] = ksdensity(well_data2(:,channel));
            
            if(optargin==2)
                plot(XI,F2,varargin{2});
            else
                plot(XI,F2,'r');
            end
            label_x = sprintf('%s level',channel_name);
            xlabel(label_x);
            ylabel('density');
            [l,n] = obj.get_well_names(well_no);
            
            switch l
                case 1
                    well_string = sprintf('well: A %d', n);
                case 2
                    well_string = sprintf('well: B %d', n);
                case 3
                    well_string = sprintf('well: C %d', n);
                case 4
                    well_string = sprintf('well: D %d', n);
                case 5
                    well_string = sprintf('well: E %d', n);
                case 6
                    well_string = sprintf('well: F %d', n);
                case 7
                    well_string = sprintf('well: G %d', n);
                case 8
                    well_string = sprintf('well: H %d', n);
            end
            
            
            
            legend(well_string,'channel density');
            
            
            hold off
            
            marker_channel_index = find(obj.marker_channels == channel);
            mean_dist = obj.well_mean_distances(well_no, marker_channel_index);
            var_dist = obj.well_var_distances(well_no, marker_channel_index);
            kurtosis_dist = obj.well_kurtosis(well_no, marker_channel_index)-obj.kurtosis(marker_channel_index);
            skew_dist = obj.well_skews(well_no, marker_channel_index)-obj.skews(marker_channel_index);
            num_peaks_after = obj.well_peaks(well_no, marker_channel_index);
            num_peaks_before = obj.peaks(marker_channel_index);
            
            [f,xi] = ksdensity(obj.data(:,channel));
            [peak_densities_before,~] = findpeaks(f);
            well_data = obj.get_well_data(well_no);
            [f2, xi2] = ksdensity(well_data(:,channel));
            [peak_densities_after,~] = findpeaks(f2);
            
            
        end
        
        function plot_well_point(obj, well_no, channel_name, mode)
            
            channel = obj.name_channel_map(channel_name);
            marker_channel_index = find(obj.marker_channels == channel);
            
            if strcmp(mode,'mean')
                
                value = obj.well_means(well_no, marker_channel_index);
                
                
                
                
                
            elseif strcmp(mode, 'var')
                
                
                value = obj.well_vars(well_no, marker_channel_index);
                
            elseif strcmp(mode, 'L2')
                
                
                
                value = obj.well_L2s(well_no, marker_channel_index);
                
            end
            
            X = [value; value];
            Y = [0; .9];
            line(X,Y);
            
        end
        
        
        function [well_corrs, plate_corrs] = plot_well_channel_channel_scatter(obj, well_no, channel1_name, channel2_name, varargin)
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            
            optargin = size(varargin,2);
            if(optargin>0)
                
                scatter(well_data(:,channel1), well_data(:, channel2), varargin{1});
            else
                scatter(well_data(:,channel1), well_data(:, channel2));
                
            end
            
            
            label_x = sprintf('%s level',channel1_name);
            label_y = sprintf('%s level',channel2_name);
            xlabel(label_x);
            ylabel(label_y);
            well_corrs = 0;
            plate_corrs = 0;
            
            %well_corrs = obj.get_well_channel_corrs(well_no, channel1_name, channel2_name);
            
            %marker_channel1 = find(obj.marker_channels == channel1);
            %marker_channel2 = find(obj.marker_channels == channel2);
            %plate_corrs = obj.corrs(marker_channel1, marker_channel2);
            
            
        end
        
        
        function [channel_data] = get_channel_data(obj, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            channel_data = obj.data(:, channel);
            
        end
        
        function [well_data] = get_well_data(obj, well_no, varargin)
            
            well_indices = find(obj.wells == well_no);
            optargin = size(varargin,2);
            
            if(optargin==1)
                
                channel_name = varargin{1};
                channel = obj.name_channel_map(channel_name);
                well_data = obj.data(well_indices, channel);
            else
                well_data = obj.data(well_indices, :);
            end
            
            
        end
        
        function plot_well_channels(obj, well_no, channel1_name, channel2_name)
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            scatter(well_data(:,channel1), well_data(:,channel2));
            
        end
        
        function well_mean = get_well_channel_mean(obj, well_no, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            well_mean = mean(well_data(:,channel))
            
        end
        
        function well_corrs = get_well_channel_corrs(obj, well_no, channel1_name, channel2_name)
            
            
            channel1 = obj.name_channel_map(channel1_name);
            channel2 = obj.name_channel_map(channel2_name);
            
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            R = corrcoef(well_data(:, channel1), well_data(:, channel2));
            well_corrs = R(1,2);
            
        end
        
        function well_var = get_well_channel_var(obj, well_no, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,:);
            well_var = var(well_data(:,channel))
            
        end
        
        function [sorted_corrs, sorted_strings] = get_well_channel_corrs_sort(obj, well_no, channel_name, computed_channels)
            
            channel = obj.name_channel_map(channel_name);
            computed_channels=transpose([obj.cell_length_col obj.dna_channels obj.marker_channels]);
            sorted_corrs = zeros(length(computed_channels),2);
            sorted_corrs(:,1) = computed_channels;
            for i=1:length(computed_channels)
                
                sorted_corrs(i,2) = obj.get_well_channel_corrs(well_no, channel, computed_channels(i));
                
            end
            
            sorted_corrs = sortrows(sorted_corrs,2);
            sorted_strings = [];
            for i=1:length(computed_channels)
                sorted_strings = [sorted_strings sprintf('%s %g \n', char(channel_headings(sorted_corrs(i,1))),sorted_corrs(i,2))];
            end
            
            
        end
        
        function no_events = num_events(obj)
            
            [no_events, ~] = size(obj.data);
        end
        
        function ws = well_size(obj, well_no)
            
            ws = size(find(obj.wells == well_no));
            
        end
        
        function obj = compute_well_names(obj, Letter, Number, swap)
            %I hate this function and the irritating
            %pattern-less way that this is assigned
            
            
            obj.well_2_names=zeros(96,2);
            obj.names_2_well = zeros(8,12);
            
            %this is just making me physically angry
            %FIGURE OUT WHICH PLATES ARE SWAPPED
            swapped_letter = [5 6 7 8 1 2 3 4];
            
            for i=1:96
                
                if(swap)
                    l = swapped_letter(Letter(i));
                else
                    l = Letter(i);
                end
                
                n = Number(i);
                obj.well_2_names(i,1) = l;
                obj.well_2_names(i,2) = n;
                obj.names_2_well(l,n) = i;
            end
            
            
            
        end
        
        function well_no = get_well_no(obj, letter, number)
            
            well_no = obj.names_2_well(letter, number);
        end
        
        
        function [letter, number ] = get_well_names(obj, well_no)
            
            letter = obj.well_2_names(well_no, 1);
            number = obj.well_2_names(well_no, 2);
            
        end
        
        
        function obj=compute_well_distances(obj)
            
            well_distances = zeros(96,2);
            
            for i=1:96
                
                [dist, ~] = obj.compute_distance_from_control(i);
                well_distances(i,1) = i;
                well_distances(i,2) = dist;
                
            end
            
            well_distances = sortrows(well_distances,-2);
            obj.well_distances = well_distances;
            
        end
        
        
        function well_counts = compute_well_counts(obj)
            
            well_counts = zeros(96,1);
            
            for i=1:96
                
                well_indices = find(obj.wells==i);
                well_counts(i) = length(well_indices);
            end
            obj.well_counts = well_counts;
            
        end
        
        
        
        
        function obj = compute_all_densities_L2s(obj, varargin)
            
            num_marker_channels = length(obj.marker_channels);
            
            obj.well_L2s = zeros(96,num_marker_channels);
            
            
            obj.well_densities = zeros(96,num_marker_channels,500);
            
            obj.channel_density_xvals = zeros(num_marker_channels, 500);
            obj.channel_densities = zeros(num_marker_channels, 500);
            
            optargin = size(varargin, 2);
            
            if(optargin == 0)
                
                for j = 1:num_marker_channels
                    
                    mc = obj.marker_channels(j);
                    [P, X] = ksdensity(obj.data(:,mc), 'npoints', 500);
                    obj.channel_density_xvals(j,:) = X;
                    obj.channel_densities(j,:) = P;
                    
                end
            else
                
                obj.channel_density_xvals = varargin{1};
                
                for j = 1:num_marker_channels
                    
                    mc = obj.marker_channels(j);
                    P = ksdensity(obj.data(:,mc), obj.channel_density_xvals(j,:));
                    obj.channel_densities(j,:) = P;
                    
                end
                
            end
            
            
            for i =1:96
                for j = 1:num_marker_channels
                    
                    mc = obj.marker_channels(j);
                    P = obj.compute_density(mc,i, obj.channel_density_xvals(j,:));
                    obj.well_densities(i,j,:) = P;
                    obj.well_L2s(i,j) = norm(P-obj.channel_densities(j,:));
                    %these have been computed at the appropriate
                    %locations
                end
            end
            
        end
        
        
        
        function well_L2 = compute_well_L2_norms(obj, well_no, varargin)
            
            num_markers = length(obj.marker_channels);
            well_L2 = zeros(1,num_markers);
            optargin = size(varargin,2);
            
            
            for i=1:num_markers
                
                channel_name = obj.channel_name_map{obj.marker_channels(i)};
                if(optargin == 0)
                    
                    well_L2(1,i) = obj.compute_distro_L2(channel_name, well_no);
                else
                    
                    well_L2(1,i) = obj.compute_distro_L2(channel_name, well_no, varargin{1});
                    
                end
                
            end
            
            
        end
        
        
        
        
        function well_knock_vector = compute_well_knockdowns(obj, well_no, control_no, up_down)
            
            num_markers = length(obj.marker_channels);
            well_knock_vector = zeros(1,num_markers);
            
            for i=1:num_markers
                
                if(up_down == 0)
                    
                    well_knock_vector(1,i) = obj.compute_knockdown_strength(well_no, control_no, obj.marker_channels(i));
                    
                else
                    
                    well_knock_vector(1,i) = obj.compute_knockup_strength(well_no, control_no, obj.marker_channels(i));
                    
                end
                
            end
            
            
        end
        
        function obj = compute_all_knockdowns(obj, control_no)
            
            num_markers = length(obj.marker_channels);
            obj.well_knockdowns = zeros(96,num_markers);
            %get rid of the control wells later
            %there are 8 control wells
            
            for i=1:8
                for j=1:12
                    well_no = obj.get_well_no(i,j);
                    obj.well_knockdowns(well_no, :) = obj.compute_well_knockdowns(well_no,control_no, 0);
                end
            end
            
            obj.well_knockups = zeros(96,num_markers);
            %there are 8 control wells
            
            for i=1:8
                for j=1:12
                    well_no = obj.get_well_no(i,j);
                    obj.well_knockups(well_no, :) = obj.compute_well_knockdowns(well_no,control_no, 1);
                end
            end
            
            
        end
        
        
        
        function tot_prob = compute_knockup_strength(obj, well_no, control_no, channel_name)
            
            channel = obj.name_channel_map(channel_name);
            well_indices_ctrl = find(obj.wells == control_no);
            well_data_ctrl = obj.data(well_indices_ctrl, channel);
            
            well_indices = find(obj.wells == well_no);
            well_data = obj.data(well_indices,channel);
            
            
            
            ctrl_size = length(well_data_ctrl);
            data_size = length(well_data);
            tot_prob = 0;
            
            for i=1:ctrl_size
                
                num_more = sum(well_data>well_data_ctrl(i));
                prob_more = (num_more/data_size) * (1/ctrl_size);
                %
                tot_prob = tot_prob+prob_more;
            end
            
            
            %             well_indices_ctrl = find(obj.wells == control_no);
            %             well_data_ctrl = obj.data(well_indices_ctrl,channel);
            %
            %             well_indices = find(obj.wells == well_no);
            %             well_data = obj.data(well_indices,channel);
            %
            %
            %             ctrl_size = length(well_data_ctrl);
            %             data_size = length(well_data);
            %             sort(well_data,'descend'); % sort data in ascending order
            %             sort(well_data_ctrl,'descend');
            %
            %             data_index = 1;
            %             ctrl_index = 1;
            %             tot_prob = 0;
            %             num_more = 0;
            %
            %             for i=1:ctrl_size
            %
            %                 ctrl_index = i;
            %
            %
            %                 if well_data(data_index)<well_data_ctrl(ctrl_index)
            %
            %                    prob_more = (num_more/data_size) * (1/ctrl_size);
            %                    tot_prob = tot_prob + prob_more;
            %                    continue;
            %                 end
            %
            %
            %                 while well_data(data_index)>well_data_ctrl(ctrl_index)
            %                     data_index = data_index+1;
            %                     if(data_index > data_size)
            %                         break;
            %                     end
            %                 end
            %
            %                 if data_index > data_size
            %                     break;
            %                 end
            %
            %                 if well_data(data_index) < well_data_ctrl(ctrl_index)
            %
            %                     num_more = data_index-1;
            %                 else
            %                     num_more = data_index;
            %                 end
            %
            %                 prob_more = (num_more/data_size) * (1/ctrl_size);
            %                 tot_prob = tot_prob + prob_more;
            %
            %
            %             end
            
        end
        
        function obj = compute_all_spearman_coeff(obj, varargin)
            
            
            obj.well_spearman= zeros(96,1);
            
            optargin = size(varargin,2);
            
            for i=1:96
                
                if(optargin==0)
                    obj.well_spearman(i) = obj.compute_well_spearman_coeff(i);
                else
                    obj.well_spearman(i) = obj.compute_well_spearman_coeff(i,varargin{1});
                end
                
            end
            
            
        end
        
        function sp_coeff = compute_well_spearman_coeff(obj, well_no, varargin)
            
            
            well_vector = transpose(obj.compute_well_means(well_no));
            
            optargin = size(varargin, 2);
            
            if(optargin==0)
                
                control_vector = transpose(obj.means);
            else
                control_vector = transpose(obj.compute_well_means(varargin{1}));
            end
            
            [r,t,p] = spear(well_vector, control_vector);
            sp_coeff = r;
            
        end
        
        
        function L2_diff = compute_distro_L2(obj, channel_name, well_no, varargin)
            
            channel = obj.name_channel_map(channel_name);
            
            if(size(varargin,2) == 0)
                
                
                non_edge_data = obj.get_nonedge_data(channel_name);
                [P_ctrl,XI]=ksdensity(non_edge_data);
                %actually get the non-edge wells here
                %average of all wells
                
                
            else
                control_no = varargin{1};
                ctrl_data = obj.get_well_data(control_no, channel_name);
                [P_ctrl,XI]=ksdensity(ctrl_data);
                
            end
            
            
            well_data = obj.get_well_data(well_no, channel_name);
            P_well = ksdensity(well_data, XI);
            
            L2_diff = norm(P_ctrl-P_well);
            
        end
        
        
        function [P_dist, XI] = compute_density(obj, channel, well_no, varargin)
            
            channel
            channel_name = obj.channel_name_map{channel};
            well_data = obj.get_well_data(well_no, channel_name);
            
            if(size(varargin,2) == 0)
                
                [P_dist, XI] = ksdensity(well_data, 'npoints', 500);
                
            else
                
                XI = varargin{1};
                P_dist = ksdensity(well_data, XI);
                
            end
            
            
            
        end
        
        function non_edge_data = get_nonedge_data(obj, channel)
            
            
            non_edge_data=[];
            
            for i=1:8
                for j=1:11
                    well_no = obj.get_well_no(i,j);
                    %non_edge_data = [non_edge_data; obj.get_well_data(well_no,channelname)];
                    well_data = obj.get_well_data(well_no);
                    non_edge_data = [non_edge_data; well_data(:,channel)];
                end
            end
            
            
        end
        
        function edge_data = get_edge_data(obj, channelname)
            
            sprintf('getting edge data')
            
            edge_data=[];
            
            edge_wells = [];
            
            %top and bottom edges
            for i=2:11
                well_no = obj.get_well_no(1,i);
                edge_data = [edge_data; obj.get_well_data(well_no,channelname)];
                edge_wells = [edge_wells well_no];
                well_no = obj.get_well_no(8,i);
                edge_data = [edge_data; obj.get_well_data(well_no,channelname)];
                edge_wells = [edge_wells well_no];
            end
            
            %side edges, not the side with control wells
            for i = 1:8
                
                well_no = obj.get_well_no(i,1);
                edge_data = [edge_data; obj.get_well_data(well_no,channelname)];
                edge_wells = [edge_wells well_no];
                
            end
            
            edge_wells = sort(edge_wells)
            
        end
        
        function tot_prob = compute_knockdown_strength(obj, channel_name, well_no, varargin)
            
            channel = obj.name_channel_map(channel_name);
            
            if(size(varargin,2)>0)
                
                control_no = varargin{1};
                well_data_ctrl = obj.get_well_data(control_no, channel_name);
                
            else
                
                well_data_ctrl = obj.data(:, channel);
                
            end
            
            
            well_data = obj.get_well_data(well_no, channel_name);
            
            
            
            ctrl_size = length(well_data_ctrl);
            data_size = length(well_data);
            tot_prob = 0;
            
            
            for i=1:ctrl_size
                
                num_less = sum((well_data<well_data_ctrl(i)));
                prob_less = (num_less/data_size) * (1/ctrl_size);
                %
                tot_prob = tot_prob + prob_less;
            end
            
            
            
            %             sort(well_data); % sort data in ascending order
            %             sort(well_data_ctrl);
            %
            %             data_index = 1;
            %             ctrl_index = 1;
            
            %             num_less = 0;
            %
            %
            %             for i=1:ctrl_size
            %
            %                 ctrl_index = i;
            %
            %
            %                 if well_data(data_index)>well_data_ctrl(ctrl_index)
            %
            %                    prob_less = (num_less/data_size) * (1/ctrl_size);
            %                    tot_prob = tot_prob + prob_less;
            %                    continue;
            %                 end
            %
            %
            %                 while well_data(data_index)<well_data_ctrl(ctrl_index)
            %                     data_index = data_index+1;
            %                     if(data_index > data_size)
            %                         break;
            %                     end
            %                 end
            %
            %                 if data_index > data_size
            %                     break;
            %                 end
            %
            %                 if well_data(data_index) > well_data_ctrl(ctrl_index)
            %
            %                     num_less = data_index-1;
            %                 else
            %                     num_less = data_index;
            %                 end
            %
            %                 prob_less = (num_less/data_size) * (1/ctrl_size);
            %                 tot_prob = tot_prob + prob_less;
            %
            %
            %             end
            
        end
        
        function obj = perform_bead_gate(obj, gate_channels, bead_threshold, nonbead_threshold)
            
            [num_events,~] = size(obj.data)
            bead_threshold
            nonbead_threshold
            num_markers = length(obj.marker_channels);
            num_dna = length(obj.dna_channels);
            obj.beads = ones(1,num_events);
            
            
            %repeat the same thing for dna channels
            
            for i=1:num_events
                
                for j=1:num_dna
                    
                    dna_channel = obj.dna_channels(j);
                    if(obj.data(i,dna_channel)>nonbead_threshold)
                        obj.beads(i) = 0;
                        break;
                    end
                    
                end
                
            end
            
            
            %bead channels can overlap with marker channels
            for i=1:num_events
                
                
                %skip the ones that are already considered as non beads
                if(obj.beads(i) == 0)
                    continue;
                end
                
                for j=1:length(gate_channels)
                    
                    bead_channel = gate_channels(j);
                    if(obj.data(i,bead_channel)<bead_threshold)
                        obj.beads(i) = 0;
                        break;
                    end
                end
                
            end
            
            %have to check that these guys are beads actually
            num_beads = length(find(obj.beads==1))
            num_non_beads = length(find(obj.beads==0))
            
        end
        
        function obj = reorder_events_by_channel(obj, channel_data)
            
            % channel =  obj.name_channel_map(channel_name);
            
            %   channel_data = obj.data(:,channel);
            [ ~, new_order] = sort(channel_data);
            obj.data = obj.data(new_order,:);
            
        end
        
        function [cdata_fracs] = split_object_into_fractions(obj, channel_data, num_fracs)
            
            cdata = obj.reorder_events_by_channel(channel_data);
            
            num_cells = floor(size(obj.data,1)/num_fracs);
            cdata_fracs = cell(1,num_fracs);
            
            for i=1:num_fracs
                
                cdata_fracs{i} = cdata;
                start_index = ((i-1)*num_cells)+1;
                finish_index = (i)*num_cells;
                cdata_fracs{i}.data = cdata.data(start_index:finish_index,:);
            end
            
        end
        
        
        function obj = structure_learn_CTBN(obj, channel_names, wanderlust_channels, max_parents, start_inds)
            
            %discretize channels
            obj = discretize_channels_GMM(channel_names);
            %compute and bin trajectories
            % I have not set the starting indices
            
            start_inds = 1:54032;
            %should be how many 0 hour values there are, just hard coded it for now, put it in the argument
            obj = obj.run_and_average_wanderlust(channel_names, 1000, 10, 100, start_inds);
            %not sure what start inds
            
            obj = obj.bin_wanderlust(100, 100);
            
            obj.ctbn_parent_map = containers.Map();
            
            for i=1:length(channel_names)
                
                best_parent_channels = obj.compute_best_parents(channel_names{i}, channel_names, max_parents);
                obj.ctbn_parent_map{channel_names{i}} = best_parent_channels;
                
            end
            
            
        end
        
        
        function best_parent_channels = compute_best_parents(obj, channel_name, channel_names, max_parents)
            
            best_parent_score = 0;
            
            best_parent_channels = [];
            %test up to n-parent combinations
            for j=1:max_parents
                
                C = nchoosek(1:length(channel_names),j);
                
                
                for k = 1:size(C,2)
                    current_combination = C(k,:)
                    
                    parent_channels = channel_names(current_combination);
                    
                    marginal_likelihood = obj.marginal_likelihood_transitions(channel_name, parent_channels);
                    marginal_likelihood = marginal_likelihood + obj.marginal_likelihood_duration(channel_name, parent_channels);
                    complexity_penalty = log(j); % or some function of it.
                    parent_score = marginal_likelihood + complexity_penalty;
                    
                    if(parent_score>best_parent_score)
                        
                        best_parent_score = parent_score;
                        best_parent_channels = parent_channels;
                    end
                    
                end
            end
        end
        
        function mlikelihood_transitions = marginal_likelihood_transitions(obj, channel_name, parent_channels)
            %computes the marginal likelihood of a particular choice of
            %parents and child
            parent_channels
            mlikelihood_transitions = 1;
            parent_channel_levels = zeros(1,length(parent_channels));
            
            for i = 1:length(parent_channels)
                
                channel = obj.name_channel_map(parent_channels{i});
                parent_channel_levels(i) = max(obj.discretized_data(:,channel));
                
            end
            
            parent_channel_levels
            
            parent_value_instantiations=zeros(prod(parent_channel_levels),length(parent_channels));
            
            
            for i=1:length(parent_channels)
                
                
                ones_vector_before = [1];
                ones_vector_after = [1];
                if(i>1)
                    ones_vector_before = ones(prod(parent_channel_levels(1:i-1)),1);
                end
                if(i<length(parent_channels))
                    
                    ones_vector_after = ones(prod(parent_channel_levels(i+1:end)),1);
                end
                run = transpose(1:parent_channel_levels(i));
                parent_value_instantiations(:,i) = kron(ones_vector_before, kron(run,ones_vector_after))
                
            end
            
            
            %the parent value instantiations stores all of the discrete
            %instantiations of the parents
            
            child_channel = obj.name_channel_map(channel_name);
            child_channel_levels = max(obj.discretized_data(:,child_channel));
            
            
            for i=1:length(parent_value_instantiations)
                
                parent_channel_values = parent_value_instantiations(i,:);
                
                for j=1:child_channel_levels
                    
                    prior_transitions_out = obj.count_transitions_out_prior(channel_name, j, prod(parent_channel_levels))
                    prod_term_numerator = gammaln(prior_transitions_out)
                    num_transitions_out = obj.count_transitions_out(channel_name, parent_channels, j, parent_channel_values)
                    prod_term_denominator = gammaln(prior_transitions_out+num_transitions_out)
                    prod_term = prod_term_numerator-prod_term_denominator
                    mlikelihood_transitions = mlikelihood_transitions + prod_term;
                    
                    %prod_term is specific to a parent value
                    %combination and a child value
                    prod_term
                end
                
                
            end
            
            for i=1:length(parent_value_instantiations)
                
                parent_channel_values = parent_value_instantiations(i,:);
                
                for j=1:child_channel_levels
                    
                    for k=1:child_channel_levels
                        
                        if (j==k)
                            continue;
                        end
                        
                        prior_transitions = obj.count_transitions_prior(channel_name, j, k, prod(parent_channel_levels));
                        
                        num_transitions = obj.count_transitions(channel_name, parent_channels, j, k, parent_channel_values);
                        
                        prod_term_numerator = gammaln(prior_transitions + num_transitions);
                        prod_term_denominator = gammaln(prior_transitions);
                        prod_term = prod_term_numerator-prod_term_denominator;
                        mlikelihood_transitions = mlikelihood_transitions + prod_term;
                    end
                    
                end
                
                
            end
            
            mlikelihood_transitions
            
            
        end
        
        function mlikelihood_duration = marginal_likelihood_duration(obj, channel_name, parent_channels)
            
            parent_channel_levels = zeros(1,length(parent_channels));
            
            
            
            
            for i = 1:length(parent_channels)
                
                channel = obj.name_channel_map(parent_channels{i});
                parent_channel_levels(i) = max(obj.discretized_data(:,channel));
                
            end
            
            parent_value_instantiations=zeros(prod(parent_channel_levels),length(parent_channels));
            
            
            for i=1:length(parent_channels)
                
                
                ones_vector_before = [1];
                ones_vector_after = [1];
                if(i>1)
                    ones_vector_before = ones(prod(parent_channel_levels(1:i-1)),1);
                end
                if(i<length(parent_channels))
                    ones_vector_after = ones(prod(parent_channel_levels(i+1:end)),1);
                end
                
                run = transpose(1:parent_channel_levels(i));
                parent_value_instantiations(:,i) = kron(ones_vector_before, kron(run, ones_vector_after));
                
            end
            
            child_channel = obj.name_channel_map(channel_name);
            child_channel_levels = max(obj.discretized_data(:,child_channel));
            
            mlikelihood_duration = 1;
            for i=1:length(parent_value_instantiations)
                
                parent_channel_values = parent_value_instantiations(i,:);
                
                for j=1:child_channel_levels
                    
                    prior_transitions_out = obj.count_transitions_out_prior(channel_name, j, prod(parent_channel_levels));
                    num_transitions_out = obj.count_transitions_out(channel_name, parent_channels, j, parent_channel_values);
                    prior_duration = obj.total_duration_prior(channel_name, j, prod(parent_channel_levels));
                    duration = obj.total_duration(channel_name, parent_channels, j, parent_channel_values);
                    
                    prod_term_numerator = gammaln(prior_transitions_out+num_transitions_out+1)+log(prior_duration^(prior_transitions_out+1));
                    prod_term_denominator = gammaln(prior_transitions_out+1)+ log(((prior_duration+duration)^(prior_transitions_out+num_transitions_out+1)));
                    prod_term = prod_term_numerator - prod_term_denominator;
                    mlikelihood_duration = mlikelihood_duration + prod_term;
                    
                    
                end
                
                
            end
            
        end
        
        function duration_prior = total_duration_prior(obj, channel_name, channel_value, num_parent_instantiations)
            
            duration = 0;
            channel = obj.name_channel_map(channel_name);
            
            [~, num_trajectories] = size(obj.trajectories);
            
            for i=1:num_trajectories
                
                current_trajectory = obj.trajectories(:,i);
                current_trajectory_times = obj.trajectory_times(:,i);
                current_trajectory_channel_values = obj.discretized_data(current_trajectory,channel);
                
                
                
                child_indices = find(current_trajectory_channel_values == channel_value);
                
                time_start = current_trajectory_times(child_indices(1));
                
                for j=2:length(child_indices)
                    
                    current_index = child_indices(j);
                    prev_index = child_indices(j-1);
                    
                    if((current_index-prev_index > 1)|(j==length(child_indices)))
                        %if the run stopped or we're at the end
                        
                        time_end = current_trajectory_times(prev_index);
                        
                        duration = duration + (time_end - time_start);
                        time_start = current_trajectory_times(current_index);
                        
                    end
                end
                
            end
            
            duration_prior = duration/num_parent_instantiations;
            
        end
        
        function duration = total_duration(obj, channel_name, parent_channel_names, channel_value, parent_channel_values)
            
            duration = 0;
            [~, num_trajectories] = size(obj.trajectories);
            channel = obj.name_channel_map(channel_name);
            for i=1:length(parent_channel_values)
                
                parent_channels(i) = obj.name_channel_map(parent_channel_names{i});
            end
            
            for i=1:num_trajectories
                
                current_trajectory = obj.trajectories(:,i);
                current_trajectory_times = obj.trajectory_times(:,i);
                current_trajectory_channel_values = obj.discretized_data(current_trajectory,channel);
                current_trajectory_parent_values = obj.discretized_data(current_trajectory, parent_channels);
                
                %select the indices where the parents are equal to the
                %parent values
                parent_indices = [];
                for j=1:length(current_trajectory_parent_values)
                    if(isequal(current_trajectory_parent_values(j,:), parent_channel_values))
                        parent_indices = [parent_indices j];
                    end
                end
                
                child_trajectory_channel_values = current_trajectory_channel_values(parent_indices);
                
                child_indices = find(child_trajectory_channel_values == channel_value);
                
                absolute_child_indices = parent_indices(child_indices);
                
                %of those select the indices where the channel is equal
                %to the channel value.
                
                %how do I add up the total value in this thingamagig?
                %so there should be series of runs that are
                %consecutive, and wherever they end is probably the
                %time that ends it
                
                %pick consecutive runs
                
                time_start = current_trajectory_times(absolute_child_indices(1));
                
                for j=2:length(absolute_child_indices)
                    
                    current_index = absolute_child_indices(j);
                    prev_index = absolute_child_indices(j-1);
                    
                    if((current_index-prev_index > 1)|(j==length(absolute_child_indices)))
                        %if the run stopped or we're at the end
                        
                        time_end = current_trajectory_times(prev_index);
                        
                        duration = duration + (time_end - time_start);
                        time_start = current_trajectory_times(current_index);
                        
                    end
                end
                
            end
            
        end
        
        
        function num_transitions_out_prior = count_transitions_out_prior(obj, channel_name, channel_value, num_parent_instantiations)
            
            num_transitions_out = 0;
            
            [~, num_trajectories] = size(obj.trajectories);
            channel = obj.name_channel_map(channel_name);
            
            for i=1:num_trajectories
                
                
                current_trajectory = obj.trajectories(:,i);
                current_trajectory_times = obj.trajectory_times(:,i);
                current_trajectory_channel_values = obj.discretized_data(current_trajectory,channel);
                
                
                
                
                for j=2:length(current_trajectory_channel_values)
                    
                    current_index = j;
                    prev_index = j-1;
                    
                    if(current_index-prev_index == 1)
                        %so parents have remained constant
                        
                        current_value = current_trajectory_channel_values(current_index);
                        prev_value = current_trajectory_channel_values(prev_index);
                        
                        if((prev_value == channel_value)&(current_value~=channel_value))
                            num_transitions_out = num_transitions_out + 1;
                        end
                        
                    end
                    
                    
                end
                
            end
            
            num_transitions_out_prior = num_transitions_out/num_parent_instantiations;
            
        end
        
        
        function num_transitions_out = count_transitions_out(obj, channel_name, parent_channel_names, channel_value, parent_channel_values)
            
            %again find the parent channel and the child channel values
            num_transitions_out = 0;
            channel = obj.name_channel_map(channel_name);
            
            for i=1:length(parent_channel_names)
                
                parent_channels(i) = obj.name_channel_map(parent_channel_names{i});
            end
            
            [~, num_trajectories] = size(obj.trajectories);
            for i=1:num_trajectories
                
                current_trajectory = obj.trajectories(:,i);
                current_trajectory_times = obj.trajectory_times(:,i);
                current_trajectory_channel_values = obj.discretized_data(current_trajectory,channel);
                current_trajectory_parent_values = obj.discretized_data(current_trajectory, parent_channels);
                
                %select the indices where the parents are equal to the
                %parent values
                parent_indices = [];
                for j=1:length(current_trajectory_parent_values)
                    if(isequal(current_trajectory_parent_values(j,:), parent_channel_values))
                        parent_indices = [parent_indices j];
                    end
                end
                
                
                for j=2:length(parent_indices)
                    
                    current_index = parent_indices(j);
                    prev_index = parent_indices(j-1);
                    
                    if(current_index-prev_index == 1)
                        %so parents have remained constant
                        
                        current_value = current_trajectory_channel_values(current_index);
                        prev_value = current_trajectory_channel_values(prev_index);
                        
                        if((prev_value == channel_value)&(current_value~=channel_value))
                            num_transitions_out = num_transitions_out + 1;
                        end
                        
                    end
                    
                    if(current_index-prev_index > 1)
                        
                        num_transitions_out = num_transitions_out+1;
                        
                    end
                end
                
            end
            
            
        end
        
        function num_transitions_prior = count_transitions_prior(obj, channel_name, channel_value_from, channel_value_to, num_parent_instantiations)
            
            
            num_transitions = 0;
            channel = obj.name_channel_map(channel_name);
            [~, num_trajectories] = size(obj.trajectories);
            for i=1:num_trajectories
                
                current_trajectory = obj.trajectories(:,i);
                current_trajectory_times = obj.trajectory_times(:,i);
                current_trajectory_channel_values = obj.discretized_data(current_trajectory,channel);
                
                
                
                for j=2:length(current_trajectory_channel_values)
                    
                    current_index =j;
                    prev_index = j-1;
                    
                    if(current_index-prev_index == 1)
                        %u has stayed steady
                        
                        val1 = current_trajectory_channel_values(prev_index);
                        val2 = current_trajectory_channel_values(current_index);
                        if((val1==channel_value_from)&(val2==channel_value_to))
                            
                            num_transitions = num_transitions+1;
                            
                        end
                        
                    end
                end
                
            end
            
            num_transitions_prior = num_transitions/num_parent_instantiations;
            
        end
        
        
        function num_transitions = count_transitions(obj, channel_name, parent_channel_names, channel_value_from, channel_value_to, parent_channel_values)
            
            num_transitions = 0;
            channel = obj.name_channel_map(channel_name);
            for i=1:length(parent_channel_names)
                parent_channels(i) = obj.name_channel_map(parent_channel_names{i});
            end
            
            [~, num_trajectories] = size(obj.trajectories);
            for i=1:num_trajectories
                
                current_trajectory = obj.trajectories(:,i);
                current_trajectory_times = obj.trajectory_times(:,i);
                current_trajectory_channel_values = obj.discretized_data(current_trajectory,channel);
                current_trajectory_parent_values = obj.discretized_data(current_trajectory, parent_channels);
                
                %select the indices where the parents are equal to the
                %parent values
                parent_indices = [];
                for j=1:length(current_trajectory_parent_values)
                    if(isequal(current_trajectory_parent_values(j,:), parent_channel_values))
                        parent_indices = [parent_indices j];
                    end
                end
                
                
                
                
                for j=2:length(parent_indices)
                    
                    current_index = parent_indices(j);
                    prev_index = parent_indices(j-1);
                    
                    if(current_index-prev_index == 1)
                        %u has stayed steady
                        
                        val1 = current_trajectory_channel_values(prev_index);
                        val2 = current_trajectory_channel_values(current_index);
                        if((val1==channel_value_from)&(val2==channel_value_to))
                            
                            num_transitions = num_transitions+1;
                            
                        end
                        
                    end
                end
                
            end
            
        end
        
        function obj = discretize_channels_GMM(obj, channel_names)
            
            [rows, columns] = size(obj.data);
            obj.discretized_data = zeros(rows,columns);
            
            %copy over teh discretization files? thats prolly easiest.
            bw=1/3;
            peak_sens=20;
            shoulder_sens=50;
            peak_sep_thr=0.1;
            shoulder_sep_thr=0.1;
            
            obj = obj.findLevels(channel_names,bw,peak_sens,shoulder_sens,peak_sep_thr,shoulder_sep_thr);%find modes / protein levels
            obj = obj.localGMMpartition(channel_names);%calculate partitions with local GMM fits between pairs of modes/levels
            
            for i=1:length(channel_names)
                
                channel = obj.name_channel_map(channel_names{i});
                channel_partitions = obj.channel_levels{channel}{4};
                prev_partition = min(obj.data(:,channel));
                num_partitions = length(channel_partitions);
                channel_data = obj.data(:,channel);
                channel_partitions(num_partitions+1) = max(channel_data);
                
                for j=1:length(channel_partitions)
                    
                    current_partition = channel_partitions(j);
                    level_indices = find((prev_partition<channel_data)&(channel_data<current_partition));
                    obj.discretized_data(level_indices,channel) = j;
                    prev_partition = current_partition;
                end
                
                
            end
            
            
            
        end
        
        
        function obj=findLevels(obj,sel_channels,bw,peak_sens,shoulder_sens,peak_sep_thr,shoulder_sep_thr)
            
            n_x=100;
            high_protein_level=6;
            
            %%%arrange data
            val_time_points=unique([min(obj.time_points),max(obj.time_points)]);
            n_time_points=length(val_time_points);
            output=cell(0);
            
            for channel_i=1:length(sel_channels)
                
                channel_name=sel_channels{channel_i};
                
                tempPeaks=[];
                peakDens=[];
                for i=1:n_time_points
                    temp_vals=obj.data(obj.time_points==val_time_points(i),obj.name_channel_map(channel_name));
                    [temp_dens,temp_x]=ksdensity(temp_vals,linspace(min(temp_vals),max(temp_vals),n_x),'bandwidth',bw);
                    [~,tmploc]=findpeaks(temp_dens,'MINPEAKHEIGHT',max(temp_dens)/peak_sens,'MINPEAKDISTANCE',ceil(n_x*high_protein_level*peak_sep_thr/max(temp_x)));
                    tempPeaks=[tempPeaks,temp_x(tmploc)];
                    peakDens=[peakDens,temp_dens(tmploc)];
                end
                [tempPeaks,tmpord]=sort(tempPeaks);
                if(length(tempPeaks)>1)
                    peakDens=peakDens(tmpord);
                    peak_sep=diff(tempPeaks);
                    groups=peak_sep<high_protein_level*peak_sep_thr;
                    groupStarts=find(diff(groups)==1)+1;
                    if(groups(1)==1)
                        groupStarts=[1,groupStarts];
                    end
                    groupEnds=find(diff(groups)==-1)+1;
                    if(groups(end)==1)
                        groupEnds=[groupEnds,length(tempPeaks)];
                    end
                    groupSingletons=setdiff([find(~groups),length(tempPeaks)],groupEnds);
                    groups=[groupStarts,groupSingletons;groupEnds,groupSingletons];
                    [~,ord]=sort(groups(2,:));
                    groups=groups(:,ord);
                    num_levels=size(groups,2);
                    channel_levels=nan(1,num_levels);
                    for i=1:num_levels
                        tempLocs=groups(1,i):groups(2,i);
                        [~,tempLoc]=max(peakDens(tempLocs));
                        channel_levels(i)=tempPeaks(tempLocs(tempLoc));
                    end
                    channel_levels=sort(channel_levels);
                else
                    channel_levels=tempPeaks;
                    num_levels=length(tempPeaks);
                end
                
                x_shoulder=[];
                x_shoulder_mag=[];
                for j=1:n_time_points
                    temp_vals=obj.data(obj.time_points==val_time_points(j),obj.name_channel_map(channel_name));
                    [temp_dens,temp_x]=ksdensity(temp_vals,linspace(min(temp_vals),max(temp_vals),n_x),'bandwidth',bw);
                    temp_dens_norm=temp_dens*mean(diff(temp_x));
                    deriv2=diff(diff(temp_dens));
                    [mags,loc]=findpeaks(-deriv2,'MINPEAKHEIGHT',max(abs(deriv2))/shoulder_sens,'MINPEAKDISTANCE',ceil(n_x*high_protein_level*shoulder_sep_thr/max(temp_x)));
                    x_shoulder_temp=temp_x(loc+1);
                    x_shoulder_dens=temp_dens(loc+1);
                    keep_loc_j=x_shoulder_dens>max(temp_dens)/peak_sens;
                    x_shoulder=[x_shoulder,x_shoulder_temp(keep_loc_j)];
                    x_shoulder_mag=[x_shoulder_mag,mags(keep_loc_j)];
                end
                keepLoc=all(abs(ones(num_levels,1)*x_shoulder-channel_levels'*ones(1,length(x_shoulder)))>high_protein_level*shoulder_sep_thr,1);
                x_shoulder=x_shoulder(keepLoc);
                x_shoulder_mag=x_shoulder_mag(keepLoc);
                [x_shoulder,tmpord]=sort(x_shoulder);
                if(length(x_shoulder)>1)
                    x_shoulder_mag=x_shoulder_mag(tmpord);
                    shoulder_sep=diff(x_shoulder);
                    groups=shoulder_sep<high_protein_level*shoulder_sep_thr;
                    groupStarts=find(diff(groups)==1)+1;
                    if(groups(1)==1)
                        groupStarts=[1,groupStarts];
                    end
                    groupEnds=find(diff(groups)==-1)+1;
                    if(groups(end)==1)
                        groupEnds=[groupEnds,length(x_shoulder)];
                    end
                    groupSingletons=setdiff([find(~groups),length(x_shoulder)],groupEnds);
                    groups=[groupStarts,groupSingletons;groupEnds,groupSingletons];
                    [~,ord]=sort(groups(2,:));
                    groups=groups(:,ord);
                    num_shoulders=size(groups,2);
                    shoulders=nan(1,num_shoulders);
                    for i=1:num_shoulders
                        tempLocs=groups(1,i):groups(2,i);
                        [~,tempLoc]=max(x_shoulder_mag(tempLocs));
                        shoulders(i)=x_shoulder(tempLocs(tempLoc));
                    end
                    shoulders=sort(shoulders);
                else
                    shoulders=x_shoulder;
                end
                if(~isempty(shoulders))
                    channel_levels=[channel_levels,shoulders];
                    channel_levels=sort(channel_levels);
                    num_levels=length(channel_levels);
                end
                
                if(~any(channel_levels<high_protein_level*peak_sep_thr))
                    channel_levels=[0,channel_levels];
                    channel_levels=sort(channel_levels);
                    num_levels=length(channel_levels);
                end
                
                output{obj.name_channel_map(channel_name)}{1}=channel_levels;
                output{obj.name_channel_map(channel_name)}{2}=channel_name;
                
            end
            
            obj.channel_levels=output;
            
        end
        
        function obj=midPartitionSet(obj,sel_channels)
            
            for channel_i=1:length(sel_channels)
                
                channel_name=sel_channels{channel_i};
                channel_loc=obj.name_channel_map(channel_name);
                
                mu=obj.channel_levels{channel_loc}{1};
                K=length(mu);
                
                partitions=[];
                if(K>1)
                    for k=1:K-1
                        partitions(k)=mean(mu(k:k+1));
                    end
                end
                
                obj.channel_levels{channel_loc}{4}=partitions;
                
            end
            
        end
        
        function obj=GMMvarFit(obj,sel_channels)
            
            sprintf('GMM variance fitting')
            
            max_iter=1000;
            tol=1e-6;%tolerance for negative log-likelihood (stop when negative log-likelihood changes by less than this amount)
            
            %%%arrange data
            val_time_points=unique([min(obj.time_points),max(obj.time_points)]);;
            n_time_points=length(val_time_points);
            data_size=nan(n_time_points,1);
            for i=1:n_time_points
                data_size(i)=sum(obj.time_points==val_time_points(i));
            end
            N=sum(data_size);
            
            for channel_i=1:length(sel_channels)
                
                channel_name=sel_channels{channel_i}
                channel_loc=obj.name_channel_map(channel_name);
                
                temp_data=cell(0);
                for i=1:n_time_points
                    temp_data{i}=obj.data(obj.time_points==val_time_points(i),obj.name_channel_map(channel_name));
                end
                data_max=max(obj.data(:,channel_loc));
                data_min=min(obj.data(:,channel_loc));
                data_var=var(obj.data(:,channel_loc));
                
                mu=obj.channel_levels{channel_loc}{1};
                K=length(mu);
                
                %%%initialise parameters
                alpha=(1/K)*ones(n_time_points,K);
                variances=data_var*ones(1,K);
                W=cell(0);
                for i=1:n_time_points
                    W{i}=nan(data_size(i),K);
                end
                Nk=nan(size(alpha));
                n_l_ell=Inf;
                i_iter=0;
                mark=0;
                n_l_ell_all=[];
                var_all=[];
                
                while(mark==0)
                    
                    i_iter=i_iter+1;
                    
                    %%%E-step
                    for i=1:n_time_points
                        for k=1:K
                            W{i}(:,k)=alpha(i,k)*pdf('norm',temp_data{i},mu(k),sqrt(variances(k)));
                        end
                        W{i}=W{i}./(sum(W{i},2)*ones(1,K));
                    end
                    
                    %%%M-step
                    for i=1:n_time_points
                        Nk(i,:)=sum(W{i},1);
                        alpha(i,:)=Nk(i,:)/data_size(i);
                    end
                    variances=zeros(1,K);
                    for i=1:n_time_points
                        variances=variances+(data_size(i)/N)*sum(W{i}.*(temp_data{i}*ones(1,K)-ones(data_size(i),1)*mu).^2,1)./Nk(i,:);
                    end
                    
                    %%%negative log-likelihood
                    n_l_ell_new=0;
                    for i=1:n_time_points
                        temp_ell=zeros(data_size(i),1);
                        for k=1:K
                            temp_ell=temp_ell+alpha(i,k)*pdf('norm',temp_data{i},mu(k),sqrt(variances(k)));
                        end
                        temp_ell=sum(log(temp_ell));
                        n_l_ell_new=n_l_ell_new-temp_ell;
                    end
                    clear temp_ell
                    mark=abs(n_l_ell-n_l_ell_new)<tol|i_iter==max_iter;
                    n_l_ell_all(i_iter)=n_l_ell_new;
                    var_all(i_iter,:)=variances;
                    n_l_ell=n_l_ell_new;
                    
                end
                
                if(i_iter==max_iter)
                    warning(['Reached maximum number of iterations = ',num2str(max_iter),' without convergence at tolerance = ',num2str(tol)])
                end
                
                [n_l_ell,loc]=min(n_l_ell_all);
                variances=var_all(loc,:);
                
                obj.channel_levels{channel_loc}{3}=variances;
                
            end
            
        end
        
        function obj=GMMpartitionSet(obj,sel_channels)
            
            %%%arrange data
            val_time_points=unique(obj.time_points);
            n_time_points=length(val_time_points);
            data_size=nan(n_time_points,1);
            for i=1:n_time_points
                data_size(i)=sum(obj.time_points==val_time_points(i));
            end
            N=sum(data_size);
            
            for channel_i=1:length(sel_channels)
                
                channel_name=sel_channels{channel_i};
                channel_loc=obj.name_channel_map(channel_name);
                
                mu=obj.channel_levels{channel_loc}{1};
                K=length(mu);
                variance=obj.channel_levels{channel_loc}{3};
                
                partitions=[];
                if(K>1)
                    for k=1:K-1
                        x_div=ceil(1000*min(mu(k:k+1)))/1000:1/1000:floor(1000*max(mu(k:k+1)))/1000;
                        probs=[pdf('norm',x_div,mu(k),sqrt(variance(k)));pdf('norm',x_div,mu(k+1),sqrt(variance(k+1)))];
                        [~,tmploc]=min(abs(diff(probs,[],1)));
                        partitions(k)=x_div(tmploc);
                        %                     plot(x_div,probs')
                        %                     hold on
                        %                     plot(ones(1,2)*partition(k),[0,max(probs(:))],'k')
                        %                     hold off
                    end
                end
                
                obj.channel_levels{channel_loc}{4}=partitions;
                
            end
            
        end
        
        function obj=localGMMpartition(obj,sel_channels)
            
            sprintf('local GMM fitting')
            
            max_iter=1000;
            tol=1e-3;%tolerance for negative log-likelihood (stop when negative log-likelihood changes by less than this amount)
            
            %%%arrange data
            val_time_points=unique([min(obj.time_points),max(obj.time_points)]);;
            n_time_points=length(val_time_points);
            
            for channel_i=1:length(sel_channels)
                
                channel_name=sel_channels{channel_i}
                channel_loc=obj.name_channel_map(channel_name);
                
                mu_all=obj.channel_levels{channel_loc}{1};
                n_comp=length(mu_all);
                K=2;
                partitions=[];
                
                if(n_comp>1)
                    
                    temp_data=cell(0);
                    for i=1:n_time_points
                        temp_data{i}=obj.data(obj.time_points==val_time_points(i),obj.name_channel_map(channel_name));
                    end
                    
                    for i_comp=1:n_comp-1
                        
                        mu_temp=mu_all(i_comp:i_comp+1);
                        fit_data=cell(0);
                        data_size=nan(n_time_points,1);
                        variances=zeros(1,2);
                        for i=1:n_time_points
                            fit_data{i}=temp_data{i}(temp_data{i}>mu_temp(1)&temp_data{i}<mu_temp(2));
                            data_size(i)=length(fit_data{i});
                            variances=variances+sum((fit_data{i}*ones(1,2)-ones(data_size(i),1)*mu_temp).^2,1);
                        end
                        N=sum(data_size);
                        variances=variances/N;
                        
                        %%%initialise parameters
                        alpha=(1/K)*ones(n_time_points,K);
                        W=cell(0);
                        for i=1:n_time_points
                            W{i}=nan(data_size(i),K);
                        end
                        Nk=nan(size(alpha));
                        n_l_ell=Inf;
                        i_iter=0;
                        mark=0;
                        n_l_ell_all=[];
                        var_all=[];
                        
                        while(mark==0)
                            
                            i_iter=i_iter+1;
                            
                            %%%E-step
                            for i=1:n_time_points
                                for k=1:K
                                    W{i}(:,k)=alpha(i,k)*pdf('norm',fit_data{i},mu_temp(k),sqrt(variances(k)));
                                end
                                W{i}=W{i}./(sum(W{i},2)*ones(1,K));
                            end
                            
                            %%%M-step
                            for i=1:n_time_points
                                Nk(i,:)=sum(W{i},1);
                                alpha(i,:)=Nk(i,:)/data_size(i);
                            end
                            variances=zeros(1,K);
                            for i=1:n_time_points
                                variances=variances+(data_size(i)/N)*sum(W{i}.*(fit_data{i}*ones(1,K)-ones(data_size(i),1)*mu_temp).^2,1)./Nk(i,:);
                            end
                            
                            %%%negative log-likelihood
                            n_l_ell_new=0;
                            for i=1:n_time_points
                                temp_ell=zeros(data_size(i),1);
                                for k=1:K
                                    temp_ell=temp_ell+alpha(i,k)*pdf('norm',fit_data{i},mu_temp(k),sqrt(variances(k)));
                                end
                                temp_ell=sum(log(temp_ell));
                                n_l_ell_new=n_l_ell_new-temp_ell;
                            end
                            clear temp_ell
                            mark=abs(n_l_ell-n_l_ell_new)<tol|i_iter==max_iter;
                            n_l_ell_all(i_iter)=n_l_ell_new;
                            var_all(i_iter,:)=variances;
                            n_l_ell=n_l_ell_new;
                            
                        end
                        
                        if(i_iter==max_iter)
                            warning(['Reached maximum number of iterations = ',num2str(max_iter),' without convergence at tolerance = ',num2str(tol)])
                        end
                        
                        [n_l_ell,loc]=min(n_l_ell_all);
                        variances=var_all(loc,:);
                        
                        %%find partition
                        x_div=ceil(1000*min(mu_temp(1)))/1000:1/1000:floor(1000*max(mu_temp(2)))/1000;
                        probs=[pdf('norm',x_div,mu_temp(1),sqrt(variances(1)));pdf('norm',x_div,mu_temp(2),sqrt(variances(2)))];
                        [~,tmploc]=min(abs(diff(probs,[],1)));
                        partitions(i_comp)=x_div(tmploc);
                        
                    end
                    
                end
                
                obj.channel_levels{channel_loc}{4}=partitions;
                
            end
            
        end
        
        
        function obj = bin_wanderlust(obj, num_bins, num_traj)
            
            
            trajectories = zeros(num_bins, num_traj);
            %bin by time?
            trajectory_times = zeros(num_bins, num_traj);
            
            
            wmax = length(obj.wanderlust_channel);
            
            bin_width = floor((wmax)/num_bins);
            bin_ends = ((1:num_bins).*bin_width);
            bin_ends = [1 bin_ends wmax];
            
            for i=2:length(bin_ends)
                
                
                bin_indices = bin_ends(i-1):bin_ends(i);
                trajectories(i-1,:) = randsample(bin_indices,num_traj, true);
                current_indices = trajectories(i-1,:);
                trajectory_times(i-1,:) = obj.wanderlust_channel(current_indices);
                
            end
            
            %later do a more sophisticated matching if need be.
            
            obj.trajectories = trajectories;
            obj.trajectory_times = trajectory_times;
            
            %ok this stories the trajectory in order so the 1st value ahs
            %the index of the first element and its wanderlust channel has
            %the value. this is an ok way of storing it.
        end
        
        
        function obj = run_and_average_wanderlust(obj, sel_channels, n_resample, start_inds)
            
            %later make it so that we specify better starting indices, user specifies them for now.
            
            
            time_pred_all = run_wanderlust(obj, sel_channels, start_inds, n_resample);
            
            
            wanderlust_channel = sum(time_pred_all,2)./n_resample;
            wanderlust_channel = wanderlust_channel./max(wanderlust_channel);
            
            
            
            obj.wanderlust_channel = wanderlust_channel;
        end
        
        function [smooth_data, smooth_index ] = plot_events_vs_time(obj, channel_name, window_size, varargin)
            
            channel =  obj.name_channel_map(channel_name);
            optargin = size(varargin,2);
            num_events = size(obj.data,1);
            
            wanderlust_channel = obj.name_channel_map('cycler');
            
            if(optargin > 0)
                
                time_data = varargin{1};
            else
                time_data = obj.data(:,wanderlust_channel);
            end
            
            [ wanderlust_value , wanderlust_order] = sort(time_data);
            time_data = time_data(wanderlust_order);
            event_data = obj.data(wanderlust_order,channel);
            
            
            smooth_data = zeros(1,num_events-window_size);
            smooth_index = zeros(1,num_events-window_size);
            
            for i=1:num_events-window_size
                
                smooth_data(i) = median(event_data(i:i+window_size));
                %smooth_index(i) = median(event_data(i:i+window_size,1));
                %ordering_time = varargin{1};
                smooth_index(i) = median(time_data(i:i+window_size));
            end
            
            %h=plot(smooth_index, smooth_data, 'LineWidth',2.0);
            %set(h,'Color',varargin{2});
            
            plot(smooth_index, smooth_data,'LineWidth',2.0);
            
            %set(gca,'XTick',[]);
            xlabel('Wanderlust Order','FontSize',14);
            ystring = sprintf('%s level',channel_name);
            ylabel(ystring, 'FontSize',14);
            %scatter(event_data(:,1), event_data(:,channel));
            
        end
        
        function plot_bead_events(obj, channel, window_size)
            
            events = find(obj.beads ==1);
            obj.plot_events_vs_time(channel, window_size, events);
            
        end
        
        function obj = perform_tcell_gate(obj, thresh)
            
            %cd4 is 14
            %cd8 is 15
            %cd7 is 31
            
            [num_events,~] = size(obj.data);
            obj.event_gate = zeros(1,num_events);
            cd4_channel = 14;
            cd8_channel = 15;
            
            
            %done computing all non-beads
            for i=1:num_events
                
                
                
                if(obj.data(i,cd4_channel)>thresh)
                    obj.event_gate(i) = 1;
                    
                elseif(obj.data(i,cd8_channel)>thresh)
                    obj.event_gate(i) = 1;
                end
                
                
            end
            
        end
        
        function obj = perform_bcell_gate(obj, thresh)
            
            %cd20 is 18
            %cd19 is 12
            %IGM-B is 22
            
            [num_events,~] = size(obj.data);
            obj.event_gate = zeros(1,num_events);
            cd20_channel = 18;
            
            
            %done computing all non-beads
            for i=1:num_events
                
                
                if(obj.data(i,cd20_channel)>thresh)
                    obj.event_gate(i) = 1;
                end
                
                
            end
            
            
        end
        
        
        
        
        function [XT, YT] = interpolate_events(obj, total_events, channel_name, mode)
            
            channel = obj.name_channel_map(channel_name);
            [num_events,~] = size(obj.data);
            
            XI = zeros(total_events,1);
            
            for i=1:total_events
                
                XI(i) = i;
                
            end
            
            indices = obj.data(:,obj.eventnum_channel);
            XI(indices)=[];
            
            YI = interp1(obj.data(:,obj.eventnum_channel),obj.data(:,channel),XI, mode);
            
            XT = zeros(total_events, 1);
            YT = zeros(total_events, 1);
            
            for i=1:num_events
                
                xval = obj.data(i,obj.eventnum_channel);
                yval = obj.data(i, channel);
                
                XT(xval) = xval;
                YT(xval) = yval;
                
            end
            
            for i=1:length(XI)
                
                
                
                XT(XI(i)) = XI(i);
                YT(XI(i)) = YI(i);
                
            end
            
            
            
        end
        
        function well_matrix=plot_channel_mean_plate_effect(obj, channel_name)
            
            well_matrix = zeros(8,12);
            channel = obj.name_channel_map(channel_name);
            
            marker_channel = find(obj.marker_channels==channel);
            %caspase_z = zscore(cdata.well_means(:,21));
            if(size(marker_channel,2)==0)
                return;
                %not a marker channel
            end
            
            for i=1:8
                for j=1:12
                    well_no = obj.get_well_no(i,j);
                    well_matrix(i,j) = obj.well_means(well_no,marker_channel);
                end
            end
            
            imagesc(well_matrix);
            
            
        end
        
        function [sorted_wells]=plot_effects_channel_metric(obj, channel_name, metric_name, num_bins)
            
            
            channel = obj.name_channel_map(channel_name);
            %histone
            
            var_data = zeros(1,96);
            
            for i=1:96
                
                effect =  obj.well_siRNA_effects(i,channel);
                
                if(strcmp(metric_name, 'var'))
                    
                    var_data(i) = effect{1}.var_diff;
                    
                elseif(strcmp(metric_name, 'mean'))
                    
                    var_data(i) = effect{1}.mean_diff;
                    
                elseif(strcmp(metric_name, 'L2'))
                    
                    var_data(i) = effect{1}.L2;
                else
                    
                    var_data(i) = effect{1}.knockdown;
                    
                end
                
            end
            
            [X,sorted_wells] = sort(var_data, 'descend');
            
            hist(var_data, num_bins);
            
        end
        
        
        function [well_events, populations, fractions] = seperate_well_populations(obj, well_no, num_clusters, channel_name)
            
            %Gausian mixture models
            
            clear X;
            
            [~,num_barcodes] = size(obj.barcoding_channels);
            
            well_data = obj.get_well_data(well_no);
            size(well_data)
            channel = obj.name_channel_map(channel_name);
            
            well_events = zeros(1,size(well_data,1));
            populations = zeros(1,size(well_data,1));
            fractions = zeros(1,num_clusters);
            
            [~, M, V, ~] = EM_GM_fast(well_data(:,channel),num_clusters);
            [~,Inds] = sort(M);
            
            for i=1:size(well_data,1)
                
                prob = zeros(1,num_clusters);
                
                
                for j=1:num_clusters
                    
                    prob(Inds(j)) = normpdf(well_data(i,channel), M(Inds(j)), V(1,1,Inds(j)));
                    
                end
                
                
                well_events(i) = well_data(i,channel);
                [~,Ind] = max(prob);
                populations(i) = Ind(1);
                
            end
            
            for i=1:num_clusters
                
                fractions(i) = length(find(populations==i))/length(populations);
                
            end
            
        end
        
        function fraction_matrix = get_responding_populations(obj, num_clusters, responding_cluster, channel_name)
            
            fraction_matrix = zeros(8,12);
            
            for i=1:8
                for j=1:12
                    well_no = obj.get_well_no(i,j);
                    
                    [~,~,fractions] = obj.seperate_well_populations(well_no, num_clusters, channel_name);
                    fraction_matrix(i,j) = fractions(responding_cluster);
                    
                end
            end
            
        end
        
        %FUNCTIONS THAT ORIGINALLY USED TO BE IN THE BEAD NORMALIZATION CLASS
        
        
        
        function obj = compute_smooth_events_window_opt(obj, corr_threshold)
            
            current_min_corr = 0;
            window_size = 50;
            step_size = 50;
            [num_samples,~] = size(obj.bead_data);
            
            while current_min_corr<corr_threshold
                
                window_size = window_size+step_size
                
                if(window_size> num_samples)
                    sprintf('CANNOT FIND WINDOW SIZE');
                    break;
                end
                
                obj = obj.compute_smooth_events_vs_time(window_size);
                obj = obj.compute_smooth_corrs();
                obj.smooth_corrs
                
                min_corr = min(obj.smooth_corrs);
                current_min_corr = min(min_corr)
                
            end
            
            window_size
            
            
        end
        
        
        function obj = compute_smooth_events_vs_time(obj, window_size)
            
            
            
            event_data = obj.bead_data;
            
            num_events = size(event_data, 1);
            num_channels = size(event_data, 2);
            obj.smooth_data = zeros(num_events-window_size, num_channels-1);
            obj.smooth_index = zeros(num_events-window_size, 1);
            
            for j=2:num_channels
                
                for i=1:num_events-window_size
                    
                    obj.smooth_data(i,j-1) = median(event_data(i:i+window_size, j));
                    obj.smooth_index(i) = median(event_data(i:i+window_size, 1));
                    
                end
                
            end
            
            
            %irregularly sampled with respect to time,
            
            
        end
        
        function obj = compute_acquisition_rates(obj, window_size)
            
            num_events = length(obj.acquisition_times);
            obj.acquisition_rate = zeros(1, num_events-window_size);
            obj.acquisition_rate_index = zeros(1, num_events-window_size);
            half_window = window_size/2;
            
            
            for i=1:num_events-window_size
                start_time = obj.acquisition_times(i);
                end_time = obj.acquisition_times(i+window_size);
                duration = end_time - start_time;
                if(duration == 0)
                    duration = .00000001;
                end
                rate = window_size/duration;
                %number of acquisitions per millisecond
                
                obj.acquisition_rate(i) = rate;
                obj.acquisition_rate_index(i) = obj.acquisition_times(i+half_window);
            end
            
        end
        
        function plot_acquisition_rate(obj, varargin)
            
            optargin = size(varargin, 2);
            
            if(optargin == 0)
                plot(obj.acquisition_rate_index, obj.acquisition_rate , 'r');
            else
                plot(obj.acquisition_rate_index, obj.acquisition_rate , varargin{1});
            end
            
        end
        
        function plot_smooth_events(obj, channel_name, varargin)
            
            optargin = size(varargin, 2);
            channel = obj.name_channel_map(channel_name);
            if(optargin == 0)
                plot(obj.smooth_index, obj.smooth_data(:,channel));
            else
                plot(obj.smooth_index, obj.smooth_data(:, channel), varargin{1});
            end
            
        end
        
        function obj = compute_smooth_corrs(obj)
            
            num_channels = size(obj.smooth_data,2);
            obj.smooth_corrs = zeros(num_channels, num_channels);
            
            for i=1:num_channels
                for j=1:num_channels
                    obj.smooth_corrs(i,j) = corr(obj.smooth_data(:,i), obj.smooth_data(:,j));
                end
            end
            
        end
        
        function obj = compute_corrections(obj)
            
            num_data_points = size(obj.smooth_data,1);
            num_channels = size(obj.smooth_data,2);
            obj.smooth_diffs = zeros(num_data_points-1, num_channels);
            obj.smooth_data_recon = zeros(num_data_points, num_channels);
            
            obj.smooth_fluctuation_vector = zeros(num_data_points-1, 1);
            
            for i=1:num_channels
                
                obj.smooth_diffs(:,i) = diff(obj.smooth_data(:,i));
                obj.smooth_fluctuation_vector = obj.smooth_fluctuation_vector + obj.smooth_diffs(:,i);
            end
            
            obj.smooth_fluctuation_vector = obj.smooth_fluctuation_vector .* (1/num_channels);
            %average the difference trends
            
            for i=1:num_channels
                
                obj.smooth_data_recon(1,i) = obj.smooth_data(1,i);
                
                for j = 1:(num_data_points-1)
                    obj.smooth_data_recon(j+1,i) = obj.smooth_data_recon(j,i)+obj.smooth_fluctuation_vector(j);
                end
                
            end
            
        end
        
        function [correction_function] = create_correction(obj, base_level, acquisition_times)
            
            num_data_points = length(obj.acquisition_times);
            correction_function = zeros(num_data_points, 1);
            
            smooth_data_points = length(obj.smooth_index);
            
            smooth_data_recon = zeros(smooth_data_points, 1);
            smooth_data_recon(1) = base_level;
            
            for i = 1:(smooth_data_points-1)
                
                smooth_data_recon(i+1) = smooth_data_recon(i) + obj.smooth_fluctuation_vector(i);
                
            end
            
            
            
            first_zero = -1;
            last_zero = 0;
            
            for i=2:length(acquisition_times)
                
                if(acquisition_times(i) == acquisition_times(i-1))
                    
                    
                    if(first_zero == -1)
                        first_zero = i-1;
                        
                    end
                    
                else
                    
                    if(first_zero > -1)
                        
                        last_zero = i-1;
                        
                        increments = (last_zero - first_zero)+1;
                        if(increments==0)
                            
                            sprintf('WARNING! Increments is 0 \n');
                        end
                        step = 1/increments;
                        sprintf('fixing acquisition time\n');
                        
                        for j=(first_zero+1):last_zero
                            
                            acquisition_times(j) = acquisition_times(j) + step;
                            step = step+(1/increments);
                            
                        end
                        
                        first_zero = -1;
                        last_zero = 0;
                        
                    end
                    
                end
                
            end
            
            
            
            first_zero = -1;
            last_zero = 0;
            smooth_index = obj.smooth_index;
            
            for i=2:length(smooth_index)
                
                if(smooth_index(i) == smooth_index(i-1))
                    
                    
                    if(first_zero == -1)
                        first_zero = i-1;
                        
                    end
                    
                else
                    
                    if(first_zero > -1)
                        
                        last_zero = i-1;
                        
                        increments = (last_zero - first_zero)+1;
                        step = 1/increments;
                        sprintf('fixing acquisition time\n');
                        
                        for j=(first_zero+1):last_zero
                            
                            smooth_index(j) = smooth_index(j) + step;
                            step = step+(1/increments);
                            
                        end
                        
                        first_zero = -1;
                        last_zero = 0;
                        
                    end
                    
                end
                
            end
            
            %if the last two timestamps are the same the above loop doesn't
            %fix them, so adding this here -rlf 20110811
            if(first_zero > -1)
                
                last_zero = i;
                
                increments = (last_zero - first_zero)+1;
                step = 1/increments;
                sprintf('fixing acquisition time\n');
                
                for j=(first_zero+1):last_zero
                    
                    smooth_index(j) = smooth_index(j) + step;
                    step = step+(1/increments);
                    
                end
                
            end
            
            
            full_data_recon = interp1(smooth_index, smooth_data_recon, acquisition_times, 'linear', 'extrap');
            
            
            correction_function = base_level - full_data_recon;
            
            
            
            
        end
        
        function obj = correct_channels(obj, sensitivity_threshold)
            %trying to return a modified version of this object
            
            num_channels = size(obj.data,2)
            %correcting all channels besides time and cell length
            
            
            
            
            for i = 3:num_channels
                
                channel = i;
                
                %obj.channel_averages = obj.channel_averages .* (1/length(obj.acquisition_times));
                
                base_level = median(obj.data(:,channel));
                
                correction_function = obj.create_correction(base_level, obj.data(:,1));
                
                low_sensitivity_indices = find(obj.data(:,channel)<sensitivity_threshold);
                
                correction_function(low_sensitivity_indices)=0;
                
                obj.data(:,channel) = obj.data(:,channel) + correction_function;
                
            end
            
            %             num_channels = length(obj.dna_channels)
            %             %marker channels, dna channels also
            %
            %
            %             for i = 1:num_channels
            %
            %                channel = obj.dna_channels(i);
            %
            %                %obj.channel_averages = obj.channel_averages .* (1/length(obj.acquisition_times));
            %
            %                base_level = median(obj.data(:,channel));
            %
            %                correction_function = obj.create_correction(base_level, obj.data(:,1));
            %
            %                obj.data(:,channel) = obj.data(:,channel) + correction_function;
            %
            %             end
            
            
            
        end
        
        function obj = find_bead_data(obj)
            
            %obj = obj.identify_beads();
            %can add that back for automation purposes
            
            obj.acquisition_times = obj.data(:,1);
            
            num_beads = length(find(obj.beads==1))
            num_non_beads = length(find(obj.beads==0))
            
            bead_events = find(obj.beads==1);
            all_channels = [1 obj.bead_channels];
            obj.bead_data = obj.data(bead_events,all_channels);
            
            
        end
        
        function obj = bead_normalize(obj, corr_threshold, approx_percentage, sensitivity_threshold)
            
            
            obj = obj.identify_beads();
            sprintf('done identifying beads \n');
            
            obj = obj.find_bead_data();
            
            sprintf('seperated bead data \n');
            
            obj = obj.compute_smooth_events_window_opt(corr_threshold);
            
            sprintf('seperated smooth events \n');
            
            obj = obj.compute_corrections();
            
            sprintf('computed corrections \n');
            
            obj = obj.correct_channels(sensitivity_threshold);
            
            sprintf('corrected channels \n');
            
        end
        
        function obj = bead_normalize_fix_window(obj, window_size, sensitivity_threshold)
            
            
            %obj = obj.identify_beads();
            %sprintf('done identifying beads \n');
            
            
            
            obj = obj.find_bead_data();
            
            sprintf('seperated bead data \n')
            
            obj = obj.compute_smooth_events_vs_time(window_size);
            
            sprintf('seperated smooth events \n')
            
            obj = obj.compute_corrections();
            
            sprintf('computed corrections \n')
            obj = obj.correct_channels(sensitivity_threshold);
            
            sprintf('corrected channels \n')
            %obj = obj.remove_beads();
        end
        
        function [COMTY] = well_louvain_clusters(obj, well_no, k_neighbors)
            
            well_data = obj.get_well_data(well_no);
            marker_data = well_data(:, obj.marker_channels);
            [num_events, num_markers] = size(marker_data);
            distance_matrix = zeros(num_events,num_events);
            
            for i=1:num_events
                for j=1:num_events
                    
                    if(j>=i)
                        continue;
                    end
                    
                    distance_matrix(i,j) = norm(marker_data(i,:)-marker_data(j,:));
                    distance_matrix(j,i) = distance_matrix(i,j);
                end
            end
            
            
            [r,c] = size(distance_matrix);
            
            
            [~,I] = sort(distance_matrix,2,'ascend');
            
            adjacency_matrix = zeros(num_events,num_events);
            
            
            
            
            for i=1:r
                
                for j=1:k_neighbors
                    %because 1 will be itself with distance 0
                    edge_sink = I(i,j+1);
                    adjacency_matrix(i,edge_sink) = 1;
                    
                end
            end
            
            [COMTY ending] = cluster_jl(adjacency_matrix,1,0,1,1);
            
        end
        
        function obj = plate_louvain_clusters(obj)
            
            obj.well_louvain_communities = cell(1,96);
            for i=1:96
                
                [COMTY] = well_louvain_cluster(obj, well_no, k_neighbors);
                obj.well_louvain_communities{i} = COMTY;
                
            end
            
        end
        
        
        
        
    end
    
end



