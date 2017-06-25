function [points_x, points_y, point_weights, density, normalized_density, xaxis,yaxis] = pairwise_visualize_matrix_min(data, varargin)            

            X = data(:, 1);
            Y = data(:, 2);
            total_cells = size(data,1);
            
            
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
                
            end
            
            %cutoff
            
            
            
            
            
            
            if(fix_limits == 0)
                
                [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[minxval minyval],[maxxval maxyval]);
                
            else
                
                [bandwidth,density,Grid_X,Grid_Y]=kde2d([X Y],num_slices,[minx miny],[maxx maxy]);
                
                
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
            matrix_to_plot = [smoothed_normalized_density side_bar];
            top_bar = [top_bar corner];
            matrix_to_plot = [matrix_to_plot; top_bar];
            
            xaxis_to_plot = [xaxis xaxis_side_bar];
            yaxis_to_plot = [yaxis; yaxis_top_bar];
            
            
            
            
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
                    %imagesc(xaxis, yaxis, smoothed_normalized_density);
                    imagesc(xaxis_to_plot,yaxis_to_plot, matrix_to_plot);
                end
                
                
                set(gca,'YDir','normal');
                %      set(gca,'XTick',[]);
                %      set(gca,'YTick',[]);
                %     set(gca, 'XTickLabel','');
                %     set(gca, 'YTickLabel','');
                
                
                %         set(gca, 'FontSize',16);
                %         xlabel(channel1_name);
                %         ylabel(channel2_name);
%                xlabel(channel1_name, 'FontSize', 16)
%                ylabel(channel2_name, 'FontSize', 16)
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